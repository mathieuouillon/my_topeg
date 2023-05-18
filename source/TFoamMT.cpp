#include <TFoamMT.h>


////////////////////////////////////////////////////////////////////////////////
/// Internal subprogram TFoam::Initialize reworked for multithreading.
/// It explores newly defined cell with help of special short MC sampling.
/// As a result, estimates of true and drive volume is defined/determined
/// Average and dispersion of the weight distribution will is found along
/// each edge and the best edge (minimum dispersion, best maximum weight)
/// is memorized for future use.
/// The optimal division point for eventual future cell division is
/// determined/recorded. Recorded are also minimum and maximum weight etc.
/// The volume estimate in all (inactive) parent cells is updated.
/// Note that links to parents and initial volume = 1/2 parent has to be
/// already defined prior to calling this routine.

void TFoamMT::Explore(TFoamCell *cell) {

    Double_t xBest = 0;
    Double_t yBest = 0;

    TFoamVect cellSize(fDim);
    TFoamVect cellPosi(fDim);

    cell->GetHcub(cellPosi, cellSize);
    cell->CalcVolume();
    Double_t intOld = cell->GetIntg();//memorize old values,
    Double_t driOld = cell->GetDriv();//will be needed for correcting parent cells


    /////////////////////////////////////////////////////
    //    Special Short MC sampling to probe cell      //
    /////////////////////////////////////////////////////
    fSum[0] = 0;
    fSum[1] = 0;
    fSum[2] = 0;
    fSum[3] = 1.0e150;
    fSum[4] = -1.0e150;
    for (int i = 0; i < fDim; i++) ((TH1D *) (*fHistEdg)[i])->Reset();// Reset histograms
    fHistWt->Reset();
    fEvMT  = 0;
    nevEff = 0.;

    for (int i = 1; i <= NbThreads; ++i) pool.push_task([&cell, i, this] { this->LocalMCLoop(cell, i); });
    pool.wait_for_tasks();

    //------------------------------------------------------------------
    //---  predefine logics of searching for the best division edge ---
    for (int k = 0; k < fDim; k++) {
        fMaskDiv[k] = 1;                      // default is all
        if (fInhiDiv[k] == 1) fMaskDiv[k] = 0;// inhibit some...
    }
    // Note that predefined division below overrule inhibition above
    Int_t kBest = -1;
    if (fOptPRD) {// quick check
        for (int k = 0; k < fDim; k++) {
            Double_t rmin = cellPosi[k];
            Double_t rmax = cellPosi[k] + cellSize[k];
            if (fXdivPRD[k] != nullptr) {
                Int_t n = (fXdivPRD[k])->GetDim();
                for (int j = 0; j < n; j++) {
                    Double_t rdiv = (*fXdivPRD[k])[j];
                    // check predefined divisions is available in this cell
                    if ((rmin + 1e-99 < rdiv) && (rdiv < rmax - 1e-99)) {
                        kBest = k;
                        xBest = (rdiv - cellPosi[k]) / cellSize[k];
                        goto ee05;
                    }
                }
            }
        }
    }
ee05:
    /////////////////////////////////////////////////////////////////////////////

    fNEffev += (Long_t) nevEff;
    Double_t nevMC   = fSum[2];
    Double_t intTrue = fSum[0] / (nevMC + 0.000001);
    Double_t intDriv = 0.;
    Double_t intPrim = 0.;

    switch (fOptDrive) {
        case 1:                                                       // VARIANCE REDUCTION
            if (kBest == -1) Varedu(fSum.data(), kBest, xBest, yBest);// determine the best edge,
            intDriv = sqrt(fSum[1] / nevMC) - intTrue;                // Foam build-up, sqrt(<w**2>) -<w>
            intPrim = sqrt(fSum[1] / nevMC);                          // MC gen. sqrt(<w**2>) =sqrt(<w>**2 +sigma**2)
            break;
        case 2:                                          // WTMAX  REDUCTION
            if (kBest == -1) Carver(kBest, xBest, yBest);// determine the best edge
            intDriv = fSum[4] - intTrue;                 // Foam build-up, wtmax-<w>
            intPrim = fSum[4];                           // MC generation, wtmax!
            break;
        default:
            Error("Explore", "Wrong fOptDrive = \n");
    }//switch
    //=================================================================================
    //hist_Neff_distrib.Fill( fLastCe/2.0+0.01, nevEff+0.01);  //
    //hist_kBest_distrib.Fill( kBest+0.50, 1.0 ); //  debug
    //hist_xBest_distrib.Fill( xBest+0.01, 1.0 ); //  debug
    //=================================================================================
    cell->SetBest(kBest);
    cell->SetXdiv(xBest);
    cell->SetIntg(intTrue);
    cell->SetDriv(intDriv);
    cell->SetPrim(intPrim);
    // correct/update integrals in all parent cells to the top of the tree
    // Double_t parIntg, parDriv;
    for (TFoamCell *parent = cell->GetPare(); parent != nullptr; parent = parent->GetPare()) {
        Double_t parIntg = parent->GetIntg();
        Double_t parDriv = parent->GetDriv();
        parent->SetIntg(parIntg + intTrue - intOld);
        parent->SetDriv(parDriv + intDriv - driOld);
    }
}


void TFoamMT::MakeAlpha() {
    if (fDim < 1) return;
    std::vector<Double_t> LocalRvec(fRNmax);
    // simply generate and load kDim uniform random numbers
    fPseRan->RndmArray(fDim, LocalRvec.data());// kDim random numbers needed
    for (Int_t k = 0; k < fDim; k++) fAlpha[k] = LocalRvec[k];
}


void TFoamMT::MakeLocalAlpha(std::vector<Double_t> &alpha) {
    if (alpha.empty()) return;
    fPseRan->RndmArray(static_cast<Int_t>(alpha.size()), alpha.data());
}

/**
 *
 * @param cell
 * @param ThNb Thread Number
 */
void TFoamMT::LocalMCLoop(TFoamCell *cell, Int_t ThNb) {
    std::vector<std::shared_ptr<TH1D>> LocalHistEdg(fDim);
    for (Int_t i = 0; i < fDim; i++) {
        TString hname   = fName + "_HistEdge_" + i + "_" + ThNb;
        TString htitle  = TString("Edge Histogram No. ") + i + "_" + ThNb;
        LocalHistEdg[i] = (std::make_shared<TH1D>(hname, htitle, fNBin, 0.0, 1.0));// Initialize histogram for each edge
        LocalHistEdg[i]->Sumw2();
    }
    TFoamVect cellSize(fDim);
    TFoamVect cellPosi(fDim);

    cell->GetHcub(cellPosi, cellSize);

    cell->CalcVolume();
    Double_t dx = cell->GetVolume();

    // ||||||||||||||||||||||||||BEGIN MC LOOP|||||||||||||||||||||||||||||
    while (fEvMT < fNSampl) {
        std::vector<Double_t> LocalAlpha(fDim);
        std::vector<Double_t> xRand(fDim);
        MakeLocalAlpha(LocalAlpha);// generate uniformly vector inside hypercube

        if (fDim > 0)
            for (Int_t j = 0; j < fDim; j++) xRand[j] = cellPosi[j] + LocalAlpha[j] * (cellSize[j]);

        Double_t wt = dx * Eval(xRand.data());

        Int_t nProj = 0;
        if (fDim > 0) {
            for (Int_t k = 0; k < fDim; k++) {
                Double_t xproj = LocalAlpha[k];
                LocalHistEdg[nProj]->Fill(xproj, wt);
                nProj++;
            }
        }

        MutexSum.lock();
        fEvMT++;
        fNCalls++;
        fSum[0] += wt;
        fSum[1] += wt * wt;
        fSum[2]++;
        if (fSum[3] > wt) fSum[3] = wt;
        if (fSum[4] < wt) fSum[4] = wt;
        nevEff = fSum[0] * fSum[0] / fSum[1];
        if (nevEff >= fNBin * fEvPerBin) {
            MutexSum.unlock();
            break;
        }
        MutexSum.unlock();
    }// ||||||||||||||||||||||||||END MC LOOP|||||||||||||||||||||||||||||
    Int_t nProj = 0;

    MutexHistEdg.lock();
    if (fDim > 0) {
        for (Int_t k = 0; k < fDim; k++) {
            dynamic_cast<TH1D *>((*fHistEdg)[nProj])->Add(LocalHistEdg[nProj].get());
            nProj++;
        }
    }
    MutexHistEdg.unlock();
}
