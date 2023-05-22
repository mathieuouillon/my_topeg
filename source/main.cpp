#include <IO.hpp>
#include <Model131V2.hpp>
#include <Model531V2.hpp>
#include <Parser.hpp>
#include <TFoamMT.h>
#include <TH2D.h>
#include <TMath.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <draw.hpp>
#include <iostream>
#include <root2lund.hpp>

int count = 0;

extern "C" {
static std::mutex FortranModelMutex;
}
extern "C" {
void mutexlock_() { FortranModelMutex.lock(); }
}
extern "C" {
void mutexunlock_() { FortranModelMutex.unlock(); }
}

struct Histograms {
    int32_t NbBin = 100;

    std::shared_ptr<TH1F> h_p0_p = std::make_shared<TH1F>("p0_p", "", NbBin, 0, 10);
    std::shared_ptr<TH1F> h_p1_p = std::make_shared<TH1F>("p1_p", "", NbBin, 0, 10);
    std::shared_ptr<TH1F> h_p2_p = std::make_shared<TH1F>("p2_p", "", NbBin, 0, 1);
    std::shared_ptr<TH1F> h_p3_p = std::make_shared<TH1F>("p3_p", "", NbBin, 0, 1);

    std::shared_ptr<TH1F> h_p0_theta = std::make_shared<TH1F>("p0_theta", "", NbBin, -10, 180);
    std::shared_ptr<TH1F> h_p1_theta = std::make_shared<TH1F>("p1_theta", "", NbBin, -1, 20);
    std::shared_ptr<TH1F> h_p2_theta = std::make_shared<TH1F>("p2_theta", "", NbBin, -10, 180);
    std::shared_ptr<TH1F> h_p3_theta = std::make_shared<TH1F>("p3_theta", "", NbBin, -10, 180);

    std::shared_ptr<TH1F> h_p0_phi = std::make_shared<TH1F>("p0_phi", "", NbBin, -10, 370);
    std::shared_ptr<TH1F> h_p1_phi = std::make_shared<TH1F>("p1_phi", "", NbBin, -10, 370);
    std::shared_ptr<TH1F> h_p2_phi = std::make_shared<TH1F>("p2_phi", "", NbBin, -10, 370);
    std::shared_ptr<TH1F> h_p3_phi = std::make_shared<TH1F>("p3_phi", "", NbBin, -10, 370);

    std::shared_ptr<TH1F> h_Q2    = std::make_shared<TH1F>("Q2", "", NbBin, 0, 10);
    std::shared_ptr<TH1F> h_W2    = std::make_shared<TH1F>("W2", "", NbBin, 0, 40);
    std::shared_ptr<TH1F> h_Gamnu = std::make_shared<TH1F>("Gamnu", "", NbBin, 0, 20);
    std::shared_ptr<TH1F> h_Xbj   = std::make_shared<TH1F>("Xbj", "", NbBin, 0, 1);
    std::shared_ptr<TH1F> h_y     = std::make_shared<TH1F>("y", "", NbBin, 0, 1);
    std::shared_ptr<TH1F> h_t     = std::make_shared<TH1F>("t", "", NbBin, 0, 3);
    std::shared_ptr<TH1F> h_log_t = std::make_shared<TH1F>("log_t", "", NbBin, 0, 3);
    std::shared_ptr<TH1F> h_phih  = std::make_shared<TH1F>("phih", "", NbBin, -10, 360);

    std::shared_ptr<TH1F> h_theta_g = std::make_shared<TH1F>("theta_g", "", NbBin, -1, 20);
};

auto setROOTOption() -> void {
    ROOT::EnableThreadSafety();
    gStyle->SetOptStat("emr");
    gStyle->SetOptFit();
    gStyle->SetNumberContours(255);
    gStyle->SetImageScaling(3.);

    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetHistLineWidth(2);
    gStyle->SetFuncWidth(2);
    gStyle->SetGridWidth(2);
    gStyle->SetLineStyleString(2, "[12 12]");// postscript dashes
    gStyle->SetMarkerStyle(kFullCircle);
    gStyle->SetMarkerSize(1.0);

    // put tick marks on top and RHS of plots
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetNdivisions(505, "x");
    gStyle->SetNdivisions(510, "y");
}
auto readConfig() -> ParserOutput {
    std::string OptFileName = "../config/config.txt";
    Parser parser(OptFileName);
    auto output      = parser.get<std::string>("output");
    auto EBene       = parser.get<double>("EBene");
    auto PHmom       = parser.get<double>("PHmom");
    auto target      = parser.get<int>("target");
    auto process     = parser.get<int>("process");
    auto model       = parser.get<int>("model");
    auto Nb_event    = parser.get<int>("Nb_event");
    auto NB_cells    = parser.get<int>("NB_cells");
    auto samples     = parser.get<int>("Samples_per_cells");
    auto seed        = parser.get<int>("seed");
    auto y_min       = parser.get<double>("y_min");
    auto y_max       = parser.get<double>("y_max");
    auto Q2_min      = parser.get<double>("Q2_min");
    auto Q2_max      = parser.get<double>("Q2_max");
    auto W2_min      = parser.get<double>("W2_min");
    auto theta_e_max = parser.get<double>("theta_e_max");
    auto t_min       = parser.get<double>("t_min");
    auto t_range     = parser.get<double>("t_range");
    auto e_helicity  = parser.get<double>("e_helicity");

    // TODO Add a safety function to check parameters
    std::cout << "These are the input parameters:    " << std::endl;
    std::cout << "Output file name                   " << output << std::endl;
    std::cout << "Electron beam energy (GeV)         " << EBene << std::endl;
    std::cout << "Nuclear beam momentum (GeV)        " << PHmom << std::endl;
    std::cout << "Target (1-8)                       " << target << std::endl;
    std::cout << "Process (1-4)                      " << process << std::endl;
    std::cout << "Model number (1-4)                 " << model << std::endl;
    std::cout << "Nb of events to generate           " << Nb_event << std::endl;
    std::cout << "Nb of TFoam cells                  " << NB_cells << std::endl;
    std::cout << "Nb of samples per cell             " << samples << std::endl;
    std::cout << "Seed                               " << seed << std::endl;
    std::cout << "y min                              " << y_min << std::endl;
    std::cout << "y Max                              " << y_max << std::endl;
    std::cout << "Q2 min (GeV^2)                     " << Q2_min << std::endl;
    std::cout << "Q2 Max (GeV^2)                     " << Q2_max << std::endl;
    std::cout << "W2 min (GeV^2)                     " << W2_min << std::endl;
    std::cout << "theta Max electron (rad)           " << theta_e_max << std::endl;
    std::cout << "t min (GeV^2)                      " << t_min << std::endl;
    std::cout << "t Range (GeV^2)                    " << t_range << std::endl;
    std::cout << "Electron beam helicity             " << e_helicity << std::endl;
    std::cout << "                                   " << std::endl;

    return {output, EBene, PHmom, target, process, model, Nb_event, NB_cells, samples, seed, y_min, y_max, Q2_min, Q2_max, W2_min, theta_e_max, t_min, t_range, e_helicity};
}
auto readModelSara() -> void {
    Histograms histograms{};

    std::unique_ptr<TFile> myFile(TFile::Open("../saraModelEventGenerated.root"));
    auto tree = std::unique_ptr<TTree>(myFile->Get<TTree>("TOPEG"));

    Int_t ievent, PIDlBeam, Nb_part;
    Double_t ElBeam, EhBeam;
    Double_t Q2, W, Nu, XBj, y, t, phih;
    Double_t part_px[10], part_py[10], part_pz[10], part_E[10];
    Double_t theta_g;
    tree->SetBranchAddress("ievent", &ievent);
    tree->SetBranchAddress("PIDlBeam", &PIDlBeam);
    tree->SetBranchAddress("ElBeam", &ElBeam);
    tree->SetBranchAddress("EhBeam", &EhBeam);
    tree->SetBranchAddress("Nb_part", &Nb_part);
    tree->SetBranchAddress("part_px", part_px);
    tree->SetBranchAddress("part_py", part_py);
    tree->SetBranchAddress("part_pz", part_pz);
    tree->SetBranchAddress("part_e", part_E);

    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("W", &W);
    tree->SetBranchAddress("Gamnu", &Nu);
    tree->SetBranchAddress("Xbj", &XBj);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("t", &t);
    tree->SetBranchAddress("phih", &phih);

    tree->SetBranchAddress("theta_g", &theta_g);

    double pi = 3.14159;

    for (int32_t iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {
        tree->GetEntry(iEntry);

        ROOT::Math::PxPyPzEVector v0{part_px[0], part_py[0], part_pz[0], part_E[0]};
        ROOT::Math::PxPyPzEVector v1{part_px[1], part_py[1], part_pz[1], part_E[1]};
        ROOT::Math::PxPyPzEVector v2{part_px[2], part_py[2], part_pz[2], part_E[2]};
        ROOT::Math::PxPyPzEVector v3{part_px[3], part_py[3], part_pz[3], part_E[3]};

        ROOT::Math::PxPyPzEVector k1{0,0,10.4,10.4};
        ROOT::Math::PxPyPzEVector q1 = k1 - v0;
        std::cout << "W2 = " << (v2 + q1).M2() << std::endl;

        histograms.h_p0_p->Fill(v0.P());
        histograms.h_p0_theta->Fill(v0.Theta() * 180. / 3.14159);
        histograms.h_p0_phi->Fill(v0.Phi() * 180. / pi + 180);

        histograms.h_p1_p->Fill(v1.P());
        histograms.h_p1_theta->Fill(v1.Theta() * 180. / pi);
        histograms.h_p1_phi->Fill(v1.Phi() * 180. / pi + 180);

        histograms.h_p2_p->Fill(v2.P());
        histograms.h_p2_theta->Fill(v2.Theta() * 180. / pi);
        histograms.h_p2_phi->Fill(v2.Phi() * 180. / pi + 180);

        histograms.h_p3_p->Fill(v3.P());
        histograms.h_p3_theta->Fill(v3.Theta() * 180. / pi);
        histograms.h_p3_phi->Fill(v3.Phi() * 180. / pi + 180);

        histograms.h_Q2->Fill(Q2);
        histograms.h_W2->Fill((v2 + q1).M2());
        histograms.h_Gamnu->Fill(Nu);
        histograms.h_Xbj->Fill(XBj);
        histograms.h_y->Fill(y);
        histograms.h_t->Fill(t);
        histograms.h_log_t->Fill(t);
        histograms.h_phih->Fill(phih);

        histograms.h_theta_g->Fill(theta_g);
    }


    const std::string path = "../data/plotsModelSara/";

    Draw::draw(histograms.h_p0_p, path, "p_{e} [GeV/c]");
    Draw::draw(histograms.h_p0_theta, path, "#theta_{e} [Deg.]");
    Draw::draw(histograms.h_p0_phi, path, "#phi_{e} [Deg.]");

    Draw::draw(histograms.h_p1_p, path, "p_{#gamma} [GeV/c]");
    Draw::draw(histograms.h_p1_theta, path, "#theta_{#gamma} [Deg.]");
    Draw::draw(histograms.h_p1_phi, path, "#phi_{#gamma} [Deg.]");

    Draw::draw(histograms.h_p2_p, path, "p_{N} [GeV/c]");
    Draw::draw(histograms.h_p2_theta, path, "#theta_{N} [Deg.]");
    Draw::draw(histograms.h_p2_phi, path, "#phi_{N} [Deg.]");

    Draw::draw(histograms.h_p3_p, path, "p_{s} [GeV/c]");
    Draw::draw(histograms.h_p3_theta, path, "#theta_{s} [Deg.]");
    Draw::draw(histograms.h_p3_phi, path, "#phi_{s} [Deg.]");

    Draw::draw(histograms.h_Q2, path, "Q^{2} [GeV^{2}/c^{2}]");
    Draw::draw(histograms.h_W2, path, "W^{2} [GeV]");
    Draw::draw(histograms.h_Gamnu, path, "#nu [GeV]");
    Draw::draw(histograms.h_Xbj, path, "x_{Bj}");
    Draw::draw(histograms.h_y, path, "y");
    Draw::draw(histograms.h_t, path, "t [GeV^{2}/c^{2}]");
    Draw::draw(histograms.h_log_t, path, "t [GeV^{2}/c^{2}]", {.logx = true});
    Draw::draw(histograms.h_phih, path, "#phi_{h} [Deg.]");

    Draw::draw(histograms.h_theta_g, path, "#theta_{#gamma} [Deg.]");
}
auto multiThread() -> void {
    ParserOutput parserOutput = readConfig();
    parserOutput.model        = 5;
    parserOutput.target       = 3;
    parserOutput.process      = 1;
    std::uint16_t modelNb     = parserOutput.model * 100 + parserOutput.target * 10 + parserOutput.process;
    std::uint16_t nbThread    = 20;
    std::cout << "Model = " << modelNb << std::endl;

    InitialsConditions initialsConditions(parserOutput.EBene, parserOutput.e_helicity, parserOutput.PHmom);
    KinematicsRange kinematicsRange(parserOutput.y_min, parserOutput.y_max, parserOutput.Q2_min, parserOutput.Q2_max, parserOutput.W2_min, parserOutput.theta_e_max, parserOutput.t_min, parserOutput.t_range);

    // std::shared_ptr<TOPEG::Model_EventGeneratorFOAM> rho = std::make_shared<TOPEG::Model531_v2>(initialsConditions, kinematicsRange);
    std::shared_ptr<TOPEG::Model_EventGeneratorFOAM> rho = std::make_shared<TOPEG::Model131V2_EventGeneratorFOAM>(
        parserOutput.EBene, parserOutput.PHmom, parserOutput.y_min, parserOutput.y_max, parserOutput.Q2_min, parserOutput.Q2_max, 
        parserOutput.W2_min, parserOutput.theta_e_max, parserOutput.t_min, parserOutput.t_range, parserOutput.e_helicity);

    const std::string gridName = "../grid_SaraModel.root";
    const std::string fileName = "../saraModelEventGenerated.root";
    auto GridFile              = std::make_unique<TFile>(gridName.c_str(), "READ");

    auto PseRan = std::make_unique<TRandom3>(parserOutput.seed);
    auto FoamX  = std::make_unique<TFoamMT>("FoamX", nbThread);// Create Simulator
    if (GridFile->IsOpen()) {
        std::cout << "Grid file opened successfully." << std::endl;
        FoamX = std::unique_ptr<TFoamMT>(reinterpret_cast<TFoamMT *>(GridFile->Get("FoamData")));
        FoamX->SetRho(rho.get());// Set 4-dim distribution, included below
        FoamX->ResetPseRan(PseRan.get());
    } else {
        std::cout << "Grid file not opened, starting initialize... " << std::endl;
        FoamX->SetkDim(rho->getDimNb());
        FoamX->SetChat(1);
        FoamX->SetnCells(1e7);  // No. of cells, default=2000
        FoamX->SetnSampl(300);  // No. of MC events in the cell MC exploration d=200
        FoamX->SetnBin(8);      // No. of bins in edge-histogram in cell exploration d=8
        FoamX->SetOptRej(1);    // Wted events for OptRej=0; wt=1 for OptRej=1 (default)
        FoamX->SetOptDrive(2);  // Maximum weight reduction, =1 for variance reduction d=2 maximum weight optimization
        FoamX->SetEvPerBin(25); // Maximum number of the effective wt=1 events/bin
        FoamX->SetMaxWtRej(1.1);// Maximum weight used to get w=1 MC events d=1.1
        FoamX->SetRho(rho.get());
        FoamX->SetPseRan(PseRan.get());
        FoamX->Initialize();

        GridFile = std::make_unique<TFile>(gridName.c_str(), "NEW");
        if (GridFile->IsOpen()) std::cout << "Grid file is created." << std::endl;
        FoamX->Write("FoamData");
    }


    // Generate event -----------------------
    std::vector<double> mcEvent = {};
    mcEvent.resize(rho->getDimNb());
    auto io = std::make_unique<IO>(fileName.c_str(), parserOutput.seed, rho);
    for (long loop = 0; loop < 10000000; ++loop) {
        FoamX->MakeEvent();
        FoamX->GetMCvect(mcEvent.data());
        rho->Transform(mcEvent);
        io->fill(modelNb, mcEvent);
        int counter = io->getCounter();
        if (counter % 1000 == 0 && counter != 0) std::cout << "nb event = " << counter << std::endl;
        if (counter == 10000) break;
    }
    io->write();
}

auto main() -> int {
    auto start = std::chrono::high_resolution_clock::now();

    setROOTOption();

    multiThread();
    readModelSara();

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " [s]" << std::endl;

    return EXIT_SUCCESS;
}
