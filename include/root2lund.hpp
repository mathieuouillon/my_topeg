#ifndef ROOT_TEST_ROOT2LUND_HPP
#define ROOT_TEST_ROOT2LUND_HPP

#include <Math/Vector4D.h>
#include <TFile.h>
#include <TTree.h>
#include <fstream>
#include <memory>

auto root2lund() -> void {
    std::unique_ptr<TFile> myFile(TFile::Open("../testSaraModelMultiThread.root"));
    auto tree = std::unique_ptr<TTree>(myFile->Get<TTree>("TOPEG"));

    int32_t nb = 0;
    std::ofstream outFile;
    outFile.open("../data/lund_file/bonus12events_" + std::to_string(nb) + ".lund");

    Int_t ievent, PIDlBeam, Nb_part;
    Double_t ElBeam, EhBeam;
    Double_t part_px[11], part_py[11], part_pz[11], part_E[11];
    tree->SetBranchAddress("ievent", &ievent);
    tree->SetBranchAddress("PIDlBeam", &PIDlBeam);
    tree->SetBranchAddress("ElBeam", &ElBeam);
    tree->SetBranchAddress("EhBeam", &EhBeam);
    tree->SetBranchAddress("Nb_part", &Nb_part);
    tree->SetBranchAddress("part_px", part_px);
    tree->SetBranchAddress("part_py", part_py);
    tree->SetBranchAddress("part_pz", part_pz);
    tree->SetBranchAddress("part_e", part_E);

    for (int32_t iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {

        if (iEntry % 10000 == 0 and iEntry != 0) {
            outFile.close();
            std::cout << "entry = " << iEntry << std::endl;
            nb++;
            outFile.open("../lund_file/bonus12events_" + std::to_string(nb) + ".lund");
        }



        tree->GetEntry(iEntry);

        ROOT::Math::PxPyPzEVector v0{part_px[0], part_py[0], part_pz[0], part_E[0]};// Electron 11
        ROOT::Math::PxPyPzEVector v1{part_px[1], part_py[1], part_pz[1], part_E[1]};// Photon 22
        ROOT::Math::PxPyPzEVector v2{part_px[2], part_py[2], part_pz[2], part_E[2]};// Neutron 2112
        ROOT::Math::PxPyPzEVector v3{part_px[3], part_py[3], part_pz[3], part_E[3]};// Proton 2212

        outFile << "4 \t 2 \t 1 \t 0 \t -1 \t 11 \t " + std::to_string(ElBeam) + "\t 2112 \t 1 \t 1\n";
        outFile << "1 \t 0 \t 1 \t 11 \t 0 \t 0 \t " + std::to_string(v0.px()) + "\t " + std::to_string(v0.py()) + "\t " + std::to_string(v0.pz()) + "\t " + std::to_string(v0.E()) + "\t 0.000511 \t 0 \t 0 \t 0\n";
        outFile << "2 \t 0 \t 1 \t 22 \t 0 \t 0 \t " + std::to_string(v1.px()) + "\t " + std::to_string(v1.py()) + "\t " + std::to_string(v1.pz()) + "\t " + std::to_string(v1.E()) + "\t 0 \t 0 \t 0 \t 0\n";
        outFile << "3 \t 0 \t 1 \t 2112 \t 0 \t 0 \t " + std::to_string(v2.px()) + "\t " + std::to_string(v2.py()) + "\t " + std::to_string(v2.pz()) + "\t " + std::to_string(v2.E()) + "\t 0.939 \t 0 \t 0 \t 0\n";
        outFile << "4 \t 0 \t 1 \t 2212 \t 0 \t 0 \t " + std::to_string(v3.px()) + "\t " + std::to_string(v3.py()) + "\t " + std::to_string(v3.pz()) + "\t " + std::to_string(v3.E()) + "\t 0.938 \t 0 \t 0 \t 0\n";
    }

    outFile.close();
}

#endif//ROOT_TEST_ROOT2LUND_HPP
