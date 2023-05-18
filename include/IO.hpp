#ifndef ROOT_TEST_IO_HPP
#define ROOT_TEST_IO_HPP

#include "Constant.hpp"
#include "Model.hpp"
#include <Math/AxisAngle.h>
#include <Math/Boost.h>
#include <Math/Rotation3D.h>
#include <Math/RotationX.h>
#include <Math/RotationY.h>
#include <Math/RotationZ.h>
#include <Math/RotationZYX.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <TFile.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TTree.h>
#include <iostream>
#include <utility>

class IO {
private:
    int counter = 0;

    std::unique_ptr<TFile> file;
    std::unique_ptr<TTree> tree;

    // Variables to be recorded in root
    Int_t ievent = 0, PIDlBeam = 0, PIDhBeam = 0;
    Int_t Nb = 0, PID[10] = {};
    Double_t ElBeam = 0, EhBeam = 0;
    Double_t Q2 = 0, W = 0, Nu = 0, XBj = 0, y = 0, t = 0, phih = 0;
    Double_t E[10] = {}, px[10] = {}, py[10] = {}, pz[10] = {};
    Double_t theta_g = 0;

    std::shared_ptr<TOPEG::Model_EventGeneratorFOAM> model;


    std::unique_ptr<TRandom> random = std::make_unique<TRandom3>();

    auto branching() -> void {
        tree->Branch("ievent", &ievent, "ievent/I");

        tree->Branch("ElBeam", &ElBeam, "ElBeam/D");
        tree->Branch("PIDlBeam", &PIDlBeam, "PIDlBeam/I");
        tree->Branch("EhBeam", &EhBeam, "EhBeam/D");
        tree->Branch("PIDhBeam", &PIDhBeam, "PIDhBeam/I");

        tree->Branch("Q2", &Q2, "Q2/D");
        tree->Branch("W", &W, "W/D");
        tree->Branch("Gamnu", &Nu, "Gamnu/D");
        tree->Branch("Xbj", &XBj, "Xbj/D");
        tree->Branch("y", &y, "y/D");
        tree->Branch("t", &t, "t/D");
        tree->Branch("phih", &phih, "phih/D");

        tree->Branch("Nb_part", &Nb, "Nb_part/I");
        tree->Branch("part_id", PID, "part_id[Nb_part]/I");
        tree->Branch("part_px", px, "part_px[Nb_part]/D");
        tree->Branch("part_py", py, "part_py[Nb_part]/D");
        tree->Branch("part_pz", pz, "part_pz[Nb_part]/D");
        tree->Branch("part_e", E, "part_e[Nb_part]/D");

        tree->Branch("theta_g", &theta_g, "theta_g/D");
    }
    auto initVars() -> void {
        PIDlBeam = 0;
        PIDhBeam = 0;
        ElBeam   = 0.;
        EhBeam   = 0.;
        Q2       = 0.;
        W        = 0.;
        Nu       = 0.;
        XBj      = 0.;
        y        = 0.;
        t        = 0.;
        phih     = 0.;

        Nb = 0;
        for (int i = 0; i < 10; i++) {
            PID[i] = 0;
            E[i]   = 0.;
            px[i]  = 0.;
            py[i]  = 0.;
            pz[i]  = 0.;
        }

        theta_g = 0;
    }
    auto fill4DOld(const std::vector<double> &MCvect) -> void {
        initVars();
        ievent++;


        PIDlBeam      = 11;
        PIDhBeam      = 1000020040;
        ElBeam        = model->GetEBene();
        EhBeam        = model->GetHBene();
        double PhBeam = model->GetHBmom();
        double mas    = model->Getmass();
        double s      = model->GetS();

        y    = MCvect.at(0);
        Q2   = MCvect.at(1);
        t    = MCvect.at(2);
        phih = MCvect.at(3);

        // Calculate the electron observables
        double ScaEne = (1 - y) * ElBeam + PhBeam * Q2 / (2 * ElBeam * (EhBeam + PhBeam));
        double ThetaE = 2 * std::asin(std::sqrt(Q2 / 4 / ElBeam / ScaEne));

        double xA = Q2 / y / (s - mas * mas);// This is xA
        double mN = 0.938;
        XBj       = xA * mas / mN;// This is xB
        Nu        = y / 2 / mas * (s - mas * mas);
        W         = std::sqrt(mas * mas - Q2 + y * (s - mas * mas));

        // Calculate the electron kinematic
        // TODO define this in class
        double PhiEle = random->Rndm() * 2 * M_PI;

        double Ek = ScaEne;
        double kx = ScaEne * std::sin(ThetaE) * std::cos(PhiEle);
        double ky = ScaEne * std::sin(ThetaE) * std::sin(PhiEle);
        double kz = -ScaEne * std::cos(ThetaE);

        // Calculate virtual photon
        double Eq = ElBeam - Ek;
        double qx = -kx;
        double qy = -ky;
        double qz = -ElBeam - kz;
        double Pq = std::sqrt(qx * qx + qy * qy + qz * qz);

        // Rotate into COM around x

        double thx = std::atan(qy / (qz + PhBeam));

        double kyp = ky * std::cos(thx) - kz * std::sin(thx);
        double kzp = kz * std::cos(thx) + ky * std::sin(thx);

        double qyp = qy * std::cos(thx) - qz * std::sin(thx);
        double qzp = qz * std::cos(thx) + qy * std::sin(thx);

        double pyp = -PhBeam * std::sin(thx);
        double pzp = PhBeam * std::cos(thx);

        // Rotate into COM around y

        double thy = std::atan(-qx / (qzp + pzp));

        double kxpp = kx * std::cos(thy) + kzp * std::sin(thy);
        double kzpp = kzp * std::cos(thy) - kx * std::sin(thy);

        double qxpp = qx * std::cos(thy) + qzp * std::sin(thy);
        double qzpp = qzp * std::cos(thy) - qx * std::sin(thy);

        double pxpp = pzp * std::sin(thy);
        double pzpp = pzp * std::cos(thy);

        // Boost into COM in z
        double bet = (qzpp + pzpp) / (Eq + EhBeam);
        double gam = 1 / std::sqrt(1 - bet * bet);

        double kzb = gam * (kzpp - bet * Ek);
        double Ekb = gam * (Ek - bet * kzpp);
        double Pkb = std::sqrt(kxpp * kxpp + kyp * kyp + kzb * kzb);

        double qzb = gam * (qzpp - bet * Eq);
        double Eqb = gam * (Eq - bet * qzpp);
        double Pqb = std::sqrt(qxpp * qxpp + qyp * qyp + qzb * qzb);

        double pzb = gam * (pzpp - bet * EhBeam);
        double Epb = gam * (EhBeam - bet * pzpp);
        double Ppb = std::sqrt(pxpp * pxpp + pyp * pyp + pzb * pzb);

        // Calculate photon kinematic base
        double Eqpb = ((Eqb + Epb) * (Eqb + Epb) - mas * mas) / (Eqb + Epb) / 2;
        double Thqb = std::acos(Eqb / Pqb + (Q2 - t) / (2 * Eqpb * Pqb));

        double tverif0 = Q2 + 2 * (Eqb * Eqpb - Pqb * Eqpb * std::cos(Thqb));

        // Longitudinal unit vector
        double lx = qxpp / Pqb;
        double ly = qyp / Pqb;
        double lz = qzb / Pqb;

        // Calculate unit vector for phi_h = 0
        double ux = qyp * (kxpp * qyp - kyp * qxpp) - qzb * (kzb * qxpp - kxpp * qzb);
        double uy = qzb * (kyp * qzb - kzb * qyp) - qxpp * (kxpp * qyp - kyp * qxpp);
        double uz = qxpp * (kzb * qxpp - kxpp * qzb) - qyp * (kyp * qzb - kzb * qyp);
        double Pu = std::sqrt(ux * ux + uy * uy + uz * uz);
        ux /= Pu;
        uy /= Pu;
        uz /= Pu;

        // Rotate unit vector by phi_h around q boosted
        double cph = std::cos(phih / 180 * M_PI);
        double sph = std::sin(phih / 180 * M_PI);

        double tx = ux * (cph + lx * lx * (1 - cph)) + uy * (lx * ly * (1 - cph) - lz * sph) + uz * (lx * lz * (1 - cph) + ly * sph);
        double ty = ux * (lx * ly * (1 - cph) + lz * sph) + uy * (cph + ly * ly * (1 - cph)) + uz * (ly * lz * (1 - cph) - lx * sph);
        double tz = ux * (lx * lz * (1 - cph) - ly * sph) + uy * (ly * lz * (1 - cph) + lx * sph) + uz * (cph + lz * lz * (1 - cph));

        // Boosted photon kine
        double cth = std::cos(Thqb);
        double sth = std::sin(Thqb);

        double bjkx = (cth * lx + sth * tx) * Eqpb;
        double bjky = (cth * ly + sth * ty) * Eqpb;
        double bjkz = (cth * lz + sth * tz) * Eqpb;

        double tverif1 = Q2 + 2 * (Eqb * Eqpb - qxpp * bjkx - qyp * bjky - qzb * bjkz);

        // Now boost back in z
        double jzpp = gam * (bjkz + bet * Eqpb);
        double Ejpp = gam * (Eqpb + bet * bjkz);

        // Rotate back around y
        double jxp = bjkx * std::cos(thy) - jzpp * std::sin(thy);
        double jzp = jzpp * std::cos(thy) + bjkx * std::sin(thy);

        // Rotate back around x
        double jx = jxp;
        double jy = bjky * std::cos(thx) + jzp * std::sin(thx);
        double jz = jzp * std::cos(thx) - bjky * std::sin(thx);

        double tverif2 = -(Eq - Ejpp) * (Eq - Ejpp) + (qx - jx) * (qx - jx) + (qy - jy) * (qy - jy) + (qz - jz) * (qz - jz);

        // Calculate Helium
        double Ppx = qx - jx;
        double Ppy = qy - jy;
        double Ppz = qz + PhBeam - jz;
        double Ppp = std::sqrt(Ppx * Ppx + Ppy * Ppy + Ppz * Ppz);
        double Epp = std::sqrt(Ppx * Ppx + Ppy * Ppy + Ppz * Ppz + mas * mas);

        // Set particles kinematic for recording
        Nb     = 3;
        PID[0] = 11;
        px[0]  = kx;
        py[0]  = ky;
        pz[0]  = kz;
        E[0]   = Ek;

        PID[1] = 22;
        px[1]  = jx;
        py[1]  = jy;
        pz[1]  = jz;
        E[1]   = Ejpp;

        PID[2] = 1000020040;
        px[2]  = Ppx;
        py[2]  = Ppy;
        pz[2]  = Ppz;
        E[2]   = Epp;

        tree->Fill();
    }
    auto fill4DNew(const std::vector<double> &MCvect) -> void {
        initVars();
        ievent++;

        PIDlBeam      = 11;
        PIDhBeam      = 1000020040;
        ElBeam        = model->GetEBene();
        EhBeam        = model->GetHBene();
        double PhBeam = model->GetHBmom();
        double mass   = model->Getmass();
        double s      = model->GetS();

        y    = MCvect.at(0);
        Q2   = MCvect.at(1);
        t    = MCvect.at(2);
        phih = MCvect.at(3);

        double xA = Q2 / y / (s - mass * mass);// This is xA
        double mN = 0.938;
        XBj       = xA * mass / mN;// This is xB
        Nu        = y / 2 / mass * (s - mass * mass);
        W         = std::sqrt(mass * mass - Q2 + y * (s - mass * mass));

        // Initialization -----------------------------------------------------------------------------
        ROOT::Math::PxPyPzEVector k1{0, 0, -ElBeam, ElBeam};
        ROOT::Math::PxPyPzEVector p1{0, 0, PhBeam, EhBeam};
        // --------------------------------------------------------------------------------------------

        // Compute scattering electron kinematic ------------------------------------------------------
        ROOT::Math::Boost boost_eN(k1.BoostToCM(p1));
        ROOT::Math::PxPyPzEVector k1b = boost_eN(k1);
        ROOT::Math::PxPyPzEVector p1b = boost_eN(p1);

        // In center of mass of (e, N) system
        double E_k2     = (1 - y) * k1b.E() + Q2 / (2 * (p1b.E() + k1b.E()));
        double theta_k2 = std::acos(1 - Q2 / (2 * k1b.E() * E_k2));
        double phi_k2   = random->Rndm() * 2 * M_PI;
        ROOT::Math::PxPyPzEVector k2b{E_k2 * std::sin(theta_k2) * std::cos(phi_k2), E_k2 * std::sin(theta_k2) * std::sin(phi_k2), -E_k2 * std::cos(theta_k2), E_k2};

        // Back in lab frame
        ROOT::Math::PxPyPzEVector k2 = boost_eN.Inverse()(k2b);
        // --------------------------------------------------------------------------------------------

        // Compute outgoing photon kinematic ----------------------------------------------------------
        // Calculate virtual photon in lab frame
        ROOT::Math::PxPyPzEVector q1 = k1 - k2;

        ROOT::Math::RotationX rotationX{std::atan(q1.y() / (q1.z() + PhBeam))};
        ROOT::Math::PxPyPzEVector k1x = rotationX(k1);
        ROOT::Math::PxPyPzEVector q1x = rotationX(q1);
        ROOT::Math::PxPyPzEVector p1x = rotationX(p1);

        ROOT::Math::RotationY rotationY{std::atan(-q1x.x() / (q1x.z() + p1x.z()))};
        ROOT::Math::PxPyPzEVector k1y = rotationY(k1x);
        ROOT::Math::PxPyPzEVector q1y = rotationY(q1x);
        ROOT::Math::PxPyPzEVector p1y = rotationY(p1x);

        // Boost in center of mass of (gamma*, n) system
        ROOT::Math::Boost boost(q1y.BoostToCM(p1y));
        ROOT::Math::PxPyPzEVector k1boost = boost(k1y);
        ROOT::Math::PxPyPzEVector q1boost = boost(q1y);
        ROOT::Math::PxPyPzEVector p1boost = boost(p1y);

        double Eq2boost      = (std::pow(q1boost.E() + p1boost.E(), 2) - mass * mass) / (2 * (q1boost.E() + p1boost.E()));
        double theta_q2boost = std::acos((Q2 - t) / (2 * Eq2boost * q1boost.R()) + q1boost.E() / q1boost.R());

        // Perpendicular vector to calcul phi
        ROOT::Math::XYZVector l = q1boost.Vect().Unit();
        ROOT::Math::XYZVector u = (q1boost.Vect().Cross(k1boost.Vect().Cross(q1boost.Vect()))).Unit();
        ROOT::Math::Rotation3D rotation{ROOT::Math::AxisAngle{l, phih / 180 * M_PI}};
        ROOT::Math::XYZVector vec_q2boost = (std::cos(theta_q2boost) * l + std::sin(theta_q2boost) * rotation(u)) * Eq2boost;
        ROOT::Math::PxPyPzEVector q2boost{vec_q2boost.x(), vec_q2boost.y(), vec_q2boost.z(), Eq2boost};

        // Boost back in lab frame
        // boost.Invert();
        // rotationY.Invert();
        // rotationX.Invert();
        ROOT::Math::PxPyPzEVector q2y = boost.Inverse()(q2boost);
        ROOT::Math::PxPyPzEVector q2x = rotationY.Inverse()(q2y);
        ROOT::Math::PxPyPzEVector q2  = rotationX.Inverse()(q2x);
        // --------------------------------------------------------------------------------------------

        // Calculate N kinematic ----------------------------------------------------------------------
        ROOT::Math::PxPyPzEVector p2 = q1 + p1 - q2;
        // --------------------------------------------------------------------------------------------

        // Set particles kinematic for recording
        Nb     = 3;
        PID[0] = 11;
        px[0]  = k2.x();
        py[0]  = k2.y();
        pz[0]  = k2.z();
        E[0]   = k2.E();

        PID[1] = 22;
        px[1]  = q2.x();
        py[1]  = q2.y();
        pz[1]  = q2.z();
        E[1]   = q2.E();

        PID[2] = 1000020040;
        px[2]  = p2.x();
        py[2]  = p2.y();
        pz[2]  = p2.z();
        E[2]   = p2.E();

        tree->Fill();
    }
    void Fill7D(const std::vector<double> &MCvect) {

        const double mn = 0.9396;
        const double mp = 0.9383;
        const double md = mn + mp - 0.00226;

        initVars();
        ievent++;
        ElBeam            = model->GetEBene();
        PIDlBeam          = 11;
        double PhBeam     = model->GetHBmom();
        EhBeam            = model->GetHBene();
        PIDhBeam          = 2212;
        double mas        = TOPEG::Constant::PROTON_MASS;
        double mass_s     = mp;
        int PID_active    = PIDhBeam;
        int PID_spectator = 2112;

        // std::cout << "ElBeam : " << ElBeam << " PhBeam : " << PhBeam << " EhBeam : " << EhBeam << std::endl;

        y    = MCvect.at(0);
        Q2   = MCvect.at(1);
        t    = MCvect.at(2);
        phih = MCvect.at(3);

        double moms   = MCvect.at(4);
        double thetas = MCvect.at(5) * M_PI / 180.;
        double phis   = MCvect.at(6) * M_PI / 180.;

        // Initialization -----------------------------------------------------------------------------
        ROOT::Math::PxPyPzEVector k1{0, 0, -ElBeam, ElBeam};
        ROOT::Math::PxPyPzEVector ps{moms * std::sin(thetas) * std::cos(phis), moms * std::sin(thetas) * std::sin(phis), moms * std::cos(thetas), std::hypot(moms, mass_s)};
        ROOT::Math::PxPyPzEVector p1{-ps.px(), -ps.py(), ps.pz(), std::sqrt(md * md - ps.E() * ps.E())};
        // --------------------------------------------------------------------------------------------

        double s      = (k1 + p1).M2();
        double xA     = Q2 / y / (s - mas * mas);
        XBj           = xA * mas / mn;
        Nu            = y / 2 / mas * (s - mas * mas);
        W             = std::sqrt(mas * mas - Q2 + y * (s - mas * mas));
        double W_test = std::sqrt(p1.E() * p1.E() - p1.P() * p1.P() - Q2 + y * (s - p1.E() * p1.E() - p1.P() * p1.P()));
        // std::cout << "W = " << W << " W_test = " << W_test << std::endl;

        // Compute scattering electron kinematic ------------------------------------------------------
        ROOT::Math::Boost boost_eN(k1.BoostToCM(p1));
        ROOT::Math::PxPyPzEVector k1b = boost_eN(k1);
        ROOT::Math::PxPyPzEVector p1b = boost_eN(p1);

        // In center of mass of (e, N) system
        double E_k2     = (1 - y) * k1b.E() + Q2 / (2 * (p1b.E() + k1b.E()));
        double theta_k2 = std::acos(1 - Q2 / (2 * k1b.E() * E_k2));
        double phi_k2   = random->Rndm() * 2 * M_PI;
        ROOT::Math::PxPyPzEVector k2b{E_k2 * std::sin(theta_k2) * std::cos(phi_k2), E_k2 * std::sin(theta_k2) * std::sin(phi_k2), -E_k2 * std::cos(theta_k2), E_k2};

        // Back in lab frame
        boost_eN.Invert();
        ROOT::Math::PxPyPzEVector k2 = boost_eN(k2b);
        // --------------------------------------------------------------------------------------------

        // Compute outgoing photon kinematic ----------------------------------------------------------
        // Calculate virtual photon in lab frame
        ROOT::Math::PxPyPzEVector q1 = k1 - k2;

        ROOT::Math::RotationX rotationX{std::atan((q1.y() + p1.y()) / (q1.z() + PhBeam))};
        ROOT::Math::PxPyPzEVector k1x = rotationX(k1);
        ROOT::Math::PxPyPzEVector q1x = rotationX(q1);
        ROOT::Math::PxPyPzEVector p1x = rotationX(p1);

        ROOT::Math::RotationY rotationY{std::atan((-q1x.x() - p1.x()) / (q1x.z() + p1x.z()))};
        ROOT::Math::PxPyPzEVector k1y = rotationY(k1x);
        ROOT::Math::PxPyPzEVector q1y = rotationY(q1x);
        ROOT::Math::PxPyPzEVector p1y = rotationY(p1x);

        // Boost in center of mass of (gamma*, n) system
        ROOT::Math::Boost boost(q1y.BoostToCM(p1y));
        ROOT::Math::PxPyPzEVector k1boost = boost(k1y);
        ROOT::Math::PxPyPzEVector q1boost = boost(q1y);
        ROOT::Math::PxPyPzEVector p1boost = boost(p1y);

        double Eq2boost      = (std::pow(q1boost.E() + p1boost.E(), 2) - mas * mas) / (2 * (q1boost.E() + p1boost.E()));
        double theta_q2boost = std::acos((Q2 - t) / (2 * Eq2boost * q1boost.R()) + q1boost.E() / q1boost.R());

        // Perpendicular vector to calcul phi
        ROOT::Math::XYZVector l = q1boost.Vect().Unit();
        ROOT::Math::XYZVector u = (q1boost.Vect().Cross(k1boost.Vect().Cross(q1boost.Vect()))).Unit();
        ROOT::Math::Rotation3D rotation{ROOT::Math::AxisAngle{l, phih / 180 * M_PI}};
        ROOT::Math::XYZVector v           = rotation(u);
        ROOT::Math::XYZVector vec_q2boost = (std::cos(theta_q2boost) * l + std::sin(theta_q2boost) * v) * Eq2boost;
        ROOT::Math::PxPyPzEVector q2boost{vec_q2boost.x(), vec_q2boost.y(), vec_q2boost.z(), Eq2boost};

        // Boost back in lab frame
        boost.Invert();
        ROOT::Math::PxPyPzEVector q2y = boost(q2boost);
        rotationY.Invert();
        ROOT::Math::PxPyPzEVector q2x = rotationY(q2y);
        rotationX.Invert();
        ROOT::Math::PxPyPzEVector q2 = rotationX(q2x);
        // --------------------------------------------------------------------------------------------

        // Calculate N kinematic ----------------------------------------------------------------------
        ROOT::Math::PxPyPzEVector p2 = q1 + p1 - q2;
        // --------------------------------------------------------------------------------------------

        // Set particles kinematic for recording
        Nb     = 4;
        PID[0] = 11;
        px[0]  = k2.x();
        py[0]  = k2.y();
        pz[0]  = -k2.z();
        E[0]   = k2.E();

        PID[1] = 22;
        px[1]  = q2.x();
        py[1]  = q2.y();
        pz[1]  = -q2.z();
        E[1]   = q2.E();

        PID[2] = PID_active;
        px[2]  = p2.x();
        py[2]  = p2.y();
        pz[2]  = -p2.z();
        E[2]   = p2.E();

        PID[3] = PID_spectator;
        px[3]  = ps.x();
        py[3]  = ps.y();
        pz[3]  = -ps.z();
        E[3]   = ps.E();

        q2.SetCoordinates(q2.x(), q2.y(), -q2.z(), q2.E());
        theta_g = q2.Theta() * 180. / 3.14159;
        // std::cout << "theta_g = " << theta_g << std::endl;
        ROOT::Math::PxPyPzEVector v1{px[1], py[1], pz[1], E[1]};
        // std::cout << "calc theta_g = " << v1.Theta() * 180. / 3.14159 << std::endl;

        if (theta_g < 4) return;
        counter++;

        tree->Fill();
    }

public:
    IO(const std::string &name, int seed, std::shared_ptr<TOPEG::Model_EventGeneratorFOAM> &model)
        : file(std::make_unique<TFile>(name.c_str(), "recreate")),
          tree(std::make_unique<TTree>("TOPEG", "Tree with simulated events from TOPEG")),
          model(model) {

        random->SetSeed(seed);
        branching();
    }

    void fill(int ModelNb, const std::vector<double> &MCvect) {
        switch (ModelNb) {
            case 161:
                fill4DOld(MCvect);
                break;
            case 131:
            case 531:
                Fill7D(MCvect);
                break;
            default:
                break;
        }
    }
    int getCounter() const { return counter; }
    auto write() -> void {
        tree->Write();
        file->Write();
    }
};


#endif//ROOT_TEST_IO_HPP
