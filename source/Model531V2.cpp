#include "../include/Model531V2.hpp"


TOPEG::Model531_v2::Model531_v2(const InitialsConditions &initialsConditions, const KinematicsRange &kinematicsRange)
    : fMassTarget(Constant::DEUTERIUM_MASS), fMassActive(Constant::NEUTRON_MASS), fMassSpectator(Constant::PROTON_MASS),
      fInitialsConditions(initialsConditions), fKinematicsRange(kinematicsRange) {
    std::cout << "Construct model 531." << std::endl;
    E_had_beam = std::hypot(initialsConditions.p_had_beam, fMassActive);
}

auto TOPEG::Model531_v2::Density(int nDim, Double_t *V) -> Double_t {
    if (nDim != fDim) throw std::invalid_argument("ERROR Nb of dimension mismatch in Model");
    std::vector<double> mc(fDim);
    for (int i = 0; i < fDim; i++) {
        if (V[i] < 0. || V[i] > 1.) throw std::invalid_argument("ERROR Variable out of bound in Model");
        mc[i] = V[i];
    }

    if (count++; count % 100000 == 0 && count != 0) {
        auto TimeP          = std::chrono::high_resolution_clock::now() - timeCount;
        auto [h, m, s, mss] = break_down_durations<std::chrono::hours, std::chrono::minutes, std::chrono::seconds, std::chrono::milliseconds>(TimeP);
        std::cout << h.count() << "h "
                  << m.count() << "min"
                  << s.count() << "s"
                  << mss.count() << "ms"
                  << "ms, for " << count << " calls." << std::endl;

        auto hrs  = std::chrono::duration_cast<std::chrono::hours>(TimeP);
        auto mins = std::chrono::duration_cast<std::chrono::minutes>(TimeP - hrs);
        auto secs = std::chrono::duration_cast<std::chrono::seconds>(TimeP - hrs - mins);
        auto ms   = std::chrono::duration_cast<std::chrono::milliseconds>(TimeP - hrs - mins - secs);
        std::cout << "Time elapsed is " << hrs.count() << "h " << mins.count() << "min " << secs.count() << "s " << ms.count() << "ms, for " << count << " calls." << std::endl;
    }

    Transform(mc);
    return Model531(mc);
}

auto TOPEG::Model531_v2::Transform(std::vector<Double_t> &mc) -> void {
    // Rescale var 4 to p spectator
    mc[4] = (mc[4] * (fKinematicsRange.p_spec_max - fKinematicsRange.p_spec_min)) + fKinematicsRange.p_spec_min;
    // Rescale var 5 to theta spectator
    // mc[5] = (mc[5] * (fKinematicsRange.theta_spec_max - fKinematicsRange.theta_spec_min)) + fKinematicsRange.theta_spec_min;
    mc[5] = mc[5] * 180.;
    // Rescale var 6 to phi spectator
    mc[6] = mc[6] * 360.;
    // Rescale var 0 to y
    mc[0] = (mc[0] * (fKinematicsRange.y_max - fKinematicsRange.y_min)) + fKinematicsRange.y_min;
    // Rescale var 1 to Q2
    mc[1] = (mc[1] * (Q2Max(mc) - fKinematicsRange.q_min)) + fKinematicsRange.q_min;
    // Rescale var 2 to t
    mc[2] = (mc[2] * fKinematicsRange.t_range) + tMin(mc);
    // Rescale var 3 to phi
    mc[3] = mc[3] * 360.;
}

auto TOPEG::Model531_v2::Model531(const std::vector<Double_t> &mc) const -> Double_t {
    double ebeam       = fInitialsConditions.E_el_beam;
    double eTarget     = 0.9396 + 0.9383 - 0.00226;
    double pol_el_beam = fInitialsConditions.pol_el_beam;

    int in    = 2;
    int igpde = 1;
    int itot  = 2;

    double Q2        = mc[1];
    double t         = -mc[2];
    double phih      = mc[3];
    double xb        = Q2 / (2 * fMassActive * fInitialsConditions.E_el_beam * mc[0]);
    double ps        = mc[4];
    double costhetas = std::cos(mc[5] / 57.2957795);
    double phis      = mc[6] / 57.2957795;
    double Dist      = 0;

    cross_boundnucleon_deut_f90_(&xb, &t, &Q2, &ebeam, &eTarget, &phih, &ps, &costhetas, &phis, &in, &igpde, &itot, &pol_el_beam, &Dist);

    return Dist;
}

auto TOPEG::Model531_v2::Q2Max(const std::vector<Double_t> &mc) const -> Double_t {
    Double_t psx    = mc[4] * std::sin(mc[5]) * std::cos(mc[6]);
    Double_t psy    = mc[4] * std::sin(mc[5]) * std::sin(mc[6]);
    Double_t psz    = mc[4] * std::cos(mc[5]);
    Double_t Es     = std::hypot(mc[4], fMassSpectator);
    Double_t phi_k2 = 0;
    ROOT::Math::PxPyPzEVector k1{0, 0, fInitialsConditions.E_el_beam, fInitialsConditions.E_el_beam};
    ROOT::Math::PxPyPzEVector v_p1{-psx, -psy, -psz, fMassTarget - Es};
    ROOT::Math::XYZVector u2{std::sin(fKinematicsRange.theta_el_scat_max) * std::cos(phi_k2), std::sin(fKinematicsRange.theta_el_scat_max) * std::sin(phi_k2), std::cos(fKinematicsRange.theta_el_scat_max)};

    Double_t s             = (k1 + v_p1).M2();
    Double_t Ek1           = fInitialsConditions.E_el_beam;
    Double_t Ep1           = v_p1.E();
    Double_t y             = mc[0];
    Double_t p1            = v_p1.P();
    Double_t costheta_p1k1 = std::cos(ROOT::Math::VectorUtil::Angle(v_p1.Vect(), k1.Vect()));
    Double_t costheta_p1u2 = std::cos(ROOT::Math::VectorUtil::Angle(v_p1.Vect(), u2));

    Double_t Ek2Max = ((1 - y) * Ek1 * (Ep1 - p1 * costheta_p1k1)) / (Ep1 - p1 * costheta_p1u2);
    Double_t Q2Max1 = 4 * Ek1 * Ek2Max * std::pow(std::sin(fKinematicsRange.theta_el_scat_max / 2), 2);
    Double_t Q2Max2 = Ep1 * Ep1 - p1 * p1 + y * (s - Ep1 * Ep1 + p1 * p1) - fKinematicsRange.W2_min;
    Double_t Q2Max3 = fKinematicsRange.xBj_max * y * (s - Ep1 * Ep1 + p1 * p1);

    double Q2Max = std::min({fKinematicsRange.q_max, Q2Max1, Q2Max2, Q2Max3}) - 0.01;
    if (Q2Max < fKinematicsRange.q_min) throw std::invalid_argument("ERROR y cut are too wide, problem in DVCSmodel::Q2Max : " + std::to_string(mc[0]) + " " + std::to_string(Q2Max1) +
                                                                    " " + std::to_string(Q2Max2) + " " + std::to_string(Q2Max3) + " " + std::to_string(fKinematicsRange.q_max));
    return Q2Max;
}

auto TOPEG::Model531_v2::tMin(const std::vector<Double_t> &mc) const -> Double_t {
    // Initialization -----------------------------------------------------------------------------
    Double_t Es     = std::hypot(mc[4], fMassSpectator);
    Double_t phi_k2 = 0;
    ROOT::Math::PxPyPzEVector k1{0, 0, -fInitialsConditions.E_el_beam, fInitialsConditions.E_el_beam};
    ROOT::Math::PxPyPzEVector p1{-(mc[4] * std::sin(mc[5]) * std::cos(mc[6])), -(mc[4] * std::sin(mc[5]) * std::sin(mc[6])), -(mc[4] * std::cos(mc[5])), fMassTarget - Es};
    Double_t y  = mc[0];
    Double_t Q2 = mc[1];
    // --------------------------------------------------------------------------------------------

    // Compute scattering electron kinematic ------------------------------------------------------
    ROOT::Math::Boost boost_eN(k1.BoostToCM(p1));
    ROOT::Math::PxPyPzEVector k1b = boost_eN(k1);
    ROOT::Math::PxPyPzEVector p1b = boost_eN(p1);

    ROOT::Math::XYZVector z(0, 0, -1);
    ROOT::Math::XYZVector v = k1b.Vect().Unit().Cross(z);
    Double_t c              = k1b.Vect().Unit().Dot(z);

    TMatrixD I(3, 3);
    I(0, 0) = 1;
    I(1, 1) = 1;
    I(2, 2) = 1;

    TMatrixD v_cross(3, 3);
    v_cross(0, 1) = -v.z();
    v_cross(0, 2) = v.y();
    v_cross(1, 0) = v.z();
    v_cross(1, 2) = -v.x();
    v_cross(2, 0) = -v.y();
    v_cross(2, 1) = v.x();
    TMatrixD R    = I + v_cross + v_cross * v_cross * (1. / (1. + c));
    ROOT::Math::Rotation3D rotation3D(R(0, 0), R(0, 1), R(0, 2), R(1, 0), R(1, 1), R(1, 2), R(2, 0), R(2, 1), R(2, 2));
    k1b = rotation3D(k1b);
    p1b = rotation3D(p1b);
    // In center of mass of (e, N) system
    Double_t E_k2     = (1 - y) * k1b.E() + Q2 / (2 * (p1b.E() + k1b.E()));
    Double_t theta_k2 = std::acos(1 - Q2 / (2 * k1b.E() * E_k2));
    ROOT::Math::PxPyPzEVector k2b{E_k2 * std::sin(theta_k2) * std::cos(phi_k2), E_k2 * std::sin(theta_k2) * std::sin(phi_k2), -E_k2 * std::cos(theta_k2), E_k2};
    k2b                          = rotation3D.Inverse()(k2b);
    ROOT::Math::PxPyPzEVector k2 = boost_eN.Inverse()(k2b);// Back in lab frame
    // --------------------------------------------------------------------------------------------

    // Compute outgoing photon kinematic ----------------------------------------------------------
    // Calculate virtual photon in lab frame
    ROOT::Math::PxPyPzEVector q1 = k1 - k2;

    // Boost in center of mass of (gamma*, n) system
    ROOT::Math::Boost boost(q1.BoostToCM(p1));
    ROOT::Math::PxPyPzEVector q1boost = boost(q1);
    ROOT::Math::PxPyPzEVector p1boost = boost(p1);

    Double_t Eq2boost = (std::pow(q1boost.E() + p1boost.E(), 2) - p1boost.E() * p1boost.E() + p1boost.P() * p1boost.P()) / (2 * (q1boost.E() + p1boost.E()));
    Double_t tmin     = Q2 + 2 * Eq2boost * (q1boost.E() - q1boost.P());
    return std::max(tmin, fKinematicsRange.t_min);
}
