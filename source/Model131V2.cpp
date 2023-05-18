#include "../include/Model131V2.hpp"

TOPEG::Model131V2_EventGeneratorFOAM::Model131V2_EventGeneratorFOAM(double E_k1, double p_p1, double y_min, double y_max, double q_min, double q_max,
                                                                    double W2_min, double theta_k2, double t_min, double t_range, double pol_k1)
    : EBene(E_k1), HBmom(p_p1), EBpol(pol_k1), ym(y_min), yM(y_max), qm(q_min), qM(q_max), Wm(W2_min), TM(theta_k2),
      tmi(t_min), tR(t_range) {
    std::cout << "Construct Model 131." << std::endl;

    mass  = Constant::NEUTRON_MASS;
    HBene = std::sqrt(HBmom * HBmom + mass * mass);
    s     = mass * mass + 2 * EBene * (HBene + HBmom);

    setKinLimits();
}

auto TOPEG::Model131V2_EventGeneratorFOAM::Transform(double *event) -> void {}

auto TOPEG::Model131V2_EventGeneratorFOAM::Transform(std::array<double, m_nDim> &mc) -> void {
    // Rescale var 0 to y
    mc[0] = (mc[0] * (yM - ym)) + ym;
    // Rescale var 1 to Q2
    double Q2M = getqM(mc);
    mc[1]      = (mc[1] * (Q2M - qm)) + qm;
    // Rescale var 2 to t
    double tmin = gettm(mc);
    mc[2]       = (mc[2] * tR) + tmin;
    // Rescale var 3 to phi_h
    mc[3] = mc[3] * 360.0;
    // Rescale var 4 to p spectator
    mc[4] = (mc[4] * (psM - psm)) + psm;
    // Rescale var 5 to theta spectator
    mc[5] = (mc[5] * (thetasM - thetasm)) + thetasm;
    // Rescale var 6 to phi spectator
    mc[6] = (mc[6] * (phisM - phism)) + phism;
}

auto TOPEG::Model131V2_EventGeneratorFOAM::Density(int nDim, Double_t *V) -> Double_t {
    std::stringstream msg;
    msg << "thread id = " << std::this_thread::get_id() << " call density function nb : " << counter << "\n";
    // std::cout << msg.str() << std::flush;

    if (nDim != m_nDim) throw std::invalid_argument(std::string("ERROR Nb of dimension mismatch in Model 131 -- dim FOAM = " + std::to_string(nDim) + " dim Model = " + std::to_string(m_nDim)));

    std::array<double, m_nDim> MC = {};
    for (int i = 0; i < m_nDim; i++) {
        if (V[i] < 0. || V[i] > 1.) throw std::invalid_argument(std::string("ERROR Variable out of bound in Model"));
        MC[i] = V[i];
    }


    Transform(MC);

    counter++;
    if (counter % 10000 == 100) {
        auto TimeP = std::chrono::high_resolution_clock::now() - TimeC;
        auto hrs   = std::chrono::duration_cast<std::chrono::hours>(TimeP);
        auto mins  = std::chrono::duration_cast<std::chrono::minutes>(TimeP - hrs);
        auto secs  = std::chrono::duration_cast<std::chrono::seconds>(TimeP - hrs - mins);
        auto ms    = std::chrono::duration_cast<std::chrono::milliseconds>(TimeP - hrs - mins - secs);
        std::cout << "Time elapsed is " << hrs.count() << "h " << mins.count() << "min " << secs.count() << "s " << ms.count() << "ms, for " << counter << " calls." << std::endl;
    }

    return Model131(MC);
}

auto TOPEG::Model131V2_EventGeneratorFOAM::deut_mom(double x) -> double {
    double A01 = 157.4, A02 = 0.234, A03 = 0.00623;
    double B01 = 1.24, B02 = 1.27, B03 = 0.220;
    double C01 = 18.3, C02 = 0, C03 = 0;

    double k = x / 0.1973;

    double n = A01 * std::exp(-B01 * k * k) / std::pow(1. + C01 * k * k, 2);
    n += A02 * std::exp(-B02 * k * k) / std::pow(1. + C02 * k * k, 2);
    n += A03 * std::exp(-B03 * k * k) / std::pow(1. + C03 * k * k, 2);

    return 3 * n * k * k;
}

auto TOPEG::Model131V2_EventGeneratorFOAM::setKinLimits() -> void {
    // Reset proper y limits if in conflict with other cuts
    double epsilon = 0.01;
    if (double Q2a = 4 * EBene * EBene * (1 - yM) * std::pow(std::tan(TM / 2.), 2); Q2a < qm) {
        yM = 1 - qm / 4 / EBene / EBene / std::pow(std::tan(TM / 2.), 2) - epsilon;
        std::cout << "ERROR y cut are too wide in DVCSmodel::SetKinLimits" << std::endl;
        std::cout << "y max is reset to " << yM << std::endl;
        if (ym > yM) exit(EXIT_FAILURE);
    }

    if (double Q2b = mass * mass + 2. * EBene * (HBene + HBmom) * ym - Wm; Q2b < qm) {
        ym = (qm - mass * mass + Wm) / (2. * EBene * (HBene + HBmom)) + epsilon;
        std::cout << "ym = " << ym << std::endl;
        std::cout << "ERROR y cut are too wide in DVCSmodel::SetKinLimits" << std::endl;
        std::cout << "y min is reset to " << ym << std::endl;
        if (ym > yM) exit(EXIT_FAILURE);
    }
}

auto TOPEG::Model131V2_EventGeneratorFOAM::getqM(const std::array<double, m_nDim> &mc) -> double {
    double epsilon = 0.01;
    double Q2a     = 4 * EBene * EBene * (1 - mc[0]) * std::pow(std::tan(TM / 2.), 2);
    double Q2b     = mass * mass + (s - mass * mass) * mc[0] - Wm;
    double Q2c     = xM * Constant::PROTON_MASS / mass * (s - mass * mass) * mc[0];

    // Select the stringent Q2 max
    double Q2M = qM;
    if (Q2M > Q2a) Q2M = Q2a;
    if (Q2M > Q2b) Q2M = Q2b;
    if (Q2M > Q2c) Q2M = Q2c;
    Q2M -= epsilon;

    if (Q2M < qm) {
        std::cout << "ERROR y cut are too wide, problem in DVCSmodel::GetqM" << std::endl;
        std::cout << mc[0] << " " << Q2a << " " << Q2b << " " << Q2c << " " << Q2M << std::endl;
        return 0.;
    }

    return Q2M;
}

auto TOPEG::Model131V2_EventGeneratorFOAM::getxa(const std::array<double, m_nDim> &mc) -> double {
    double xb = mass / Constant::PROTON_MASS * (mc[1] / ((s - mass * mass) * mc[0]));
    if (xb > xM) {
        std::cout << "WARNING: strange kinematics! Fix xb value in DVCSmodel::Getxa" << std::endl;
        std::cout << xb << " " << mc[1] << " " << mc[0] << " " << std::endl;
        xb = xM;
    }

    double xa = Constant::PROTON_MASS / mass * xb;
    return xa;
}

auto TOPEG::Model131V2_EventGeneratorFOAM::gettm(const std::array<double, m_nDim> &mc) -> double {
    double xa   = getxa(mc);
    double csi  = xa / (2. - xa);
    double tmin = 4. * mass * mass * csi * csi / (1. - csi * csi);
    return std::max(tmin, tmi);
}

auto TOPEG::Model131V2_EventGeneratorFOAM::Model131(const std::array<double, m_nDim> &mc) -> double {

    double q      = mc[1];
    double xb     = q / (2 * Constant::NUCLEON_MASS * EBene * mc[0]);
    double t      = mc[2];
    double f      = mc[3] / 57.2957795;
    double ps     = mc[4];
    double thetas = mc[5] / 57.2957795;
    double crosssect = 0;

    return crosssect;
}
