#ifndef ROOT_TEST_MODEL131V2_HPP
#define ROOT_TEST_MODEL131V2_HPP

#include <Constant.hpp>
#include <Model.hpp>
#include <TFoamIntegrand.h>
#include <array>
#include <chrono>
#include <complex>
#include <fmt/core.h>
#include <iostream>
#include <thread>

namespace TOPEG {
    class Model131V2_EventGeneratorFOAM : public Model_EventGeneratorFOAM {
    private:
        static constexpr int m_nDim = 7;
        // Parameters from table 4.1 in Moh.thesis
        const double alp  = 2.5;
        const double q0   = 1.0;
        const double d    = 0.4;
        const double xc   = 0.2;
        const double c    = 0.2;
        const double b    = 11.0;
        const double beta = 12.0;

        const double xM = 1.5;// TODO : in the config file
        double EBene = 0, HBene = 0, HBmom = 0, EBpol = 0, s = 0, mass = 0;
        double ym = 0, yM = 0, qm = 0, qM = 0, Wm = 0, TM = 0, tmi = 0, tR = 0;
        double psM = 1, psm = 0, thetasM = 180, thetasm = 0, phisM = 360, phism = 0;

        std::atomic<int> counter = 0;

        std::chrono::time_point<std::chrono::high_resolution_clock> TimeC = std::chrono::high_resolution_clock::now();

        std::array<double, m_nDim> MCV = {};
        std::array<double, m_nDim> ExV = {};

        template<class... Durations, class DurationIn>
        std::tuple<Durations...> break_down_durations(DurationIn d) {
            std::tuple<Durations...> output;
            using discard = int[];
            (void) discard{0, (void(((std::get<Durations>(output) = std::chrono::duration_cast<Durations>(d)), (d -= std::chrono::duration_cast<DurationIn>(std::get<Durations>(output))))), 0)...};
            return output;
        }

        auto getqM(const std::array<double, m_nDim> &mc) -> double;
        auto getxa(const std::array<double, m_nDim> &mc) -> double;
        auto gettm(const std::array<double, m_nDim> &mc) -> double;
        auto setKinLimits() -> void;

        auto deut_mom(double x) -> double;

    public:
        Model131V2_EventGeneratorFOAM(double b1, double b2, double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9);
        auto Transform(std::vector<double> &mc) -> void override {}
        auto Transform() -> void override {}
        auto Transform(double *event) -> void override;
        auto Transform(std::array<double, m_nDim> &mc) -> void;
        auto Density(int nDim, Double_t *V) -> Double_t override;
        auto Model131(const std::array<double, m_nDim> &mc) -> double;

        auto getDimNb() -> int override { return m_nDim; }
        auto GetEBene() -> double override { return EBene; }
        auto GetHBmom() -> double override { return HBmom; };
        auto GetHBene() -> double override { return HBene; };
        auto Getmass() -> double override { return mass; };
        auto GetS() -> double override { return s; };
        auto GetDim() -> int override { return m_nDim; };
    };
}// namespace TOPEG
#endif//ROOT_TEST_MODEL131_HPP
