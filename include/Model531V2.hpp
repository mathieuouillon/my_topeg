#ifndef ROOT_TEST_MODEL531V2_HPP
#define ROOT_TEST_MODEL531V2_HPP

#include <Constant.hpp>
#include <Math/AxisAngle.h>
#include <Math/Boost.h>
#include <Math/Rotation3D.h>
#include <Math/RotationX.h>
#include <Math/RotationY.h>
#include <Math/RotationZ.h>
#include <Math/RotationZYX.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <Model.hpp>
#include <TFoamIntegrand.h>
#include <TMatrixD.h>
#include <array>
#include <chrono>
#include <complex>
#include <fmt/core.h>
#include <iostream>
#include <memory>
#include <mutex>
#include <thread>
#include <tuple>
#include <vector>


extern "C" {
void cross_boundnucleon_deut_f90_(double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, double *, double *);
}

namespace TOPEG {

    class Model531_v2 : public Model_EventGeneratorFOAM {
    private:
        const Int_t fDim = 7;

        std::atomic<Int_t> count = 0;

        Double_t fMassTarget    = 0;
        Double_t fMassActive    = 0;
        Double_t fMassSpectator = 0;

        Double_t E_had_beam = 0;

        InitialsConditions fInitialsConditions;
        KinematicsRange fKinematicsRange;

        std::chrono::time_point<std::chrono::high_resolution_clock> timeCount = std::chrono::high_resolution_clock::now();

        // Functions
        [[nodiscard]] auto Model531(const std::vector<Double_t> &mc) const -> Double_t;
        [[nodiscard]] auto Q2Max(const std::vector<Double_t> &mc) const -> Double_t;
        [[nodiscard]] auto tMin(const std::vector<Double_t> &mc) const -> Double_t;

        template<class... Durations, class DurationIn>
        std::tuple<Durations...> break_down_durations(DurationIn d) {
            std::tuple<Durations...> output;
            using discard = int[];
            (void) discard{0, (void(((std::get<Durations>(output) = std::chrono::duration_cast<Durations>(d)), (d -= std::chrono::duration_cast<DurationIn>(std::get<Durations>(output))))), 0)...};
            return output;
        }

    public:
        Model531_v2(const InitialsConditions &initialsConditions, const KinematicsRange &kinematicsRange);
        auto Density(int nDim, Double_t *V) -> Double_t override;
        auto Transform(double *event) -> void override {}
        auto Transform() -> void override {}
        auto Transform(std::vector<double> &mc) -> void override;


        auto getDimNb() -> int override { return fDim; }
        auto GetEBene() -> double override { return fInitialsConditions.E_el_beam; }
        auto GetHBmom() -> double override { return fInitialsConditions.p_had_beam; };
        auto GetHBene() -> double override { return E_had_beam; };
        auto Getmass() -> double override { return fMassTarget; };
        auto GetS() -> double override { return 0; };
        auto GetDim() -> int override { return fDim; };
    };

}// namespace TOPEG

#endif//ROOT_TEST_MODEL531V2_HPP
