#ifndef ROOT_TEST_MODELT_HPP
#define ROOT_TEST_MODELT_HPP

#include <TFoamIntegrand.h>
#include <array>

struct InitialsConditions {
    Double_t E_el_beam   = 0;
    Double_t pol_el_beam = 0;
    Double_t p_had_beam  = 0;

    InitialsConditions(Double_t tE_el_beam, Double_t tpol_el_beam, Double_t tp_had_beam) : E_el_beam(tE_el_beam), pol_el_beam(tpol_el_beam), p_had_beam(tp_had_beam) {}
};

struct KinematicsRange {
    Double_t y_min             = 0.;
    Double_t y_max             = 0.;
    Double_t q_min             = 0.;
    Double_t q_max             = 0.;
    Double_t W2_min            = 0.;
    Double_t theta_el_scat_max = 0.;
    Double_t t_min             = 0.;
    Double_t t_range           = 0.;
    Double_t xBj_max           = 1.5;
    Double_t p_spec_max        = 0.25;
    Double_t p_spec_min        = 0.04;
    Double_t theta_spec_max    = 160.;
    Double_t theta_spec_min    = 20.;

    KinematicsRange(Double_t y_min, Double_t y_max, Double_t q_min, Double_t q_max, Double_t W2_min, Double_t theta_el_scat_max, Double_t t_min, Double_t t_range)
        : y_min(y_min), y_max(y_max), q_min(q_min), q_max(q_max), W2_min(W2_min), theta_el_scat_max(theta_el_scat_max), t_min(t_min), t_range(t_range) {}
};

namespace TOPEG {

    class Model_EventGeneratorFOAM : public TFoamIntegrand {
    private:
        static constexpr int m_dim = 0;

    public:
        virtual auto Transform() -> void                        = 0;
        virtual auto Transform(double *event) -> void           = 0;
        virtual auto Transform(std::vector<double> &mc) -> void = 0;
        virtual auto getDimNb() -> int                          = 0;
        virtual auto GetEBene() -> double                       = 0;
        virtual auto GetHBmom() -> double                       = 0;
        virtual auto GetHBene() -> double                       = 0;
        virtual auto Getmass() -> double                        = 0;
        virtual auto GetS() -> double                           = 0;
        virtual auto GetDim() -> int                            = 0;
    };

}// namespace TOPEG

#endif//ROOT_TEST_MODELT_HPP
