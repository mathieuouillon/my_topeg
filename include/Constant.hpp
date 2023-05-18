#ifndef ROOT_TEST_CONSTANT_HPP
#define ROOT_TEST_CONSTANT_HPP

namespace TOPEG {

    struct Constant {
        inline static constexpr double PROTON_MASS    = 0.9382720882;                    // in GeV
        inline static constexpr double NEUTRON_MASS   = 0.9395654205;                    // in GeV
        inline static constexpr double NUCLEON_MASS   = (PROTON_MASS + NEUTRON_MASS) / 2;// in GeV
        inline static constexpr double DEUTERIUM_MASS = 1.87561294257;                   // in GeV
    };

}// namespace TOPEG

#endif//ROOT_TEST_CONSTANT_HPP
