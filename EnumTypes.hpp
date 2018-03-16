#ifndef ZSVM_ENUM_TYPES_HPP_INCLUDED
#define ZSVM_ENUM_TYPES_HPP_INCLUDED

#include <map> // for std::map
#include <string> // for std::string

namespace zsvm {

    struct ExchangeStatistics {

        enum class Type {
            BOSON, FERMION, DISTINGUISHABLE
        }; // enum class Type

        const Type type;

        explicit ExchangeStatistics(Type type) noexcept;

        static const std::map<std::string, Type> MAP;

    }; // struct ExchangeStatistics

    struct DispersionRelation {

        enum class Type {
            RADIAL_POWER_LAW
        }; // enum class Type

        const Type type;
        const double strength;
        const double exponent;
        const std::string carrier;

        explicit DispersionRelation(
                Type type, double strength, double exponent,
                const std::string &carrier);

        static const std::map<std::string, Type> MAP;

    }; // struct DispersionRelation

    struct ConfiningPotential {

        enum class Type {
            RADIAL_POWER_LAW
        }; // enum class Type

        const Type type;
        const double strength;
        const double exponent;
        const std::string carrier;

        explicit ConfiningPotential(
                Type type, double strength, double exponent,
                const std::string &carrier);

        static const std::map<std::string, Type> MAP;

    }; // struct ConfiningPotential

    struct PairwisePotential {

        enum class Type {
            RADIAL_POWER_LAW
        }; // enum class Type

        const Type type;
        const double strength;
        const double exponent;
        const std::string carrier;

        explicit PairwisePotential(
                Type type, double strength, double exponent,
                const std::string &carrier);

        static const std::map<std::string, Type> MAP;

    }; // struct PairwisePotential

} // namespace zsvm

#endif // ZSVM_ENUM_TYPES_HPP_INCLUDED
