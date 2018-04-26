#ifndef ZSVM_ENUM_TYPES_HPP_INCLUDED
#define ZSVM_ENUM_TYPES_HPP_INCLUDED

// C++ standard library headers
#include <map>
#include <string>

namespace zsvm {

    struct ExchangeStatistics {

        enum class Type {
            BOSON, FERMION, DISTINGUISHABLE
        }; // enum class Type

        const Type type;

        explicit ExchangeStatistics(Type type) noexcept
                : type(type) {}

        static const std::map<std::string, Type> &map() {
            static const std::map<std::string, Type> MAP{
                    {"boson",           Type::BOSON},
                    {"fermion",         Type::FERMION},
                    {"none",            Type::DISTINGUISHABLE},
                    {"distinguishable", Type::DISTINGUISHABLE},
            };
            return MAP;
        };

    }; // struct ExchangeStatistics

    template <typename T>
    struct DispersionRelation {

        enum class Type {
            RADIAL_POWER_LAW
        }; // enum class Type

        const Type type;
        const T strength;
        const T exponent;
        const std::string carrier;

        explicit DispersionRelation(
                Type type, T strength, T exponent, std::string carrier)
                : type(type), strength(strength),
                  exponent(exponent), carrier(std::move(carrier)) {}

        static const std::map<std::string, Type> &map() {
            static const std::map<std::string, Type> MAP{
                    {"radial_power_law", Type::RADIAL_POWER_LAW},
            };
            return MAP;
        };

    }; // struct DispersionRelation

    template <typename T>
    struct ConfiningPotential {

        enum class Type {
            RADIAL_POWER_LAW
        }; // enum class Type

        const Type type;
        const T strength;
        const T exponent;
        const std::string carrier;

        explicit ConfiningPotential(
                Type type, T strength, T exponent, std::string carrier)
                : type(type), strength(strength),
                  exponent(exponent), carrier(std::move(carrier)) {}

        static const std::map<std::string, Type> &map() {
            static const std::map<std::string, Type> MAP{
                    {"radial_power_law", Type::RADIAL_POWER_LAW},
            };
            return MAP;
        };

    }; // struct ConfiningPotential

    template <typename T>
    struct PairwisePotential {

        enum class Type {
            RADIAL_POWER_LAW
        }; // enum class Type

        const Type type;
        const T strength;
        const T exponent;
        const std::string carrier;

        explicit PairwisePotential(
                Type type, T strength, T exponent, std::string carrier)
                : type(type), strength(strength),
                  exponent(exponent), carrier(std::move(carrier)) {}

        static const std::map<std::string, Type> &map() {
            static const std::map<std::string, Type> MAP{
                    {"radial_power_law", Type::RADIAL_POWER_LAW},
            };
            return MAP;
        };

    }; // struct PairwisePotential

} // namespace zsvm

#endif // ZSVM_ENUM_TYPES_HPP_INCLUDED
