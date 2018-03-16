#include "EnumTypes.hpp"


zsvm::ExchangeStatistics::ExchangeStatistics(Type type) noexcept
        : type(type) {}


const std::map<std::string, zsvm::ExchangeStatistics::Type>
        zsvm::ExchangeStatistics::MAP = {
        {"boson",           Type::BOSON},
        {"fermion",         Type::FERMION},
        {"none",            Type::DISTINGUISHABLE},
        {"distinguishable", Type::DISTINGUISHABLE},
};


zsvm::DispersionRelation::DispersionRelation(
        Type type, double strength, double exponent,
        const std::string &carrier)
        : type(type), strength(strength),
          exponent(exponent), carrier(carrier) {}


const std::map<std::string, zsvm::DispersionRelation::Type>
        zsvm::DispersionRelation::MAP = {
        {"radial_power_law", Type::RADIAL_POWER_LAW},
};


zsvm::ConfiningPotential::ConfiningPotential(
        Type type, double strength, double exponent,
        const std::string &carrier)
        : type(type), strength(strength),
          exponent(exponent), carrier(carrier) {}


const std::map<std::string, zsvm::ConfiningPotential::Type>
        zsvm::ConfiningPotential::MAP = {
        {"radial_power_law", Type::RADIAL_POWER_LAW},
};


zsvm::PairwisePotential::PairwisePotential(
        Type type, double strength, double exponent,
        const std::string &carrier)
        : type(type), strength(strength),
          exponent(exponent), carrier(carrier) {}


const std::map<std::string, zsvm::PairwisePotential::Type>
        zsvm::PairwisePotential::MAP = {
        {"radial_power_law", Type::RADIAL_POWER_LAW},
};
