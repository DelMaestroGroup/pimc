#ifndef FACTORY_POTENTIAL_H 
#define FACTORY_POTENTIAL_H

#include <functional>
#include <map>
#include <memory>
#include <string>
#include "setup.h"

class PotentialFactory {
public:
    enum class Type {
        Interaction,
        External
    };

    static PotentialFactory& instance() {
        static PotentialFactory factory;
        return factory;
    }

    template<Type T>
    void registerCreator(const std::string& name, std::function<PotentialBase*()> creator) {
        if constexpr (T == Type::Interaction) {
            creatorsInteraction[name] = creator;
        } else {
            creatorsExternal[name] = creator;
        }
    }

    template<Type T>
    PotentialBase* create(const std::string& name) const {
        const auto& creators = (T == Type::Interaction) ? creatorsInteraction : creatorsExternal;
        auto it = creators.find(name);
        if (it != creators.end()) {
            return it->second();
        }
        return nullptr;
    }

    template<Type T>
    std::vector<std::string> getNames() const {
        std::vector<std::string> names;
        const auto& creators = (T == Type::Interaction) ? creatorsInteraction : creatorsExternal;
        for (const auto& pair : creators) {
            names.push_back(pair.first);
        }
        return names;
    }

private:
    PotentialFactory() = default;
    std::map<std::string, std::function<PotentialBase*()>> creatorsInteraction;
    std::map<std::string, std::function<PotentialBase*()>> creatorsExternal;
};

#define GET_SETUP() auto& setup = Setup::instance();
#define NO_SETUP() do {} while(0);

#define REGISTER_POTENTIAL(NAME, TYPE, POTENTIAL_TYPE, WITH_PARAMS, ...) \
    struct TYPE ## _ ## POTENTIAL_TYPE ## _Potential_Register { \
        TYPE ## _ ## POTENTIAL_TYPE ## _Potential_Register() { \
            PotentialFactory::instance().registerCreator<PotentialFactory::Type::POTENTIAL_TYPE>(NAME, []() { \
                WITH_PARAMS \
                return new TYPE(__VA_ARGS__); \
            }); \
        } \
    }; \
    static TYPE ## _ ## POTENTIAL_TYPE ## _Potential_Register global_##TYPE##_##POTENTIAL_TYPE##_potential_register;

#define REGISTER_INTERACTION_POTENTIAL(NAME, TYPE, WITH_PARAMS, ...) \
    REGISTER_POTENTIAL(NAME, TYPE, Interaction, WITH_PARAMS, __VA_ARGS__)

#define REGISTER_EXTERNAL_POTENTIAL(NAME, TYPE, WITH_PARAMS, ...) \
    REGISTER_POTENTIAL(NAME, TYPE, External, WITH_PARAMS, __VA_ARGS__)

#endif
