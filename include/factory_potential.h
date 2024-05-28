#ifndef FACTORY_POTENTIAL_H 
#define FACTORY_POTENTIAL_H

#include <functional>
#include <map>
#include <memory>
#include <string>
#include "setup.h"

class PotentialFactory {
public:
    static PotentialFactory& instance() {
        static PotentialFactory factory;
        return factory;
    }

    void registerCreator(const std::string& name, std::function<PotentialBase*()> creator) {
        creators[name] = creator;
    }

    PotentialBase* create(const std::string& name) const {
        auto it = creators.find(name);
        if (it != creators.end()) {
            return it->second();
        }
        return nullptr;
    }

private:
    PotentialFactory() = default;
    std::map<std::string, std::function<PotentialBase*()>> creators;
};

#define GET_SETUP() auto& setup = Setup::instance();
#define NO_SETUP() do {} while(0);

#define REGISTER_INTERNAL_POTENTIAL(NAME, TYPE, WITH_PARAMS, ...) \
    struct TYPE ## _Internal_Potential_Register { \
        TYPE ## _Internal_Potential_Register() { \
            PotentialFactory::instance().registerCreator(NAME, []() { \
                WITH_PARAMS \
                return new TYPE(__VA_ARGS__); \
            }); \
        } \
    }; \
    static TYPE ## _Internal_Potential_Register global_##TYPE##_internal_potential_register;

#define REGISTER_EXTERNAL_POTENTIAL(NAME, TYPE, WITH_PARAMS, ...) \
    struct TYPE ## _External_Potential_Register { \
        TYPE ## _External_Potential_Register() { \
            PotentialFactory::instance().registerCreator(NAME, []() { \
                WITH_PARAMS \
                return new TYPE(__VA_ARGS__); \
            }); \
        } \
    }; \
    static TYPE ## _External_Potential_Register global_##TYPE##_external_potential_register;

#endif
