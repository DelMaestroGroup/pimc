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
        Internal,
        External
    };

    static PotentialFactory& instance() {
        static PotentialFactory factory;
        return factory;
    }

    template<Type T>
    void registerCreator(const std::string& name, std::function<PotentialBase*()> creator) {
        if constexpr (T == Type::Internal) {
            creatorsInternal[name] = creator;
        } else {
            creatorsExternal[name] = creator;
        }
    }

    template<Type T>
    PotentialBase* create(const std::string& name) const {
        const auto& creators = (T == Type::Internal) ? creatorsInternal : creatorsExternal;
        auto it = creators.find(name);
        if (it != creators.end()) {
            return it->second();
        }
        return nullptr;
    }

private:
    PotentialFactory() = default;
    std::map<std::string, std::function<PotentialBase*()>> creatorsInternal;
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

#define REGISTER_INTERNAL_POTENTIAL(NAME, TYPE, WITH_PARAMS, ...) \
    REGISTER_POTENTIAL(NAME, TYPE, Internal, WITH_PARAMS, __VA_ARGS__)

#define REGISTER_EXTERNAL_POTENTIAL(NAME, TYPE, WITH_PARAMS, ...) \
    REGISTER_POTENTIAL(NAME, TYPE, External, WITH_PARAMS, __VA_ARGS__)

#endif
