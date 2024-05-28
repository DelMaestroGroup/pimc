#ifndef FACTORY_POTENTIAL_H 
#define FACTORY_POTENTIAL_H

#include <functional>
#include <map>
#include <memory>
#include <string>

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

#define REGISTER_INTERNAL_POTENTIAL(NAME, TYPE, ...) \
    const std::string TYPE::name = NAME; \
    struct TYPE ## _Register { \
        TYPE ## _Register() { \
            PotentialFactory::instance().registerCreator(NAME, []() { \
                auto& setup = Setup::instance(); \
                return new TYPE(__VA_ARGS__); \
            }); \
        } \
    }; \
    static TYPE ## _Register global_##TYPE##_register;

#define REGISTER_EXTERNAL_POTENTIAL(NAME, TYPE, ...) \
    const std::string TYPE::name = NAME; \
    struct TYPE ## _Register { \
        TYPE ## _Register() { \
            PotentialFactory::instance().registerCreator(NAME, []() { \
                auto& setup = Setup::instance(); \
                return new TYPE(__VA_ARGS__); \
            }); \
        } \
    }; \
    static TYPE ## _Register global_##TYPE##_register;

#endif
