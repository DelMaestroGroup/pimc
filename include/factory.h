/**
 * @file factory.h
 * @author Adrian Del Maestro
 * @date 06.13.2016
 *
 * @brief Factory class definitions.
 */

#ifndef FACTORY_H 
#define FACTORY_H

/** 
 * An abstract factory class which creates new object instances based on a std::string
 * descripter and constructor signature.
 *
 * We use the factory method design pattern 
 * (http://en.wikipedia.org/wiki/Factory_(object-oriented_programming)) to
 * create a new instance of an estimator object.
 *
 * Design borrowed heavily from: 
 * @see http://www.gamedev.net/page/resources/_/technical/general-programming/creating-a-generic-object-factory-r2097
 */

template <typename CtorSignature> class Factory;

template<class BaseType, class ...ParamType>
class Factory<BaseType (ParamType...)>
{
    protected:
        /* The constructor signature (*) refers to a function pointer (the constructor) */
        using CreateObjectFunc = BaseType (*)(ParamType...);

    public:

        /** The names of all objects instatiated by the factory */
        std::vector<std::string> getNames() const {
            std::vector<std::string> names;
            for(auto const& createIter: _create)
                names.push_back(createIter.first);
            return names;
        }

        /** Singleton access */
        Factory<BaseType (ParamType...)> *getInstance()
        {
            static Factory<BaseType (ParamType...)> fact;
            return &fact;
        
        }

        /** Overload () to return a singleton instance */
        const Factory<BaseType (ParamType...)> * operator() () const { return getInstance();}
        Factory<BaseType (ParamType...)> * operator() () { return getInstance();}

        /** Return an instantiated object with a given name */
        BaseType Create(std::string name, ParamType ...param) {
            typename std::map<std::string,CreateObjectFunc>::const_iterator objItr = _create.find(name);
            if (objItr != _create.end()) {
                return (objItr->second)(param...);
                /* auto created = (objItr->second)(param...); */
                /* created->name1 = name; */
                /* return created; */
            }
            return nullptr;
        }

        /** Register the derived type with a descriptive name in the std::map */
        template<class DerivedType>
            bool Register(std::string name)
            {
                _create[name] = &createObj<DerivedType>;
                return true;
            }

    protected:
        std::map<std::string,CreateObjectFunc> _create;        // The name->constructor std::map

        /* Create the new object */
        template <class DerivedType>
            static BaseType createObj(ParamType ...param) 
            { return new DerivedType(param...);}

};

class EstimatorBase;
class MoveBase;
class ActionBase;
class WaveFunctionBase;
class Path;
class MTRand;
class LookupTable;

/* Typedefs used for actually creating factories */
typedef Factory<EstimatorBase* (Path &, ActionBase *, MTRand &, double)> EstimatorFactory;
typedef Factory<EstimatorBase* (Path &, Path &, ActionBase *, ActionBase *, MTRand &, double)> MultiEstimatorFactory;
typedef Factory<MoveBase* (Path &, ActionBase *, MTRand &)> MoveFactory;
typedef Factory<WaveFunctionBase* (const Path &, LookupTable &)> WaveFunctionFactory;

/* template<typename BaseType, class DerivedType, class ...ParamType> */
/* BaseType CreateObject(ParamType ...param) */
/* { */
/*    return new DerivedType(param...); */
/* } */

#endif
