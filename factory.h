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
 * An abstract factory class which creates new object instances based on a string
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
        /* The constructor signature */
        using CreateObjectFunc = BaseType (*)(ParamType...);

    public:
        /** The names of all objects instatiated by the factory */
        static vector<string> names;

        vector<string> getNames() const {return names;};

        /** Singleton access */
        Factory<BaseType (ParamType...)> *Instance()
        {
            static Factory<BaseType (ParamType...)> fact;
            return &fact;
        
        }

        /** Overload () to return a singleton instance */
        const Factory<BaseType (ParamType...)> * operator() () const { return Instance();}
        Factory<BaseType (ParamType...)> * operator() () { return Instance();}

        /** Return an instantiated object with a given name */
        BaseType Create(string name, ParamType ...param) {
            typename map<string,CreateObjectFunc>::const_iterator objItr = _create.find(name);
            if (objItr != _create.end()) {
                return (objItr->second)(param...);
                /* auto created = (objItr->second)(param...); */
                /* created->name1 = name; */
                /* return created; */
            }
            return nullptr;
        }

        /** Register the derived type with a descriptive name in the map */
        template<class DerivedType>
            bool Register(string name)
            {
                _create[name] = &createObj<DerivedType>;

                // Add the name to the list
                // NB: this had been commented out with an external init
                // function but this seems to be working now.
                names.push_back(name);
                return true;
            }

    protected:
        map<string,CreateObjectFunc> _create;        // The name->constructor map

        /* Create the new object */
        template <class DerivedType>
            static BaseType createObj(ParamType ...param) 
            { return new DerivedType(param...);}

};

template<class BaseType, class ...ParamType>
vector<string> Factory<BaseType (ParamType...)>::names;


class EstimatorBase;
class MoveBase;
class ActionBase;
class Path;
class MTRand;

/* Typedefs used for actually creating factories */
typedef Factory<EstimatorBase* (Path &, ActionBase *, MTRand &, double)> EstimatorFactory;
typedef Factory<EstimatorBase* (Path &, Path &, ActionBase *, ActionBase *, MTRand &, double)> MultiEstimatorFactory;
typedef Factory<MoveBase* (Path &, ActionBase *, MTRand &)> MoveFactory;

/* template<typename BaseType, class DerivedType, class ...ParamType> */
/* BaseType CreateObject(ParamType ...param) */
/* { */
/*    return new DerivedType(param...); */
/* } */

#endif
