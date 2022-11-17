#ifndef UTOPIA_MODELS_MULAN_MA_HH_ORGANISM_FACTORY
#define UTOPIA_MODELS_MULAN_MA_HH_ORGANISM_FACTORY

#include "helper.hh"

/**
 * @brief Throw an error if config does not match organismFactory parameters
 * @details If e.g. config value is not in specified range
 * 
 */
template <typename INTERFACE, template<int ...> typename >
struct OrganismFactoryError
{
    template<typename ...Ts>
    static inline INTERFACE* create(const std::vector<int>& /*c*/, Ts&&... /*ps*/)
    {
        throw std::runtime_error ("Could not create organism");
    }
};


/**
 * @brief Create any possible specialization of T according to config
 * @details See specialization. 
 *          WARNING!!! Compile time scales exponential with 
 *          #template parameters of T
 * 
 * @tparam INTERFACE The base class of template T
 * @tparam template <int ...CFGPs>T derived from INTERFACE 
 */
template < typename INTERFACE, template<int ...> typename T, typename ... >
struct OrganismFactory {};

/**
 * @brief Specialization: Throw Error if parameter reaches 0
 * @details If no parameter matches the value given in runtime config
 * 
 * @tparam INTERFACE The base class of template T
 * @tparam template <int ...CFGPs>T derived from INTERFACE 
 * @tparam As 
 * @tparam Bs 
 */
template<typename INTERFACE, template<int ...> typename T, int ...As, int ...Bs>
struct OrganismFactory<INTERFACE, T, IntTuple< 0, As... >, IntTuple< Bs...>> : OrganismFactoryError<INTERFACE, T> {};

/**
 * @brief Specialization: Exit condition
 * @details If the first tuple is empty, every parameter 
 *          has been fixed by the runtime config => Success       
 *
 * @tparam INTERFACE The base class of template T
 * @tparam template <int ...CFGPs>T derived from INTERFACE 
 * @tparam Bs The fixed parameters equal to runtime config
 */
template <typename INTERFACE, template<int ...> typename T, int ...Bs >
struct OrganismFactory <INTERFACE, T, IntTuple<>, IntTuple< Bs... > >  {

    template<typename ...Ts>
    [[maybe_unused]] static inline INTERFACE* create(const std::vector<int>& /*c*/, Ts&&... ps ) {
        return new T <Bs...>(ps...);
    }

};


/**
 * @brief Specialization: Build Factory recursively from runtime config
 * @details Lockpicking:  Iterate over all the parameters and fix them
 * 
 * @tparam INTERFACE The base class of template T
 * @tparam template <int ...CFGPs>T derived from INTERFACE 
 * @tparam A Value to fix
 * @tparam As Fixed values
 * @tparam Bs Remaining values
 */
template <typename INTERFACE, template<int .../*CFGPs*/> typename T, int A, int ...As, int ...Bs >
struct OrganismFactory < INTERFACE, T, IntTuple< A, As... >, IntTuple< Bs...> >  {

    /**
     * @brief 
     * 
     * @tparam Ts 
     * @param c config vector of ints
     * @param ps parameters to pass to constructor
     * @return INTERFACE* 
     */
    template<typename ...Ts>
    [[maybe_unused]] static inline INTERFACE* create(const std::vector<int>& c, Ts&&... ps) {

        if (c[sizeof...(Bs)] == A - 1){
            return OrganismFactory <INTERFACE, T, IntTuple< As... >, IntTuple< Bs..., A - 1> >::create(c, ps...);
        } else {

            return OrganismFactory <INTERFACE, T, IntTuple< A - 1, As... >, IntTuple< Bs...> >::create(c, ps...);
        }
    }

};


#endif // UTOPIA_MODELS_MULAN_MA_HH_ORGANISM_FACTORY