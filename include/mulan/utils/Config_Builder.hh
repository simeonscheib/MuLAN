#ifndef UTOPIA_MODELS_MULAN_MA_HH_CONFIG_BUILDER
#define UTOPIA_MODELS_MULAN_MA_HH_CONFIG_BUILDER

#include <utopia/core/model.hh>
#include "helper.hh"

using namespace Utopia::Models::MuLAN_MA;

/**
 * @brief Throw an error if config does not match Factory parameters
 * @details If e.g. config value is not in specified range
 * 
 */
struct FactoryError
{
    template<typename ParentType>
    [[maybe_unused]] static inline void create(const std::vector<int>& /*c*/, ParentType& /*pp*/)
    {
        throw std::runtime_error ("Could not create domain");
    }
};


/**
 * @brief Create any possible model to choose from according to config
 * @details See specialization. 
 *          WARNING!!! Compile time scales exponential with 
 *          #template parameters of domain
 * 
 * @tparam  
 */
template < typename ... >
struct Factory {};


/**
 * @brief Specialization: Throw Error if parameter reaches 0
 * @details If no parameter matches the value given in runtime config
 * 
 * @tparam As 
 * @tparam Bs 
 */
template<typename RUNNER, int ...As, int ...Bs>
struct Factory< RUNNER, IntTuple< 0, As... >, IntTuple< Bs...>> : FactoryError {};

/**
 * @brief Specialization: Exit condition
 * @details If the first tuple is empty every parameter 
 *          has been fixed by the runtime config => Success       
 * 
 * @tparam Bs The fixed parameters equal to runtime config
 */
template <typename RUNNER, int ...Bs >
struct Factory <RUNNER, IntTuple<>, IntTuple< Bs... > >  {

    template< typename ParentType>
    [[maybe_unused]] static inline void create(const std::vector<int>& /*c*/, ParentType& pp ) {
        ((std::cout << Bs),...);
        RUNNER::template run<ParentType, Bs...>(pp);
    }

};

/**
 * @brief Specialization: Build Factory recursively from runtime config
 * @details Lockpicking:  Iterate over all the parameters and fix them
 * 
 * @tparam A Value to fix
 * @tparam As Fixed values
 * @tparam Bs Remaining values
 */
template <typename RUNNER, int A, int ...As, int ...Bs >
struct Factory < RUNNER, IntTuple< A, As... >, IntTuple< Bs...> >  {

    template<typename ParentType>
    [[maybe_unused]] static inline void create(const std::vector<int>& c, ParentType& pp ) {

        if (c[sizeof...(Bs)] == A - 1){
            Factory < RUNNER, IntTuple< As... >, IntTuple< Bs..., A - 1> >::create(c, pp);
        } else {
            if (A == 0) {
                return;
            }

            Factory < RUNNER, IntTuple< A - 1, As... >, IntTuple< Bs...> >::create(c, pp);
        }
    }

};

/**
 * @brief Wrapper to call Factory and generate Model
 * 
 * @tparam ParentType Pseudo Parent of model
 * @tparam CFG_MAP Map config to ints
 */
template<typename ParentType, typename CFG_MAP>
class ModelBuilder {

private:

    template< template<typename...> typename MAP, typename KEY, typename VAL>
    inline KEY find_key(const MAP<KEY, int>& map, VAL val) {
        KEY key = map.begin()->first;
        for (auto &i : map) {
            if (i.second == val) {
                key = i.first;
                return key;
            }
        }
        std::cout << "Could not find key for " << val << std::endl;
        return "";
    }

    std::vector<int> config_v;

public:

    /**
     * @brief Construct a new Model Builder object 
     * 
     * @param pp Pseudo parent
     * @param map map config to ints
     * @param name Model Name
     */
    ModelBuilder(ParentType& pp, const CFG_MAP& map, const std::string& name) {

        auto config = pp.get_cfg()[name];

        for (auto c_tuple : map){
            auto cfg_val = std::get<0>(c_tuple);
            auto cfg_map = std::get<1>(c_tuple);
            auto keyword = Utopia::get_as<std::string>(cfg_val, config);
            if (cfg_map.count(keyword) == 1){
                std::cout << "Using " << cfg_val << ": " << keyword << std::endl;
            } else {
                std::cout << "No option: " << keyword << std::endl;
                std::cout << "Using Default: " << find_key(cfg_map, 0) << std::endl;
            }
            config_v.push_back(cfg_map[keyword]);
        }

    };

    /**
     * @brief Destroy the Model Builder object
     * 
     */
    ~ModelBuilder() = default;

    /**
     * @brief Call the factory with config
     * 
     * @tparam Ts The numbers of possible values the domain can receive
     * @param pp Pseudo parent
     */
    template<typename RUNNER, int ...Ts>
    void build(ParentType& pp)
    {
        Factory< RUNNER, IntTuple< Ts... >, IntTuple<> >::create(this->config_v, pp);
    }
};

#endif // UTOPIA_MODELS_MULAN_MA_HH_CONFIG_BUILDER
