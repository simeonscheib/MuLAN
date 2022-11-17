#ifndef UTOPIA_MODELS_MULAN_MA_HH_ORG_MANAGER
#define UTOPIA_MODELS_MULAN_MA_HH_ORG_MANAGER

#include "Organism_Factory.hh"

/**
 * @brief The OrganismManager manages the species in a given domain
 * @details It allows to register species templates at compile time
 *          and chooses a specialization at run time from configuration
 *          
 * 
 * @tparam DOM_T    Domain type
 * @tparam FUN      Wrapper that holds the signature of the "add" function of DOM_T
 *                  and returns a function pointer to the "add" function via operator()
 * @tparam ORG_Ts   The species templates 
 *                  <class DOM, class PSPACE, int... CFGs>
 *                  that can be added to the domain
 */
template <typename DOM_T, template <typename> typename FUN, template<int ...> class ...ORG_Ts>
class OrganismManager {
private:

    //using FUNC = std::_Mem_fn<int& (DOM_T::*)(double, double)>;
    using FUNC = typename FUN<DOM_T>::type;
    //using FUNC = std::_Mem_fn<FUN>;
    DOM_T& dom;

    /**
     * @brief Hold the number of species registered with the manager
     * 
     */
    static const std::size_t num_species = sizeof...(ORG_Ts);

public:

    /**
     * @brief Type of Domain
     */
    using dom_t = DOM_T;

    /**
     * @brief Construct a new Organism Manager object
     * 
     * @param dom   The domain
     */
    explicit OrganismManager(DOM_T& dom) : dom(dom) {};

    /**
     * @brief   Wrapper class to build OrganismAdder for each 
     *          organism template
     * 
     * @tparam ORG_T    species template
     */
    template < template<int ...> class ORG_T >
    class Add {
    public:

        /**
         * @brief Interface for OrganismAdder
         * 
         */
        class OrganismAdder {
        public:

            /**
             * @brief   Add an organism of template type ORG_T
             * @details The type of the added organism depends on the 
             *          run time configuration provided at organism registration
             *          via `OrganismManager::register_`
             * 
             * @param domain    Reference to the domain
             * @param mass      Initial mass
             * @param pspace    Parameter space
             * @return FUNC     Function that creates the new vertex
             *                  or returns existing vertex with same parameter space
             */
            virtual inline FUNC add() = 0;
        
            /**
             * @brief Get the organism count of specialized species
             * 
             * @return int& Reference to the static counter
             */
            virtual int& get_spec_count() = 0;
        };


        template <int ...Ns>
        class OrganismAdderSpec : public OrganismAdder {
        public:
            /**
             * @brief The specialized species type
             * 
             */
            using organism_t = ORG_T<Ns...>;

            using f_t = FUN<DOM_T>;

            f_t f;

            virtual inline FUNC add(){

                return this->f.template operator()<organism_t>();
            };

            virtual int& get_spec_count(){
                return organism_t::spec_count;
            }
        };

    };

protected:
    /**
     * @brief Tuple of all registered species templates
     * @details Holds a pointer to an Organism Adder for each species
     * 
     */
    std::tuple<typename Add< ORG_Ts>::OrganismAdder*...> adders;
public:

    /**
     * @brief Get the number of species registered
     * 
     * @return size_t 
     */
    static constexpr size_t get_number_of_species() {
        return num_species;
    }

    /**
     * @brief   Finalize registration process of a species by providing a
     *          run time configuration
     * @details Choose a specialization of the species template according to
     *          configuration
     * 
     * @tparam ORG_T    The species template
     * @tparam MAXs     The numbers of config options
     * @param config_v  Vector of ints according to configuration
     * @return auto     Return pointer to the OrganismAdder
     */
    template <template<int ...> class ORG_T, int ...MAXs>
    auto register_ (std::vector<int> config_v) {

        // Choose OrganismAdder according to provided species
        auto& a = std::get<typename Add< ORG_T>::OrganismAdder*>(adders);

        // Check if the species is already registered
        if (a == nullptr) {

            a = OrganismFactory<typename Add< ORG_T>::OrganismAdder, 
                                                    Add< ORG_T>::template OrganismAdderSpec,  
                                                    IntTuple< MAXs... >, 
                                                    IntTuple<> >::create(config_v);
        } else {
            std::cout << "Could not register organism. Already registered" << std::endl;
        }

        return a;
    };

    /**
     * @brief Add an organism of species ORG_T using the previous registration
     * 
     * @tparam ORG_T    Species
     * @param mass      Initial mass
     * @param pspace    Parameter space
     * @return auto     return the newly created vertex or existing
     */
    template <template<int ...> class ORG_T, typename ...Ts>
    inline auto add_(Ts... ts) {

        auto& a = std::get<typename Add< ORG_T>::OrganismAdder*>(adders);

        if (a == nullptr) {
            std::cout << "Could not add organism. Not registered" << std::endl;
        }

        //return a->add(this->dom, ts...);
        auto f = a->add();

        return f(this->dom, ts...);

    };

};

#endif // UTOPIA_MODELS_MULAN_MA_HH_ORG_MANAGER