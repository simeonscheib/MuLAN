#include <boost/exception/detail/shared_ptr.hpp>
#include <boost/unordered_map.hpp>

#ifndef UTOPIA_MODELS_ORGANISM_BASE_HH
#define UTOPIA_MODELS_ORGANISM_BASE_HH


namespace Utopia::Models::MuLAN_MA {

    /**
     * @brief Organism Base
     * @details Base class for all kinds of organisms
     * 
     * @tparam ORG_T Child Type
     * @tparam DOM_T Domain Type
     * @tparam PSPACE_T Parameter Space Type
     * @tparam INTEGRATOR Integrator Type
     */
    template <typename ORG_T, typename DOM_T, typename PSPACE_T, typename INTEGRATOR>
    class organism_base {
    protected:
        /**
         * @brief Domain
         * @details Back-pointer to domain allows usage of domain parameters in DEs
         * 
         */
        [[maybe_unused]] DOM_T* dom;

        // Every organism has a change in Biomass during a time differential
        // virtual CURRENCY dxdt() = 0;

    public:
        using org_ptr = boost::shared_ptr<ORG_T>;
        using edge_cont = typename DOM_T::edge_cont;

        using pspace_t = PSPACE_T;
        using currency = typename INTEGRATOR::currency;



        const static int type_id = 0;

        /**
         * @brief Count number of species
         * 
         */
        [[maybe_unused]] static int spec_count;

        /**
         * @brief Parameter Space
         * @details Every organism has its parameter space
         * 
         */
        PSPACE_T parameters;

        /**
         * @brief The Integrator 
         * @details Manages Integration
         * 
         */
        INTEGRATOR integrator;


        /**
         * @brief Construct a new organism base object
         * 
         * @param d Domain
         * @param org This Child Object
         */
        explicit organism_base(DOM_T* d, ORG_T& org) : integrator(org){
            this->dom = d;
            

        };

        /**
         * @brief Destroy the organism base object
         * 
         */
        ~organism_base() = default;

        /**
         * @brief Get the value
         * @details E.g. Biomass of the organism
         * 
         * @return currency& 
         */
        virtual currency& get_value() {
            return this->integrator.get_value();
        };

        /**
         * @brief Set the value
         * 
         * @param val 
         */
        virtual void set_value(const currency& val) {
            return this->integrator.set_value(val);
        };


        /**
         * @brief Calculate the interaction of two organisms
         * 
         * @param org2 The second organism
         * @return edge_cont 
         */
        virtual edge_cont calc_interaction_coeff(const org_ptr org2) const = 0;

        /**
         * @brief Add the interaction to the sum vector
         * 
         * @param val 
         * @param organism_2 
         */

        template<int N>
        void add_edge_cont(const edge_cont& val, org_ptr organism_2) {
            static_cast<ORG_T&>(*this).template add_edge_cont<N>(val, organism_2);
        };

        /**
         * @brief Change count for organism
         * 
         * @param c 
         */
        virtual void count(int c) = 0;

    };




} // namespace Utopia::Models::MuLAN_MA

#endif //UTOPIA_MODELS_ORGANISM_BASE_HH
