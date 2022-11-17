#ifndef UTOPIA_MODELS_PRODUCER_HH
#define UTOPIA_MODELS_PRODUCER_HH

#include "../organism.hh"
#include "consumer.hh"

namespace Utopia {
namespace Models {
namespace MuLAN_MA {

    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    class consumer;

    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    class primary_producer: public organism<DOM_T, PSPACE_T>{

    public:

        using pp_ptr [[maybe_unused]] = boost::shared_ptr<primary_producer>;
        using org_ptr = typename organism<DOM_T, PSPACE_T>::org_ptr;
        using currency = typename DOM_T::integrator_t::currency;

        // Interaction between producers CFGs[0]
        // 0: No interaction
        // 1: Interaction like Gause Lotka Volterra, Gaussian interaction kernel
        static constexpr int pp_interaction_type = get_from_pack<int, 0, CFGs...>();

        double dxdt(const currency& x, const double tpdt) override {
            double logistic_term;

            if constexpr ( DOM_T::env_func == 2 || DOM_T::env_func == 3 ) {
                this->parameters.params[1] = this->dom->S(this->parameters.trait[0], tpdt);
            } else {

            }

            if constexpr ( pp_interaction_type == 0 ) {
                logistic_term = x/this->parameters.params[1];
            } else if constexpr ( pp_interaction_type == 1 ) {
                logistic_term = this->integrator[1]/this->parameters.params[1];
            } else {
                throw std::invalid_argument("No valid 'pp_interaction'. Check config.");
            }
            double ret = x * ( this->parameters.params[0] * (1 - logistic_term) - this->integrator[0]);

            this->integrator.clear_sum();

            return ret;

        }

        // Producers ID: 0
        const static int type_id = 0;

        // Store number of species
        static int spec_count;

        /**
         * @brief Construct a producer with access to  domain d
         *
         * @param d Domain
         */
        explicit primary_producer(DOM_T* d) : organism<DOM_T, PSPACE_T>(d) {
            int sum_size;

            // Adapt sum_size
            if constexpr ( pp_interaction_type == 0 ){
                sum_size = 1;
            } else if constexpr ( pp_interaction_type == 1 ){
                sum_size = 2;
            } else {
                throw std::invalid_argument("No valid 'pp_interaction'. Check config.");
            }
            this->integrator.resize(sum_size);
        };

        /**
         * @brief Construct producer with access to  domain d and initial mass
         *
         * @param d             Domain
         * @param initial_mass  initial biomass
         */
        explicit primary_producer(DOM_T* d, double initial_mass) : primary_producer(d) {

            this->set_value(initial_mass);

        };

        /**
         * @brief destruct producer
         */
        ~primary_producer() = default;


        /**
         * @brief Calculate the interaction edge of this and org2
         *
         * @param org2  Organism
         * @return
         */
        typename DOM_T::edge_cont calc_interaction_coeff(const org_ptr org2) const override {
            if (org2->parameters.type == 0) {
                if constexpr (pp_interaction_type == 1 ) {

                    if (org2.get() == this) {
                        return typename DOM_T::edge_cont{1.0};
                    }

                } else {
                    return typename DOM_T::edge_cont{0.0};
                }
            }
            return typename DOM_T::edge_cont{0.0};
        }


        void add_edge_cont1(const typename DOM_T::edge_cont& val, org_ptr organism_2) override {
            if constexpr ( pp_interaction_type == 1 ) {
                const currency &mass = organism_2->get_value();
                organism_2->integrator[1] += val[0] * mass;
            } else {
                return;
            }
        }

        void add_edge_cont2(const typename DOM_T::edge_cont& /*val*/, const org_ptr /*organism_2*/) override {

            return;
        }

        // void set_parameters();

        virtual void count(int c) override{

            primary_producer<DOM_T, PSPACE_T, CFGs...>::spec_count += c;
        }

    };

    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    int primary_producer<DOM_T, PSPACE_T, CFGs...>::spec_count = 0;


} // namespace MuLAN_MA
} // namespace Models
} // namespace Utopia

#endif //UTOPIA_MODELS_PRODUCER_HH