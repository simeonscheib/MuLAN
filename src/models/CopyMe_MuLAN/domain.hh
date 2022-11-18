//
// Created by Simeon Scheib on 2019-01-31.
//
#include "./MuLAN_base/_domain.hh"
#include "parameter_space.hh"
#include "vertex_wrapper.hh"
#include "organism.hh"
#include "integrators/rkck.hh"

#ifndef UTOPIA_MY_DOMAIN_HH
#define UTOPIA_MY_DOMAIN_HH

namespace Utopia::Models::MuLAN_MA {

    // TODO: Modify template. Define the class put on the nodes and edges as well as parameter/trait
    /**
     * @brief Domain Interface
     * @details Interface for domain class after filtering configs
     * 
     */
    template <typename CURRENCY, std::size_t SUM_SIZE, int INTEGRATOR, int OPTION_1>
    class domain_interface : public domain_base<domain_interface<CURRENCY, SUM_SIZE, INTEGRATOR, OPTION_1>,
                                        // The currency of the system. e.g. double
                                        CURRENCY,
                                        // Custom parameter space per orgaqnism
                                        parameter_space,
                                        // Vertex Type
                                        vertex_wrapper<
                                            organism<
                                                    domain_interface<
                                                            CURRENCY,
                                                            SUM_SIZE,
                                                            INTEGRATOR,
                                                            OPTION_1>,
                                                    parameter_space>
                                        >,
                                        // edge type
                                        edge_wrapper<
                                            // In this case an array with a single weight
                                            std::array<double, 1> 
                                        >,
                                        // How often to iterate over edges. Use to calculate nested sums
                                        SUM_SIZE
                                    > {
    public:

        /**
         * @brief Define the chosen Integrator
         * 
         */
        using integrator_t = std::conditional_t< INTEGRATOR == 0,
            rkck<CURRENCY, organism<domain_interface<CURRENCY, SUM_SIZE, INTEGRATOR, OPTION_1>, parameter_space> >,
            euler<CURRENCY, organism<domain_interface<CURRENCY, SUM_SIZE, INTEGRATOR, OPTION_1>, parameter_space> > >;

        using currency = typename integrator_t::currency;


        /**
         * @brief Construct domain_interface
         *
         */
        domain_interface() {std::cout << SUM_SIZE << std::endl;};//= default;
        ~domain_interface() = default;


        void set_error(const currency& error){

            if constexpr(INTEGRATOR == 0){
                integrator_t::w_error = error;
            } else {
                std::cout << "Integrator has no error estimate!" << std::endl;
            }
        }

       // TODO: Add members here. You can use virtual functions
    };



    /**
     * @brief Domain Class
     * @details Specify organism, vertex and edge type
     * 
     */
    template <  typename CURRENCY,
                int SUM_SIZE,
                int INTEGRATOR,
                int OPTION_1,
                int OPTION_2,
                int ...Ns>
    class domain : public domain_interface<CURRENCY, SUM_SIZE, INTEGRATOR, OPTION_1> {

    public:

        /**
         * @brief Currency of the domain
         * @details e.g. Biomass
         * 
         */
        using currency = CURRENCY;

        /**
         * @brief Base Type of domain
         * @details e.g. domain_interface. This is very important and should always be set. If you
         * do not use
         * 
         */
        using base_t [[maybe_unused]] = domain_interface<CURRENCY, SUM_SIZE, INTEGRATOR, OPTION_1>;

    };


} // namespace Utopia::Models::MuLAN_MA


#endif //UTOPIA_MY_DOMAIN_HH
