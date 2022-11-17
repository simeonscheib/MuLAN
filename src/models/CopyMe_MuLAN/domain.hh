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


    /**
     * @brief Domain Interface
     * @details Interface for domain class after filtering configs
     * 
     */
    template <typename CURRENCY, int INTEGRATOR, std::size_t SUM_SIZE>
    class domain_interface : public domain_base<domain_interface<CURRENCY, INTEGRATOR, SUM_SIZE>,
                                        CURRENCY, 
                                        parameter_space, 
                                        vertex_wrapper<
                                            organism<
                                                    domain_interface<
                                                            CURRENCY,
                                                            INTEGRATOR,
                                                            SUM_SIZE,
                                                            ENV_FUNC>,
                                                    parameter_space>
                                        >, 
                                        edge_wrapper<
                                            std::array<double, 1> 
                                        >,
                                        SUM_SIZE
                                    > {
    public:

        /**
         * @brief Define the chosen Integrator
         * 
         */
        using integrator_t = std::conditional_t< INTEGRATOR == 0, 
            rkck<CURRENCY, organism<domain_interface<CURRENCY, INTEGRATOR, SUM_SIZE>, parameter_space> >,
            euler<CURRENCY, organism<domain_interface<CURRENCY, INTEGRATOR, SUM_SIZE, ENV_FUNC>, parameter_space> > >;

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

        virtual double S(double, double) = 0;
    };



    /**
     * @brief Domain Class
     * @details Specify organism, vertex and edge type
     * 
     */
    template <  typename CURRENCY, 
                int INTEGRATOR,
                int SUM_SIZE,
                int ...Ns>
    class domain : public domain_interface<CURRENCY, INTEGRATOR, SUM_SIZE> {

    public:

        /**
         * @brief Currency of the domain
         * @details e.g. Biomass
         * 
         */
        using currency [[maybe_unused]] = CURRENCY;

        /**
         * @brief Base Type of domain
         * @details e.g. domain_interface
         * 
         */
        using base_t [[maybe_unused]] = domain_interface<CURRENCY, INTEGRATOR, SUM_SIZE>;

    };


} // namespace Utopia::Models::MuLAN_MA


#endif //UTOPIA_MY_DOMAIN_HH
