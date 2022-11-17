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
    template <typename CURRENCY, int INTEGRATOR, std::size_t SUM_SIZE, int ENV_FUNC>
    class domain_interface : public domain_base<domain_interface<CURRENCY, INTEGRATOR, SUM_SIZE, ENV_FUNC>,
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
            rkck<CURRENCY, organism<domain_interface<CURRENCY, INTEGRATOR, SUM_SIZE, ENV_FUNC>, parameter_space> >,
            euler<CURRENCY, organism<domain_interface<CURRENCY, INTEGRATOR, SUM_SIZE, ENV_FUNC>, parameter_space> > >;

        using currency = typename integrator_t::currency;

        /**
         * @brief The type of env function to be used encoded as int
         * @details Needed by MuLAN_MA Model class
         *
         */
        static constexpr int  env_func = ENV_FUNC;

        /**
         * @brief Parameter vector of the domain
         * @details Holds Domain Parameters
         *
         */
        std::vector<double> params;

        /**
         * @brief Construct domain_interface
         *
         */
        domain_interface() {std::cout << SUM_SIZE << std::endl;};//= default;
        ~domain_interface() = default;

        // Buffer size
        int bsize = 0;

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
                int STEP_FUNC,
                int SUM_SIZE,
                int ENV_FUNC,
                int ...Ns>
    class domain : public domain_interface<CURRENCY, INTEGRATOR, SUM_SIZE, ENV_FUNC> {

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
        using base_t [[maybe_unused]] = domain_interface<CURRENCY, INTEGRATOR, SUM_SIZE, ENV_FUNC>;

        /**
         * @brief The type of step function to be used encoded as int
         * @details Needed by MuLAN_MA Model class
         * 
         */
        static constexpr int  step_func = STEP_FUNC;

        /**
         * @brief Carrying Capacity function
         * 
         * @param z Niche Position
         * @param dt_loc timestep to add to current time
         * @return double 
         */
        virtual double S(double z, double dt_loc){

            if constexpr (ENV_FUNC == 0) {
                double tmp = z / this->params[1];
                return this->params[0] * exp(-1.0 * 0.5 * tmp * tmp);
            } else if constexpr (ENV_FUNC == 1){
                double tmp = z / this->params[1];
                double tmp2 = this->params[0] * exp(-1.0 * 0.5 * tmp * tmp)
                        + this->params[2] * cos(z/this->params[3]);
                return std::max(tmp2, 0.0);
            } else if constexpr (ENV_FUNC == 2){
                double tmp = z / this->params[1];
                double tmp2 = this->params[0] * exp(-1.0 * 0.5 * tmp * tmp)
                       + this->params[2] * cos(z/this->params[3])
                       * sin(this->params[4] * (this->time + dt_loc));
                return std::max(tmp2, 0.0);
            } else if constexpr (ENV_FUNC == 3){
                double tmp = z / this->params[1];
                double tmp2 =  this->params[0] * exp(-1.0 * 0.5 * tmp * tmp)
                         + this->params[2] * cos(z/this->params[3]) * cos(z/this->params[3])
                         * sin(this->params[4] * (this->time + dt_loc))
                         * sin(this->params[4] * (this->time + dt_loc));
                return std::max(tmp2, 0.0);
            } else {
                throw std::invalid_argument("No valid 'env_func' function. Check config.");
            }
        }

    };


} // namespace Utopia::Models::MuLAN_MA


#endif //UTOPIA_MY_DOMAIN_HH
