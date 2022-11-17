#ifndef UTOPIA_MODELS_INTEGRATORS_RKCK
#define UTOPIA_MODELS_INTEGRATORS_RKCK

#include <cmath>

namespace Utopia::Models::MuLAN_MA {


    /**
     * @brief Integration Scheme for Euler
     * @details Holds everything needed for integration. Perform Time step by performing step() 'steps' times
     *
     * @tparam CURRENCY The currency (e.g. double)
     * @tparam ORG_T
     */
    template<typename CURRENCY, typename ORG_T>
    class euler {

    public:
        using currency = CURRENCY;

    private:
        /**
         * @brief Backreference to the organism
         *
         */
        ORG_T& org;

        /**
         * @brief Hold values
         * @details xs[0] holds value/biomass of the organism. Others temp values for substeps
         *
         */
        currency x;

        /**
         * @brief Store here every sum in the differential equations
         * @details If the k values depend on sums (coupled diff. eq.) perform summing between substeps
         *
         */
        std::vector<currency> sum;

    public:

        /**
         * @brief Global number of steps for the Integrator
         *
         */
        inline static const int steps = 1;

        /**
         * @brief Store change of the previous step
         *
         */
        currency last_dxdt = 0.0;


        /**
         * @brief Store previous step-size
         *
         */
        currency last_dt = 0.0;

        /**
         * @brief Construct a new euler object
         *
         * @param org
         */
        explicit euler(ORG_T& org) : org(org){

        };

        /**
         * @brief Get the value
         *
         * @return currency&
         */
        currency& get_value() {
            return this->x;
        }

        /**
         * @brief Get the value
         *
         * @return const currency&
         */
        const currency& get_value() const {
            return this->x;
        }

        /**
         * @brief Get Error
         *
         * @return currency&
         */
        currency get_error() const {
            return 0.5 * this->last_dxdt * this->last_dt * this->last_dt;
        }

        /**
         * @brief Set the value
         *
         * @param val
         */
        void set_value(const currency& val) {
            this->x = val;
        }

        /**
         * @brief Access sum
         *
         * @tparam T
         * @param i
         * @return currency&
         */
        template <typename T>
        currency& operator[] (const T& i) {
            static_assert(std::is_integral_v<T>, "i has to be integral");
            return this->sum[i];

        }

        template <typename T>
        const currency& operator[] (const T& i) const {
            static_assert(std::is_integral_v<T>, "i has to be integral");
            return this->sum[i];
        }

        /**
         * @brief Clear the sum
         *
         */
        void clear_sum(){
            std::fill(this->sum.begin(), this->sum.end(), 0.0);
        }

        /**
         * @brief Step
         * @details Iterate over steps by calling step() 6 times. Between calls sums should be calculated
         *
         * @param dt
         */
        void step(double dt){
            this->last_dt = dt;
            this->last_dxdt = org.dxdt(this->x, 0.0);
            this->x = this->x + this->last_dxdt * dt;
        }

        /**
         * @brief resize sum vector
         *
         * @param size
         */
        void resize(const int& size){
            this->sum.resize(size);
        }

        /**
         * @brief Calculate new step size individual level
         * @details Use this to calculate new step size based on individual error
         *
         * @param dt
         * @param new_dt
         */
        virtual void calc_new_step_size_s(double /*dt*/, double& /*new_dt*/) {

        }

        /**
         * @brief Calculate new step size on global level
         * @details Use this to calculate new step size based on global error. Calculate e.g. sums over errors in calc_new_step_size_s()
         *
         * @param dt
         * @param new_dt
         */
        static void calc_new_step_size(double /*dt*/, double& /*new_dt*/) {

        }

    };

    /**
     * @brief Integration Scheme for Runge-Kutta Cash-Karp
     * @details Holds everything needed for integration. Perform Time step by performing step() 'steps' times
     *
     * @tparam CURRENCY The currency (e.g. double)
     * @tparam ORG_T
     */
    template<typename CURRENCY, typename ORG_T>
    class rkck {

    public:
        using currency = CURRENCY;

    private:
        /**
         * @brief Backreference to the organism
         *
         */
        ORG_T& org;

        /**
         * @brief Error
         * @details Estimated by scheme
         *
         */
        currency error;

        /**
         * @brief Count the step number
         * @details RKCK performs 6 substeps
         *
         */
        int stepnum = 0;


        /**
         * @brief Safety Factor
         * @details Limits growth of timestep
         *
         */
        double beta = 0.95;

        /**
         * @brief Integration Scheme
         * @details Perform Substeps according to RKCK definition
         *
         * @tparam N
         * @param dt
         */
        template<int N>
        void substep(double dt){
            if constexpr (N == 0) {
                this->ks[0] = org.dxdt(this->xs[0], 0.0);
                this->xs[1] = this->xs[0] + this->ks[0] * dt / 5.0;
            } else if constexpr (N == 1) {
                this->ks[1] = org.dxdt(this->xs[1], dt / 5.0);
                this->xs[2] = this->xs[0] + this->ks[0] * (3.0/ 40.0) + this->ks[1] * dt * (9.0/ 40.0);
            } else if constexpr (N == 2) {
                this->ks[2] = org.dxdt(this->xs[2], (3.0 / 10.0) * dt);
                this->xs[3] = this->xs[0] + this->ks[0] * (3.0 / 10.0) * dt - this->ks[1] * (9.0/10.0) * dt  + this->ks[2] * (6.0 / 5.0) * dt;
            } else if constexpr (N == 3) {
                this->ks[3] = org.dxdt(this->xs[3], (3.0 / 5.0) * dt);
                this->xs[4] = this->xs[0] - this->ks[0] * (11.0 / 54.0) * dt + this->ks[1] * (5.0/2.0) * dt  - this->ks[2] * (70.0 / 27.0) * dt + this->ks[3] * (35.0 / 27.0) * dt ;
            } else if constexpr (N == 4) {
                this->ks[4] = org.dxdt(this->xs[4], dt);
                this->xs[5] = this->xs[0] + this->ks[0] * (1631.0 / 55296.0) * dt + this->ks[1] * (175.0/512.0) * dt  + this->ks[2] * (575.0 / 13824.0) * dt + this->ks[3] * (44275.0 / 110592.0) * dt + this->ks[4] * (253.0 / 4096.0) * dt ;
            } else if constexpr (N == 5) {
                this->ks[5] = org.dxdt(this->xs[5], (7.0 / 8.0) * dt);
                this->last_dxdt = (this->ks[0] * (37.0/378.0) + this->ks[2] * (250.0/621.0) + this->ks[3] * (125.0/594.0) + this->ks[5] * (512.0/1771.0));
                this->xs[0] = this->xs[0] + this->last_dxdt * dt;
                this->error = this->xs[0] - (this->xs[0] + ( (2825.0 / 27648.0) * this->ks[0] + (18575.0 / 48384.0) * this->ks[2] + (13525.0 / 55296.0) * this->ks[3] + (277.0 / 14336.0) * this->ks[4] + (1.0 / 4.0) * this->ks[5]) * dt);   
            }
        };

        /**
         * @brief Hold values
         * @details xs[0] holds value/biomass of the organism. Others temp values for substeps
         *
         */
        std::array<currency,6> xs;

        /**
         * @brief Hold k values
         * @details Needed for substeps
         *
         */
        std::array<currency,6> ks;

        /**
         * @brief Store here every sum in the differential equations
         * @details If the k values depend on sums (coupled diff. eq.) perform summing between substeps
         *
         */
        std::vector<currency> sum;

    public:

        /**
         * @brief Error Threshold
         * @details Estimate time step with error and w_error
         *
         */
        inline static currency w_error = 0.1;

        /**
         * @brief Global number of steps for the Integrator
         *
         */
        inline static const int steps = 6;
        [[maybe_unused]] inline static const int error_step = 5;

        /**
         * @brief Store change of the previous step
         *
         */
        currency last_dxdt = 0.0;

        /**
         * @brief Construct a new rkck object
         *
         * @param org
         */
        explicit rkck(ORG_T& org) : org(org){

        };

        /**
         * @brief Get the value
         *
         * @return currency&
         */
        currency& get_value() {
            return this->xs[this->stepnum];
        }

        /**
         * @brief Get the value
         *
         * @return const currency&
         */
        const currency& get_value() const {
            return this->xs[this->stepnum];
        }

        /**
         * @brief Get Error
         *
         * @return const currency&
         */
        const currency& get_error() const {
            return this->error;
        }

        /**
         * @brief Set the value
         *
         * @param val
         */
        void set_value(const currency& val) {
            this->xs[this->stepnum] = val;
        }

        /**
         * @brief Set the safety value
         *
         * @param val
         */
        void set_beta(const double& val) {
            this->beta = val;
        }

        /**
         * @brief Access sum
         *
         * @tparam T
         * @param i
         * @return currency&
         */
        template <typename T>
        currency& operator[] (const T& i) {
            static_assert(std::is_integral_v<T>, "i has to be integral");
            return this->sum[i];

        }

        template <typename T>
        const currency& operator[] (const T& i) const {
            static_assert(std::is_integral_v<T>, "i has to be integral");
            return this->sum[i];
        }

        /**
         * @brief Clear the sum
         *
         */
        void clear_sum(){
            std::fill(this->sum.begin(), this->sum.end(), 0.0);
        }

        /**
         * @brief Step
         * @details Iterate over steps by calling step() 6 times. Between calls sums should be calculated
         *
         * @param dt
         */
        void step(double dt){
            if (this->stepnum == 0) {
                this->substep<0>(dt);
            } else if (this->stepnum == 1) {
                this->substep<1>(dt);
            }else if (this->stepnum == 2) {
                this->substep<2>(dt);
            }else if (this->stepnum == 3) {
                this->substep<3>(dt);
            }else if (this->stepnum == 4) {
                this->substep<4>(dt);
            }else if (this->stepnum == 5) {
                this->substep<5>(dt);
                this->stepnum = -1;
            }
            this->stepnum++;
        }

        /**
         * @brief resize sum vector
         *
         * @param size
         */
        void resize(const int& size){
            this->sum.resize(size);
        }

        /**
         * @brief Calculate new step size individual level
         * @details Use this to calculate new step size based on individual error
         *
         * @param dt
         * @param new_dt
         */
        virtual void calc_new_step_size_s(double dt, double& new_dt) {
            double err_q = fabs(this->w_error / this->error);
            if ( err_q < 1.0 ) {
                new_dt = std::min(new_dt, this->beta * dt * pow(err_q, 0.2));
            } else {
                new_dt = std::min(new_dt, this->beta * dt * pow(err_q, 0.25));
            }

        }

        /**
         * @brief Calculate new step size on global level
         * @details Use this to calculate new step size based on global error. Calculate e.g. sums over errors in calc_new_step_size_s()
         *
         * @param dt
         * @param new_dt
         */
        static void calc_new_step_size(double /*dt*/, double& /*new_dt*/) {

        }

    };

} // namespace Utopia::Models::MuLAN_MA

#endif //UTOPIA_MODELS_INTEGRATORS_RKCK
