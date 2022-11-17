#include "consumer.hh"

#ifndef UTOPIA_CONSUMER_IMPL_HH
#define UTOPIA_CONSUMER_IMPL_HH


namespace Utopia::Models::MuLAN_MA {


    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    consumer<DOM_T, PSPACE_T, CFGs...>::consumer(DOM_T* d) : consumer::organism(d), buffer(d->bsize + 1){

        this->integrator.resize(sum_size);

        for(int i = 0; i < d->bsize; i++) {

            buffer.push_back(0.0);

        }

    }

    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    consumer<DOM_T, PSPACE_T, CFGs...>::consumer(DOM_T *d, double initial_mass) : consumer(d) {

        this->integrator.set_value(initial_mass);

    }

    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    typename DOM_T::edge_cont consumer<DOM_T, PSPACE_T, CFGs...>::calc_interaction_coeff(const org_ptr org2) const {

        if (org2->parameters.type == 0) {
            //typename primary_producer<DOM_T, PSPACE_T>::pp_ptr tmp = boost::dynamic_pointer_cast<primary_producer<DOM_T, PSPACE_T> >(org2);

            double c = this->parameters.params[3];
            double y =  this->parameters.trait[1];
            double k = (org2->parameters.trait[0] - this->parameters.trait[0])/y;
            double tmp = exp(-c * y)/(sqrt(2 * M_PI) * y) * exp(-0.5 * k * k);

            return typename DOM_T::edge_cont{tmp};
            //return typename DOM_T::edge_cont{this->calc_z_coeff(*tmp)};

        } else if (org2->parameters.type == 1) {


            // cons_ptr tmp = boost::dynamic_pointer_cast<consumer>(org2);

            return typename DOM_T::edge_cont{0.0};

        } else {

            assert(org2->parameters.type == 1);

        }

        return typename DOM_T::edge_cont{0.0};
    }

    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    double consumer<DOM_T, PSPACE_T, CFGs...>::get_influx(){
        double b = this->parameters.params[1];
        return b * this->integrator[0];
    }

    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    double consumer<DOM_T, PSPACE_T, CFGs...>::dxdt(const currency& x, [[maybe_unused]] const double tpdt) {

        double R = this->parameters.params[0];
        double b = this->parameters.params[1];
        double m = this->parameters.params[2];

        double ret = x * R * (b * this->integrator[0] - m);



        this->integrator.clear_sum();

        if (ret >= 0.0){
            this->buffer.push_back(ret);

            ret = this->buffer[0];

            this->buffer.pop_front();

            return ret;
        } else {

            this->buffer.push_back(0.0);

            ret += this->buffer[0];

            this->buffer.pop_front();

            return ret;

        }
    }

    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    void consumer<DOM_T, PSPACE_T, CFGs...>::add_edge_cont1(const typename DOM_T::edge_cont& val, org_ptr organism_2) {

        if constexpr ( response_func == 1 ) {
            const currency &mass = organism_2->get_value();
            const currency &x = this->get_value();

            this->integrator[0] += mass * val[0];
            organism_2->integrator[0] += x * val[0];
        } else if constexpr ( response_func == 3 ) {
            const currency &mass = organism_2->get_value();

            this->integrator[1] += mass * val[0];
        } else if constexpr ( response_func == 4 ) {
            const currency &mass = organism_2->get_value();

            // this->integrator[1] += this->parameters.params[4] * mass * val[0] *  mass * val[0]
            //                         + this->parameters.params[5] * mass * val[0];

            this->integrator[1] += mass * val[0];
        } else {
            throw std::invalid_argument("No valid 'response_step' function. Check config.");
        }
    }
    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    void consumer<DOM_T, PSPACE_T, CFGs...>::add_edge_cont2(const typename DOM_T::edge_cont& val, const org_ptr organism_2) {

        if constexpr ( response_func == 1 ) {

        } else if constexpr ( response_func == 3 ) {
            const currency &mass = organism_2->get_value();
            const currency &x = this->get_value();

            auto tmp = val[0]/(1 + this->integrator[1] * this->parameters.params[4]);

            this->integrator[0] += mass * tmp;
            organism_2->integrator[0] += x * tmp;
        } else if constexpr ( response_func == 4 ) {
            const currency &mass = organism_2->get_value();
            const currency &x = this->get_value();

            auto tmp = val[0]/(1 + this->integrator[1]);

            this->integrator[0] += mass * mass * tmp;
            organism_2->integrator[0] += x * mass * tmp;
        } else {
            throw std::invalid_argument("No valid 'response_step' function. Check config.");
        }
    }


} // namespace Utopia::Models::MuLAN_MA

#endif //UTOPIA_CONSUMER_IMPL_HH
