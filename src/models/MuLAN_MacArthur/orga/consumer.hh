#ifndef UTOPIA_MODELS_CONSUMER_HH
#define UTOPIA_MODELS_CONSUMER_HH

#include <boost/circular_buffer.hpp>
#include "../organism.hh"
#include "primary_producer.hh"
#include "utils/helper.hh"

namespace Utopia::Models::MuLAN_MA {

    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    class primary_producer;

    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    class consumer: public organism<DOM_T, PSPACE_T>{
    protected:

        //double G();

        boost::circular_buffer<double> buffer;

    public:


        using org_ptr = typename organism<DOM_T, PSPACE_T>::org_ptr;
        using currency = typename DOM_T::integrator_t::currency;

        static constexpr int response_func = get_from_pack<int, 0, CFGs...>();

        inline static const int sum_size = 2;

        double dxdt(const currency& x, [[maybe_unused]] const double tpdt) override;

        const static int type_id = 1;
        static int spec_count;

        explicit consumer(DOM_T* d);
        explicit consumer(DOM_T* d, double initial_mass);

        ~consumer() = default;

        typename DOM_T::edge_cont calc_interaction_coeff(org_ptr org2) const override;


        void add_edge_cont1(const typename DOM_T::edge_cont& val, org_ptr organism_2) override;
        void add_edge_cont2(const typename DOM_T::edge_cont& val, org_ptr organism_2) override;

        double get_niche_width() override {

            return this->parameters.trait[1];

        }

        double get_influx() override;

        void count(int c) override {

            consumer<DOM_T, PSPACE_T, CFGs...>::spec_count += c;
        }

    };

    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    int consumer<DOM_T, PSPACE_T, CFGs...>::spec_count = 0;


} // namespace Utopia::Models::MuLAN_MA

#endif //UTOPIA_MODELS_CONSUMER_HH