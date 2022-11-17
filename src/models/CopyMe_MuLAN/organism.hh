#include "MuLAN_base/_organism.hh"
#include "integrators/rkck.hh"

#ifndef UTOPIA_MODELS_MY_ORGANISM
#define UTOPIA_MODELS_MY_ORGANISM

namespace Utopia::Models::MuLAN_MA {

    /**
     * @brief Organism
     * @details Concretize the organism
     *
     * @tparam DOM_T
     * @tparam PSPACE_T
     */
    template <typename DOM_T, typename PSPACE_T>
    class organism : public organism_base<organism<DOM_T, PSPACE_T>, DOM_T, PSPACE_T, typename DOM_T::integrator_t> {

    public:

        // using typename organism_base<organism<DOM_T, PSPACE_T>, DOM_T, PSPACE_T, double>::currency;
        using currency = typename DOM_T::integrator_t::currency;
        using edge_cont = typename DOM_T::edge_cont;
        using org_ptr = typename organism<DOM_T, PSPACE_T>::org_ptr;

        /**
         * @brief Construct a new organism object
         *
         * @param d
         */
        explicit organism(DOM_T* d) : organism::organism_base(d, *this) {

        };


        // TODO: Add virtual functions that should be available for vertices


        /**
         * @brief Calculate the change rate
         * @details May depend on parameters, current value x and the integrator (e.g. sums, ..)
         *
         * @param x Value
         * @return currency
         */
        virtual currency dxdt(const currency& x, const double tpdt) = 0;


        // TODO: Add more functions if you need more than 2 nested sums
        /**
         * @brief Interface to add Edges to the sum vector of an organism pair
         * @details Call virtual functions e.g. add_edge_contX of organisms
         *
         * @tparam N Depth of sums e.g. if sum is sum of sums
         * @param val Edge value (scalar or vector)
         * @param organism_2 The second organism
         */
        template<int N>
        void add_edge_cont(const edge_cont& val, org_ptr organism_2) {
            if constexpr ( N == 0 ) {
                this->add_edge_cont1(val, organism_2);
            } if constexpr ( N == 1 ) {
                this->add_edge_cont2(val, organism_2);
            } else {
                return;
            }
        };

        virtual void add_edge_cont1(const edge_cont& val, org_ptr organism_2) = 0;
        virtual void add_edge_cont2(const edge_cont& val, org_ptr organism_2) = 0;

    };


} // namespace Utopia::Models::MuLAN_MA

#endif //UTOPIA_MODELS_MY_ORGANISM
