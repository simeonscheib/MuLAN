#ifndef UTOPIA_MODELS_SPECIES_CLASS_2_HH
#define UTOPIA_MODELS_SPECIES_CLASS_2_HH

#include "../organism.hh"

// Include all other species classes
#include "species_class_1.hh"

namespace Utopia {
namespace Models {
namespace MuLAN_MA {

    // TODO: introduce all the other species
    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    class species_class_1;

    // Inherit from MuLAN organism
    // TODO: change name
    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    class species_class_2: public organism<DOM_T, PSPACE_T>{

    public:

        using pp_ptr [[maybe_unused]] = boost::shared_ptr<species_class_2>;
        using org_ptr = typename organism<DOM_T, PSPACE_T>::org_ptr;
        using currency = typename DOM_T::integrator_t::currency;

        // TODO: get compile time parameters from pack
        static constexpr int option_1 = get_from_pack<int, 0, CFGs...>();


        double dxdt(const currency& x, const double tpdt) override {

            // TODO: Add differential equation here

        }

        // TODO: give it a layer id
        const static int type_id = 1;

        // Store number of species
        static int spec_count;

        /**
         * @brief Construct a producer with access to  domain d
         *
         * @param d Domain
         */
        explicit species_class_2(DOM_T* d) : organism<DOM_T, PSPACE_T>(d) {
            // TODO: construcion
        }; // TODO: rename

        /**
         * @brief Construct producer with access to  domain d and initial mass
         *
         * @param d             Domain
         * @param initial_mass  initial biomass
         */
        explicit species_class_2(DOM_T* d, double initial_mass) : species_class_2(d) {

            this->set_value(initial_mass);

        }; // TODO: change name

        /**
         * @brief destruct producer
         */
        ~species_class_2() = default; // TODO: change name


        /**
         * @brief Calculate the interaction edge of this and org2
         *
         * @param org2  Organism
         * @return
         */
        typename DOM_T::edge_cont calc_interaction_coeff(const org_ptr org2) const override {

            // TODO: calculate interaction edge
            return typename DOM_T::edge_cont{0.0};
        }


        void add_edge_cont1(const typename DOM_T::edge_cont& val, org_ptr organism_2) override {

            // TODO: calculate most inner sums

        }

        void add_edge_cont2(const typename DOM_T::edge_cont& /*val*/, const org_ptr /*organism_2*/) override {

            // TODO: calculate second most inner sums
        }

        // void set_parameters();

        virtual void count(int c) override{

            species_class_2<DOM_T, PSPACE_T, CFGs...>::spec_count += c; // TODO: change name
        }

    };

    // TODO: change name
    template <typename DOM_T, typename PSPACE_T, int ...CFGs>
    int species_class_2<DOM_T, PSPACE_T, CFGs...>::spec_count = 0;


} // namespace MuLAN_MA
} // namespace Models
} // namespace Utopia

#endif //UTOPIA_MODELS_SPECIES_CLASS_2_HH