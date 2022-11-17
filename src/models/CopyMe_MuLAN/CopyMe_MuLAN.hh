#ifndef UTOPIAMODELSMULAN_COPYME_MULAN_HH
#define UTOPIAMODELSMULAN_COPYME_MULAN_HH

// Utopia includes
#include <utopia/core/model.hh>
#include <utopia/core/types.hh>

#include <functional>

// MuLAN includes
#include "utils/MuLAN_MA_utils.hh"
#include "utils/Organism_Factory.hh"
#include "utils/Organism_Manager.hh"

#include "domain.hh"

#include "orga/species_class_1.hh"
#include "orga/species_class_2.hh"

namespace Utopia::Models::MuLAN_MA {

    /**
     * @brief Wrapper for adding nodes to the graph
     * @details Define a function of teh domain to use
     *
     * @tparam DOM_T
     */
    template <typename DOM_T>
    struct DomainAddFunctor {
        /*
         * @brief The type of function that takes a domain and necessary parameters to add an Organism
         */
        using type = std::function<typename DOM_T::vertex_t& (DOM_T&, typename DOM_T::currency, typename DOM_T::pspace_t)>;

        /*
         * @brief   Return a function that calls the correct adding function for an Organism ORG_T into a domain
         *
         * @return Lambda   Function of type "type"
         */
        template <typename ORG_T>
        auto operator() (){
            // Change me if you use a custom domain
            return [](DOM_T& dom, typename DOM_T::currency a, typename DOM_T::pspace_t b) -> typename DOM_T::vertex_t& {
                return dom.template add_<ORG_T>(a, b);
            };
        }
    };

    /// Typehelper to define types of MuLAN model
    using MuLAN_MAModelTypes = ModelTypes<>;

    /**
     * @brief Multilayered Adaptive Network Model
     * @details The Model Class
     *
     * @tparam DOM_T
     */
    template <typename DOM_T>
    class CopyMe_MuLANModel:
            public Model<CopyMe_MuLANModel<DOM_T>, MuLAN_MAModelTypes>
    {
    public:
        /// The base model type
        using Base = Model<CopyMe_MuLANModel, MuLAN_MAModelTypes>;

        /// Data type that holds the configuration
        using Config [[maybe_unused]] = typename Base::Config;

        /// Data type of the group to write model data to, holding datasets
        using DataGroup [[maybe_unused]] = typename Base::DataGroup;

        /// Data type for a dataset
        using DataSet [[maybe_unused]] = typename Base::DataSet;

        /// Data type of the shared RNG
        using RNG [[maybe_unused]] = typename Base::RNG;


        // TODO: Specializations your species classes here. Add more if necessary
        template<int ...CFGs>
        using species_class_1_t = species_class_1<typename DOM_T::base_t, typename DOM_T::pspace_t, CFGs...>;

        template<int ...CFGs>
        using species_class_2_t = species_class_2<typename DOM_T::base_t, typename DOM_T::pspace_t, CFGs...>;


        /**
         * @brief Organism Manager Type
         * @details Register species templates that should be available
         *
         */
        using organism_manager_t = OrganismManager<DOM_T, DomainAddFunctor, species_class_1_t, species_class_2_t>;
        // TODO:  Add all species classes to the  manger



    private:
        // Base members: _time, _name, _cfg, _hdfgrp, _rng, _monitor

        // -- Members of this model -- //

        /**
         * @brief Domain where the Network forms
         *
         */
        DOM_T _dom;

        /**
         * @brief Organism Manager adds Species to the domain
         *
         */
        organism_manager_t _org_mngr;

        /**
         * @brief Constant number of registered species templates
         *
         */
        static constexpr std::size_t number_species = organism_manager_t::get_number_of_species();

    public:
        /**
         * @brief Construct a new MuLAN_MAModel object
         *
         * @tparam ParentModel
         * @param name Name of this model instance
         * @param parent The parent model this model instance resides in
         */
        template<class ParentModel>
        CopyMe_MuLANModel (const std::string name,
                       ParentModel &parent)
                :
                // Initialize first via base model
                Base(name, parent),

                // Initialize OrganismManager with domain
                _org_mngr(_dom)

        {

            // TODO: rename
            std::vector<int> species_class_1_conf_v;

            // TODO: Fill in configuration mappings for compile time specialization
            std::map<std::string, int> interaction_map_species_class_1 = {
                    {"A", 0},
                    {"B", 1}
            };

            species_class_1_conf_v.push_back(interaction_map_species_class_1[get_as<std::string>("species_class_1_p1", this->_cfg)]);

            // Register species_class_1_p1
            // TODO: change number of options
            this->_org_mngr.template register_<species_class_1_t, 2>(species_class_1_conf_v);

            std::vector<int> species_class_2_conf_v;

            // Response Function predator prey
            std::map<std::string, int> species_class_2_map = {
                    {"A", 0},
                    {"B", 1},
                    {"C", 2}
            };

            species_class_2_conf_v.push_back(species_class_2_map[get_as<std::string>("species_class_2_p1", this->_cfg)]);


            // Initializes the first species
            this->initialize_species();

            // save number of species classes
            this->_hdfgrp->add_attribute("num_species_classes", this->_org_mngr.get_number_of_species());


        }

        /**
         * @brief Initialize the system according to config
         *
         */
        void initialize_species(){

            initialize_species_class_1();
            initialize_species_class_2();

        }

        void initialize_species_class_1() {

            // TODO: Initialize species_class_1 here
            // TODO: Change the names and arguments according to previous configuration
            // Use this to add a species to the domain
            this->_org_mngr.template add_<species_class_2_t>(mass, ps);

            this->_log->debug("species_class_1 initialized");

        }

        void initialize_species_class_2() {

            // TODO: Initialize species_class_2 here
            // TODO: Change the names and arguments according to previous configuration
            // Use this to add a species to the domain
            this->_org_mngr.template add_<species_class_2_t>(mass, ps);

            this->_log->debug("species_class_2 initialized");

        }


        // Runtime functions ......................................................

        /* perform one timestep
         * One time step in this context is dt2
         *
        */
        void perform_step [[maybe_unused]]()
        {

            // TODO: Add code to perform a step. You can use the MuLAN domain integrator

            // mutate
            this->mutation();

        }

        void mutation () {

            // TODO: Add code for mutation subroutine

        }


        /// Monitor model information

        void monitor [[maybe_unused]] ()
        {
            this->_monitor.set_entry("Number producers: ", *(this->pp_spec_count));
            this->_monitor.set_entry("Number consumers: ", *(this->cons_spec_count));
        }


        /// Write data
        void write_data [[maybe_unused]] ()
        {

            // TODO: Save graph properties e.g.
            // TODO: you have to change the return value of the lambdas accordingly
            save_graph_vertex_value(this->_dom.graph, this->_hdfgrp, std::to_string(this->_time), "_types",
                                    [](auto& vw){return static_cast<char>(vw.get_type());});

            if ( get_as<bool>("safe_graph",  this->_cfg) ) {
                this->_log->debug("Saving Graph");
                save_graph_vertex_value(this->_dom.graph, this->_hdfgrp, std::to_string(this->_time), "_graph_vertex",
                                        [](auto& vw){return vw.id;});
                save_graph_edges(this->_dom.graph, this->_hdfgrp, std::to_string(this->_time), "_graph_edges");
                save_graph_edge_value(this->_dom.graph, this->_hdfgrp, std::to_string(this->_time), "_graph_edge_v",
                                      [](auto& vw){return vw.interaction_vector[0];});
            }

            /*
            save_container(Add container here, this->_hdfgrp, std::to_string(this->_time), "name_here",
                               [](auto &it) { return it; });
            */

        }


        // Getters and setters ....................................................
        // Add getters and setters here to interface with other model

    };

}

#endif //UTOPIAMODELSMULAN_COPYME_MULAN_HH
