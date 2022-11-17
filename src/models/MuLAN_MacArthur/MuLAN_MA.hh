#ifndef UTOPIA_MODELS_MULAN_MA_HH
#define UTOPIA_MODELS_MULAN_MA_HH


// Utopia-related includes
#include <utopia/core/model.hh>
#include <utopia/core/types.hh>

/*#include <dune/utopia/base.hh>
#include <dune/utopia/core/setup.hh>
#include <dune/utopia/core/model.hh>
#include <dune/utopia/core/apply.hh>
#include <dune/utopia/core/types.hh>
*/

#include <functional>

#include "domain.hh"
#include "utils/MuLAN_MA_utils.hh"
#include "utils/Organism_Factory.hh"
#include "utils/Organism_Manager.hh"

#include "orga/primary_producer.hh"
#include "orga/consumer_impl.hh"

//#define TWO_SUMS

namespace Utopia::Models::MuLAN_MA {

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
            return [](DOM_T& dom, typename DOM_T::currency a, typename DOM_T::pspace_t b) -> typename DOM_T::vertex_t& {return dom.template add_<ORG_T>(a, b);};
        }
    };

    /// Typehelper to define types of MuLAN model
    using MuLAN_MAModelTypes = ModelTypes<>;


    /// The MuLAN Model

    /**
     * @brief Multilayered Adaptive Network Model
     * @details The Model Class
     * 
     * @tparam DOM_T 
     */
    template <typename DOM_T>
    class MuLAN_MAModel:
        public Model<MuLAN_MAModel<DOM_T>, MuLAN_MAModelTypes>
    {
    public:
        /// The base model type
        using Base = Model<MuLAN_MAModel, MuLAN_MAModelTypes>;

        /// Data type that holds the configuration
        using Config [[maybe_unused]] = typename Base::Config;

        /// Data type of the group to write model data to, holding datasets
        using DataGroup [[maybe_unused]] = typename Base::DataGroup;

        /// Data type for a dataset
        using DataSet [[maybe_unused]] = typename Base::DataSet;

        /// Data type of the shared RNG
        using RNG [[maybe_unused]] = typename Base::RNG;

        using trait_t = typename DOM_T::pspace_t::trait_t;
        using trait_t_base = typename DOM_T::pspace_t::trait_t::value_type;

        template<int ...CFGs>
        using primary_producer_t = primary_producer<typename DOM_T::base_t, typename DOM_T::pspace_t, CFGs...>;

        template<int ...CFGs>
        using consumer_t = consumer<typename DOM_T::base_t, typename DOM_T::pspace_t, CFGs...>;

        /**
         * @brief Organism Manager Type
         * @details Register species templates that should be available 
         * 
         */
        using organism_manager_t = OrganismManager<DOM_T, DomainAddFunctor, primary_producer_t, consumer_t>;

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

        /**
         * @brief Add Producer by config
         * 
         */
        typename organism_manager_t::template Add<primary_producer_t>::OrganismAdder* Add_Producer;

        /**
         * @brief Reference to static organism counter (producer)
         * 
         */
        int* pp_spec_count = nullptr;

        /**
         * @brief Add Consumer by config
         * 
         */
        typename organism_manager_t::template Add<consumer_t>::OrganismAdder* Add_Consumer;

        /**
         * @brief Reference to static organism counter (consumer)
         * 
         */
        int* cons_spec_count = nullptr;

        /**
         * @brief Time step for Population Dynamics
         * 
         */
        double _dt = 0.01;

        /**
         * @brief Time step for mutation
         * 
         */
        double _dt2 = 0.1;

        /**
         * @brief Time step minimum
         *
         */
        double _min_dt = 0.0;

        /**
         * @brief Which Traits are mutable
         * @details Entry i corresponds to trait entry i
         * 
         */
        std::array<bool, number_species> mutate;


        // TODO: Set by config
        /**
         * @brief Trait change per mutation
         * 
         */
        trait_t change = {3, 1};

        /**
         * @brief Rate of mutation
         * 
         */
        double mutation_rate{};

        // Normal distribution
        //std::normal_distribution<double> dist;
        std::uniform_real_distribution<double> dist;

        std::vector<double> dts;



    public:
        /**
         * @brief Construct a new MuLAN_MAModel object
         * 
         * @tparam ParentModel 
         * @param name Name of this model instance
         * @param parent The parent model this model instance resides in
         */
        template<class ParentModel>
        MuLAN_MAModel (const std::string name,
                     ParentModel &parent)
        :
            // Initialize first via base model
            Base(name, parent),

            // Initialize OrganismManager with domain
            _org_mngr(_dom),

            // initialize time steps
            _dt(get_as<double>("dt", this->_cfg)),
            _dt2(get_as<double>("dt2", this->_cfg)),

            // initialize normal distribution
            dist(0.0, 1.0)
        {
            // initialize domain parameters
            this->_dom.params = get_as<std::vector<double> >("dom_params", this->_cfg);

            this->_hdfgrp->add_attribute("dt2", this->_dt2);

            if (this->_cfg["S"]) {
                this->_log->debug("Setting S_0. Overrides domain parameter 1");
                this->_dom.params[0] = get_as<double>("S", this->_cfg);
            }

            if (this->_cfg["min_dt"]) {
                this->_min_dt = get_as<double>("min_dt", this->_cfg);
            }

            // TODO: fix for adaptive step
            this->_dom.bsize = get_as<int>("delay", this->_cfg);

            // Set interaction tolerance
            this->_dom.set_interaction_tolerance(get_as<double>("interaction_tolerance", this->_cfg));

            // Set Biomass threshold
            this->_dom.set_bm_threshold(get_as<double>("bm_threshold", this->_cfg));

            std::vector<int> producer_conf_v;

            // Interaction between producers
            std::map<std::string, int> interaction_map = {
                    {"none", 0},
                    {"LV", 1}
            };

            producer_conf_v.push_back(interaction_map[get_as<std::string>("pp_interaction", this->_cfg)]);

            // Register the producer type
            this->Add_Producer = this->_org_mngr.template register_<primary_producer_t, 2>(producer_conf_v);
            this->pp_spec_count = &this->Add_Producer->get_spec_count();

            std::vector<int> consumer_conf_v;

            // Response Function predator prey
            std::map<std::string, int> response_func_map = {
                    {"none", 0},
                    {"linear", 1},
                    {"linear_cut", 2},
                    {"type_2", 3},
                    {"type_3", 4}
            };

            consumer_conf_v.push_back(response_func_map[get_as<std::string>("response_func", this->_cfg)]);

            // Register the consumer type
            this->Add_Consumer = this->_org_mngr.template register_<consumer_t, 5>(consumer_conf_v);
            this->cons_spec_count = &this->Add_Consumer->get_spec_count();
            
            // set mutation on/off
            this->mutate = get_as<std::array<bool, number_species>>("mutation", this->_cfg);

            // set mutation rate
            this->mutation_rate = get_as<double>("mutation_rate", this->_cfg);

            // Set error for rkck
            this->_dom.set_error(get_as<double>("error", this->_cfg));

            // Initializes the first species
            this->initialize_species();

            // save number of species classes
            this->_hdfgrp->add_attribute("num_species_classes", this->_org_mngr.get_number_of_species());

            // Write initial state
            // this->write_data();
            for (size_t i = 0; i < this->mutate.size(); i++){
                this->_log->debug("Mutate Trait {1}: {2}", i, this->mutate[i]);
            }

        }

        /**
         * @brief Initialize the system according to config 
         * 
         */
        void initialize_species(){
            /*
            if (this->_cfg["init_graph"]){

                auto path = get_as<std::string>("init_graph", this->_cfg);
                std::ifstream file(path);

                CSVRow row;

                while(file >> row)
                {
                    double type = to_<double>(row[0]);
                    double trait = to_<double>(row[1]);
                    double width = to_<double>(row[2]);
                    double mass = to_<double>(row[3]);
                    std::cout << mass << std::endl;
                }
            } else {

            }
            */
                initialize_pp();
                initialize_consumer();

        }

        void initialize_pp() {
            double pp_trait = 0;
            double pp_init_mass = get_as<double>("pp_init_mass",  this->_cfg);
            std::vector<double> pp_params = {get_as<double>("r", this->_cfg), this->_dom.S(pp_trait, 0.0)};
            // this->_dom.template add_<primary_producer>(pp_init_mass, {0, {pp_trait, 0}, pp_params});

            for (auto i : this->_cfg["init_pp"]) {

                // Gaussian initialization
                if ( get_as<std::string>("type", i) == "gaussian" ) {
                    typename DOM_T::pspace_t ps = {0, {pp_trait, 0}, pp_params};
                    this->_org_mngr.template add_<primary_producer_t>(pp_init_mass, ps);
                    pp_trait = pp_trait + 1.0;

                    trait_t_base pp_trait_max = std::numeric_limits<trait_t_base>::max();

                    if (this->_cfg["pp_trait_max"]) {
                        pp_trait_max = get_as<trait_t_base>("pp_trait_max", this->_cfg);
                    }


                    while (this->_dom.S(pp_trait, 0.0) > this->_dom.bm_threshold && pp_trait < pp_trait_max) {
                        pp_params[1] = this->_dom.S(pp_trait, 0.0);
                        //this->_dom.template add_<primary_producer>(pp_init_mass, {0, {pp_trait, 0}, pp_params});
                        typename DOM_T::pspace_t ps = {0, {pp_trait, 0}, pp_params};
                        this->_org_mngr.template add_<primary_producer_t>(pp_init_mass, ps);
                        //this->_dom.template add_<primary_producer>(pp_init_mass, {0, {-pp_trait, 0}, pp_params});
                        typename DOM_T::pspace_t ps2 = {0,{-pp_trait, 0}, pp_params};
                        this->_org_mngr.template add_<primary_producer_t>(pp_init_mass, ps2);

                        pp_trait = pp_trait + 1.0;
                    }
                }

                // Initialize range start - end with step
                if ( get_as<std::string>("type", i) == "range" ) {
                    double start = get_as<double>("start", i);
                    double end = get_as<double>("end", i);
                    double step = get_as<double>("step", i);

                    if ( end < start ) {
                        std::swap(start, end);
                    }

                    while ( start < end ) {
                        pp_params[1] = this->_dom.S(start, 0.0);
                        typename DOM_T::pspace_t ps = {0, {start, 0}, pp_params};
                        this->_org_mngr.template add_<primary_producer_t>(pp_init_mass, ps);
                        start += step;
                    }

                }

                // Initialize a single producer
                if ( get_as<std::string>("type", i) == "single" ) {
                    double trait = get_as<double>("trait", i);
                    pp_params[1] = this->_dom.S(trait, 0.0);
                    double init_mass;
                    if (i["mass"]) {
                       init_mass = get_as<double>("mass", i);
                    } else {
                        init_mass = pp_init_mass;
                    }
                    typename DOM_T::pspace_t ps = {0, {trait, 0}, pp_params};
                    this->_org_mngr.template add_<primary_producer_t>(init_mass, ps);

                }

            }
            this->_log->debug("Producers initialized");

        }

        void initialize_consumer() {

            // initialize species
            if (this->_cfg["init_species"]) {

                for (auto i : this->_cfg["init_species"])
                {
                    double mass;
                    //int type_v;
                    trait_t trait;
                    std::vector<double> params;

                    trait = get_as<trait_t>("trait", i);

                    if (i["mass"]){
                        mass = get_as<double>("mass", i);
                    } else {
                        mass = 10.0;
                    }

                    // Get params
                    if (i["params"]){
                        params = get_as<std::vector<double> >("params", i);
                    } else {
                        params = {1};
                    }

                    // Use global nieche pos
                    if (this->_cfg["fixed_p"]) {
                        this->_log->debug("Setting fixed p. Overrides trait 1");
                        trait[0] = get_as<double>("fixed_p", this->_cfg);
                        this->mutate[0] = false;
                    }

                    // Use global nieche width
                    if (this->_cfg["fixed_y"]) {
                        this->_log->debug("Setting fixed y. Overrides trait 2");
                        trait[1] = get_as<double>("fixed_y", this->_cfg);
                        this->mutate[1] = false;
                    }

                    // Use global handling rate width
                    if (this->_cfg["fixed_h"]) {
                        this->_log->debug("Setting fixed h. Overrides param 5");
                        params[4] = get_as<double>("fixed_h", this->_cfg);
                    }

                    // Use global cost width
                    if (this->_cfg["fixed_c"]) {
                        this->_log->debug("Setting fixed c. Overrides param 4");
                        params[3] = get_as<double>("fixed_c", this->_cfg);
                    }

                    // Use global mortality
                    if (this->_cfg["m"]) {
                        this->_log->debug("Setting m. Overrides consumer parameter 3");
                        params[2] = get_as<double>("m", this->_cfg);
                    }

                    // Use global food effectiveness
                    if (this->_cfg["b"]) {
                        this->_log->debug("Setting b. Overrides consumer parameter 2");
                        params[1] = get_as<double>("b", this->_cfg);
                    }

                    //this->_dom.template add_<consumer>(mass, {1, trait, params});
                    typename DOM_T::pspace_t ps = {1, trait, params};
                    this->_org_mngr.template add_<consumer_t>(mass, ps);
            

                }
            }

            this->_log->debug("Consumers initialized");
            
        }


        // Runtime functions ......................................................

        /* perform one timestep
         * One time step in this context is dt2
         *
        */
        void perform_step [[maybe_unused]]()
        {


            // With time step log

            // first perform dt1 steps for the system to develop
            double sum = 0.0;
            while (sum <= this->_dt2) {

                sum += this->_dt;

                if constexpr (DOM_T::step_func == 0) {

                } else if constexpr (DOM_T::step_func == 1) {

                    this->dts.push_back(this->_dt);
                } else {
                    throw std::invalid_argument("No valid 'perform_step' function. Check config.");
                }

                this->_dt = this->_dom.step(this->_dt);

                if (sum + this->_dt > this->_dt2) {
                    this->_dt = this->_dt2 - sum;
                }

                if (this->_dt < this->_min_dt)  {
                    this->_dt = this->_min_dt;
                }

            }

            // Mutate system
            this->mutation();

        }

        void mutation () {

            // TODO: Use multiple mutation methods

            // at dt2 mutate system
            double f;           // random number
            trait_t trait;      // trait value
            double mass;        // mass value

            trait_t_base trait_change;

            // iterate over all species
            for (auto v : this->_dom.get_vertices() ) {

                auto tmp = this->_dom[v];

                int type_val = tmp.get_type();

                // only active ones
                if (!tmp.active || type_val == 0 ){
                    continue;
                }

                // set values
                trait = tmp.get_trait();
                mass = tmp.get_mass();

                bool has_changed = false;
                // get the species type

                // check if random number exceeds value
                // depends on organism mass
                for (size_t j = 0; j < mutate.size(); ++j) {


                    // TODO: no need to iterate if mutation is off
                    // mutation on/off
                    if (!this->mutate[j]){
                        continue;
                    }
                    // get random number
                    f = this->dist(*this->_rng);

                    if ( f < this->mutation_rate && mass > this->_dom.bm_threshold * 5 ) {

                        if ( f < this->mutation_rate / 2.0 ) {
                            trait_change = -this->change[j];
                        } else {
                            trait_change = this->change[j];
                        }

                        // mutate the trait
                        trait[j] = trait[j] + trait_change;

                        if ( trait_change != 0 ) {
                            has_changed = true;
                        }

                    }
                }
                // depending on type add new species to network
                // initiate with one percent of parent species

                if ( !has_changed ) {
                    continue;
                }
                if ( type_val == 1 ) {

                    //auto tmp_vw = this->_dom.template add_<consumer>(0.01 * mass, {type_val, trait,
                    //                                                                 {tmp.org->parameters.params}});
                    typename DOM_T::pspace_t ps = {type_val, trait, {tmp.org->parameters.params}};
                    auto tmp_vw = this->_org_mngr.template add_<consumer_t>(0.01 * mass, ps);
                }

                // energy conservation in mutation process
                tmp.set_mass(mass*0.99);


            }

            // Revive consumers
            if (*(this->cons_spec_count) <= 0) {
                initialize_consumer();
            }

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

            // call save_dom
            // save_dom(this->_dom, this->_hdfgrp, std::to_string(this->_time));
            save_graph_vertex_value(this->_dom.graph, this->_hdfgrp, std::to_string(this->_time), "_types",
                                    [](auto& vw){return static_cast<char>(vw.get_type());});

            save_graph_vertex_value(this->_dom.graph, this->_hdfgrp, std::to_string(this->_time), "_traits",
                                    [](auto& vw){return vw.get_trait()[0];});

            save_graph_vertex_value(this->_dom.graph, this->_hdfgrp, std::to_string(this->_time), "_masses",
                                    [](auto& vw){return std::max(vw.get_mass(), 0.0);});

            save_graph_vertex_value(this->_dom.graph, this->_hdfgrp, std::to_string(this->_time), "_niche_w",
                                    [](auto& vw){return vw.org->get_niche_width();});

            if (get_as<bool>("safe_params",  this->_cfg)) {
                save_graph_vertex_value(this->_dom.graph, this->_hdfgrp, std::to_string(this->_time), "_params",
                                        [](auto &vw) { return vw.org->parameters.params; });
            }

            if ( get_as<bool>("safe_graph",  this->_cfg) ) {
                this->_log->debug("Saving Graph");
                save_graph_vertex_value(this->_dom.graph, this->_hdfgrp, std::to_string(this->_time), "_graph_vertex",
                                        [](auto& vw){return vw.id;});
                save_graph_edges(this->_dom.graph, this->_hdfgrp, std::to_string(this->_time), "_graph_edges");
                save_graph_edge_value(this->_dom.graph, this->_hdfgrp, std::to_string(this->_time), "_graph_edge_v",
                                      [](auto& vw){return vw.interaction_vector[0];});
            }

            save_graph_vertex_value(this->_dom.graph, this->_hdfgrp, std::to_string(this->_time), "_influx",
                                    [](auto& vw){return vw.org->get_influx();});

            if constexpr (DOM_T::step_func == 1) {
                save_graph_vertex_value(this->_dom.graph, this->_hdfgrp, std::to_string(this->_time), "_error",
                                        [](auto &vw) { return vw.org->integrator.get_error(); });


                save_container(this->dts, this->_hdfgrp, std::to_string(this->_time), "dts",
                               [](auto &it) { return it; });
            } else {

            }

            this->dts.clear();

        }


        // Getters and setters ....................................................
        // Add getters and setters here to interface with other model

    };

} // namespace Utopia::Models::MuLAN_MA

#endif // UTOPIA_MODELS_MULAN_MA_HH
