#include <boost/graph/graphviz.hpp>
#include "_organism.hh"
#include "_edge_wrapper.hh"

#ifndef UTOPIA_MODELS_DOMAIN_BASE_HH
#define UTOPIA_MODELS_DOMAIN_BASE_HH


namespace Utopia::Models::MuLAN_MA {


/**
 * @brief Domain Class
 * @details Describes an environment or part of an environment where species networks can form depending on attributes of the domain.
 * 
 * @tparam DOM_T Child Type
 * @tparam PSPACE_T Parameter Space Type of the Entities in the Graph
 * @tparam VERTEX_T Vertex Type holding the Entities
 * @tparam EDGE_T Edge Type
 */
template <typename DOM_T,
            typename CURRENCY,
            typename PSPACE_T,
            typename VERTEX_T,
            typename EDGE_T,
            std::size_t SUM_SIZE>
class domain_base {

protected:
    /**
     * @brief Counter to give new species a new id
     * 
     */
    int oid = 0;

    template<int N>
    void calculate_sums();

    template<int N, int MAX>
    void calculate_all_sums_iter();

    /*
     * @brief time of the simulation
     */
    double time = 0.0;

public:

    // set aliases
    using edge_t = EDGE_T;
    using vertex_t = VERTEX_T;
    using graph_t = boost::adjacency_list<boost::listS, boost::listS, boost::directedS, VERTEX_T, EDGE_T>;
    using vertex_desc_t = typename boost::graph_traits<graph_t>::vertex_descriptor;
    using edge_desc_t = typename boost::graph_traits<graph_t>::edge_descriptor;
    using vertex_iter = typename boost::graph_traits<graph_t>::vertex_iterator;
    using out_edge_iter = typename boost::graph_traits<graph_t>::out_edge_iterator;
    using in_edge_iter = typename boost::graph_traits<graph_t>::in_edge_iterator;
    using edge_cont = typename EDGE_T::edge_cont;

    using pspace_t = PSPACE_T;
    using currency = CURRENCY;

    using hash_t = typename PSPACE_T::hash_value;

    // map parameters space to vertex to avoid multiple species with same properties
    using parameter_space_map = boost::unordered_map<PSPACE_T, vertex_desc_t, hash_t>;

    /**
     * @brief Construct a new domain base object
     * 
     */
    domain_base();


    /**
     * @brief Interaction Threshold
     * @details Defines the threshold to which interactions are considered
     * 
     */
    double interaction_tolerance = 1.0e-3;

    /**
     * @brief Biomass Threshold
     * @details Defines the threshold an organism needs to survive
     * 
     */
    double bm_threshold = 0.05;

    /**
     * @brief Food Network
     * 
     */
    graph_t graph;

    /**
     * @brief Parameter Space Map
     * @details Map parameter spaces to vertex descriptor
     * 
     */
    parameter_space_map psm;

    /**
     * @brief Add a species of type ORG_T
     * 
     * @tparam ORG_T Species type
     * @param mass Initial Mass
     * @param ps Parameter Space
     * @return VERTEX_T& 
     */
    template <template<class, class> class ORG_T>
    VERTEX_T& add_(double mass, PSPACE_T ps);

    /**
     * @brief Add a species of type ORG_T
     * 
     * @tparam ORG_T Species type
     * @param mass Initial Mass
     * @param ps Parameter Space
     * @return VERTEX_T& 
     */
    template <typename ORG_T>
    VERTEX_T& add_(currency mass, PSPACE_T ps);

    /**
     * @brief Calculate sums depending on graph
     */
    inline void calculate_all_sums();

    /**
     * @brief Perform one time_step of dt
     * 
     * @param dt Time Step
     * @return double New Time Step (adaptive Time step)
     */
    double step(double dt);

    /**
     * @brief Set the Interaction Tolerance
     * 
     * @param new_tolerance 
     */
    void set_interaction_tolerance(double new_tolerance);

    /**
     * @brief Set the Biomass Threshold
     * 
     * @param new_threshold 
     */
    void set_bm_threshold(double new_threshold);

    /**
     * @brief Access vertex v
     * 
     * @param v 
     * @return VERTEX_T& 
     */
    VERTEX_T& operator[](vertex_desc_t v);

    /**
     * @brief Access vertex v via its parameters
     *
     * @param v
     * @return VERTEX_T&
     */
    VERTEX_T& operator[](pspace_t v);

    /*
     * @brief get time of the simulation
     */
    double get_time() const {
        return this->time;
    }

    /*
     * @brief Set the simulation time
     * @details For special porposes _domain allows changing simulation time
     *          Note that changing time is ill advised under the temporal prime directive
     */
    void set_time(double t) {
        this->time = std::move(t);
    }

    /*
     * @brief Simplify access to vertex iterator_range
     *
     * @return boost::iterator_range<vertex_iter>
     */
    inline auto get_vertices() const {
        return boost::make_iterator_range(vertices(this->graph));
    }

    /*
     * @brief Simplify access to out_edge iterator_range
     *
     * @param v Vertex to get out_edges for
     * @return boost::iterator_range<out_edge_iter>
     */
    inline auto get_out_edges(const vertex_desc_t& v) const {
        return boost::make_iterator_range(out_edges(v, this->graph));
    }

    /*
     * @brief Simplify access to in_edge iterator_range
     *
     * @param v Vertex to get in_edges for
     * @return boost::iterator_range<in_edge_iter>
     */
    inline auto get_in_edges(const vertex_desc_t& v) const {
        return boost::make_iterator_range(in_edges(v, this->graph));
    }

    /*
     * @brief Simplify access to out_degree
     *
     * @param v Vertex to get out_degree for
     * @return int Degree
     */
    inline auto get_out_degree(const vertex_desc_t& v) const {
        return out_degree(v, this->graph);
    }

    /*
     * @brief Simplify access to in_degree
     *
     * @param v Vertex to get in_degree for
     * @return int Degree
     */
    inline auto get_in_degree(const vertex_desc_t& v) const {
        return in_degree(v, this->graph);
    }

protected:

    /**
     * @brief Add edges to given vertex
     * 
     * @param reference The Vertex
     */
    void add_edges(vertex_desc_t reference);

};

template <typename DOM_T, typename CURRENCY, typename PSPACE_T, typename VERTEX_T, typename EDGE_T, std::size_t SUM_SIZE>
domain_base<DOM_T, CURRENCY, PSPACE_T, VERTEX_T, EDGE_T, SUM_SIZE>::domain_base() = default;

// Take an initial mass and a parameter space to add species
template <typename DOM_T, typename CURRENCY, typename PSPACE_T, typename VERTEX_T, typename EDGE_T, std::size_t SUM_SIZE>
template <typename ORG_T>
VERTEX_T& domain_base<DOM_T, CURRENCY, PSPACE_T, VERTEX_T, EDGE_T, SUM_SIZE>::add_(currency mass, PSPACE_T ps) {

    ps.type = ORG_T::type_id;

    // insert pair in map
    auto iter_ins_pair = this->psm.insert(std::make_pair(ps, vertex_desc_t()));

    // If the species does not exist ...
    if (iter_ins_pair.second){

        // ... add it
        vertex_desc_t v = add_vertex(this->graph);

        iter_ins_pair.first->second = v;

        VERTEX_T& tmp_vertex = this->graph[v];

        // And set
        tmp_vertex.id = oid;
        tmp_vertex.org = boost::make_shared<ORG_T>(static_cast<DOM_T*>(this), mass);

        tmp_vertex.org->parameters = ps;

        // calculate the edges from and to the new species
        this->add_edges(v);

        // increment id
        oid++;

        tmp_vertex.org->count(+1);

        return tmp_vertex;

    } else {

        // ... add the biomass to the existing species
        vertex_desc_t v_tmp = iter_ins_pair.first->second;

        VERTEX_T& tmp_vertex = this->graph[v_tmp];

        tmp_vertex.set_mass(tmp_vertex.get_mass() + mass);

        if (!tmp_vertex.active){
            tmp_vertex.active = true;
            tmp_vertex.org->count(1);
        }

        return tmp_vertex;

    }
}

template <typename DOM_T, typename CURRENCY, typename PSPACE_T, typename VERTEX_T, typename EDGE_T, std::size_t SUM_SIZE>
template <template<typename, typename> typename ORG_T>
VERTEX_T& domain_base<DOM_T, CURRENCY, PSPACE_T, VERTEX_T, EDGE_T, SUM_SIZE>::add_(double mass, PSPACE_T ps) {

    return this->add_<ORG_T<DOM_T, PSPACE_T> >(mass, ps);

}

// method to add edges from and to a species
template <typename DOM_T, typename CURRENCY, typename PSPACE_T, typename VERTEX_T, typename EDGE_T, std::size_t SUM_SIZE>
void domain_base<DOM_T, CURRENCY, PSPACE_T, VERTEX_T, EDGE_T, SUM_SIZE>::add_edges(vertex_desc_t reference) {


    edge_cont z1, z2;

    // iterate over all organisms
    for (auto v : this->get_vertices() ) {

        // call the interaction coefficient method of the species with every other species
        z1 = this->graph[reference].org->calc_interaction_coeff(this->graph[v].org);

        // accepts the interaction if the value exceeds a threshold
        for (auto const &i: z1) {
            if (i > this->interaction_tolerance) {

                // add the edge
                add_edge(reference, v, {z1}, this->graph);
                break;
            }
        }

        // do the same the other way round (asymmetrical)
        z2 = this->graph[v].org->calc_interaction_coeff(this->graph[reference].org);

        // avoid two edges to the organism itself
        if ( reference != v ){
            for (auto const &i: z2) {
                if ( i > this->interaction_tolerance ) {

                    // add edge
                    add_edge(v, reference, {z2}, this->graph);
                    break;

                }
            }
        }
    }

}

template <typename DOM_T, typename CURRENCY, typename PSPACE_T, typename VERTEX_T, typename EDGE_T, std::size_t SUM_SIZE>
template<int N>
void domain_base<DOM_T, CURRENCY, PSPACE_T, VERTEX_T, EDGE_T, SUM_SIZE>::calculate_sums() {

    vertex_desc_t vertex1;

    for (auto v: this->get_vertices() ) {

        VERTEX_T &vertex2 = this->graph[v];

        if (!vertex2.active) {
            continue;
        }

        // iterate over all edges
        for (auto e : this->get_out_edges(v) ) {

            vertex1 = target(e, this->graph);

            if (!this->graph[vertex1].active) {
                continue;
            }

            // Add the edge values to the vertex
            vertex2.org->template add_edge_cont<N>(this->graph[e].interaction_vector, this->graph[vertex1].org);

        }
    }
}

template <typename DOM_T, typename CURRENCY, typename PSPACE_T, typename VERTEX_T, typename EDGE_T, std::size_t SUM_SIZE>
template<int N, int MAX>
void domain_base<DOM_T, CURRENCY, PSPACE_T, VERTEX_T, EDGE_T, SUM_SIZE>::calculate_all_sums_iter() {
    if constexpr ( N == MAX ) {
        return;
    } else {
        this->template calculate_sums<N>();
        this->calculate_all_sums_iter<N + 1, MAX>();
    }
}

template <typename DOM_T, typename CURRENCY, typename PSPACE_T, typename VERTEX_T, typename EDGE_T, std::size_t SUM_SIZE>
inline void domain_base<DOM_T, CURRENCY, PSPACE_T, VERTEX_T, EDGE_T, SUM_SIZE>::calculate_all_sums() {
    this->calculate_all_sums_iter<0, SUM_SIZE>();
}

    // perform one timestep of the whole domain
template <typename DOM_T, typename CURRENCY, typename PSPACE_T, typename VERTEX_T, typename EDGE_T, std::size_t SUM_SIZE>
double domain_base<DOM_T, CURRENCY, PSPACE_T, VERTEX_T, EDGE_T, SUM_SIZE>::step(double dt) {

    // Increment Time
    this->time += dt;

    double new_dt = std::numeric_limits<double>::max();

    for (size_t i = 0; i < DOM_T::integrator_t::steps; i++)
    {
    
    
        /*
        * In the first step the edges are added to the 'sum' vector of a vertex
        * This allows to calculate sums over all links e.g. sum all biomass to share carrying capacity
        * How exactly the edge is added is defined for each species class.
        */

        this->calculate_all_sums();
    
        if (i == DOM_T::integrator_t::steps - 1) {


            // perform a step for each organism
            for (auto v : this->get_vertices() ) {
                VERTEX_T& vertex2 = this->graph[v];

                vertex2.org->integrator.step(dt);

                vertex2.org->integrator.calc_new_step_size_s(dt, new_dt);

                if (vertex2.get_mass() < this->bm_threshold && vertex2.active){
                    vertex2.set_mass(0.0);
                    vertex2.active = false;
                    vertex2.org->count(-1);
                }
            }
            DOM_T::integrator_t::calc_new_step_size(dt, new_dt);
        } else  {
            // perform a step for each organism
            for (auto v : this->get_vertices() ) {

                VERTEX_T& vertex2 = this->graph[v];
                vertex2.org->integrator.step(dt);

            }
        }
    
    }

    // If timestep has not been adapted use the old one
    if (new_dt == std::numeric_limits<double>::max()){
        return dt; 
    } else {
        return new_dt;

    }
}

// Getter & Setter


template <typename DOM_T, typename CURRENCY, typename PSPACE_T, typename VERTEX_T, typename EDGE_T, std::size_t SUM_SIZE>
VERTEX_T& domain_base<DOM_T, CURRENCY, PSPACE_T, VERTEX_T, EDGE_T, SUM_SIZE>::operator[](const vertex_desc_t v) {

    return this->graph[v];
}

template <typename DOM_T, typename CURRENCY, typename PSPACE_T, typename VERTEX_T, typename EDGE_T, std::size_t SUM_SIZE>
VERTEX_T& domain_base<DOM_T, CURRENCY, PSPACE_T, VERTEX_T, EDGE_T, SUM_SIZE>::operator[](const pspace_t v) {

    return this->graph[this->psm[v]];
}

template <typename DOM_T, typename CURRENCY, typename PSPACE_T, typename VERTEX_T, typename EDGE_T, std::size_t SUM_SIZE>
void domain_base<DOM_T, CURRENCY, PSPACE_T, VERTEX_T, EDGE_T, SUM_SIZE>::set_interaction_tolerance(double new_tolerance) {

    this->interaction_tolerance = new_tolerance;

}

template <typename DOM_T, typename CURRENCY, typename PSPACE_T, typename VERTEX_T, typename EDGE_T, std::size_t SUM_SIZE>
void domain_base<DOM_T, CURRENCY, PSPACE_T, VERTEX_T, EDGE_T, SUM_SIZE>::set_bm_threshold(double new_threshold) {

    this->bm_threshold = new_threshold;

}


} // namespace Utopia::Models::MuLAN_MA

#endif //UTOPIA_MODELS_DOMAIN_BASE_HH
