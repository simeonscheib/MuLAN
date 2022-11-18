#ifndef UTOPIA_MODELS_MULAN_UTILS_HH
#define UTOPIA_MODELS_MULAN_UTILS_HH

//#include <dune/utopia/core/logging.hh>
//#include <dune/utopia/data_io/hdfgroup.hh>
#include <functional>


namespace Utopia {
namespace Models {
namespace MuLAN_MA {

    template <typename DOM_T>
    [[maybe_unused]] std::shared_ptr<Utopia::DataIO::HDFGroup> save_dom(DOM_T &dom,
                                         const std::shared_ptr<Utopia::DataIO::HDFGroup>& parent_grp,
                                         const std::string& name)
    {
        // Collect some information on the graph
        auto num_vertices = boost::num_vertices(dom.graph);


        // Create the group for the graph and store metadata in its attributes
        auto grp = parent_grp->open_group(name);

        // save the number of species
        grp->add_attribute("num_species", num_vertices);

        // Initialize datasets to store
        // types
        auto dset_types = grp->open_dataset("_types", {num_vertices});
        // traits
        auto dset_traits = grp->open_dataset("_traits", {num_vertices});
        // biomasses
        auto dset_masses = grp->open_dataset("_masses", {num_vertices});
        // last change rate dx/dt
        auto dset_nw = grp->open_dataset("_niche_w", {num_vertices});


        // Save
        auto [v, v_end] = boost::vertices(dom.graph);

        dset_types->write(v, v_end,
                       [&](auto vd){return static_cast<char>(dom.graph[vd].get_type());}
        );

        dset_traits->write(v, v_end,
                          [&](auto vd){return dom.graph[vd].get_trait()[0];}
        );

        dset_masses->write(v, v_end,
                          [&](auto vd){return std::max(dom.graph[vd].get_mass(), 0.0);}
        );

        dset_nw->write(v, v_end,
                           [&](auto vd){return dom.graph[vd].org->get_niche_width();}
        );

        auto log = spdlog::get("data_io");
        log->debug("Domain saved.");

        // Return the newly created group
        return grp;
    }


    /** This function writes the result of a function applied to all vertices of a boost::graph into a HDFGroup.
    *
    * @tparam GraphType
    * @tparam FUNC_T
    *
    * @param g The graph to save
    * @param parent_grp The parent HDFGroup the graph should be stored in
    * @param name The name the newly created graph group should have
    * @param dset_name The name of the property to be saved
    * @param get_value The functional applied to all vertices of the graph in order to extract the property
    *
    * @return std::shared_ptr<HDFGroup> The newly created graph group
    */

    template <typename GraphType, typename FUNC_T>
    std::shared_ptr<Utopia::DataIO::HDFGroup> save_graph_vertex_value(GraphType &g,
                                                       const std::shared_ptr<Utopia::DataIO::HDFGroup>& parent_grp,
                                                       const std::string& name,
                                                       const std::string& dset_name,
                                                       FUNC_T&& get_value)
    {
        // Collect some information on the graph
        auto num_vertices = boost::num_vertices(g);

        auto log = spdlog::get("data_io");

        // Create the group for the graph and store metadata in its attributes
        auto grp = parent_grp->open_group(name);

        // save the number of species
        grp->add_attribute("num_species", num_vertices);

        // Initialize datasets to store
        auto dset = grp->open_dataset(dset_name, {num_vertices});

        // Save
        auto [v, v_end] = boost::vertices(g);

        if (num_vertices != 0) {
            dset->write(v, v_end,
                        [&get_value, &g](auto &vd) { return get_value(g[vd]); }
            );
        }

        // debug info
        log->debug("Graph {} ({} property) saved.", name, dset_name, num_vertices);

        // Return the newly created group
        return grp;
    }


    /** This function writes the result of a function applied to all edges of a boost::graph into a HDFGroup.
    *
    * @tparam GraphType
    * @tparam FUNC_T
    *
    * @param g The graph to save
    * @param parent_grp The parent HDFGroup the graph should be stored in
    * @param name The name the newly created graph group should have
    * @param dset_name The name of the property to be saved
    * @param get_value The functional applied to all vertices of the graph in order to extract the property
    *
    * @return std::shared_ptr<HDFGroup> The newly created graph group
    */

    template <typename GraphType>
    std::shared_ptr<Utopia::DataIO::HDFGroup> save_graph_edges(GraphType &g,
                                                                      const std::shared_ptr<Utopia::DataIO::HDFGroup>& parent_grp,
                                                                      const std::string& name,
                                                                      const std::string& dset_name
                                                                      )
    {
        // Collect some information on the graph
        std::size_t num_edges = boost::num_edges(g);

        auto log = spdlog::get("data_io");

        // Create the group for the graph and store metadata in its attributes
        auto grp = parent_grp->open_group(name);

        // save the number of species
        grp->add_attribute("num_edges", num_edges);

        // Initialize datasets to store
        auto dset = grp->open_dataset(dset_name, {2, num_edges});

        // Save
        auto [e, e_end] = boost::edges(g);

        if (num_edges != 0) {
            dset->write(e, e_end,
                        [&g](auto vd) { return g[source(vd, g)].id; }
            );

            dset->write(e, e_end,
                        [&g](auto vd) { return g[target(vd, g)].id; }
            );
        }

        // debug info
        log->debug("Graph {} ({} property) saved.", name, dset_name, num_edges);

        // Return the newly created group
        return grp;
    }


    /** This function writes the result of a function applied to all edges of a boost::graph into a HDFGroup.
    *
    * @tparam GraphType
    * @tparam FUNC_T
    *
    * @param g The graph to save
    * @param parent_grp The parent HDFGroup the graph should be stored in
    * @param name The name the newly created graph group should have
    * @param dset_name The name of the property to be saved
    * @param get_value The functional applied to all vertices of the graph in order to extract the property
    *
    * @return std::shared_ptr<HDFGroup> The newly created graph group
    */

    template <typename GraphType, typename FUNC_T>
    std::shared_ptr<Utopia::DataIO::HDFGroup> save_graph_edge_value(GraphType &g,
                                                               const std::shared_ptr<Utopia::DataIO::HDFGroup>& parent_grp,
                                                               const std::string& name,
                                                               const std::string& dset_name,
                                                               FUNC_T&& get_value)
    {
        // Collect some information on the graph
        auto num_edges = boost::num_edges(g);

        auto log = spdlog::get("data_io");

        // Create the group for the graph and store metadata in its attributes
        auto grp = parent_grp->open_group(name);

        // save the number of species
        grp->add_attribute("num_edges", num_edges);

        // Initialize datasets to store
        auto dset = grp->open_dataset(dset_name, {num_edges});

        // Save
        auto [e, e_end] = boost::edges(g);

        if (num_edges != 0) {
            dset->write(e, e_end,
                        [&get_value, &g](auto vd) { return get_value(g[vd]); }
            );
        }

        // debug info
        log->debug("Graph {} ({} property) saved.", name, dset_name, num_edges);

        // Return the newly created group
        return grp;
    }

    template <typename ContainerType, typename FUNC_T>
    std::shared_ptr<Utopia::DataIO::HDFGroup> save_container(ContainerType &container,
                                                       const std::shared_ptr<Utopia::DataIO::HDFGroup>& parent_grp,
                                                       const std::string& name,
                                                       const std::string& dset_name,
                                                       FUNC_T&& get_value)
    {

        auto log = spdlog::get("data_io");

        // Create the group for the graph and store metadata in its attributes
        auto grp = parent_grp->open_group(name);

        auto size = container.size();

        // Initialize datasets to store
        auto dset = grp->open_dataset(dset_name, {size});

        // Save
        auto v = container.begin();
        auto v_end = container.end();

        dset->write(v, v_end, get_value);

        // debug info
        log->debug("Container {} ({} property) saved.", name, dset_name);

        // Return the newly created group
        return grp;
    }


} // namespace MuLAN_MA
} // namespace Models
} // namespace Utopia

#endif // UTOPIA_MODELS_MULAN_UTILS_HH
