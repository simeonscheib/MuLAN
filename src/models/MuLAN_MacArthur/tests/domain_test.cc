#define BOOST_TEST_MODULE domain_test
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "../domain.hh"

#include "../orga/primary_producer.hh"
#include "../orga/consumer_impl.hh"

namespace Utopia {
namespace Models {
namespace MuLAN_MA {

using Doms = boost::mpl::list<domain<double, 0, 0, 1, 0> >;

BOOST_AUTO_TEST_CASE_TEMPLATE (necessary_funcs, Dom, Doms)
{
    Dom dom;

    // Test setter for Error
    dom.set_error(0.2);
    BOOST_TEST( Dom::integrator_t::w_error == 0.2 );

    dom.set_interaction_tolerance(1.0e-2);
    BOOST_TEST( dom.interaction_tolerance == 1.0e-2 );

    dom.params = std::vector<double>({100.0, 100.0});
    BOOST_TEST( dom.S(0.0, 0.0) == 100.0 );

    std::vector<double> pp_params = {10, dom.S(0.0, 0.0)};
    typename Dom::pspace_t ps = {0, {0, 0}, pp_params};

    auto pp = dom.template add_<primary_producer<typename Dom::base_t, typename Dom::pspace_t, 0> >(5, ps);

    BOOST_TEST( num_vertices(dom.graph) == 1 );
    BOOST_TEST( num_edges(dom.graph) == 0 );

    BOOST_TEST( pp.get_mass() == 5 );
    BOOST_TEST( pp.org->parameters == ps );

    typename Dom::pspace_t ps_c1 = {1, {0, 2}, {0.7, 0.8, 5.0, 0.0}};
    auto c1 = dom.template add_<consumer<typename Dom::base_t, typename Dom::pspace_t, 1> >(5, ps_c1);
    BOOST_TEST( c1.get_mass() == 5 );
    BOOST_TEST( c1.org->parameters == ps_c1 );

    typename Dom::pspace_t ps_c2 = {1, {50, 2}, {0.7, 0.8, 5.0, 0.0}};
    auto c2 = dom.template add_<consumer<typename Dom::base_t, typename Dom::pspace_t, 0> >(5, ps_c2);

    // Three vertices
    BOOST_TEST(num_vertices(dom.graph) == 3);

    // Only one consumer is close enough to the producer to form an edge
    BOOST_TEST(num_edges(dom.graph) == 1);

    // Test the graph building
    for (auto v : dom.get_vertices() ) {
        for (auto e: dom.get_out_edges(v) ) {

            // No self edges
            BOOST_TEST( v != target(e, dom.graph));

            // Check if the edge connects the correct vertices
            BOOST_TEST( dom.graph[v] == c1 );
            BOOST_TEST( dom.graph[target(e, dom.graph)] == pp );

            // Does the edge have the correct value defined by
            BOOST_TEST( dom.graph[v].org->calc_interaction_coeff(dom.graph[target(e, dom.graph)].org)
                        == dom.graph[e].interaction_vector );

            // Is there only one edge
            BOOST_TEST( dom.get_out_degree(v) == 1);
        }
    }

    // Propagate the system

    double dt = 0.01;
    int iters = 1000;

    for (int i = 0; i < iters; i++) {
        dt = dom.step(dt);
    }

    auto m = dom[ps].get_mass();
    auto m1 = dom[ps_c1].get_mass();
    auto m2 = dom[ps_c2].get_mass();

    // This should hold for this particular system for t-> infinity

    // The producer is eaten by consumer 1
    // and cant grow over 50% for these params
    BOOST_TEST( m < 50 );

    // Comsumer 1 has food and can sustain biomass
    BOOST_TEST( m1 > 0 );

    // Consumer 2 has no food and dies out
    BOOST_TEST( m2 == 0 );

}


} // namespace MuLAN_MA
} // namespace Models
} // namespace Utopia