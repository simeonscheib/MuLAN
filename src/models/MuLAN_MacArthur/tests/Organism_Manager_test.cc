#define BOOST_TEST_MODULE integrator_test
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/mpl/list.hpp>

#include "utils/Organism_Manager.hh"
#include "../domain.hh"

#include "../orga/primary_producer.hh"
#include "../orga/consumer_impl.hh"

namespace Utopia {
namespace Models {
namespace MuLAN_MA {

    template <typename VERTEX_T, typename CURRENCY, typename PSPACE_T>
    class dummy_domain_base {
        int counter = -1;

    public:

        using vertex_t = VERTEX_T;
        using currency = CURRENCY;
        using pspace_t = PSPACE_T;


        std::vector<std::pair<double, double> > v;
        int cfg = 0;

        template <typename ORG_T>
        vertex_t& add_(currency mass, pspace_t pspace) {

            std::cout << "add_ function called" << std::endl;

            v.push_back(std::make_pair(mass, pspace));

            ORG_T o;
            this->cfg = o.cfg;

            counter++;
            return this->counter;

        }
    };

    template <typename VERTEX_T>
    class dummy_domain : public dummy_domain_base<VERTEX_T, double, double> {
    public:
        using base_t = dummy_domain<VERTEX_T>;
    };

    template <typename VERTEX_T, typename CURRENCY, typename PSPACE_T, typename T>
    class dummy_domain2 : public dummy_domain_base<VERTEX_T, CURRENCY, PSPACE_T> {
    public:
        using base_t = dummy_domain_base<VERTEX_T,  CURRENCY, PSPACE_T>;
    };

    template <class A, class B, int...CFGs>
    class dummy_org {
    public:
        static constexpr int cfg = get_from_pack<int, 0, CFGs...>();
        static int spec_count;
    };

    template <typename DOM_T>
    struct wrapper {

        using type = std::function<typename DOM_T::vertex_t& (DOM_T&, typename DOM_T::currency, typename DOM_T::pspace_t)>;

        template <typename ORG_T>
        auto operator() (){

            return [](DOM_T& dom,
                    typename DOM_T::currency a,
                    typename DOM_T::pspace_t b) -> typename DOM_T::vertex_t& {
                        return dom.template add_<ORG_T>(a, b);
                    };
        }
    };

    template <typename A, typename B, int ...CFGs>
    [[maybe_unused]] int dummy_org<A, B, CFGs...>::spec_count = 0;

    template <int... Ns>
    using dummy_org_t = dummy_org<typename dummy_domain<int>::currency,
                                    typename  dummy_domain<int>::pspace_t,
                                    Ns...>;



    using Mngrs = boost::mpl::list<OrganismManager<dummy_domain<int>, wrapper , dummy_org_t>,
                            OrganismManager<dummy_domain2<int, int, int, int>, wrapper , dummy_org_t>
                >;

    // Test functions with dummies
    BOOST_AUTO_TEST_CASE_TEMPLATE (add, mngr, Mngrs)
    {

        // Test if number of species templates is computed correctly
        BOOST_TEST(mngr::get_number_of_species() == 1);
        typename mngr::dom_t d;
        mngr m(d);

        [[maybe_unused]] auto a = m.template register_<dummy_org_t, 4>({2});

        boost::test_tools::output_test_stream output;
        {
            cout_redirect guard( output.rdbuf( ) );

            m.template add_<dummy_org_t>(2.0, 2.0);
        }

        // Test if arguments are passed correctly
        BOOST_TEST((d.v[0].first == 2.0 && d.v[0].second ==  2.0));

        // Test if template parameters are set correctly
        BOOST_TEST(d.cfg == 2);

        // Test if correct function called
        BOOST_TEST(output.is_equal("add_ function called\n"));

        // Test if return value is passed correctly
        BOOST_TEST(m.template add_<dummy_org_t>(1.0, 1.0) == 1);
        BOOST_TEST(m.template add_<dummy_org_t>(2.0, 2.0) == 2);

    }

    template <int... Ns>
    using primary_producer_t = primary_producer<domain<double, 0, 1, 1, 0>::base_t,
                                                    domain<double, 0, 1, 1, 0>::pspace_t,
                                                    Ns...>;

    template <int... Ns>
    using consumer_t = consumer<domain<double, 0, 1, 1, 0>::base_t,
                                domain<double, 0, 1, 1, 0>::pspace_t,
                                Ns...>;

    using Mngrs2 = boost::mpl::list<OrganismManager<domain<double, 0, 1, 1, 0>, wrapper , primary_producer_t, consumer_t>
    >;

    // Test with "real" domain and organisms
    BOOST_AUTO_TEST_CASE_TEMPLATE (real_dom, mngr, Mngrs2)
    {
        // Test if number of species templates is computed correctly
        BOOST_TEST(mngr::get_number_of_species() == 2);

        typename mngr::dom_t d;
        mngr m(d);

        auto a = m.template register_<primary_producer_t, 2>({1});
        auto b = m.template register_<consumer_t, 2>({1});

        typename mngr::dom_t::pspace_t ps = {0, {1, 0}, {1,2,3}};
        typename mngr::dom_t::pspace_t ps2 = {0, {1, 1}, {1,2,3}};

        auto v1 = m.template add_<primary_producer_t>(2.0, ps);
        // Test if return value is correct and can access mass
        BOOST_TEST(v1.get_mass() == 2.0);

        // Same test
        for (auto v : boost::make_iterator_range(vertices(d.graph))) {
            // Test if arguments are passed correctly
            BOOST_TEST((d[v].get_mass() == 2.0));
        }

        // Same Test (access via parameter_space)
        BOOST_TEST((d[ps].get_mass() == 2.0));

        // Add a second pp
        m.template add_<primary_producer_t>(3.0, ps2);

        // Check if correct number of species has been added
        BOOST_TEST(a->get_spec_count() == 2);
        BOOST_TEST(b->get_spec_count() == 0);

    }

}
}
}