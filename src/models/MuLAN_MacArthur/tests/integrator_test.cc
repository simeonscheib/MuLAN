#define BOOST_TEST_MODULE integrator_test
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "../integrators/rkck.hh"
#include "../organism.hh"


namespace Utopia {
namespace Models {
namespace MuLAN_MA {

template <typename DOM_T>
class dummy_organism {

public:
    using currency = typename DOM_T::integrator_t::currency;

    const static int type_id = 0;

    /**
     * @brief Count number of species
     * 
     */
    static int spec_count;

    /**
     * @brief The Integrator 
     * @details Manages Integration
     * 
     */
    typename DOM_T::integrator_t integrator;

    /**
     * @brief Construct a new organism base object
     * 
     * @param d Domain
     * @param org This Child Object
     */
    explicit dummy_organism() : integrator(*this){

    };

    /**
     * @brief Destroy the organism base object
     * 
     */
    ~dummy_organism() = default;

    /**
     * @brief Get the value
     * @details E.g. Biomass of the organism
     * 
     * @return currency& 
     */
    virtual currency& get_value() {
        return this->integrator.get_value();
    };

    /**
     * @brief Set the value
     * 
     * @param val 
     */
    virtual void set_value(const currency& val) {
        return this->integrator.set_value(val);
    };

    /**
     * @brief Calculate the change rate
     * @details May depend on parameters, current value x and the integrator (e.g. sums, ..)
     *
     * @param x Value
     * @return currency
     */
    virtual currency dxdt(const currency& x, const double tpdt ) = 0;

};

/**
 * @brief dummy organism with simple logistic growth dgl
 * 
 * @tparam DOM_T 
 */
template <typename DOM_T>
class dummy_organism_1 : public dummy_organism<DOM_T> {
public:
    using currency = typename DOM_T::integrator_t::currency;
    /**
     * @brief Calculate the change rate
     * @details May depend on parameters, current value x and the integrator (e.g. sums, ..)
     *
     * @param x Value
     * @return currency
     */
    virtual currency dxdt(const currency& x, const double tpdt) {

        currency logistic_term = (x)/100.0;

        return x * (1 - logistic_term);
    };

};


/**
 * @brief Dummy domain for euler
 * 
 */
class dummy_domain_euler {
public:
    using org_t = dummy_organism_1<dummy_domain_euler>;
    using integrator_t = euler<double, dummy_organism<dummy_domain_euler> >;
    static constexpr typename integrator_t::currency error = 1.0;

    org_t orga;

    dummy_domain_euler () {
        orga.set_value(1.0);
    };
};

/**
 * @brief Dummy domain for rkck
 * 
 */
class dummy_domain_rkck {
public:
    using org_t = dummy_organism_1<dummy_domain_rkck>;
    using integrator_t = rkck<double, dummy_organism<dummy_domain_rkck> >;
    static constexpr typename integrator_t::currency& error = integrator_t::w_error;

    org_t orga;

    dummy_domain_rkck() {
        this->orga.set_value(1.0);
        this->orga.integrator.set_beta(0.95);
    }
};

/**
 * @brief Exact Solution of DGL
 * 
 * @param t 
 * @return double 
 */
double exact_value(double t) {

    return ( 100.0 )/( 1 + exp(-t) * 99.0 );
}


using AllDoms = boost::mpl::list<dummy_domain_euler, dummy_domain_rkck>;

// Test the integrators 
BOOST_AUTO_TEST_CASE_TEMPLATE (case1, ThisDom, AllDoms)
{
    ThisDom dom;

    auto last_value = dom.orga.get_value();
    
    double dt = 0.01;

    double t = 0.0;

    // integrate 1000 steps
    for (size_t steps = 0; steps < 1000; steps++) {

        double new_dt = std::numeric_limits<double>::max();
        
        for (size_t i = 0; i < ThisDom::integrator_t::steps; i++)
        {
        
            if (i == ThisDom::integrator_t::steps - 1) {

                dom.orga.integrator.step(dt);

                dom.orga.integrator.calc_new_step_size_s(dt, new_dt);
                
                ThisDom::integrator_t::calc_new_step_size(dt, new_dt);

            } else  {

                dom.orga.integrator.step(dt);

            }
        
        }

        // Test if error from analytical solution is smaller than specified error
        BOOST_TEST( exact_value(t) - dom.orga.get_value() <= ThisDom::error );

        t += dt;

        if (new_dt != std::numeric_limits<double>::max()){
            dt = new_dt;
        }

        // Test if solution is monotonically increasing
        BOOST_TEST(dom.orga.get_value() >= last_value - ThisDom::error);

        // Cross-check error estimate with specification
        BOOST_TEST(dom.orga.integrator.get_error() <= ThisDom::error);

        last_value = dom.orga.get_value();
    }
}


} // namespace MuLAN_MA
} // namespace Models
} // namespace Utopia