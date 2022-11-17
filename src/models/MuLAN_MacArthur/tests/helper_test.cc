#define BOOST_TEST_MODULE helper_test
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "utils/helper.hh"


namespace Utopia {
namespace Models {
namespace MuLAN_MA {

template<int N, int... pack>
void test_get_from_int_pack() 
{
    if constexpr ( N > 0 ) {
        constexpr int i = get_from_pack<int, N, pack...>();
        BOOST_TEST( i == N );
        test_get_from_int_pack<N - 1, pack...>();
    }
}

BOOST_AUTO_TEST_CASE (get_from_pack)
{
    test_get_from_int_pack<5, 0, 1, 2, 3, 4, 5>();
}


} // namespace MuLAN_MA
} // namespace Models
} // namespace Utopia