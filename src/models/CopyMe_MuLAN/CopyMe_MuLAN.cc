#include <iostream>
#include <utopia/core/model.hh>

#include "CopyMe_MuLAN.hh"
#include "utils/Config_Builder.hh"

using namespace Utopia::Models::MuLAN_MA;

/**
 * @brief Wrapper to run a model
 * @details Pass this object to the ModelBuilder together with a run config
 *
 * @tparam Ts Types to pass to the domain class
 */
template < typename ...Ts>
class Run {
public:

    /**
     * @brief Run the MuLAN_MAModel with domain<Ts..., Bs...>
     *
     * @tparam ParentType
     * @tparam Bs       ints to pass to the domain
     * @param pp        Parent model
     */
    template<typename ParentType, int ...Bs>
    static inline void run(ParentType& pp){
        MuLAN_MAModel<domain<Ts..., Bs...> > model("MuLAN_MA", pp);
        model.run();
    }
};

int main (int /*argc*/, char** argv)
{
    // a map of config parameters to integer
    using cfg_map = std::list<std::tuple<std::string, std::map<std::string, int> > >;

    try {

        /// Initialize the PseudoParent from config file path
        Utopia::PseudoParent pp(argv[1]);

        /// map config entries to ints for model building
        // TODO: Add your own configuration options that are evaluated ad compile time
        const cfg_map config_map {
                {"config_parameter_1", {
                                       {"A1", 0},
                                       {"A2", 1}
                               }
                },
                {"config_parameter_2", {
                                       {"B1", 0},
                                       {"B2", 1}
                               }
                }
        };

        /// ModelBuilder Object
        auto Builder = ModelBuilder(pp, config_map, "MuLAN_MA");

        /// Run the Model with Run::run() config from pp
        /// integers behind the runner class denote the number of config parameters
        // TODO: Adapt number of parameter for every configuration option
        Builder.build<Run<double>, 2,2>(pp);

        return 0;
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Exception occurred!" << std::endl;
        return 1;
    }
}