#include <iostream>

#include <utopia/core/model.hh>
// #include <utopia/core/types.hh>
#include "MuLAN_MA.hh"

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
    using cfg_map = std::list<std::tuple<std::string, std::map<std::string, int> > >;
    try {
        /// Dune::MPIHelper::instance(argc, argv);

        /// Initialize the PseudoParent from config file path
        Utopia::PseudoParent pp(argv[1]);

        /// map config entries to ints for model building
        const cfg_map config_map {
            {"integrator", {
                    {"rkck", 0},
                    {"euler", 1}
                }
            },
            {"step_func", {
                    {"normal", 0},
                    {"normal_with_step", 1}
                }
            },
            {"response_func", {
                    {"none", 0},
                    {"linear", 1},
                    {"linear_cut", 1},
                    {"type_2", 2},
                    {"type_3", 2}
                }
            },
            {"env_func", {
                    {"gaussian", 0},
                    {"gausscos", 1},
                    {"gausscostd", 2},
                    {"gausscos2td", 3}
                }
            }
        };

        /// ModelBuilder Object
        auto Builder = ModelBuilder(pp, config_map, "MuLAN_MA");

        /// Run the Model with Run::run() config from pp
        /// integers behind the runner class denote the number of config parameters
        Builder.build<Run<double>, 2,2,3,4>(pp);

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
