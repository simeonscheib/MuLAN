
#ifndef UTOPIA_MODELS_EDGE_WRAPPER_BASE
#define UTOPIA_MODELS_EDGE_WRAPPER_BASE

namespace Utopia { namespace Models { namespace MuLAN_MA {


    /**
     * @brief Edge Wrapper stores edge weights
     *
     * @tparam EDGE_CONT
     */
    template <typename EDGE_CONT>
    struct edge_wrapper {

        using edge_cont = EDGE_CONT;

        EDGE_CONT interaction_vector;
    };



} } } // namespace Utopia::Models::MuLAN_MA

#endif //UTOPIA_MODELS_EDGE_WRAPPER_BASE
