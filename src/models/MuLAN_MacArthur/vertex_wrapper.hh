#include "./MuLAN_base/_vertex_wrapper.hh"

#ifndef UTOPIA_MY_VERTEX_WRAPPER_HH
#define UTOPIA_MY_VERTEX_WRAPPER_HH

namespace Utopia::Models::MuLAN_MA {

    /**
     * @brief Vertex Wrapper
     * @details Concretize Vertex Wrapper Base
     *
     * @tparam ORG_T
     */
    template <typename ORG_T>
    struct vertex_wrapper : vertex_wrapper_base<ORG_T> {

        /**
         * @brief Construct a new vertex wrapper object
         *
         */
        vertex_wrapper() = default;

        /**
         * @brief Construct a new vertex wrapper object
         *
         * @param id
         */
        explicit vertex_wrapper(int id);

        /**
         * @brief Destroy the vertex wrapper object
         *
         */
        ~vertex_wrapper() = default;

    };


    template <typename ORG_T>
    vertex_wrapper<ORG_T>::vertex_wrapper(int id) : vertex_wrapper_base<ORG_T>(id) {

    }




} // namespace Utopia::Models::MuLAN_MA


#endif //UTOPIA_MY_VERTEX_WRAPPER_HH
