//
// Created by Simeon Scheib on 2019-01-31.
//

#ifndef UTOPIA_PARAMETER_SPACE_BASE_HH
#define UTOPIA_PARAMETER_SPACE_BASE_HH

// Every organism can be identified by its parameter space
// Two organisms with identical parameter space are equal
namespace Utopia::Models::MuLAN_MA {

    /**
     * @brief Parameter Space Base
     *
     * @tparam TRAIT_T Trait Type
     * @tparam PARAM_T Parameter Type
     */
    template <typename TRAIT_T, typename PARAM_T>
    struct parameter_space_base {
        /**
         * @brief Type/Layer ID
         * @details Integer to identify layer
         *
         */
        int type = 0;

        using trait_t [[maybe_unused]] = TRAIT_T;

        /**
         * @brief Trait
         * @details Holds unique parameters of the organism
         *
         */
        TRAIT_T trait;

        /**
         * @brief Parameters
         * @details Holds general parameters like e.g. growth rate
         *
         */
        PARAM_T params;

        /**
         * @brief Compare 2 Parameter Space Objects
         * @details Parameter Space Objects are equal if Type, parameter and trait are equal
         *
         * @param rhs
         * @return true
         * @return false
         */
        bool operator==(const parameter_space_base& rhs) const;

        /**
         * @brief Calculate hash for a parameter_space object
         *
         * @tparam PSPACE_T
         */
        template <typename PSPACE_T>
        struct hash_value_base{
            virtual std::size_t operator()(const PSPACE_T& a) const = 0;
        };
    };

    template <typename TRAIT_T, typename PARAM_T>
    bool parameter_space_base<TRAIT_T, PARAM_T>::operator==(const parameter_space_base<TRAIT_T, PARAM_T>& rhs) const {

        return (this->type == rhs.type) && (this->trait == rhs.trait) && (this->params == rhs.params);

    }

} // namespace Utopia::Models::MuLAN_MA

#endif //UTOPIA_PARAMETER_SPACE_BASE_HH
