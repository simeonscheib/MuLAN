//
// Created by Simeon Scheib on 2019-01-31.
//

#ifndef UTOPIA_VERTEX_WRAPPER_BASE_HH
#define UTOPIA_VERTEX_WRAPPER_BASE_HH

namespace Utopia::Models::MuLAN_MA{

    /**
     * @brief Vertex_wrapper is the object stored on the graph
     * @details Holds a pointer to an organism of type ORG_T
     *
     * @tparam ORG_T
     */
    template <typename ORG_T>
    struct vertex_wrapper_base {

        using pspace_t = typename ORG_T::pspace_t;

        using trait_t = typename pspace_t::trait_t;

        /**
         * @brief Construct a new vertex wrapper base object
         *
         */
        vertex_wrapper_base() = default;

        /**
         * @brief Construct a new vertex wrapper base object
         *
         * @param id
         */
        explicit vertex_wrapper_base(int id) : id(id) {};

        /**
         * @brief Destroy the vertex wrapper base object
         *
         */
        ~vertex_wrapper_base() = default;

        /**
         * @brief ID of the vertex
         *
         */
        int id{};

        /**
         * @brief Is the vertex active?
         *
         */
        bool active = true;

        /**
         * @brief Pointer to organism
         *
         */
        typename ORG_T::org_ptr org;

        /**
         * @brief Get the type of the organism
         *
         * @return const int&
         */
        [[nodiscard]] const int& get_type() const;

        /**
         * @brief Get the id of the organism
         *
         * @return const int&
         */
        [[maybe_unused]] [[nodiscard]] const int& get_id() const;

        /**
         * @brief Get the value/biomass of the organism
         *
         * @return const ORG_T::currency&
         */
        const typename ORG_T::currency& get_mass() const;

        /**
         * @brief Get the trait of the organism
         *
         * @return const trait_t&
         */
        const trait_t& get_trait() const;

        /**
         * @brief Set the trait of teh organism
         *
         * @param new_trait
         */
        [[maybe_unused]] void set_trait(const trait_t& new_trait);

        /**
         * @brief Set the value/biomass of the organism
         *
         * @param new_mass
         */
        void set_mass(const typename ORG_T::currency& new_mass);

        bool operator==(const vertex_wrapper_base v2) const {
            return this->id == v2.id && this->org == v2.org && this->active == v2.active;
        }

        friend std::ostream& operator<< (std::ostream& stream, const vertex_wrapper_base v2) {
            stream << "[ Id: " << v2.id << " ";
            stream << "Organism: " << v2.org->parameters<< " ";
            stream << "Active: " << v2.active<< " ] ";
            return stream;
        }

    };


    // Getter & Setter

    template <typename ORG_T>
    const int& vertex_wrapper_base<ORG_T>::get_type() const {
        return this->org->parameters.type;
    }

    template <typename ORG_T>
    const typename ORG_T::currency& vertex_wrapper_base<ORG_T>::get_mass() const {
        return this->org->get_value();
    }

    template <typename ORG_T>
    const typename vertex_wrapper_base<ORG_T>::trait_t& vertex_wrapper_base<ORG_T>::get_trait() const {
        return this->org->parameters.trait;
    }

    template <typename ORG_T>
    [[maybe_unused]] void vertex_wrapper_base<ORG_T>::set_trait(const trait_t& new_trait) {

        this->org->parameters.trait = new_trait;

    }

    template <typename ORG_T>
    void vertex_wrapper_base<ORG_T>::set_mass( const typename ORG_T::currency& new_mass) {

        this->org->set_value(new_mass);

    }

    template<typename ORG_T>
    [[maybe_unused]] const int &vertex_wrapper_base<ORG_T>::get_id() const {
        return this->id;
    }

} // namespace Utopia::Models::MuLAN_MA

#endif //UTOPIA_VERTEX_WRAPPER_BASE_HH
