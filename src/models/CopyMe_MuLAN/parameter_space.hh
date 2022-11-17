//
// Created by Simeon Scheib on 2019-01-31.
//
#include "./MuLAN_base/_parameter_space.hh"

#ifndef UTOPIA_MY_PARAMETER_SPACE_HH
#define UTOPIA_MY_PARAMETER_SPACE_HH

namespace Utopia::Models::MuLAN_MA {


    // TODO: Substitute array and vector by trait space and parameter container
    /**
     * @brief Parameter Space
     * @details Holds parameters and traits
     *
     */
    struct parameter_space : parameter_space_base<std::array<double, 2>, std::vector<double> >  {

        /**
         * @brief Hash Value
         * @details Calculate an hash value for a parameter_space object to store in inordered_map
         */
        struct [[maybe_unused]] hash_value : hash_value_base<parameter_space>{

            /**
             * @brief Calculate the hash
             * @details     Combine the hashes of all members of the parameter_space
             *              To get an unique hash for each object
             *
             * @param a     The object to calculate the hash for
             * @return      The hash value
             */
            std::size_t operator()(const parameter_space& a) const override {
                std::size_t seed = boost::hash_range(a.trait.begin(), a.trait.end());
                boost::hash_combine(seed, boost::hash_range(a.params.begin(), a.params.end()));
                boost::hash_combine(seed, a.type);
                return seed;
            }
        };

        /**
         * @brief Stream operator
         * @details Allow usage of stream operator to print the parameter_space
         *
         * @param stream
         * @param ps
         * @return stream object
         */
        friend std::ostream& operator<< (std::ostream& stream, const parameter_space& ps) {
            stream << "[ Type: " << ps.type << " ";

            stream << "Trait: ";
            for (auto const &i : ps.trait) {
                stream << i << " ";
            }

            stream << "Parameters: ";
            for (auto const &i : ps.params) {
                stream << i << " ";
            }
            stream << "] ";
            return stream;
        }

    };

    // calculate hash for a parameter_space object


} // namespace Utopia::Models::MuLAN_MA


#endif //UTOPIA_MY_PARAMETER_SPACE_HH
