#include <charconv>

#ifndef UTOPIA_MODELS_MULAN_MA_HH_HELPER
#define UTOPIA_MODELS_MULAN_MA_HH_HELPER

/**
 * @brief Hold arbitrary number of rvalues
 * 
 * @tparam Ns 
 */
template<int ... Ns>
struct IntTuple {
};

/**
 * @brief Get the ith Entry in a parameter pack of type T
 *
 * @tparam T        Type
 * @tparam i        index
 * @tparam Ints     Parameter pack of type T
 * @return          The ith entry
 */
template<class T,std::size_t i, T... Ints>
constexpr T get_from_pack() {
    constexpr T arr[] = {Ints...};
    return arr[i];
}

template <typename T>
T to_(std::string_view s) {
    if (T v; std::from_chars(s.begin(), s.end(), v).ec == std::errc{}) {
        return v;
    } else {
        throw std::runtime_error("Can't convert " + std::string(s));
    }
}


template <std::size_t I, typename T>
struct indexed {
    using type = T;
};

template <typename Is, typename ...Ts>
struct indexer;

template <std::size_t ...Is, typename ...Ts>
struct indexer<std::index_sequence<Is...>, Ts...>
    : indexed<Is, Ts>...
{ };

template <std::size_t I, typename T>
static indexed<I, T> select(indexed<I, T>);


/**
 * @brief Get the nth element of a parameter pack of types
 *
 * @tparam I        index
 * @tparam Ts       Types parameter pack
 */
template <std::size_t I, typename ...Ts>
using nth_element = typename decltype(select<I>(
    indexer<std::index_sequence_for<Ts...>, Ts...>{}
))::type;


/**
 * @brief Redirect the cout sstream to another buffer
 * @details The output is reset once the object is destroyed
 *
 * @param  new_buff The new buffer
 */
struct cout_redirect {

    private:
        std::streambuf* old_buff;

    public:
        cout_redirect( std::streambuf* new_buff )
                : old_buff(std::cout.rdbuf(new_buff)){}

        ~cout_redirect( ) {
            std::cout.rdbuf(old_buff);
        }
};

class CSVRow
{
public:
    std::string_view operator[](std::size_t index) const
    {
        return std::string_view(&m_line[m_data[index] + 1], m_data[index + 1] -  (m_data[index] + 1));
    }
    std::size_t size() const
    {
        return m_data.size() - 1;
    }
    void readNextRow(std::istream& str)
    {
        std::getline(str, m_line);

        m_data.clear();
        m_data.emplace_back(-1);
        std::string::size_type pos = 0;
        while((pos = m_line.find(',', pos)) != std::string::npos)
        {
            m_data.emplace_back(pos);
            ++pos;
        }
        // This checks for a trailing comma with no data after it.
        pos   = m_line.size();
        m_data.emplace_back(pos);
    }
private:
    std::string         m_line;
    std::vector<int>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}

#endif // UTOPIA_MODELS_MULAN_MA_HH_HELPER