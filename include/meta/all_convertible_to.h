// This is a recursive template that checks whether all types are convertible to T

#if __cplusplus < 201703L
// Disable this for C++17 or later, use fold expressions
#include <type_traits>

template<typename T, typename... Rest>
struct all_convertible_to;

template<typename T>
struct all_convertible_to<T> : std::true_type {};

template<typename T, typename First, typename... Rest>
struct all_convertible_to<T, First, Rest...>
    : std::conditional<
        std::is_convertible<First, T>::value,
        all_convertible_to<T, Rest...>,
        std::false_type
    >::type
{};

#endif
