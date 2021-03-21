#pragma once

#include "base_func.hpp"


namespace _LIBMAN {

template<typename T, T VALUE>
struct static_parametr {

};
template<typename T, T VALUE>
struct static_value : static_parametr<T, VALUE> {
    static const T value = VALUE;
    operator T () const {
        return VALUE;
    }
    static_value(int = 0) {}
};
template <int N, int M, int L>
struct index {};

template <typename T, int N, int M, int L>
void zeroize_helper(T** const data, index<N, M, L>)
{
    zeroize_helper(data, index<N/2, M, L>());
    zeroize_helper(data, index<N/2, M + N/2, L>());
    zeroize_helper(data, index<N/2, M + N/2, M + N + L/2>());
}
template <typename T, int M, int L>
void zeroize_helper(T** const data, index<1, M, L>) {
    zeroize_helper(data, index<0, M/2, L>());
    zeroize_helper(data, index<0, M/2, M + L/2>());
}
template <typename T, int M>
void zeroize_helper(T** const data, index<1, M, 0>) {
    data[0][M] = T(1);
}

template <typename T, int N, int M>
void zeroize(T (&data)[N][M])
{
    zeroize_helper(data, index<N, M, 0>());
}


}
