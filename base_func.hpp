#pragma once

#include <boost/mpl/bool.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/proto/proto.hpp>
#include <boost/multi_array.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include <type_traits>
#include <array>
#include <complex>
#include <random>
#include <ctime>
#include <chrono>
#include <thread>


namespace mpl = boost::mpl;
namespace proto = boost::proto;

using float50 = boost::multiprecision::cpp_bin_float_50;

static constexpr int NZ {36}; //Кол наблюдаемых точек

static constexpr int step = 4;
static constexpr int freq = 4;

static constexpr int length = step * freq;
static const double edg_cube = .3;

static constexpr int sampleCountX = 64;
static constexpr int sampleCountZ = 64;
static constexpr int heightMapGridStepX = 4;
static constexpr int heightMapGridStepZ = 4;
static const float sampleMin = -16.0f;
static const float sampleMax = 16.0f;

///Factorial

template<size_t N>
struct factorial {
    static constexpr size_t value = N*factorial<N - 1>::value;
};
template<>
struct factorial<0> {
    static constexpr size_t value = 1;
};

///
///Вычисление степени
///
namespace my {
template<class T, int N>
    struct helper
    {
        static constexpr T pow(const T x){
            return helper<T, N-1>::pow(x) * x;
        }
    };

    template<class T>
    struct helper<T, 1>
    {
        static constexpr T pow(const T x){
            return x;
        }
    };

    template<class T>
    struct helper<T, 0>
    {
        static constexpr T pow(const T x){
            return T(1);
        }
    };
    template<int N, class T>
    T constexpr pow(T const x)
    {
        return helper<T, N>::pow(x);
    }
}

// is vector for multiply vfnrbx on vector
template <typename F, size_t N, typename Vector>
struct is_vector : boost::mpl::false_ {};
template <typename F, size_t N>
struct is_vector<F, N, boost::array<F, N>> : boost::mpl::true_ {
    static constexpr size_t size = N;
};
template <typename F, size_t N>
struct is_vector<F, N, std::array<F, N>> : boost::mpl::true_ {
    static constexpr size_t size = N;
};


/*!
 * struct meta_loop
 */
template <size_t N, size_t I, class Closure>
typename std::enable_if_t<(I == N)> is_meta_loop(Closure&) {}

template <size_t N, size_t I, class Closure>
typename std::enable_if_t<(I < N)> is_meta_loop(Closure& closure) {
    closure.template apply<I>();
    is_meta_loop<N, I + 1>(closure);
}
template <size_t N, class Closure>
void meta_loop(Closure& closure) {
    is_meta_loop<N, 0>(closure);
}
template <size_t N, class Closure>
void meta_loopUV(Closure& closure) {
    is_meta_loop<N, 1>(closure);
}
template <size_t N, size_t K, class Closure>
void meta_loop_KN(Closure& closure) {
    is_meta_loop<N, K>(closure);
}
///++
///
/*!
 * struct meta_loop_inv
 */
template <int N, int I, class Closure>
typename std::enable_if_t<(I < 0)> is_meta_loop_inv(Closure&) {}

template <int N, int I, class Closure>
typename std::enable_if_t<(I >= 0)> is_meta_loop_inv(Closure& closure) {
    closure.template apply<I>();
    is_meta_loop_inv<0, I - 1>(closure);
}
template <int N, class Closure>
void meta_loop_inv(Closure& closure) {
    is_meta_loop_inv<0, N>(closure);
}

///++

/////Calculate Binom

template<size_t N, size_t K>
struct BC {
    static constexpr size_t value = factorial<N>::value / factorial<K>::value / factorial<N - K>::value;
};
/*!
 * struct abstract_sum
 */
template<class Closure>
struct abstract_sum_closures {
    typedef typename Closure::value_type value_type;
    abstract_sum_closures(Closure &closure) :  closure(closure), result(value_type(0)){}

    template<unsigned I>
    void apply(){
        result += closure.template value<I>();
    }
    Closure &closure;
    value_type result;
};

template<size_t N, class Closure>
typename Closure::value_type abstract_sums(Closure &closure) {
    abstract_sum_closures<Closure> my_closure(closure);
    meta_loop<N>(my_closure);
    return my_closure.result;
}

/*!
 * struct abstract_subtract
 */
template<class Closure>
struct abstract_subtract_closures {
    typedef typename Closure::value_type value_type;
    abstract_subtract_closures(Closure &closure) :  closure(closure), result(value_type()){}

    template<unsigned I>
    void apply(){
        result -= closure.template value<I>();
    }
    Closure &closure;
    value_type result;
};

template<unsigned N, class Closure>
typename Closure::value_type abstract_subtract(Closure &closure) {
    abstract_subtract_closures<Closure> my_closure(closure);
    meta_loop<N>(my_closure);
    return my_closure.result;
}
/*!
 * struct abstract_mult
 */
template<class Closure>
struct abstract_multiple_closures {
    using value_type = typename Closure::value_type;
    abstract_multiple_closures(Closure &closure) : closure(closure), result(value_type(1)){}
    template<size_t I>
    void apply(){
        result *= closure.template value<I>();
    }
    Closure &closure;
    value_type result;
};
template<size_t K, class Closure>
typename Closure::value_type abstract_multiple(Closure &closure) {
    abstract_multiple_closures<Closure> my_closure(closure);
    meta_loop<K>(my_closure);
    return my_closure.result;
}

/*!
 * struct abstract_divide
 */
template<class Closure>
struct abstract_divide_closures {
    typedef typename Closure::value_type value_type;
    abstract_divide_closures(Closure &closure) :  closure(closure), result(value_type(1)){}

    template<unsigned I>
    void apply(){
        result /= closure.template value<I>();
    }
    Closure &closure;
    value_type result;
};

template<unsigned N, class Closure>
typename Closure::value_type abstract_divide(Closure &closure) {
    abstract_subtract_closures<Closure> my_closure(closure);
    meta_loop<N>(my_closure);
    return my_closure.result;
}
///meta func
///
// is vector for multiply vfnrbx on vector
template <size_t N, typename Array>
struct is_arrayd : boost::mpl::false_ {};
template <size_t N>
struct is_arrayd<N, std::array<double, N>> : boost::mpl::true_ {
    static constexpr size_t size = N;
};

template<typename T> struct is_double : boost::mpl::false_ {};
template<> struct is_double<double> : boost::mpl::true_ {};

template<typename T> struct is_complexd : boost::mpl::false_ {};
template<> struct is_complexd<std::complex<double>> : boost::mpl::true_ {};

///Gen random lon lat
template<size_t N, typename T, typename Array, typename = std::enable_if_t<std::is_same_v<Array, std::array<T, N>>>>
inline void gen_rand_array(Array& A, T min, T max, int n = 0) {
    auto start = std::chrono::system_clock::now();
    std::this_thread::sleep_for(std::chrono::milliseconds(n));
    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::mt19937 gen{static_cast<std::uint32_t>(end_time)};
    std::uniform_real_distribution<> dist{min, max};
    for(size_t i = 0; i < N; ++i)
            A[i] = dist(gen);
}
