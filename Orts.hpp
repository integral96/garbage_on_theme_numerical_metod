#pragma once

#include <boost/noncopyable.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered_set.hpp>

#include <Eigen/Dense>

#include "base_func.hpp"

struct XYZ {
    float x;
    float y;
    float z;
    XYZ(float x, float y, float z) : x(x), y(y), z(z) {}
    inline bool operator == (const XYZ& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};
struct Hasher
{
  inline std::size_t operator()(const XYZ& k) const {
      size_t seed{};
      boost::hash_combine(seed, boost::hash_value(k.x));
      boost::hash_combine(seed, boost::hash_value(k.y));
      boost::hash_combine(seed, boost::hash_value(k.z));
      return seed;
  }
};

///Function
static auto func_([](auto x, auto y) { return my::pow<2>(x) + my::pow<2>(y); });
///Realisaciia
///
namespace _TPS {

template<size_t N, typename T>
class coordinates_set  : boost::noncopyable {
public:
    using value_type = T;
    using type_array = std::array<value_type, N>;
    using type_matrix = Eigen::Matrix<value_type, N, N>;
    using type_pair  = std::array<std::pair<value_type, value_type>, N>;
    using type_XYZ  = boost::unordered_set<XYZ, Hasher>;
private:
    const value_type X_range_min, X_range_max, Y_range_min, Y_range_max;
    type_array X;
    type_array Y;
    type_matrix Z;
    type_pair XY;
    type_XYZ  XYZ_;
    int n{};
    void init_X() {
        gen_rand_array<N>(X, X_range_min, X_range_max, n);
        std::sort(X.begin(), X.end());
//        for(size_t i = 0; i < N; ++i ) {
//            X[i] = (X_range_min + X_range_max)/2 +
//                    (X_range_max - X_range_min)/2*std::cos((2*i + 1)*M_PI/(2*(N + 1)));
//        }
    }
    void init_Y() {
        gen_rand_array<N>(Y, Y_range_min, Y_range_max, n);
    }
    void init_Z() {
        for(size_t i = 0; i < N; ++i)
            for(size_t j = 0; j < N; ++j)
                Z(i, j) = func_(get_X()[i], get_Y()[j]);
    }
public:
    explicit coordinates_set(const value_type x_min, const value_type x_max, const value_type y_min, const value_type y_max, int n = 0) :
                                X_range_min(x_min), X_range_max(x_max), Y_range_min(y_min), Y_range_max(y_max), n(n) {
        init_X();
        init_Y();
        init_Z();
    }
public:

    //
    const type_array& get_X() const {
        return X;
    }
    const type_array& get_Y() const {
        return Y;
    }
    const type_matrix& get_Z() const {
        return Z;
    }

    //print
    void print() {
        std::cout << "X = {";
        for(const auto& x : X) {
            std::cout << x << ", ";
        }
        std::cout << "}" << std::endl;
        std::cout << "Y = {";
        for(const auto& y : Y) {
            std::cout << y << ", ";
        }
        std::cout << "}" << std::endl;
        std::cout << "Z = {";
            std::cout << "(" << Z << "), ";
        std::cout << "}" << std::endl;
    }
};
}
