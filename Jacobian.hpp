#pragma once

#include <boost/math/differentiation/autodiff.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/adaptor/strided.hpp>


#include "base_func.hpp"


namespace diff = boost::math::differentiation;



template <typename X, typename Y>
struct diff_polynom {
    template<size_t J>
    struct iner {
        iner(const diff_polynom& thiss) : thiss(thiss) {}
        template<size_t I>
        diff::promote<X, Y> func() const {
          return my::pow<I>(thiss.x_)*my::pow<J>(thiss.y_);
        }
    private:
        const diff_polynom& thiss;
    };
    diff_polynom(const X& x, const Y& y) : x_(x), y_(y) {}
private:
    const X& x_;
    const Y& y_;
};

///Func for polynom
template<size_t I, size_t J, typename X, typename Y>
diff::promote<X, Y> func(const X& x, const Y& y) {
    return my::pow<I>(x)*my::pow<J>(y);
}
///Func for Z function
template<typename X, typename Y>
diff::promote<X, Y> func_Z(const X& x, const Y& y) {
//    return y + .005*y*y + .005*y*y*y +
//            .02*x + .008*x*y + .005*x*y*y + .0001*x*y*y*y +
//           .005*x*x + .0001*x*x*y + .00002*x*x*y*y + .00003*x*x*y*y*y +
//           .00005*x*x*x + .000003*x*x*x*y + .0000008*x*x*x*y*y + .0000005*x*x*x*y*y*y;

    return my::pow<2>(x) + my::pow<2>(y);
}


///Func for polynom
template <size_t I, size_t J, size_t Nx, size_t Ny>
constexpr auto DERIV_(double x_, double y_) {
        auto const variables = diff::make_ftuple<double, Nx, Ny>(x_, y_);
        auto const& x = std::get<0>(variables);  // Up to Nw derivatives at w=12
        auto const& y = std::get<1>(variables);
        auto const result = func<I, J>(x, y);
        return result.derivative(Nx, Ny);
    }

///Func for Z function
template <size_t Nx, size_t Ny>
constexpr auto DERIV_FUNC_Z(double x_, double y_) {
        auto const variables = diff::make_ftuple<double, Nx, Ny>(x_, y_);
        auto const& x = std::get<0>(variables);  // Up to Nw derivatives at w=12
        auto const& y = std::get<1>(variables);
        auto const result = func_Z(x, y);
        return result.derivative(Nx, Ny);
    }


//////JAKOBIAN
///
template<size_t PX, size_t PY, size_t DX, size_t DY>
struct base_ {
    static constexpr size_t px = PX;
    static constexpr size_t py = PY;
    static constexpr size_t dx = DX;
    static constexpr size_t dy = DY;
};

template<size_t I>
using type_base  =
typename mpl::if_c<(I >= 0 && I < 4),
typename mpl::vector< base_<I, 0, 0, 0>, base_<I, 0, 1, 0>, base_<I, 0, 0, 1>, base_<I, 0, 1, 1> >::type,
typename mpl::if_c<(I >= 4 && I < 8),
typename mpl::vector< base_<(I - 4),  1, 0, 0>, base_<(I - 4), 1, 1, 0>, base_<(I - 4), 1, 0, 1>, base_<(I - 4), 1, 1, 1> >::type,
typename mpl::if_c<(I >= 8 && I < 12),
typename mpl::vector< base_<(I - 8),  2, 0, 0>, base_<(I - 8), 2, 1, 0>, base_<(I - 8), 2, 0, 1>, base_<(I - 8), 2, 1, 1> >::type,
typename mpl::vector< base_<(I - 12), 3, 0, 0>, base_<(I - 12), 3, 1, 0>, base_<(I - 12), 3, 0, 1>, base_<(I - 12), 3, 1, 1> >::type
                >::type>::type>::type;

/////////////////////
template<size_t I, typename Matrix>
struct JACOBIAN_I {
    using type_I = typename type_base<I>::type;

    using value_type = typename Matrix::value_type;
    JACOBIAN_I(const value_type x0, const value_type x1, const value_type y0, const value_type y1, Matrix& JACOB) :
        x0(x0), x1(x1), y0(y0), y1(y1), JACOB(JACOB) {
    }
    template<size_t J>
    void apply() const {
        using type_J = typename mpl::at_c<type_I, J/4>::type;
        if constexpr((J == 0) || (J == 4) || (J == 8) || (J == 12)) {
            JACOB(I, J) = DERIV_<type_J::px, type_J::py, type_J::dx, type_J::dy>(x0, y0);
        }
        if constexpr((J == 1) || (J == 5) || (J == 9) || (J == 13)) {
            JACOB(I, J) = DERIV_<type_J::px, type_J::py, type_J::dx, type_J::dy>(x1, y0);
        }
        if constexpr((J == 2) || (J == 6) || (J == 10) || (J == 14)) {
            JACOB(I, J) = DERIV_<type_J::px, type_J::py, type_J::dx, type_J::dy>(x0, y1);
        }
        if constexpr((J == 3) || (J == 7) || (J == 11) || (J == 15)) {
            JACOB(I, J) = DERIV_<type_J::px, type_J::py, type_J::dx, type_J::dy>(x1, y1);
        }
    }
private:
    const value_type x0;
    const value_type x1;
    const value_type y0;
    const value_type y1;
    Matrix& JACOB;
};
template<typename Matrix>
struct JACOBIAN_J {
    using value_type = typename Matrix::value_type;
    JACOBIAN_J(const value_type x0, const value_type x1, const value_type y0, const value_type y1, Matrix& JACOB) :
        x0(x0), x1(x1), y0(y0), y1(y1), JACOB(JACOB) {
    }
    template<size_t I>
    void apply() const {
        JACOBIAN_I<I, Matrix> closure(x0, x1, y0, y1, JACOB);
        meta_loop<16>(closure);
    }
private:
    const value_type x0;
    const value_type x1;
    const value_type y0;
    const value_type y1;
    Matrix& JACOB;
};
template<class Matrix>
inline void INIT_JACOBIAN(typename Matrix::value_type x0, typename Matrix::value_type x1,
                          typename Matrix::value_type y0, typename Matrix::value_type y1, Matrix &matrix) {
    JACOBIAN_J<Matrix> closure(x0, x1, y0, y1, matrix);
    meta_loop<16>(closure);
}
//////JAKOBIAN FINISH
