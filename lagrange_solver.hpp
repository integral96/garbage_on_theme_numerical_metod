#pragma once

#include "Orts.hpp"
#include "lagrange_polynom.hpp"



namespace _LAG {

template<size_t N, size_t M>
using RESULT_LAG = RESULT_M<N, std::array<double, N>,
        std::array<double, M>>;

///Proto Lagrange

struct LagrangeGrammar : proto::or_<
    proto::terminal< proto::_ >,
    proto::plus< LagrangeGrammar, LagrangeGrammar>,
    proto::minus< LagrangeGrammar, LagrangeGrammar>
> {};

template<typename Expr> struct LagrangeExpr;
struct LagrangeDomain
    : proto::domain<proto::generator<LagrangeExpr>, LagrangeGrammar> {};

struct SubscriptCntxt : proto::callable_context<const SubscriptCntxt> {
        typedef std::pair<double, double> result_type;

        int index;
        SubscriptCntxt(int index_) :  index(index_) {}

        // polinom element
        template<typename Lagrange>
        std::pair<double, double> operator()(proto::tag::terminal, const Lagrange& pol) const {
            return pol.get_RES()[index];
        }
        template<typename E1, typename E2>
        std::pair<double, double> operator()(proto::tag::plus, const E1& e1, const E2& e2) const {
            return std::make_pair((proto::eval(e1, *this).first + proto::eval(e2, *this).first),
                                  (proto::eval(e1, *this).second + proto::eval(e2, *this).second));
        }

        template<typename E1, typename E2>
        std::pair<double, double> operator()(proto::tag::minus, const E1& e1, const E2& e2) const {
            return std::make_pair((proto::eval(e1, *this).first - proto::eval(e2, *this).first),
                                  (proto::eval(e1, *this).second - proto::eval(e2, *this).second));
        }
};

template<typename Expr>
struct LagrangeExpr : proto::extends<Expr, LagrangeExpr<Expr>, LagrangeDomain> {
        explicit LagrangeExpr(const Expr& e)
            : proto::extends<Expr, LagrangeExpr<Expr>, LagrangeDomain>(e) {
        }
        typename proto::result_of::eval< Expr, SubscriptCntxt>::type
        operator [](int i) const {
            const SubscriptCntxt ctx(i);
            return proto::eval(*this, ctx);
        }
};

template<size_t N, size_t M, typename T>
class LAGRANG_SET : public _TPS::coordinates_set<N, T>

{
    using value_type = double;
    using coordinates = _TPS::coordinates_set<N, T>;
    using type_array = std::array<value_type, M>;
    value_type shift{};
    std::array<value_type, M> x_orig{};
    std::array<value_type, M> y_orig{};

    value_type x_min{};
    value_type x_max{};
    value_type y_min{};
    value_type y_max{};
    value_type step_hM{};
    value_type step_hN{};

    constexpr void result() {
        RESULT_LAG<N, M> closure(x_orig, shift, y_orig, coordinates::get_X(), coordinates::get_Y());
        meta_loop<M>(closure);
    }
public:
    LAGRANG_SET (value_type x_min, value_type x_max, value_type y_min, value_type y_max, value_type shift, int n = 0) :
        coordinates(x_min, x_max, y_min, y_max, n),
        shift(shift), x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max) {
        step_hM = (x_max - x_min)/(value_type(M - 1));
        step_hN = (y_max - y_min)/(value_type(M - 1));
        for(size_t i = 0; i < M; ++i ) {
            x_orig[i] = x_min + step_hM*i;
        }
        result();
    }
    const value_type& get_step() const {
        return step_hN;
    }
    const value_type& get_x_min() const {
        return x_min;
    }
    const value_type& get_y_min() const {
        return y_min;
    }
    const auto& get_RES() const {
        std::array<std::pair<double, double>, M> tmp;
        for(size_t i  = 0; i < M; ++i ) {
            tmp[i] = std::make_pair(x_orig[i], y_orig[i]);
        }
        return tmp;
    }
    const auto& get_x_orig(size_t i) const {
        return x_orig[i];
    }
    const auto& get_y_orig(size_t j) const {
        return y_orig[j];
    }
    std::pair<double, double>& operator[](int i) { return this->get_RES()[i]; }
    const std::pair<double, double>& operator[](int i) const { return this->get_RES()[i]; }

    template<typename Expr>
    LAGRANG_SET& operator = ( const Expr& expr ) {
        for(int i = 0; i < M; ++i) {

            const SubscriptCntxt ctx(i);
            this->get_RES()[i] = proto::eval(proto::as_expr<LagrangeDomain>(expr), ctx);
        }
        return *this;
    }

};


template<typename>
struct IsLagrange : mpl::false_ {};
template<>
struct IsLagrange<LAGRANG_SET<15, 128, float50>> : mpl::true_  {};

namespace LagrangeOps {
    // This defines all the overloads to make expressions involving
    // Vector objects to build expression templates.
    BOOST_PROTO_DEFINE_OPERATORS(IsLagrange, LagrangeDomain)
}
}
