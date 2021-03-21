#pragma once

#include "base_func.hpp"

#include <memory>
#include <vector>
#include <utility>
#include <Eigen/Dense>
#include <boost/range/algorithm.hpp>
#include <boost/range/adaptor/strided.hpp>
#include <boost/fusion/container.hpp>
#include <boost/fusion/tuple.hpp>
#include <boost/fusion/include/pair.hpp>
#include <boost/fusion/include/for_each.hpp>

#include "Jacobian.hpp"

namespace fusion = boost::fusion;

namespace _DIRIHLET {

template<typename T>
struct float50TypeOf { using type = float50; };
template<template<typename> class T,typename R>
struct float50TypeOf<T<R>> { using type = R; };


template<class Derived, typename T = typename float50TypeOf<Derived>::type>
class SolverBase
{
  public:
    using VectorType = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    using map_type = fusion::map<std::pair<const VectorType, char>,
                                 std::pair<const VectorType, char>,
                                 std::pair<const VectorType, char>>;

    using MatrixType = typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorView = fusion::set<fusion::vector<const VectorType, mpl::int_<0>, Eigen::InnerStride<Eigen::Dynamic>>>;
    using MatrixView = fusion::set<fusion::vector<const MatrixType, mpl::int_<0>, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>>;

  protected:

    std::vector<T> xData, yData, zData;
    std::unique_ptr<map_type> maView;

    std::unique_ptr<VectorView> X,Y;
    std::unique_ptr<MatrixView> Z;

  private:

    friend Derived;

    SolverBase(){ }


  private:

    template<typename F>
    struct has_Solver
    {
      private:
        using yes = mpl::true_;
        using no  = mpl::false_;


        template<typename U> static auto test(int) -> decltype(std::declval<U>().init_ZNACH(), yes());
        template<typename  > static no   test(...);

      public:
        static constexpr bool value = std::is_same_v<decltype(test<F>(0)), yes>;
    };

    template<typename F>
    typename std::enable_if_t<has_Solver<F>::value>
    callSetupSolver() {
        static_cast<Derived*>(this)->setupSolver();
    }

    template<typename F>
    typename std::enable_if_t<!has_Solver<F>::value>
    callSetupSolver(){ }


};

}
