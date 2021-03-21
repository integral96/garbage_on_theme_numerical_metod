#pragma once


#include <boost/coroutine2/all.hpp>
#include <boost/range/algorithm.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>

#include <iostream>
#include <thread>
#include <future>
#include <functional>
#include <mutex>

#include <Eigen/Dense>
#include <Eigen/LU>

#include "Orts.hpp"
#include "Jacobian.hpp"
#include "task_dirihlet_base.hpp"


namespace _TPS {

template<typename T>
struct bicubic_solver : public coordinates_set<NZ, T>, public _DIRIHLET::SolverBase<bicubic_solver<T>> {
private:
    using value_type = T;
    using Orts = coordinates_set<NZ, T>;
    using ColVector4 = Eigen::Matrix<value_type, 4,  1>;
    using RowVector4 = Eigen::Matrix<value_type, 1,  4>;

    using matrix16   = Eigen::Matrix<value_type, 16, 16>;
    using vector16   = Eigen::Matrix<value_type, 16, 1>;

    using Vector16Array = boost::multi_array<vector16, 2>;
    Vector16Array RES_A;
    vector16 FUNC;

    matrix16 JACOBIAN;
    matrix16 INV_JCOB;

    const value_type x_min;
    const value_type x_max;
    const value_type y_min;
    const value_type y_max;
public:
    bicubic_solver(const value_type x_min, const value_type x_max, const value_type y_min, const value_type y_max) :
        x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max),
        Orts(x_min, x_max, y_min, y_max), RES_A(boost::extents[NZ][NZ])  {
    }
    bicubic_solver(const bicubic_solver& solver) : x_min(solver.x_min), x_max(solver.x_max),
        y_min(solver.y_min), y_max(solver.y_max),
        Orts(solver.x_min, solver.x_max, solver.y_min, solver.y_max), RES_A(solver.RES_A)  {
    }
    void print() const {
        std::cout << "FUNC = \n"     << FUNC << std::endl;
        std::cout << "JACOBIAN = \n" << JACOBIAN << std::endl;
        std::cout << "INV_JCOB = \n" << INV_JCOB << std::endl;
//        for ( int i = 0; i < NZ - 1; ++i)
//            for( int j = 0; j < NZ - 1; j++)
//                std::cout << RES_A[i][j] << " ";
//        std::cout << std::endl;
    }

    int get_x_index_to_left_of(value_type x) const {
      auto xrng = std::make_pair(Orts::get_X().data(), Orts::get_X().data() + Orts::get_X().size());
      return boost::lower_bound( xrng, x) - boost::begin(xrng) - 1;
    }

    int get_x_index_to_right_of(value_type x) const {
      return this->get_x_index_to_right_of(x) + 1;
    }


    int get_y_index_below(value_type y) const {
      auto yrng = std::make_pair(Orts::get_Y().data(), Orts::get_Y().data() + Orts::get_Y().size());
      return boost::lower_bound( yrng, y) - boost::begin(yrng) - 1;
    }

    int get_y_index_above(value_type y) const
    {
      return this->get_y_index_below(y) + 1;
    }

    auto bucubic(const value_type x, const value_type y) {

        int i = this->get_x_index_to_left_of(x);
        int j = this->get_y_index_below(y);

        if(i < 0)
          i = 0;
        if(j < 0)
          j = 0;

        init_ZNACH();

        value_type xL = Orts::get_X()[i + 1] - Orts::get_X()[i];
        value_type yL = Orts::get_Y()[j + 1] - Orts::get_Y()[j];


        RowVector4 vx;
        vx[0] = 1;
        vx[1] = (x - Orts::get_X()[i])/xL;
        vx[2] = vx[1] * vx[1];
        vx[3] = vx[2] * vx[1];

        ColVector4 vy;
        vy[0] = 1;
        vy[1] = (y - Orts::get_Y()[j])/yL;
        vy[2] = vy[1] * vy[1];
        vy[3] = vy[2] * vy[1];

        value_type tmp{};
        for(size_t k1 = 0; k1 < 4; ++k1) {
            value_type tmp1{};
            for(size_t k2 = 0; k2 < 4; ++k2) {
                tmp1 += RES_A[i][j](k1*4 + k2)*vx[k1]*vy[k2];
            }
            tmp += tmp1;
        }
        if(std::abs(tmp - func_(Orts::get_X()[i], Orts::get_Y()[j])) < 20.)
            return tmp;

    }
    inline void init_ZNACH() {
        std::mutex mutex;
        tbb::parallel_for( tbb::blocked_range2d<int>(0, NZ - 1, 0, NZ - 1),
                           [this, &mutex](tbb::blocked_range2d<int> r)
        {
            for ( int i = r.rows().begin(); i < r.rows().end(); ++i) {
                    for( int j = r.cols().begin(); j < r.cols().end(); j++) {
                        mutex.lock();

                        value_type xL = Orts::get_X()[i + 1] - Orts::get_X()[i];
                        value_type yL = Orts::get_Y()[j + 1] - Orts::get_Y()[j];

                        auto dx0 = (Orts::get_X()[std::min(i + 1, NZ - 1)] - Orts::get_X()[std::max(i - 1, 0)])/xL;
                        auto dx1 = (Orts::get_X()[std::min(i + 2, NZ - 1)] - Orts::get_X()[std::max(i, 0)])/xL;

                        auto dy0 = (Orts::get_Y()[std::min(j + 1, NZ - 1)] - Orts::get_Y()[std::max(j - 1,0)])/yL;
                        auto dy1 = (Orts::get_Y()[std::min(j + 2, NZ - 1)] - Orts::get_Y()[std::max(j, 0)])/yL;

                        auto X0 = Orts::get_X()[i];
                        auto X1 = Orts::get_X()[i + 1];
                        auto Y0 = Orts::get_Y()[j];
                        auto Y1 = Orts::get_Y()[j + 1];
                        mutex.unlock();

                        FUNC << DERIV_FUNC_Z<0, 0>(X0, Y0), DERIV_FUNC_Z<0, 0>(X1, Y0),
                                DERIV_FUNC_Z<0, 0>(X0, Y1), DERIV_FUNC_Z<0, 0>(X1, Y1),

                                dx0*DERIV_FUNC_Z<1, 0>(X0, Y0), dx1*DERIV_FUNC_Z<1, 0>(X1, Y0),
                                dx0*DERIV_FUNC_Z<1, 0>(X0, Y1), dx1*DERIV_FUNC_Z<1, 0>(X1, Y1),

                                dy0*DERIV_FUNC_Z<0, 1>(X0, Y0), dy1*DERIV_FUNC_Z<0, 1>(X1, Y0),
                                dy0*DERIV_FUNC_Z<0, 1>(X0, Y1), dy1*DERIV_FUNC_Z<0, 1>(X1, Y1),

                                dx0*dy0*DERIV_FUNC_Z<1, 1>(X0, Y0), dx1*dy1*DERIV_FUNC_Z<1, 1>(X1, Y0),
                                dx0*dy0*DERIV_FUNC_Z<1, 1>(X0, Y1), dx1*dy1*DERIV_FUNC_Z<1, 1>(X1, Y1);


                        INIT_JACOBIAN(X0, X1, Y0, Y1, JACOBIAN);
                        mutex.lock();
                        INV_JCOB =  JACOBIAN.inverse();
                        mutex.unlock();

                        RES_A[i][j] = INV_JCOB*FUNC;

                    }
                }

        });

    }

};
}
