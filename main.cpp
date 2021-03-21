#include <iostream>
#include "main_viz.hpp"

#include "Bicubic.hpp"
#include "lagrange_solver.hpp"
#include "libman_solver.hpp"
#include "task_dirihlet_base.hpp"

struct print_type
    {
        template <class T>
        void operator() (T) const
        {
            auto const & ti = BOOST_CORE_TYPEID(T);
            std::cout << boost::core::demangled_name( ti ) << std::endl;
        }
    };

int main(int argc, char **argv)
{
//    viz(argc, argv);

//    mpl::for_each<type_base<2>::type>(print_type());

//    _LAG::LAGRANG_SET<NZ, sampleCountX> lag1(sampleMin, sampleMax, sampleMin, sampleMax, 0);
//    _LAG::LAGRANG_SET<NZ, sampleCountX> lag2(sampleMin, sampleMax, sampleMin, sampleMax, 0, 1000);

      float stepX = (sampleMax - sampleMin) / float(sampleCountX - 1);
      float stepZ = (sampleMax - sampleMin) / float(sampleCountZ - 1);

      std::cout << "OUT => " << std::endl;
      tbb::parallel_for( tbb::blocked_range2d<int>(0, sampleCountX, 0, sampleCountZ),
                         [stepX, stepZ](tbb::blocked_range2d<int> r)
      {
                _TPS::bicubic_solver solver(sampleMin, sampleMax, sampleMin, sampleMax);
          for ( int i = r.rows().begin(); i < r.rows().end(); ++i) {
              float z = qMin(sampleMax, (i * stepZ + sampleMin));
//              std::cout << lag1.get_x_orig(i) << std::endl;
                  for( int j = r.cols().begin(); j < r.cols().end(); j++) {
                      float x = qMin(sampleMax, (j * stepX + sampleMin));
                      std::cout << solver.bucubic(x, z) << std::endl;
                  }
          }
      });
//    double test[8][8];
//    _LIBMAN::zeroize(test);
//    std::cout << test[0][2] << std::endl;
      return 0;
}
