#if __cplusplus == 201103L

#include <iostream>

#include "mtk.h"
#include "mtk_dense_matrix.h"
#include "mtk_1d_grid_node.h"
#include "mtk_2d_uniform_cartesian_staggered_grid.h"
#include "mtk_roots.h"
int main () {

  std::cout << "Testing MTK 2D Grid" << std::endl;

  mtk::MTK_2DUniformStaggeredGrid *GRID = new mtk::MTK_2DUniformStaggeredGrid(0,1,12, 0,1,12);
  GRID->Print2DArray();
    
  return 0;

}
 

#else
#include <iostream>
int main () {
  std::cout << "This code HAS to be compiled to support C++11." << std::endl;
  std::cout << "On shell: sudo apt-get install build-essentials" << std::endl;
  std::cout << "Exiting..." << std::endl;
}
#endif

