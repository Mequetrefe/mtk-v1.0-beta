#if __cplusplus == 201103L

#include <iostream>

#include "mtk.h"
#include "mtk_2d_lap.h"
#include "mtk_dense_matrix.h"


int main () {

  std::cout << "Testing MTK 2D lap " << std::endl;
  
  int Order_of_Accuracy=2;
  int NumCellsX=6;
  int NumCellsY=6;
//  int nn = NumCellsX+2;
//  int mm = NumCellsY+1;

  mtk::MTK_Real Tau=1e-6;

  mtk::MTK_DenseMatrix *ll2d;
 
//  mtk::MTK_Real *ll2dR;
  
  mtk::MTK_2DLap *lap = new mtk::MTK_2DLap(Order_of_Accuracy,Tau);

  lap = lap->Construct2DLap();
 
  
  if (lap == nullptr) {
    return -1;
  }


  ll2d=lap->ReturnAsMatrix(NumCellsX, 0.0,1.0, NumCellsY,0.0,1.0);
  cout << *ll2d << endl;

  

  free(ll2d);

  
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
