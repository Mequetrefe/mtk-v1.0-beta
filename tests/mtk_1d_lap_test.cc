#if __cplusplus == 201103L

#include <iostream>

#include "mtk.h"
#include "mtk_1d_lap.h"
#include "mtk_dense_matrix.h"


int main () {

  std::cout << "Testing MTK 1D Laplacian" << std::endl;
  
  int Order_of_Accuracy=2;
  int NumCells=3*Order_of_Accuracy;
  
  
  
  mtk::MTK_Real Tau=1e-6;
  
  mtk::MTK_DenseMatrix *ll;
  mtk::MTK_1DLap *lap = new mtk::MTK_1DLap(Order_of_Accuracy,Tau);

  lap = lap->Construct1DLap();
  
  if (lap == nullptr) {
    return -1;
  }
  
  
  
  ll=lap->ReturnAsMatrix(NumCells, 0.0,1.0);


  
  cout << "Laplacian as a dense Matrix"<<endl;
  std::cout << *ll<<std::endl;


  
  
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
