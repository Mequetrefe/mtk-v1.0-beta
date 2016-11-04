#if __cplusplus == 201103L

#include <iostream>

#include "mtk.h"
#include "mtk_2d_grad.h"
#include "mtk_dense_matrix.h"


int main () {

  std::cout << "Testing MTK 2D Grad" << std::endl;
  
  int Order_of_Accuracy=10;
  int NumCellsX=Order_of_Accuracy*3;
  int NumCellsY=Order_of_Accuracy*3;

  mtk::MTK_Real Tau=1e-6;

  mtk::MTK_DenseMatrix *gg2d;
  mtk::MTK_Real *gg2dR;
  
  mtk::MTK_2DGrad *grad = new mtk::MTK_2DGrad(Order_of_Accuracy,Tau);

  grad = grad->Construct2DGrad();
  
  
  if (grad == nullptr) {
    return -1;
  }


  gg2d=grad->ReturnAsMatrix(NumCellsX, 0.0, 1.0, NumCellsY, 0.0, 1.0);
   std::cout << "Grad 2D as MTK_DenseMatrix" << std::endl;

   std::cout << *gg2d << std::endl;
   gg2dR=grad->ReturnAsReal(NumCellsX, 0.0, 1.0,NumCellsY, 0.0, 1.0);
   
   
 

  free(gg2d);
  free(gg2dR);
  
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
