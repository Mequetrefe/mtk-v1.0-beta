#if __cplusplus == 201103L

#include <iostream>

#include "mtk.h"
#include "mtk_2d_div.h"
#include "mtk_dense_matrix.h"


int main () {

  std::cout << "Testing MTK 2D Divergence" << std::endl;
  
  int Order_of_Accuracy=8;
  int NumCellsX=Order_of_Accuracy*3;
  int NumCellsY=Order_of_Accuracy*3;

  
  mtk::MTK_Real Tau=1e-6;
  mtk::MTK_Real *dd2dR;
  mtk::MTK_DenseMatrix *dd2d;
  
  mtk::MTK_2DDiv *div = new mtk::MTK_2DDiv(Order_of_Accuracy,Tau);

  div = div->Construct2DDiv();
  
  if (div == nullptr) {
    return -1;
  }


  dd2d=div->ReturnAsMatrix(NumCellsX, 0.0,1.0, NumCellsY, 0.0, 1.0);
  std::cout << "Div 2D as MTK_DenseMatrix" << std::endl;
  std::cout << *dd2d << std::endl;
  
  dd2dR=div->ReturnAsReal(NumCellsX, 0.0, 1.0, NumCellsY, 0.0, 1.0);

  std::cout << "Div 2D as double (MTK_Real)" << std::endl;

  for(auto ii=0; ii<((NumCellsX+1)*NumCellsY+(NumCellsY+1)*NumCellsX);ii++)
  {
      for(auto jj=0;jj<(NumCellsX+2)*(NumCellsY+2); jj++)
      {
          std::cout <<dd2dR[ii*(NumCellsX+2)*(NumCellsY+2)+jj]<<" ";
  
      }
      std::cout<<std::endl;
  }
  free (dd2dR);
  free (dd2d);
  
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
