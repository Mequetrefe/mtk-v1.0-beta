#if __cplusplus == 201103L

#include <iostream>

#include "mtk.h"
#include "mtk_2d_lap.h"
#include "mtk_dense_matrix.h"
#include "mtk_2d_bc_desc.h"

mtk::MTK_Real f1(mtk::MTK_Real xx, mtk::MTK_Real yy)
{
    return 1.0;
}

mtk::MTK_Real f2(mtk::MTK_Real xx, mtk::MTK_Real yy)
{
    return 4.0;
}

mtk::MTK_Real f3(mtk::MTK_Real xx, mtk::MTK_Real yy)
{
    return 2.0;
}

mtk::MTK_Real f4(mtk::MTK_Real xx, mtk::MTK_Real yy)
{
    return 3.0;
}

int main () {

  std::cout << "Testing MTK 2D BC descriptor" << std::endl;
  
  int Order_of_Accuracy=2;
  int NumCellsX=6;
  int NumCellsY=14;
  mtk::MTK_DenseMatrix *ll2d;
  mtk::MTK_2DBCDesc *BC_Desc;
  mtk::MTK_Real RHS[(NumCellsX+2)*(NumCellsY+2)];

  std::cout << "Creating a MTK 2D Grid" << std::endl;


  mtk::MTK_Real Tau=1e-6;

  mtk::MTK_2DLap *lap = new mtk::MTK_2DLap(Order_of_Accuracy,Tau, NumCellsX, 0.0,1.0, NumCellsY, 0.0,1.0);
  
  mtk::MTK_2DUniformStaggeredGrid *GRID = new mtk::MTK_2DUniformStaggeredGrid(0,1,NumCellsX, 0,1,NumCellsY);


  lap = lap->Construct2DLap();
   
  if (lap == nullptr) {
    return -1;
  }

  ll2d=lap->ReturnAsMatrix(NumCellsX, 0.0,1.0, NumCellsY, 0.0,1.0);
  
  BC_Desc = new mtk::MTK_2DBCDesc(lap, GRID,f1, f2, f3, f4);

  BC_Desc->MTK_2DBCDirichlet(ll2d);     

    //The Laplacian Matrix will be returned with BC applied 
  std::cout << "Back in main() " << std::endl;
  std::cout << *ll2d << std::endl;

  BC_Desc->MTK_2DBoundaryRHS(RHS);     
  for (auto ii=0; ii<(NumCellsX+2)*(NumCellsY+2);ii++)
  {
      std::cout << RHS[ii]<<" " ;
  }
  std::cout<<std::endl;
 
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
