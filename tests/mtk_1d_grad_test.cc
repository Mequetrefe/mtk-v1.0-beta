#if __cplusplus == 201103L

#include <iostream>

#include "mtk.h"
#include "mtk_1d_grad.h"
#include "mtk_dense_matrix.h"


int main () {

  std::cout << "Testing MTK 1D Gradients" << std::endl;
  
  int Order_of_Accuracy=10;
  int NumCells=3*Order_of_Accuracy;

  
  mtk::MTK_Real Tau=1e-6    ;

  mtk::MTK_DenseMatrix *gg;

  mtk::MTK_Real *ww;
  mtk::MTK_Real *ss;
  mtk::MTK_Real *ee;
  mtk::MTK_1DGrad *grad = new mtk::MTK_1DGrad(Order_of_Accuracy,Tau);

  grad = grad->Construct1DGrad();
  
  if (grad == nullptr) {
    return -1;
  }


  gg=grad->ReturnAsMatrix(NumCells, 0.0,1.0);
    
  ss = grad->ReturnStencil();
  ww = grad->ReturnWeights();
  ee = grad->ReturnExtraRows(); 
  std::cout << *gg << std::endl;
  
  cout << "Weights" << endl;
  for(auto ii=0;ii<Order_of_Accuracy;ii++)
    cout << ww[ii]<< endl;
  
  cout << endl;
  
  cout << "Stencil"<<endl;
  for(auto ii=0;ii<Order_of_Accuracy;ii++)
    cout << ss[ii]<< endl;

  cout << endl;
  
  cout << "Mimetic Boundaries"<<endl;
  for(auto ii=0;ii<Order_of_Accuracy/2*3*Order_of_Accuracy/2;ii++)
    cout << ee[ii]<< endl;

 
  
  free(ss);
  free(ww);
  free(ee);
  free(gg);
  
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
