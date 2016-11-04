#if __cplusplus == 201103L

#include <iostream>

#include "mtk.h"
#include "mtk_1d_div.h"
#include "mtk_dense_matrix.h"


int main () {

  std::cout << "Testing MTK 1D Divergences" << std::endl;
  
  int order = 10;

  mtk::MTK_1DDiv *div = new  mtk::MTK_1DDiv(order,1e-6);

  div = div->Construct1DDiv();

  if (div == nullptr) {
    return -1;
  }

  mtk::MTK_DenseMatrix* dd;
  mtk::MTK_Real *ss;


  dd = div->ReturnAsMatrix(24,0.0,1.0);
    
  ss = div->ReturnStencil();

  

  std::cout << *dd << std::endl;
  
if(order >2)
{    
  mtk::MTK_Real *ww;
  ww = div->ReturnWeights();
  cout << "Weights" << endl;
  for(auto ii=0;ii<order;ii++)
    cout << ww[ii]<< endl;
free (ww);
}  
  cout<<endl;
  
  cout << "Stencil"<<endl;
   for(auto ii=0;ii<order;ii++)
    cout << ss[ii]<< endl;
free (ss);
  cout<<endl;
  

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
