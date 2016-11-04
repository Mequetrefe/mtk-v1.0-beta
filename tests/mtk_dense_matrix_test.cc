#include <iostream>

#include "mtk.h"
#include "mtk_dense_matrix.h"

int main () {
  
  std::cout << "Testing MTK_DenseMatrix class." << std::endl;
  
  mtk::MTK_DenseMatrix *xx;
  mtk::MTK_DenseMatrix *yy;
  mtk::MTK_DenseMatrix *zz;
  mtk::MTK_DenseMatrix *pp;
  int rr1,rr2,cc1,cc2;
  mtk::MTK_Real* rr;
  
  rr1=6;
  cc1=5;
  
  rr2=5;
  cc2=6;
  
  xx=new mtk::MTK_DenseMatrix(rr1,cc1);
  yy=new mtk::MTK_DenseMatrix(rr2,cc2);
  zz=new mtk::MTK_DenseMatrix(rr1*rr2,cc1*cc2);
  pp=new mtk::MTK_DenseMatrix(rr1,cc2);
  rr=new mtk::MTK_Real[rr1*cc1];
   
  xx->SetValue(0,0,1);
  xx->SetValue(0,1,2);
  xx->SetValue(0,2,3);
  xx->SetValue(0,3,4);
  xx->SetValue(0,4,5);
  xx->SetValue(1,0,6);
  xx->SetValue(1,1,7);
  xx->SetValue(1,2,8);
  xx->SetValue(1,3,9);
  xx->SetValue(1,4,10);
  xx->SetValue(2,0,11);
  xx->SetValue(2,1,12);
  xx->SetValue(2,2,13);
  xx->SetValue(2,3,14);
  xx->SetValue(2,4,15);
  xx->SetValue(3,0,16);
  xx->SetValue(3,1,17);
  xx->SetValue(3,2,18);
  xx->SetValue(3,3,19);
  xx->SetValue(3,4,20);
  xx->SetValue(4,0,21);
  xx->SetValue(4,1,22);
  xx->SetValue(4,2,23);
  xx->SetValue(4,3,24);
  xx->SetValue(4,4,25);
  xx->SetValue(5,0,26);
  xx->SetValue(5,1,27);
  xx->SetValue(5,2,28);
  xx->SetValue(5,3,29);
  xx->SetValue(5,4,30);
  

  yy->SetValue(0,0,1);
  yy->SetValue(0,1,2);
  yy->SetValue(0,2,3);
  yy->SetValue(0,3,4);
  yy->SetValue(0,4,5);
  yy->SetValue(0,5,6);
  yy->SetValue(1,0,7);
  yy->SetValue(1,1,8);
  yy->SetValue(1,2,9);
  yy->SetValue(1,3,10);
  yy->SetValue(1,4,11);
  yy->SetValue(1,5,12);
  yy->SetValue(2,0,13);
  yy->SetValue(2,1,14);
  yy->SetValue(2,2,15);
  yy->SetValue(2,3,16);
  yy->SetValue(2,4,17);
  yy->SetValue(2,5,18);
  yy->SetValue(3,0,19);
  yy->SetValue(3,1,20);
  yy->SetValue(3,2,21);
  yy->SetValue(3,3,22);
  yy->SetValue(3,4,23);
  yy->SetValue(3,5,24);
  yy->SetValue(4,0,25);
  yy->SetValue(4,1,26);
  yy->SetValue(4,2,27);
  yy->SetValue(4,3,28);
  yy->SetValue(4,4,29);
  yy->SetValue(4,5,30);
  
  
  for (auto ii=0;ii<rr1;ii++)
  {
      for(auto jj=0;jj<cc1;jj++)
      {
          for(auto kk=0;kk<rr2;kk++)
          {
              for(auto ll=0;ll<cc2;ll++)
              {
                  
                  zz->SetValue(kk+ii*rr2,ll+jj*cc2, xx->GetValue(ii,jj)*yy->GetValue(kk,ll));
              }
          }
      }
  }
  

  std::cout << "Direct Kronecker " << std::endl <<*zz << std::endl;
  zz->Kronecker(*xx,*yy, *zz);
  std::cout << "Kronecker xx*yy from MTK_DenseMatrix"<< std::endl <<*zz << std::endl;
  zz->Kronecker(*yy,*xx, *zz);
  std::cout << "Kronecker yy*xx from MTK_DenseMatrix"<< std::endl <<*zz << std::endl;
 


  
    rr=xx->DenseMatrix_To_MTK_Real(*xx)  ;
  std::cout << "Matrix x as MTK_Real Using Continnum: "<<std::endl;
  for (auto ii=0; ii<rr1*cc1 ; ii++)
  {
   std::cout <<rr[ii] << " " ;
  }
  std::cout <<std::endl;
  
  
  pp=  pp->MatrixMultiplication(*xx,*yy);
  std::cout<<"Product x*y from matrix mult using dgemm"<<std::endl<< *pp << std::endl;  

  
  free(xx);
  free(yy);  

  free(zz);
  free(pp);
  return 0;
}
