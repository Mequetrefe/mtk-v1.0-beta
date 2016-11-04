/*
 * 
 Define problem to solve:

 Consider a steady diffusive reactive process of the
 form:

     -lap(u(x, p)) = s(x, p), x \in [a,b] = [0,1],
 where p \in IR^k, k >= 0 is a 1D array of k parameters.

 Function "source" defines the source term which will
 end up being the right-hand side (RHS) of the arising system
 of equations.

 We consider Robin's boundary conditions of the form:

 alpha*u(a) - beta*u'(a) = west_bndy_value
alpha*u(b) + beta*u'(b) = east_bndy_value
*/


#if __cplusplus == 201103L

#include <iostream>
#include <fstream>
#include <cmath>

#include "mtk.h"
#include "mtk_1d_lap.h"
#include "mtk_dense_matrix.h"
#include "mtk_lapack_facade.h"

void source(mtk::MTK_Real *xx,mtk::MTK_Real *yy,mtk::MTK_Real lambda, int nn)
{
    for( auto ii=0;ii<nn;ii++)
    {
        yy[ii]=-lambda*lambda*exp(lambda*xx[ii])/(exp(lambda)-1.0);
    }
}

void exact_sol(mtk::MTK_Real *xx,mtk::MTK_Real *yy,mtk::MTK_Real lambda, int nn)
{
    for( auto ii=0;ii<nn;ii++)
    {
        yy[ii]=(exp(lambda*xx[ii])-1.0)/(exp(lambda)-1.0);
    }
}


int main (int argc, char *argv[]) {

// Parematers for the discrete mimetic operators
  int Order_of_Accuracy=4;
  int NumCells=3*Order_of_Accuracy;

// Choose your mimetic threshold

  mtk::MTK_Real Tau=1e-6;
  
  
  std::cout<<argc<<" "<<argv[0]<<std::endl;
  
  Order_of_Accuracy=atoi(argv[1]);
  NumCells=atoi(argv[2]);
  Tau=atof(argv[3]);

// Number of cell centers plus the boundary conditions
// east and west.
  int nn = NumCells+2;
  int mm = NumCells+2;
  
  std::streambuf *psbuf;
  std::streambuf *backup;
  std::ofstream filestr;

  mtk::MTK_Real West = 0.0;
  mtk::MTK_Real East = 1.0;
  mtk::MTK_Real west_bndy_value = -1.0;
  mtk::MTK_Real east_bndy_value = 0.0;

  mtk::MTK_Real lambda=-1.0;
  mtk::MTK_Real alpha=-exp(lambda);
  mtk::MTK_Real beta=(exp(lambda)-1.0)/lambda;
  mtk::MTK_Real staggering_freq = (East-West)/NumCells;
// The Laplacian will be returned as matrix

  
  mtk::MTK_DenseMatrix *ll;
  mtk::MTK_DenseMatrix *llT;
  mtk::MTK_DenseMatrix *gg;
  
  
//Pointer to the object MTK_1DLap, yet to be instanced. 

  mtk::MTK_1DLap *lap = new mtk::MTK_1DLap(Order_of_Accuracy,Tau);
  mtk::MTK_1DGrad *grad = new mtk::MTK_1DGrad(Order_of_Accuracy,Tau);

// Unidimensional grid from 0 to 1:

  mtk::MTK_Real *XX;
// Relative error calculations
  mtk::MTK_Real norm_sol;
  mtk::MTK_Real norm_err;
  

  
  std::cout << "Testing MTK 1D Laplacian solving a Poisson Equation" << std::endl;
  



  XX=new mtk::MTK_Real[mm];
  if (XX == nullptr)
  {
    std::cout << "Problem getting memory for the staggered grid. "<< endl << "Exiting now..."<<endl;
    return -1;
  }
  else
  {
    XX[0]=West;
    XX[1]=West+staggering_freq/2.0;
    XX[mm-1]=East;
    XX[mm-2]=East-staggering_freq/2.0;
    for (auto ii=2;ii<mm-2;ii++)
    {
        XX[ii]=XX[ii-1]+staggering_freq;
    }
  } 
  
// Known Solution holder:

  mtk::MTK_Real *Known_Sol;
  Known_Sol=new mtk::MTK_Real[mm];
  if (Known_Sol == nullptr)
  {
    std::cout << "Problem getting memory for the known solution holder. "<< endl << "Exiting now..."<<endl;
    return -1;
  }
  else
  {
   exact_sol(XX,Known_Sol,lambda,mm);
   
  } 

// Calculated solution:

  mtk::MTK_Real *YY;
  YY=new mtk::MTK_Real[mm];
  if (YY == nullptr)
  {
    std::cout << "Problem getting memory for the calculated solution. "<< endl << "Exiting now..."<<endl;
    return -1;
  }
  else
  {
   source(XX,YY,lambda,mm);
   YY[0]=west_bndy_value;
   YY[mm-1]=east_bndy_value;
  } 

 
// Build and instance of MTK1DLap:

  lap = lap->Construct1DLap();
  grad = grad->Construct1DGrad();
  
  if (lap == nullptr) {
    std::cout << "Unable to construct the Laplacian Operator.... Exiting. "<<endl;
    return -1;
  }


  ll = lap->ReturnAsMatrix(NumCells, West, East);
  
  gg = grad->ReturnAsMatrix(NumCells, West, East);
  
  llT=new mtk::MTK_DenseMatrix(nn,mm);
  
//  std::cout << "Laplacian order " << Order_of_Accuracy<< std::endl; 

  //std::cout << *ll << std::endl ;
      
  

   for(auto ii=0;ii<nn;ii++)
  {
      
    for(auto jj=0;jj<mm;jj++)
    {
       ll->SetValue(ii,jj,-ll->GetValue(ii,jj));

       
    }

      
}

  
  std::cout << "Applying gradient as BC " << endl; 

  
    for(auto jj=0;jj<mm;jj++)
    {
  
        ll->SetValue(0,jj,-beta*gg->GetValue(00,jj) );

        ll->SetValue(nn-1,jj, beta*gg->GetValue((nn-2),jj));
       }
   

    


//  std::cout << *ll << std::endl ;

  std::cout << "Applying Diritchlet conditions"<< endl; 


  ll->SetValue(0,0,(ll->GetValue(0,0) +alpha));
  
  ll->SetValue(nn-1,mm-1,(ll->GetValue(nn-1,mm-1) +alpha));
  
  
  //std::cout << *ll << std::endl ;

//   std::cout << "Staggered Grid:" << endl;
//   for(auto ii=0;ii<mm;ii++)
//      std::cout << XX[ii] << " " << endl;
//   std::cout << "Source function:" << endl;
//   for(auto ii=0;ii<mm;ii++)
//      std::cout << YY[ii] << " " << endl;

//Transpose for dgesv
  ll->Transpose(*ll,*llT); 
// Needed parameters for dgesv
    int N=mm; 
    int nrhs = 1;
    int info = 0;
    int lda = nn;
    int ldb=mm;
    int *ipiv;
    ipiv = new int[mm];
    mtk::dgesv_(&N,&nrhs,llT->DenseMatrix_To_MTK_Real(*llT),&lda,ipiv,YY,&ldb,&info);
   if (!info){
     std::cout << endl<<"Solution Obtained!" <<endl;
   }
   else
   {
     std::cout << "Error from dgesv" << info <<endl;
     return 1;
   }
  for(auto ii=0;ii<nn;ii++)
//     std::cout << YY[ii] << " " << endl;
 
//   std::cout<<"Writing results to file \"results.txt\" in current directory."<<endl;
//   
//   filestr.open("results.txt");
//   backup=std::cout.rdbuf();
//   psbuf=filestr.rdbuf();
//   std::cout.rdbuf(psbuf);
//   
//    
//   std::cout << "Testing MTK 1D Laplacian solving a Poisson Equation" << std::endl;
//   
// 
//   std::cout << "Laplacian order " << Order_of_Accuracy<< " For " << NumCells << " cell centers."<<endl; 
//   std::cout << "Mimetic Threshold set to "<< Tau <<endl; 
//   std::cout<<"--------------------------------------------------------------------------------------"<<endl;
//   std::cout << "Staggered Grid"<<"\t"<<"        Calculated Solution" <<"\t"<<"Exact Solution"<<"\t"<<"\t"<<"Abs. Error"<<endl;
//   std::cout<<"--------------------------------------------------------------------------------------"<<endl;
//   
//   for(auto ii=0;ii<nn;ii++)
//   {
//     std::cout <<std::fixed<<std::setprecision(12)<<std::right<<XX[ii]<<"\t"<<"\t"<<YY[ii]<<"\t"<<"\t"<<Known_Sol[ii]<<"\t"<<"\t"<<abs(YY[ii]-Known_Sol[ii])<<endl;
//   }
//   std::cout<<"--------------------------------------------------------------------------------------"<<endl;
  norm_sol=0.0;
  norm_err=0.0;
  for(auto ii=0;ii<nn;ii++)
  {
    norm_sol=norm_sol+Known_Sol[ii]*Known_Sol[ii];
    norm_err=norm_err+abs(YY[ii]-Known_Sol[ii])*abs(YY[ii]-Known_Sol[ii]);
  }
  norm_sol=sqrt(norm_sol);
  norm_err=sqrt(norm_err);
  std::cout <<"Relative Error: "<< norm_err/norm_sol << endl;
/*  
  std::cout.rdbuf(backup);
  filestr.close();*/
   delete[] XX; 
   delete[] YY;
   delete[] Known_Sol;
   delete[] ipiv;
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
