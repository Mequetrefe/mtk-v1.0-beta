/*
 * 
 Define problem to solve:


 */


#if __cplusplus == 201103L

#include <iostream>
#include <fstream>
#include <cmath>

#include "mtk.h"
#include "mtk_2d_lap.h"
#include "mtk_dense_matrix.h"
#include "mtk_lapack_facade.h"

#include "mtk_2d_bc_desc.h"


/*
 * For a RHS F(x,y)=exp(-xy) * (x^4y^2 + x^2y^4 -4x^3y -4xy^3 +2x^2 +2y^2)
 * 
 *  the exact solution is
 * 
 * Lap(F(X))=x^2y^2 * exp(-xy)
 * 
 */

void source(mtk::MTK_Real *xx,mtk::MTK_Real *yy,mtk::MTK_Real *ss, int lap_size)
{
    for( auto ii=0;ii<lap_size;ii++)
    {
        ss[ii]=xx[ii]*yy[ii]*exp(-xx[ii]*xx[ii]/2.0 -yy[ii]*yy[ii]/2.0)*(xx[ii]*xx[ii]+yy[ii]*yy[ii]-6.0);
    }
}

void exact_sol(mtk::MTK_Real *xx,mtk::MTK_Real *yy,mtk::MTK_Real *ss, int lap_size)
{
    for( auto ii=0;ii<lap_size;ii++)
    {
        ss[ii]=xx[ii]*yy[ii]*exp(-xx[ii]*xx[ii]/2.0 -yy[ii]*yy[ii]/2.0);
    }
}

// This is South
mtk::MTK_Real fSouth(mtk::MTK_Real xx, mtk::MTK_Real yy)
{
 return xx*yy*exp(-xx*xx/2.0 -yy*yy/2.0);
    
}

// West boundary

mtk::MTK_Real fWest(mtk::MTK_Real xx, mtk::MTK_Real yy)
{
 return xx*yy*exp(-xx*xx/2.0 -yy*yy/2.0);
}

//East Boundary
mtk::MTK_Real fEast(mtk::MTK_Real xx, mtk::MTK_Real yy)
{
 return xx*yy*exp(-xx*xx/2.0 -yy*yy/2.0);
    
}

mtk::MTK_Real fNorth(mtk::MTK_Real xx, mtk::MTK_Real yy)
{
 return xx*yy*exp(-xx*xx/2.0 -yy*yy/2.0);
    
}

int main () {

  mtk::MTK_Real (*West)(mtk::MTK_Real, mtk::MTK_Real);
  mtk::MTK_Real (*East)(mtk::MTK_Real, mtk::MTK_Real);
  mtk::MTK_Real (*South)(mtk::MTK_Real, mtk::MTK_Real);
  mtk::MTK_Real (*North)(mtk::MTK_Real, mtk::MTK_Real);
  
  
  
  West=&fWest;
  East=&fEast;
  South=&fSouth;
  North=&fNorth;

  // The Laplacian will be returned as matrix and we need 
    // the transpose for BLAS


  
  mtk::MTK_DenseMatrix *ll;
  mtk::MTK_DenseMatrix *llT;

  // Parematers for the discrete mimetic operators
  int Order_of_Accuracy=8;
  int NumCellsX = 100;
  int NumCellsY = 100;

// Bidimensional grid :

  mtk::MTK_Real *XX;
  mtk::MTK_Real *YY;
// Choose your mimetic threshold

  mtk::MTK_Real Tau=1e-9;


// Number of cell centers plus the boundary conditions
// east and west.
  
  // mm rows (Y) by nn  columns (X)
  
  int nn = NumCellsX+2;
  int mm = NumCellsY+2;
  
  std::streambuf *psbuf;
  std::streambuf *backup;
  std::ofstream filestr;

  mtk::MTK_Real West_Bndy = 0.0;
  mtk::MTK_Real East_Bndy = 1.0;
  mtk::MTK_Real South_Bndy = 0.0;
  mtk::MTK_Real North_Bndy = 1.0;
  
  mtk::MTK_DenseMatrix *MeshX; 
  mtk::MTK_DenseMatrix *MeshY; 
 // mtk::MTK_DenseMatrix *Calculated_Solution; 
  mtk::MTK_DenseMatrix *Analytic_Solution; 
  mtk::MTK_DenseMatrix *RHS; 

  mtk::MTK_Real staggering_freqX = (East_Bndy-West_Bndy)/NumCellsX;
  mtk::MTK_Real staggering_freqY = (North_Bndy - South_Bndy)/NumCellsY;

  mtk::MTK_2DUniformStaggeredGrid *GRID = new mtk::MTK_2DUniformStaggeredGrid(0,1,NumCellsX, 0,1,NumCellsY);
  
//Pointer to the object MTK_1DLap, yet to be instanced. 

  mtk::MTK_2DLap *lap = new mtk::MTK_2DLap(Order_of_Accuracy,Tau);
mtk::MTK_2DBCDesc *BC_Desc ;
  
  std::cout << "Testing MTK 2D Laplacian solving a Poisson Equation" << std::endl;
  

// Build and instance of MTK2DLap:

  lap = lap->Construct2DLap();

  
  if (lap == nullptr) {
    std::cout << "Unable to construct the Laplacian Operator.... Exiting. "<<endl;
    return -1;
  }


  ll = lap->ReturnAsMatrix(NumCellsX, West_Bndy, East_Bndy, NumCellsY, South_Bndy, North_Bndy);


  llT=new mtk::MTK_DenseMatrix(lap->LapSize(),lap->LapSize());
  
  
  std::cout << "Laplacian order " << Order_of_Accuracy<< "Matrix Size: "; 

  std::cout << lap->LapSize() << " by "<< lap->LapSize()<<std::endl ;
  
 // Creating the Staggered Grid Values XX,YY     
  
  XX=new mtk::MTK_Real[nn];
  if (XX == nullptr)
  {
    std::cout << "Problem getting memory for the staggered grid. "<< endl << "Exiting now..."<<endl;
  }
  else
  {
    XX[0]=West_Bndy;
    XX[1]=West_Bndy+staggering_freqX/2.0;
    XX[nn-1]=East_Bndy;
    XX[nn-2]=East_Bndy-staggering_freqX/2.0;
    for (auto ii=2;ii<nn-2;ii++)
    {
        XX[ii]=XX[ii-1]+staggering_freqX;
    }
  } 

   YY=new mtk::MTK_Real[mm];
  if (YY == nullptr)
  {
    std::cout << "Problem getting memory for the staggered grid. "<< endl << "Exiting now..."<<endl;
  }
  else
  {
    YY[0]=South_Bndy;
    YY[1]=South_Bndy+staggering_freqY/2.0;
    YY[mm-1]=North_Bndy;
    YY[mm-2]=North_Bndy-staggering_freqY/2.0;
    for (auto ii=2;ii<mm-2;ii++)
    {
        YY[ii]=YY[ii-1]+staggering_freqY;
    }
  } 
 
 //Forming the Mesh as 2 Dense Matrices
 
  MeshX = new mtk::MTK_DenseMatrix(mm,nn);
  for (auto ii=0;ii<mm;ii++)
  {
      for(auto jj=0;jj<nn;jj++)
      {
          MeshX->SetValue(ii,jj,XX[jj]);
      }
  }
  std::cout << "Mesh in X direction: " << std::endl;
  std::cout << *MeshX << std::endl;
  
 
  MeshY = new mtk::MTK_DenseMatrix(mm,nn);
  for (auto ii=0;ii<mm;ii++)
  {
      for(auto jj=0;jj<nn;jj++)
      {
          MeshY->SetValue(ii,jj,YY[ii]);
      }
  }
  std::cout << "Mesh in Y direction: "<< std::endl;
  std::cout << *MeshY << std::endl;
  
  
  
//  Calculated_Solution = new mtk::MTK_DenseMatrix(mm,nn);
  Analytic_Solution = new mtk::MTK_DenseMatrix(mm,nn);
  RHS = new mtk::MTK_DenseMatrix(mm,nn);
   
  
  // Fill the Matrix Analytic_Solution with the values from the function
  
  exact_sol(MeshX->DenseMatrix_To_MTK_Real(*MeshX),MeshY->DenseMatrix_To_MTK_Real(*MeshY),Analytic_Solution->DenseMatrix_To_MTK_Real(*Analytic_Solution),lap->LapSize());
  
  // Put the RHS values inside 
  
  source(MeshX->DenseMatrix_To_MTK_Real(*MeshX),MeshY->DenseMatrix_To_MTK_Real(*MeshY),RHS->DenseMatrix_To_MTK_Real(*RHS),lap->LapSize());

  std::cout << "Source Function Before applying BC: "<< std::endl;
 // std::cout << *RHS << std::endl;
  
  
  BC_Desc = new mtk::MTK_2DBCDesc(lap, GRID,South, North, West, East);

  BC_Desc->MTK_2DBCDirichlet(ll);     

 BC_Desc->MTK_2DBoundaryRHS(RHS->DenseMatrix_To_MTK_Real(*RHS));
 
  std::cout << "Analytic Solution: "<< std::endl;
//  std::cout << *Analytic_Solution << std::endl;

  std::cout << "Source Function After applying BC "<< std::endl;
//  std::cout << *RHS << std::endl;
  
 
  
//Transpose for dgesv
  llT->Transpose(*ll,*llT);
//  cout<< *llT <<endl;
  
 
// Needed parameters for dgesv
    int N=lap->LapSize(); 
    int nrhs = 1;
    int info = 0;
    int lda = lap->LapSize();
    int ldb=lap->LapSize();
    int *ipiv;
    ipiv = new int[lap->LapSize()];
     mtk::dgesv_(&N,&nrhs,llT->DenseMatrix_To_MTK_Real(*llT),&lda,ipiv,RHS->DenseMatrix_To_MTK_Real(*RHS),&ldb,&info);
    if (!info){
      std::cout << endl<<"Solution Obtained!" <<std::endl;
    }
    else
    {
      std::cout << "Error from dgesv : " << info <<std::endl;
      return 1;
    }
  //    std::cout << *RHS <<  std::endl;
 
      
  std::cout<<"Writing results to file \"results-2d.txt\" in current directory.\b"<<endl;
  
  filestr.open("results-2D.txt");
  backup=std::cout.rdbuf();
  psbuf=filestr.rdbuf();
  std::cout.rdbuf(psbuf);
  
   
  std::cout << "Testing MTK 2D Laplacian solving a Poisson Equation" << std::endl;
  

  std::cout << "Laplacian order " << Order_of_Accuracy<< " For " << NumCellsX << " cell centers in X."<<endl; 
  std::cout << "and " << NumCellsY << " cell centers in Y."<<endl; 
  std::cout << "Mimetic Threshold set to "<< Tau <<endl; 
  std::cout<<"--------------------------------------------------------------------------------------"<<endl;
  std::cout << "\tx\t\t\ty\t       Analytic Solution\tCalculated Solution\tAbs. Error"<<std::endl;
  std::cout<<"--------------------------------------------------------------------------------------"<<endl;
  
  auto norm_sol=0.0;
  auto norm_err=0.0;
  for(auto ii=0;ii<nn;ii++)
  {
      for(auto jj =0; jj<nn; jj++)
      {
       auto xA= Analytic_Solution->DenseMatrix_To_MTK_Real(*Analytic_Solution)[ii*nn+jj];
       auto xC=RHS->DenseMatrix_To_MTK_Real(*RHS)[ii*nn+jj];
           norm_sol=norm_sol+xA*xA;
           norm_err=norm_err+(xA-xC)*(xA-xC);
    std::cout <<std::fixed<<std::setprecision(12);
    std::cout <<std::right<<XX[jj]<<"\t"<<"\t"<<YY[ii]<<"\t";
    std::cout<<"\t"<<xA<<"\t";
    std::cout<<"\t"<<xC;
    std::cout<<"\t\t"<<abs(xA-xC)<<std::endl;
      }
  }
  std::cout<<"--------------------------------------------------------------------------------------"<<endl;

  norm_sol=sqrt(norm_sol);
  norm_err=sqrt(norm_err);
  std::cout <<"Relative Error: "<< norm_err/norm_sol << endl;

  std::cout.rdbuf(backup);
  filestr.close();
   delete[] XX; 
   delete[] YY;
//   delete[] Known_Sol;
   delete[] ipiv;
  std::cout <<"Relative Error: "<< norm_err/norm_sol << endl;
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
