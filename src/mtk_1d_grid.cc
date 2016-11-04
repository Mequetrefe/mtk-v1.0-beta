/*!
 * \file mtk_grid_1D.cc
 *
 * \brief Definition of the default constructor for the
 * MTKCCSSparseMatrix class.
 *
 * This file contains the definition of the default constructor for the
 * MTKCCSSparseMatrix.
 *
 * \date: Monday, September 03, 2012
 * \version: 2012-09-03.
 * \author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
 */
 /* Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
 * State University (SDSU).
 *
 * http:www.csrc.sdsu.edu/mtk/
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification are permitted provided that the following conditions are met:
 *
 * -# Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * -# Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * -# Neither the name of the CSRC, SDSU nor the names of its contributors may
 * be used to endorse or promote products derived from this software without
 * specific prior written permission.
 * -# Modifications whether they are partial or complete; whether they are
 * additions or eliminations should be reported through an email at:
 * esanchez@sciences.sdsu.edu
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>


#include "mtk_1d_grid.h"


using namespace std;
using namespace mtk;

/*! Constructs a default MTKGrid1D instance. */

/*! Constructs a default MTKGrid1D instance. */
MTKGrid1D::MTKGrid1D(void):
  nNodes_(0),         /*!< Number of rows. */
  nDummy_(0),         /*!< Number of columns. */
  datanelm_(0),
  datanDummy_(0),
  dx_(0),
  xmin_(0),
  xmax_(0),
  Grid_(NULL),
  CBMetric_(NULL),
  NBMetric_(NULL),
  Data_(NULL),
  btype_(NullBind),
  InterpC2N_(NULL),
  InterpN2C_(NULL) {
}

/*! Constructs MTK_DenseMatrix from a 1D array dense matrix. */

/*! Constructs MTKCCSSparseMatrix from a dense matrix which is implemented
 * as a 1D array.
 * \param *array INPUT: Pointer to the beginning of the 1D array implementing
 * the dense matrix.
 * \param num_rows INPUT: Number of rows of the dense matrix.
 * \param num_cols INPUT: Number of cols of the dense matrix.
 * \param unitary_offset INPUT: States if ...
 */
MTKGrid1D::MTKGrid1D(int nDummy,
                     MTK_Real dx,
                     MTK_Real xmin,
                     MTK_Real xmax )
    : datanelm_(0),
      datanDummy_(0),
      CBMetric_(NULL),
      NBMetric_(NULL),
      Data_(NULL),
      btype_(NullBind),
      InterpC2N_(NULL),
      InterpN2C_(NULL)
      {

  /* Verify first that the requested sizes make sense: */
  if (nDummy < 0) {
    fprintf(stderr, "Incorrect parameter at line %d.\n", __LINE__ - 1);
  }
  if (dx <= 0) {
    fprintf(stderr, "Incorrect parameter at line %d.\n", __LINE__ - 1);
  }

        /*!< Number of columns. */
  nDummy_= nDummy;
  dx_ = dx;
  xmin_ = xmin;
  xmax_ = xmax;

  nNodes_ = (int)( (xmax_ - xmin_) / dx_ + 1.5 );


  Grid_ = (MTK_Real *)malloc((nNodes_+2*nDummy_)*sizeof(MTK_Real));
  if (Grid_ == (MTK_Real *) NULL) {
    fprintf(stderr, "Allocation error at line %d.\n", __LINE__ - 2);
  }
 }

  bool MTKGrid1D::SetGridValues (void) {

    int ii;

    for ( ii=0; ii<nNodes_+2*nDummy_; ii++) {
      Grid_[ii]= xmin_ - nDummy_*dx_ + ii*dx_;

    }
    return true;
  }

  MTK_Real MTKGrid1D::GetGridValue(int ii){
    return Grid_[ii];
  }

  bool MTKGrid1D::SetDataValues (int ii, MTK_Real vv){
    Data_[ii] = vv;
    return true;
  }

  MTK_Real MTKGrid1D::GetDataValue(int ii){
    return Data_[ii];
  }


  void MTKGrid1D::ComputeMetric(BindingType1D btype){


  switch(btype) {

    case CenterBound_NoBC: {
      int nelm = nNodes_-1;
      CBMetric_= (MTK_Real *)malloc(nelm*sizeof(MTK_Real));
      //printf("Now Here Really!!\n");
      if(!CBMetric_) {
        fprintf(stderr, "Binding Type error at line %d.\n",__LINE__ - 1);
      }

      int i;
      for (i=0; i < nelm; i++){
        CBMetric_[i]= 1.0 /
        (Grid_[i+nDummy_+1]-Grid_[i+nDummy_]);
      }
      break;
    }

    case NodeBound_NoBC: {
      int nelm = nNodes_;
      NBMetric_= (MTK_Real *)malloc(nelm*sizeof(MTK_Real));

      if(!NBMetric_) {
        fprintf(stderr, "Binding Type error at line %d.\n",__LINE__ - 1);
      }

      int i;
      for (i=0; i < nelm; i++){
        NBMetric_[i]= 2.0 /
        (Grid_[i+nDummy_+1]-Grid_[i+nDummy_-1]);
      }
      break;
    }

    case  NullBind: {
      fprintf(stderr, "Binding Not Set");
      break;
    }

    case GradNaturalGrid: {
      fprintf(stderr, "%d",nDummy_);
      break;
    }

    case  DivNaturalGrid: {
      fprintf(stderr, "%d",nDummy_);
      break;
    }

    case  NodeBound: {
      fprintf(stderr, "%d",nDummy_);
      break;
    }

    case  CenterBound: {
      fprintf(stderr, "%d",nDummy_);
      break;
    }

    default: {
      fprintf(stderr, "Binding Type error at line %d.\n", __LINE__ - 2);
      break;
    }
  }
}

bool MTKGrid1D::CalcInterpC2N(void){
  int nRow=nNodes_;
  int nCol=4;
  MTK_DenseMatrix *InterpC2N_ = new MTK_DenseMatrix(nRow, nCol);

  if (nDummy_<2) return false;

  int i;
  MTK_Real a,b,c,d;
  for(i=0;i<nRow;i++){
    a=((Grid_[i+nDummy_-2]+Grid_[i+nDummy_-1])*0.5-Grid_[i+nDummy_]);
    b=((Grid_[i+nDummy_-1]+Grid_[i+nDummy_  ])*0.5-Grid_[i+nDummy_]);
    c=((Grid_[i+nDummy_  ]+Grid_[i+nDummy_+1])*0.5-Grid_[i+nDummy_]);
    d=((Grid_[i+nDummy_+1]+Grid_[i+nDummy_+2])*0.5-Grid_[i+nDummy_]);
    InterpC2N_->SetValue (-b*c*d/( (a-b)*(a-c)*(a-d) ), i, 0);
    InterpC2N_->SetValue (-a*c*d/( (b-a)*(b-c)*(b-d) ), i, 1);
    InterpC2N_->SetValue (-a*b*d/( (c-a)*(c-b)*(c-d) ), i, 2);
    InterpC2N_->SetValue (-a*b*c/( (d-a)*(d-b)*(d-c) ), i, 3);
  }
//   InterpC2N_->PrintBlock();
  return true;
}

bool MTKGrid1D::CalcInterpN2C(void){
  int nRow=nNodes_-1;
  int nCol=4;
  MTK_DenseMatrix *InterpN2C_ = new MTK_DenseMatrix(nRow, nCol);

  int i;
  MTK_Real a,b,c,d;
  for(i=0;i<nRow;i++){
    a=(Grid_[i+nDummy_-1]-(Grid_[i+nDummy_]+Grid_[i+nDummy_+1])*0.5);
    b=(Grid_[i+nDummy_  ]-(Grid_[i+nDummy_]+Grid_[i+nDummy_+1])*0.5);
    c=(Grid_[i+nDummy_+1]-(Grid_[i+nDummy_]+Grid_[i+nDummy_+1])*0.5);
    d=(Grid_[i+nDummy_+2]-(Grid_[i+nDummy_]+Grid_[i+nDummy_+1])*0.5);
    InterpN2C_->SetValue (-b*c*d/( (a-b)*(a-c)*(a-d) ), i, 0);
    InterpN2C_->SetValue (-a*c*d/( (b-a)*(b-c)*(b-d) ), i, 1);
    InterpN2C_->SetValue (-a*b*d/( (c-a)*(c-b)*(c-d) ), i, 2);
    InterpN2C_->SetValue (-a*b*c/( (d-a)*(d-b)*(d-c) ), i, 3);
  }
//   InterpN2C_->PrintBlock();
  return true;
}

bool MTKGrid1D::BindDataonGrid1D(BindingType1D btype){

//   if (btype != btype_){
//     fprintf(stderr, "Incompatible Binding Type!!!\n");
//     return false;
//   }
  btype_ = btype;
  if (btype_ == NullBind){
    printf("Binding Metric Not Constructed");


    return false;
  }

  switch (btype_){

    case NullBind:{
      fprintf(stderr, "Need Bind Type");
    }
    case GradNaturalGrid:{
      datanelm_=nNodes_+1;
      datanDummy_=0;
      break;
    }
    case DivNaturalGrid:{
      fprintf(stderr, "Blah");
      break;
    }
    case NodeBound_NoBC:{
      datanelm_=nNodes_;
      datanDummy_=0;
      break;
    }
    case NodeBound:{
      datanelm_=nNodes_+2*nDummy_;
      datanDummy_=nDummy_;
      break;
    }
    case CenterBound:{
      datanelm_= nNodes_+2*nDummy_-1;
      datanDummy_=nDummy_;

      if (Data_ == NULL){
        Data_ = (MTK_Real*) calloc(datanelm_,sizeof(MTK_Real) );
        if (Data_ == NULL){
          printf("Could not allocate data!!!\n");
          return false;
        }
      }
      break;
    }
    default:{
      return(NULL);
      break;
    }
  }


  return true;
}

  void MTKGrid1D::PrintArray(){
    int NN = nNodes_ + 2*nDummy_;
    int ii;

    for (ii = 0; ii<NN; ii++){
      fprintf(stdout,"%5.5g\t",Grid_[ii]);
    }
    fprintf(stdout,"\n");
  }

  void MTKGrid1D::PrintCenterBound_NoBC(){
    int NN = nNodes_ - 1;
    int ii;

    for (ii = 0; ii<NN; ii++){
      fprintf(stdout,"%5.5g\t",CBMetric_[ii]);
    }
    fprintf(stdout,"\n");
  }

void MTKGrid1D::PrintNodeBound_NoBC(){
  int NN = nNodes_;
  int ii;

  for (ii = 0; ii<NN; ii++){
    fprintf(stdout,"%5.5g\t",NBMetric_[ii]);
  }
  fprintf(stdout,"\n");
}


MTK_Real MTKGrid1D::GetData1Df(int Index){

    int I0,I1,I2,I3;
    MTK_Real a,b,c,d;
    MTK_Real  a0,a1,a2,a3;
    MTK_Real  x,IndexFrac;

    if ( (Index<0) || (Index > nNodes_-1) ){
      fprintf(stderr, "Bad Index");
      return(0);
    }

    switch (btype_){

      case NullBind:{
        fprintf(stderr, "Need Bind Type");

      case NodeBound:
        if (floorf(Index)==Index){
          return(Data_[nDummy_+(int)Index]);
        }
        else{
          IndexFrac=Index-floorf(Index);
          I0=(int) floorf(Index-1.0f)+nDummy_;
          I1=(int) floorf(Index)+nDummy_;
          I2=(int) floorf(Index+1.0f)+nDummy_;
          I3=(int) floorf(Index+2.0f)+nDummy_;

          x=(1.0f-IndexFrac)*Grid_[I1]+IndexFrac*Grid_[I2];
          a=(Grid_[I0]-x);
          b=(Grid_[I1]-x);
          c=(Grid_[I2]-x);
          d=(Grid_[I3]-x);
          a0=-b*c*d/( (a-b)*(a-c)*(a-d) );
          a1=-a*c*d/( (b-a)*(b-c)*(b-d) );
          a2=-a*b*d/( (c-a)*(c-b)*(c-d) );
          a3=-a*b*c/( (d-a)*(d-b)*(d-c) );
          return( a0*Data_[I0]+
          a1*Data_[I1]+
          a2*Data_[I2]+
          a3*Data_[I3] );
        }
        break;
      case CenterBound:
        if ( (floorf(Index)+0.5)==Index ){
          return(Data_[nDummy_+(int) floorf(Index)]);
        }
        else{
          IndexFrac=Index-floorf(Index);
          I0=(int) (floor(Index+0.5f)-2+nDummy_);
          I1=(int) (floor(Index+0.5f)-1+nDummy_);
          I2=(int) (floor(Index+0.5f)+nDummy_);
          I3=(int) (floor(Index+0.5f)+1+nDummy_);
          x=(1.0f-IndexFrac)*Grid_[(int)
          floorf(Index)+nDummy_]+IndexFrac*Grid_[(int)
          floorf(Index+1.0f)+nDummy_];
          a=((Grid_[I0]+Grid_[I0+1])*0.5f-x);
          b=((Grid_[I1]+Grid_[I1+1])*0.5f-x);
          c=((Grid_[I2]+Grid_[I2+1])*0.5f-x);
          d=((Grid_[I3]+Grid_[I3+1])*0.5f-x);
          a0=-b*c*d/( (a-b)*(a-c)*(a-d) );
          a1=-a*c*d/( (b-a)*(b-c)*(b-d) );
          a2=-a*b*d/( (c-a)*(c-b)*(c-d) );
          a3=-a*b*c/( (d-a)*(d-b)*(d-c) );
          return( a0*Data_[I0]+
          a1*Data_[I1]+
          a2*Data_[I2]+
          a3*Data_[I3] );
        }
        break;
      default:
        fprintf(stderr, "Unsupported Binding");
            return(0);
    }
}
}
/*! Constructs MTKCCSSparseMatrix from a dense matrix which is implemented
 * as a 1D array.
 * \param *array INPUT: Pointer to the beginning of the 1D array implementing
 * the dense matrix.
 * \param num_rows INPUT: Number of rows of the dense matrix.
 * \param num_cols INPUT: Number of cols of the dense matrix.
 * \param unitary_offset INPUT: States if ...
 */
MTKGrid1D::~MTKGrid1D() {

  if (Grid_ != (MTK_Real *) NULL)
    delete [] Grid_;

  if (CBMetric_ != (MTK_Real *) NULL)
    delete [] CBMetric_;

  if (Data_ != (MTK_Real *) NULL)
    delete [] Data_;
}
