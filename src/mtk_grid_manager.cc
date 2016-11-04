#include <stdio.h>
#include <stdlib.h>

#include "mtk_enums.h"
#include "mtk_grid_1D.h"
#include "mtk_grid_manager.h"

using namespace std;
using namespace mtk;


MTKGrid1D* MTKGridManager::MultData1D (MTKNumber mult,
                                       MTKGrid1D* grid_in) {

  MTKGrid1D *Mult = new MTKGrid1D(grid_in->nDummy_,
                                             grid_in->dx_,
                                             grid_in->xmin_,
                                             grid_in->xmax_);


  Mult->btype_ = grid_in->btype_;
  Mult->nNodes_ = grid_in->nNodes_;
  Mult->datanDummy_ = grid_in->datanDummy_ ;
  Mult->datanelm_ = grid_in->datanelm_;

  int nelm = Mult->nNodes_-1;
  int ii;

  for ( ii=0; ii<Mult->nNodes_+2*Mult->nDummy_; ii++)
    Mult->Grid_[ii]= (Mult->xmin_) - (Mult->nDummy_)*(Mult->dx_) + ii*(Mult->dx_);


  Mult->CBMetric_= (MTKNumber *)malloc(nelm*sizeof(MTKNumber));

  for (int ii =0; ii<nelm; ii++)
    Mult->CBMetric_[ii] = grid_in->CBMetric_[ii];

  Mult->Data_ = (MTKNumber*) calloc( (Mult->datanelm_),sizeof(MTKNumber) );
    if (Mult->Data_ == NULL){
      printf("Could not allocate data!!!\n");
      return false;
    }

  for (ii=0;ii<Mult->nNodes_-1;ii++) {
    Mult->Data_[ii+Mult->nDummy_] = grid_in->Data_[ii+Mult->nDummy_]*mult;
    //printf("%f\n",Mult->Data_[ii]);
  }
  return Mult;
}

MTKGrid1D* MTKGridManager::CGM1DDiv( MTKGrid1D* const grid_in ){

  if (grid_in->btype_ != CenterBound && grid_in->btype_ != CenterBound_NoBC) {
    fprintf(stderr, "Incompatible btype_ Type");
  }


  if (grid_in->btype_ == NullBind){
    fprintf(stderr, "btype_ Metric Not Constructed");
  }

  if (grid_in->CBMetric_ == NULL) {
    fprintf(stderr, "No Metric");
  }

  MTKGrid1D* Div = new MTKGrid1D(grid_in->datanDummy_, grid_in->dx_, grid_in->xmin_, grid_in->xmax_);
  int nelm = grid_in->nNodes_-1;
  int ii;

  Div->btype_ = grid_in->btype_;
  Div->nNodes_ = grid_in->nNodes_;
  Div->datanDummy_ = grid_in->datanDummy_ ;
  Div->datanelm_ = grid_in->datanelm_;

  for ( ii=0; ii<Div->nNodes_+2*Div->nDummy_; ii++)
    Div->Grid_[ii]= (Div->xmin_) - (Div->nDummy_)*(Div->dx_) + ii*(Div->dx_);

  Div->CBMetric_= (MTKNumber *)malloc(nelm*sizeof(MTKNumber));

  for (int ii=0; ii<nelm; ii++)
    Div->CBMetric_[ii] = grid_in->CBMetric_[ii];

  Div->Data_ = (MTKNumber*) calloc( (Div->datanelm_),sizeof(MTKNumber) );
  if (Div->Data_ == NULL){
    printf("Could not allocate data!!!\n");
    return false;
  }

  switch (Div->btype_){

    case NullBind:{
      fprintf(stderr, "Need Bind Type");
      break;
    }
    case GradNaturalGrid:{
      fprintf(stderr, "Blah1");
      break;
    }
    case DivNaturalGrid:{
      fprintf(stderr, "Blah2");
      break;
    }
    case NodeBound_NoBC:{
      fprintf(stderr, "Blah3");
      break;
    }
    case NodeBound:{
      fprintf(stderr, "Blah4");
      break;
    }
    case CenterBound:{
      for(int i=0; i<nelm; i++) {
        Div->Data_[i+grid_in->datanDummy_]= grid_in->CBMetric_[i]*grid_in->GetData1Df(i+1)-
        grid_in->CBMetric_[i]*grid_in->GetData1Df(i);
      }
      break;
    }

    case CenterBound_NoBC:{
      fprintf(stderr, "Blah5\n");
      break;
    }
    default:{
      return(NULL);
    }
  }
  return Div;
}


MTKGrid1D* MTKGridManager::Set2Sum1D (MTKGrid1D* const grid_in1,
                                      MTKGrid1D* const grid_in2,
                                      MTKNumber c){


  //if (grid_in1->Grid_ !=grid_in2->Grid_)
    //return(NULL);

  if (grid_in1->btype_!=NodeBound && grid_in1->btype_!=CenterBound)
    return(NULL);

  if (grid_in2->btype_!=NodeBound && grid_in2->btype_!=CenterBound)
    return(NULL);

  MTKGrid1D* Output = new MTKGrid1D(grid_in1->datanDummy_, grid_in1->dx_, grid_in1->xmin_, grid_in1->xmax_);

  Output->btype_ = grid_in1->btype_;
  Output->nNodes_ = grid_in1->nNodes_;
  Output->datanDummy_ = grid_in1->datanDummy_ ;
  Output->datanelm_ = grid_in1->datanelm_;

  int nelm = Output->nNodes_-1;
  int ii;

  for ( ii=0; ii<Output->nNodes_+2*Output->nDummy_; ii++)
   Output->Grid_[ii]= (Output->xmin_) - (Output->nDummy_)*(Output->dx_) + ii*(Output->dx_);

  Output->CBMetric_= (MTKNumber *)malloc(nelm*sizeof(MTKNumber));

  for (ii=0; ii<nelm; ii++)
    Output->CBMetric_[ii] = grid_in1->CBMetric_[ii];

  Output->Data_ = (MTKNumber*) calloc( (Output->datanelm_),sizeof(MTKNumber) );

  if (Output->Data_ == NULL){
    printf("Could not allocate data!!!\n");
    return false;
  }

  MTKNumber shiftIndex;

  if (Output->btype_==NodeBound) {
    shiftIndex=0.0;
  }
  else {
    shiftIndex=0.5;

  }

  for(ii=Output->nDummy_;(ii<Output->datanelm_-Output->nDummy_);ii++) {
    Output->Data_[ii]=grid_in1->Data_[ii]+
    c*grid_in2->GetData1Df((MTKNumber)(ii-grid_in2->nDummy_) + shiftIndex);
  }
    return Output;
}

  bool MTKGridManager::Apply1DPeriodic(MTKGrid1D *grid_in){
    grid_in->SetDataValues(0,grid_in->GetDataValue(grid_in->datanelm_-grid_in->datanDummy_-2));
    grid_in->SetDataValues(1,grid_in->GetDataValue(grid_in->datanelm_-grid_in->datanDummy_-1));
    grid_in->SetDataValues(grid_in->datanelm_-2,grid_in->GetDataValue(grid_in->datanDummy_));
    grid_in->SetDataValues(grid_in->datanelm_-1,grid_in->GetDataValue(grid_in->datanDummy_+1));

    return true;

  }


