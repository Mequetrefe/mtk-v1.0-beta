/*!
  \file mtk_grid_1D.h

  \brief Includes the definition of the class MTKGrid1D.

  This file contains the definition of the class MTKCDenseMatrix.
  This class implements a sparse matrix using Compressed Column Scheme (CCS) to
  ensure storing only the non-zero elements of a given dense matrix. As a
  consequence of the selected storage scheme, no assumptions are being made in
  terms of the actual structure of the matrix.

  \date: Sunday, August 14 2011 11:42 PM
  \version: 2011-14-08-01.
  \author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
 */
 /*
  Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
  State University (SDSU).

  http://www.csrc.sdsu.edu/mtk/

  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification are permitted provided that the following conditions are met:

  -# Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.
  -# Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
  -# Neither the name of the CSRC, SDSU nor the names of its contributors may
  be used to endorse or promote products derived from this software without
  specific prior written permission.
  -# Modifications whether they are partial or complete; whether they are
  additions or eliminations should be reported through an email at:
  esanchez@sciences.sdsu.edu

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
 */

#ifndef MTK_INCLUDE_MTK_GRID_1D_H
#define MTK_INCLUDE_MTK_GRID_1D_H

#include "mtk_enums.h"
#include "mtk_dense_matrix.h"
#include "mtk_roots.h"
namespace mtk
{

//class MTK_DenseMatrix;

/*! \brief Compatible class MTKCSSparseMatrixDouble.

  This class implements the Compressed Row Scheme (CCS) format for sparse
  matrices. This is a general data type that makes no assumptions about the
  sparsity of the original matrix. This is also compatible with well known
  sparse matrices implementations such as the one in
  <a href="http://crd.lbl.gov/~xiaoye/SuperLU/">SuperLU</a>.
 */
class MTKGrid1D {

  friend class MTKGridManager;
  
public:
  /*! Default constructor.  */
  /*! This functions CONSTRUCTS a default CRS sparse matrix. */
  MTKGrid1D(void);

  /*! Constructs a CRS sparse matrix form a 1D array.  */
  /*! This functions CONSTRUCTS a CRS sparse matrix from an already predefined
   * matrix which is stored as an one dimensional array.
   * \param *array Input pointer to the array which hold the original matrix.
   * \param num_rows Number of rows of the original matrix.
   * \param num_cols Number of rows of the original matrix.
   * \return The status of the creation as a logical value.
   */
  MTKGrid1D(int nDummy, MTK_Real dx,
         MTK_Real xmin, MTK_Real xmax );

  /*! Default destructor.  */
  /*! This functions DESTRUCTS an already created dense matrix. */

  bool SetGridValues (void);

  bool SetDataValues (int ii, MTK_Real vv);

  MTK_Real GetGridValue(int ii);

  MTK_Real GetDataValue(int ii);

  MTKGrid1D* CGM1DDiv(BindingType1D btype);

  int nNodes() {

    return nNodes_;
  }

  int nDummy() {

    return nDummy_;
  }

  int datanelm() {

    return datanelm_;
  }

  int datanDummy() {

    return datanDummy_;
  }

  MTK_Real GetData1Df(int Index);

  void ComputeMetric(BindingType1D btype);

  bool CalcInterpC2N(void);

  bool CalcInterpN2C(void);

  bool BindDataonGrid1D(BindingType1D btype);

  bool Apply1DPeriodic(void);

  MTKGrid1D* MultData1D(MTK_Real);

  void PrintArray();

  void PrintCenterBound_NoBC();

  void PrintNodeBound_NoBC();


  /*// Default destructor.
  ! This functions DESTRUCTS an already created dense matrix.
  //int nRow() {

   // return nRow_;
 // }

  Default destructor.
  ! This functions DESTRUCTS an already created dense matrix. */
  //int nCol() {

    //return nCol_;
  //}

  /*! Default destructor.  */
  /*! This functions DESTRUCTS an already created dense matrix. */
  //MTK_Real GetValue(int rr, int cc);

  /*! Default destructor.  */
  /*! This functions DESTRUCTS an already created dense matrix. */
  //void PrintBlock();

  /*! Default destructor.  */
  /*! This functions DESTRUCTS an already created dense matrix. */
  ~MTKGrid1D();

private:
  int nNodes_;         /*!< Number of rows. */
  int nDummy_;         /*!< Number of columns. */
  int datanelm_;
  int datanDummy_;
  MTK_Real dx_;
  MTK_Real xmin_;
  MTK_Real xmax_;
  MTK_Real *Grid_;
  MTK_Real *CBMetric_;
  MTK_Real *NBMetric_;
  MTK_Real *Data_;
  BindingType1D btype_;
  MTK_DenseMatrix *InterpC2N_;
  MTK_DenseMatrix *InterpN2C_;
};
}
#endif
