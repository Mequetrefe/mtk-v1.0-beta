/*!
\file mtk_matrix.h

\brief Definition of the representation of a matrix in the MTK.
 
Definition of the representation for the matrices implemented in the MTK.

\date: June 20, 2014, 11:37 AM

\version: 1.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at sciences dot sdsu dot edu

\bug No known bugs.
*/

/*
Copyright (C) 2015 Computational Science Research Center (CSRC) at San Diego
State University (SDSU).

Website for the project: http://www.csrc.sdsu.edu/mtk/

All rights reserved.

Redistribution and use in source and binary forms, with or without modification
are permitted provided that the following conditions are met:

-# Redistributions of source code must retain the above copyright notice, this
list of conditions and this disclaimer.
-# Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.
-# Neither the name of the CSRC, SDSU nor the names of its contributors may be
used to endorse or promote products derived from this software without specific
prior written permission.
-# Modifications whether they are partial or complete; whether they are
additions or eliminations should be reported through an email at:

esanchez at sciences dot sdsu dot edu

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#ifndef MTK_INCLUDE_MTK_MATRIX_H
#define	MTK_INCLUDE_MTK_MATRIX_H

#include "mtk_roots.h"

namespace mtk {
  
/*!
\class MTK_Matrix
 
\ingroup c03-data_structures

\brief Definition of the representation of a matrix in the MTK.

Definition of the representation for the matrices implemented in the MTK.
*/
class MTK_Matrix {
  
  friend class MTK_LAPACKFacade;

  public:
    
    /*!
    \brief Default constructor.
    */
    MTK_Matrix();

    /*!
    \brief Destructor.
    */
    ~MTK_Matrix();

    /*!
    \brief Gets the number of rows.
     
    \return Number of rows of the matrix.
    */
    int num_rows() const;

    /*!
    \brief Gets the number of rows.
     
    \return Number of rows of the matrix.
    */    
    int num_cols() const;

    /*!
    \brief Gets the number of values.
     
    \return Number of values of the matrix.
    */
    int num_values() const;

    /*!
    \brief Gets the number of non-zero values.
     
    \return Number of non-zero values of the matrix.
    */
    int num_non_zero() const;
    
    /*!
    \brief Gets the number of zeros.
     
    \return Number of zeros of the matrix.
    */
    int num_zero() const;

    /*!
    \brief Gets the number of null values.
     
    See http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf
     
    \return Number of null values of the matrix.
    */
    int num_null() const;
    
    /*!
    \brief Gets information regarding the matrix's unitary offset.

    Let \f$ N \f$ be any generic size. Some examples relative to matrices online
    are explained using a \f$[1,N]\f$ index type, rather than a \f$[0,N - 1]\f$
    index type. When the unitary offset is turned on, the matrix will adopt a
    \f$[1,N]\f$ index type, thus matching the examples online.
     
    \return Information regarding the matrix's unitary offset
    */
    int unitary_offset() const;
    
    /*!
    \brief Gets the number of rank.
     
    \return Rank of the matrix.
    */
    int rank() const;
    
    /*!
    \brief Gets the number of unknowns the matrix is solving for.
     
    See http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf
     
    \return Number of unknowns the matrix is solving for.
    */
    int num_unknowns() const;
    
    /*!
    \brief Gets the bandwidth.
     
    \return Bandwidth of the matrix.
    */
    int bandwith() const;
    
    /*!
    \brief Gets the matrix tolerance.
     
    For some problems, we must discriminate among non-zeros and zeros. For
    this, we use this tolerance.
     
    \return Matrix tolerance of the matrix.
    */
    MTK_Real matrix_tol() const;
    
    /*!
    \brief Gets the absolute density.
     
    See http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf
     
    \return Absolute density of the matrix.
    */
    MTK_Real abs_density() const;
    
    /*!
    \brief Gets the relative density.
     
    See http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf
     
    \return Relative density of the matrix.
    */
    MTK_Real rel_density() const;
    
    /*!
    \brief Gets the Absolute sparsity.
     
    See http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf
     
    \return Absolute sparsity of the matrix.
    */
    MTK_Real abs_sparsity() const;
    
    /*!
    \brief Gets the Relative sparsity.
     
    See http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf
     
    \return Relative sparsity of the matrix.
    */
    MTK_Real rel_sparsity() const;

    /*!
    \brief Sets the number of rows of the matrix.
    
    \param[in] num_rows Number of rows.
    */
    void set_num_rows (int num_rows);
    
    /*!
    \brief Sets the number of columns of the matrix.
    
    \param[in] num_cols Number of columns.
    */
    void set_num_cols (int num_cols);
    
    /*!
    \brief Sets the number of non zeros.
    
    \param[in] num_non_zero Number of non-zeros.
    */
    void set_num_non_zero (int num_non_zero);
    
    /*!
    \brief Sets the number of zeros.
    
    \param[in] num_zero Number of zeros.
    */
    void set_num_zero (int num_zero);
    
    /*!
    \brief Sets the number of null values.
    
    See http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf
     
    \param[in] num_null Number of null values.
    */
    void set_num_null (int num_null);
    
    /*!
    \brief Sets the unitary offset.
     
    \sa num_null
    
    See http://www.csrc.sdsu.edu/research_reports/CSRCR2013-01.pdf
     
    \param[in] unitary_offset Number of null values.
    */
    void set_unitary_offset (int unitary_offset);
    
    /*!
    \brief Sets the number of unknowns.
     
    \sa num_unknowns
    
    \param[in] num_unknowns Number of unknowns to solve for.
    */
    void set_num_unknowns (int num_unknowns);

  private:
    
    /*!
    \brief Updates the attributes that must be computed as function of others.
    */
    void UpdateMatrixState();

    int num_rows_;        //!< Number of rows.
    int num_cols_;        //!< Number of columns.
    int num_values_;      //!< Number of values in matrix.
    int num_non_zero_;    //!< Number of non-zero values.
    int num_zero_;        //!< Number of zeros.
    int num_null_;        //!< Number of null values.
    int unitary_offset_;  //!< Is position 0 position 1?
    int rank_;            //!< Rank of the matrix.
    int num_unknowns_;    //!< Number of unknowns.
    int bandwith_;        //!< Bandwidth of the matrix.

    MTK_Real abs_density_;  //!< Absolute density of matrix.
    MTK_Real rel_density_;  //!< Relative density of matrix.
    MTK_Real abs_sparsity_; //!< Absolute sparsity of matrix.
    MTK_Real rel_sparsity_; //!< Relative sparsity of matrix.
};

}

#endif	// End of: MTK_INCLUDE_MTK_MATRIX_H
