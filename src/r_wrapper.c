/****************************************************************************
 * Copyright (C) 2014 by Taylor Arnold, Veeranjaneyulu Sadhanala,           *
 *                       Ryan Tibshirani                                    *
 *                                                                          *
 * This file is part of the glmgen library / package.                       *
 *                                                                          *
 *   glmgen is free software: you can redistribute it and/or modify it      *
 *   under the terms of the GNU Lesser General Public License as published  *
 *   by the Free Software Foundation, either version 2 of the License, or   *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   glmgen is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU Lesser General Public License for more details.                    *
 *                                                                          *
 *   You should have received a copy of the GNU Lesser General Public       *
 *   License along with glmgen. If not, see <http://www.gnu.org/licenses/>. *
 ****************************************************************************/

#include <Rinternals.h>
#include <string.h>
#include <stdio.h>



SEXP tf_predict_R (SEXP sBeta, SEXP sX, SEXP sN, SEXP sK, SEXP sX0, SEXP sN0,
    SEXP sNLambda, SEXP sZeroTol)
{
  /* Initialize all of the variables */
  int i;
  double * beta;
  double * x;
  double * x0;
  int n;
  int n0;
  int k;
  int nlambda;
  double zero_tol;
  double * pred;

  beta = REAL(sBeta);
  x = REAL(sX);
  x0 = REAL(sX0);
  n = asInteger(sN);
  n0 = asInteger(sN0);
  k = asInteger(sK);
  nlambda = asInteger(sNLambda);
  zero_tol = asReal(sZeroTol);

  /* Output */
  SEXP sPred;
  PROTECT(sPred = allocVector(REALSXP, n0 * nlambda));
  pred = REAL(sPred);

  for (i = 0; i < nlambda; i++)
  {
    tf_predict_gauss(beta + n*i, x, n, k, x0, n0, pred + n0*i, zero_tol);
  }

  /* Free the allocated objects for the gc and return the output as a list */
  UNPROTECT(1);
  return sPred;
}





