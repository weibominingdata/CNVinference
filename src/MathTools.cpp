#include "MathTools.h"
#include <float.h>  // for data type constant


/*
 *
 * borrowed from glm.c
 *
 */


/*
 *  glm.c
 *  
 *
 *  Created by Wei Sun on 5/11/10.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>



/**********************************************************************
 *
 * Some very simple linear algebra functions
 *
 * Originally from mla.cpp of R/CNVTools 
 * 
 
 Package: CNVtools
 Type: Package
 Title: A package to test genetic association with CNV data
 Version: 1.42.0
 Date: 2009-10-06
 Author: Chris Barnes <christopher.barnes@imperial.ac.uk> and 
 Vincent Plagnol <vincent.plagnol@cimr.cam.ac.uk>
 Maintainer: Chris Barnes <christopher.barnes@imperial.ac.uk>
 Description: This package is meant to facilitate the testing of Copy Number 
 Variant data for genetic association, typically in case-control studies.
 License: GPL-3
 Depends: survival
 biocViews: DNACopyNumber,GeneticVariability
 Packaged: 2010-04-23 00:36:41 UTC; biocbuild
 
 **********************************************************************/

const double MathTools::_PI = 3.141592653589793238462643;

const double MathTools::_LOG10_2 = 0.301029995663981195213738894724;

#ifndef __max(a,b)
#define __max(a,b)  (((a) > (b)) ? (a) : (b))
#endif

#ifndef __min(a,b)
#define __min(a,b)  (((a) < (b)) ? (a) : (b))
#endif

#define CHAR_BIT      8         /* number of bits in a char */

#define INT_MAX       2147483647    /* maximum (signed) int value */


MathTools::MathTools(void)
{
}

MathTools::~MathTools(void)
{
}


/**********************************************************************
 *
 * calculate mean value of y, return residuals (if resid!=0) or center
 *
 * the part about stratum is removed from the original code of CNVtools
 * 
 **********************************************************************/

int MathTools::wcenter(double *y, int n, double *weight, int resid, double *ynew) 
{
  int i = 0, s=0, empty = 0;
  double swt=0.0, swy=0.0, wi, epsilon=1e-8;
  
  if (weight) {
    for (i=0; i<n; i++) {
      wi   = weight[i];
      swt += wi;
      swy += wi*y[i];
    }
  }else {
    for (i=0; i<n; i++) {
      swy += y[i];
    }
    swt = (double) n;
  }
  swy /= swt;
  
  if (swt>epsilon) {
    if(resid){
      for (i=0; i<n; i++){ ynew[i] = y[i] - swy; }
    }else {
      for (i=0; i<n; i++){ ynew[i] = swy; }
    }
  }else {
    empty = 1;
  }
  
  return(empty);
}

/**********************************************************************
 *
 * Replace y by residual from (weighted) regression through the origin
 * 
 **********************************************************************/

int MathTools::wresid(double *y, int n, double *weight, double *x, 
           double *ynew, double *beta) 
{
  double swxy, swxx, wi, xi, wx;
  int i;
  
  swxy = swxx = 0.0;
  if (weight) {
    for (i=0; i<n; i++) {
      wi = weight[i];
      xi = x[i];
      wx = wi*xi;
      swxy += wx*y[i];
      swxx += wx*xi;
    }
  } else {
    for (i=0; i<n; i++) {
      xi = x[i];
      swxy += xi*y[i];
      swxx += xi*xi;
    }
  }
  
  if (swxx>0) {
    swxy /= swxx;
    *beta = swxy;
    for (i=0; i<n; i++) {
      if (weight[i] > 0.0) {
        ynew[i] = y[i] - swxy*x[i];
      } else {
        ynew[i] = y[i];
      }
    }
    return(n);
  }else{ 
    return(0);
  }
}

/**********************************************************************
 *
 * Weighted sum of squares
 * 
 **********************************************************************/

double MathTools::wssq(double *y, int n, double *weights) 
{  
  double res = 0.0, yi;
  int i;
  
  if (weights) {
    for (i=0; i<n; i++) {
      yi = y[i];
      res += weights[i]*yi*yi;
    }
  }else {
    for (i=0; i<n; i++) {
      yi = y[i];
      res += yi*yi;
    }
  }
  
  return(res);
}

/**********************************************************************
 *
 * utilit functions for glm_fit
 *
 * Originally from R/family and glm_test.cpp of R/CNVTools 
 *
 * Codes invovled strata are removed
 * 
 **********************************************************************/

/* 
 
 Variance function
 
 family:
 1    Binomial
 2    Poisson
 3    Gaussian
 4    Gamma
 5    Negagtive Binomial
 */

double MathTools::varfun(int family, double mu, double phi){
  switch (family) {
    case 1: return((mu*(1.0-mu)));  /* Binomial */
    case 2: return(mu);             /* Poisson */
    case 3: return(1.0);            /* Gaussian */
    case 4: return(mu*mu);          /* Gamma */
    case 5: return(mu + mu*mu*phi); /* Negative Binomial */
    default: return(0.0);
  }
}

/* Valid values for fitted value, mu. 
 
 If, during iteration, an invalid value is returned, the case is omitted 
 
 */

int MathTools::muvalid(int family, double mu) {
  double minb = 0.0001, maxb = 0.9999, minp = 0.0001;
  double gammaMin = 0.001;
  switch (family) {
    case 1: return(mu>minb && mu<maxb);    /* Binomial */
    case 2: return(mu>minp);               /* Poisson  */
    case 4: return(mu>gammaMin);           /* Gamma    */
    case 5: return(mu>minp);               /* Negative Binomial */
    default: return(1);                    /* Gaussian */
  }
}

/* Link function
 
 Link
 1    Logit
 2    Log
 3    Identity
 4    Inverse
 
 Note that a canonical link shares the code of the corresponding family
 so that the test for canonical link is (link==family)
 
 */

double MathTools::linkfun(int link, double mu) {
  switch (link) {
    case 1: {     /* Logit */
      if (mu == 1) return  HUGE_VAL;
      if (mu == 0) return -HUGE_VAL;
      return(log(mu/(1.0-mu)));
    }
    case 2: return(log(mu));             /* Log */
    case 3: return(mu);                  /* Identity */
    case 4: return(-1.0/mu);             /* Inverse */
    default: return 0.0;
  }
}

double MathTools::invlink(int link, double eta) {
  switch (link) {
    case 1: {  /* Logit */
      /*
       if(eta > 30){ return( 1 - 2e-16) }
       if(eta < -30){ return( 2e-16) }
       */
      if (eta ==  HUGE_VAL) return (1); 
      if (eta == -HUGE_VAL) return (0); 
      return(exp(eta)/(1.0+exp(eta))); 
    }
    case 2: return(exp(eta));                /* Log */
    case 3: return(eta);                     /* Identity */
    case 4: return(-1.0/eta);                /* Inverse */
    default: return(0.0);
  }
}

/* dlink = d eta / dmu */
double MathTools::dlink(int link, double mu) {
  switch (link) {
    case 1: return(-1.0/(mu*(1.0-mu)));    /* Logit */
    case 2: return(1.0/mu);                /* Log */
    case 3: return(1.0);                   /* Identity */
    case 4: return(-1.0/(mu*mu));          /* Inverse */
    default: return 0.0;
  }
}

/**********************************************************************
 *
 * utility functions for glm_fit
 *
 * following the implementation as R/stats/family.R 
 * 
 **********************************************************************/

void MathTools::initialize(int family, double* y, double* mu, int N, double* nTotal_binom) 
{
  int i;
  
  if (family==BINOMIAL) {         /* Binomial */
    for (i=0; i<N; i++) {
      if (y[i] <0 || nTotal_binom[i] < 0) {
        printf("negative values not allowed for the Binomial family");
      }
      if (y[i] > nTotal_binom[i]) {
        printf("# of success is larger than # of total trials in Binomial family");
      }
      mu[i] = (y[i] + 0.5)/(nTotal_binom[i]+1.0);
    }
  }else if (family==POISSON) {    /* Poisson */
    for (i=0; i<N; i++) {
      if (y[i] < 0) {
        printf("negative values not allowed for the Poisson family");
      }      
      mu[i] = y[i] + 0.1;
    }
  }else if (family==GAUSSIAN) {   /* Gaussian*/
    for (i=0; i<N; i++) {
      mu[i] = y[i];
    }
  }else if (family==GAMMA){       /* Gamma */
    for (i=0; i<N; i++) {
      if (y[i] <= 0) {
        printf("non-poistive values not allowed for the Gamma family");
      }      
      mu[i] = y[i] + 0.1;
    }
  }else if (family==NB){          /* Negagtive Binomial */
    for (i=0; i<N; i++) {
      if (y[i] < 0) {
        printf("negative values not allowed for the Negative Binomial family");
      }else if (y[i] < 0.01) {
        mu[i] = y[i] + 0.1667;
      }else{
        mu[i] = y[i];
      }
    }
  }else {
    printf("invaid family");
  }
}

/**********************************************************************
 *
 * glmFit
 *
 * Originally from glm_test.cpp of R/CNVTools 
 *
 * There are several channges:
 * 
 * (1)  For binomial distribution, allow two sets of inputs, 
 *      y and nTotal_binom, i.e. the number of succes and the total 
 *      number of trials
 * 
 * (2)  Codes invovled strata are removed
 * 
 * (3)  The way to handle the situation of no covaraite
 * 
 * (4)  The way to handel least square problem with Guassian family
 * 
 * (5)  Initial values of mu
 *
 * (6)  Add support for negative binomial distribution
 *
 
 Package: CNVtools
 Type: Package
 Title: A package to test genetic association with CNV data
 Version: 1.42.0
 Date: 2009-10-06
 Author: Chris Barnes <christopher.barnes@imperial.ac.uk> and 
 Vincent Plagnol <vincent.plagnol@cimr.cam.ac.uk>
 Maintainer: Chris Barnes <christopher.barnes@imperial.ac.uk>
 Description: This package is meant to facilitate the testing of Copy Number 
 Variant data for genetic association, typically in case-control studies.
 License: GPL-3
 Depends: survival
 biocViews: DNACopyNumber,GeneticVariability
 Packaged: 2010-04-23 00:36:41 UTC; biocbuild
 
 
 Input:
 
 family       GLM family (see below)
 link         Link function (see below)
 N            # units
 M            # X variables
 y            y-variable (N-vector)
 X            If M>0, N*M matrix of X variables
 maxit        Maximum number of iterations of IRLS algorithm
 conv         Proportional change in weighted sum of squares residuals to
              declare convergence
 init         If true (non-zero), the iteration starts from initial estimates 
              of fitted values (see below). This option has no effect if
              no iteration is required
 
 Output:
 
 rank         rank of X after regression on strata
 Xb           orthogonal basis for X space (N*rank matrix)
 fitted       fitted values 
 resid        working residuals (on linear predictor scale) (N-vector)
 weights      weights (N-vector)
 scale        scale factor (scalar)
 df_resid     residual degrees of freedom
 beta         regression coeficien of the last covariate

 Return
 
 0            convergence
 1            no convergence after maxit iterations
 
 **********************************************************************/

int MathTools::glmFit(int* familyR, int* linkR, int* dims, int* nIter,
           double *y, double *prior, double *offset, 
           double *X, double *nTotal_binom, double *convR, 
           int *rank, double *Xb, double *fitted, double *resid, 
           double *weights, double *phi, int* trace, 
           double *scale, int *df_resid, double* beta) 
{
  double epsilon = 1e-8;       /* Singularity threshold */
  int N, M, maxit, init, useOffset;
  int i = 0, j=0, Nu, dfr, irls;
  int empty = 0, invalid;
  int x_rank = 0, convg = 0, iter = 0;
  
  N = dims[0];
  M = dims[1];
  maxit = dims[2];
  init  = dims[3];
  useOffset = dims[4];
  
  int family  = *familyR;
  int link    = *linkR;
  double conv = *convR;
  
  double mu, ri, wi, pi, D, Vmu, wss, wss_last, ssx, ssr;
  double *xi, *xbi, *xbj, zmu, wsum;
  //double *z = (double *) malloc(N*sizeof(double));  //working y
  double *z = new double[N];
    
  if(family > 6 || family < 1){
    printf("family=%d, ", family);
    printf("Invalid family!\n");
  }
  
  /* Is iteration necessary? */
  irls = !((family==GAUSSIAN) && (link==IDENTITY));
  
  /* ----------------------------------------------------------------*
   * by default, initialize mu (fitted) by y itself, with neccesary 
   * modification, e.g., y + 0.1 to avoid log(0) for Poisson family
   * ----------------------------------------------------------------*/
  if (!init) {
    initialize(family, y, fitted, N, nTotal_binom);
  }
  
  /* ----------------------------------------------------------------*
   * Initialize wi (weights) and (standardized) residual 
   * In IRLS, z_i = eta_i + (y_i - mu_i)*(d_eta/d_mu)_i
   *
   * In the following code:
   * ri = (y_i - mu_i)*(d_eta/d_mu)_i
   * wi = pi*Vmu is the weight,  
   *
   * for those invlaid values, set their weight to be 0.
   * ----------------------------------------------------------------*/
  Nu      = 0;
  invalid = 0;
  wsum    = 0.0;
  
  for (i=0; i<N; i++) {
    mu = fitted[i];
    pi = prior[i];
    
    if (!muvalid(family, mu)) {
      invalid = 1;
      pi = 0.0;
    }
    
    if (!(pi)) { 
      wi = ri = 0.0; 
    }else {
      Nu ++;
      Vmu = varfun(family, mu, *phi);
      
      if (link == family) {
        ri = (y[i] - mu)/Vmu;
        wi = pi*Vmu;
      }else {
        D  = dlink(link, mu);
        ri = D*(y[i] - mu);
        wi = pi/(D*D*Vmu);
      }
    }
    
    weights[i] = wi;
    resid[i]   = ri;
	// not easy to find a reasonable epsilon
    //if (weights[i] < epsilon) weights[i] = 0.;
    
    wsum += weights[i];
  }
  
  /* ----------------------------------------------------------------*
   * If summation of all weights is too small, stop 
   * ----------------------------------------------------------------*/
  
  if (wsum < epsilon) {
    if(*trace)
      printf("  glmFit: empty model, summation of all weights are small!\n");
    return(0);
  }
  
  if(*trace > 3){
    printf("\n  glmFit: finish initialization, N=%d, M=%d, family=%d, and irls=%d\n", N, M, family, irls);
  }
  
  /* ----------------------------------------------------------------*
   * If M=0, there is only an intercept 
   * ----------------------------------------------------------------*/
  
  if (M == 0){
  /* IRLS algorithm */
  while(iter<maxit && !convg) {
		if (*trace > 5) {
		  printf("    glmFit: iteration %d: ", iter);
		}

	    
		for (i=0; i<N; i++) {
		  /**
		   * current estimate of eta + (y-mu)/gradient
		   *
		   * linkfun(link, fitted[i]) = eta
		   * resid[i] = (y-mu)/gradien
		   */
		  z[i] = linkfun(link, fitted[i]) + resid[i];
		}
	    
		if (useOffset) {
		  for (i=0; i<N; i++) z[i] -= offset[i];
		}
	            
		empty = wcenter(z, N, weights, 1, resid);  //removes the mean from z
	    
		if (empty == 1) {
		  if(*trace > 0)
			printf("  glmFit: empty model, summation of all weights are small!\n");
		  return(0);
		}
	    
		/**
		 * tries to fit the regression line (no intercept) to the residuals
		 */
	    
		xi  = X;
		xbi = Xb;
		x_rank = 0;
	    
		for (i=0; i<M; i++, xi+=N) {
		  wcenter(xi, N, weights, 1, xbi);
		  ssx = wssq(xbi, N, weights);
		  xbj = Xb;
	      
		  for (j=0; j<x_rank; j++, xbj+=N) wresid(xbi, N, weights, xbj, xbi, beta);
	      
		  ssr = wssq(xbi, N, weights);
	      
		  // printf("i=%d, ssx=%.2e, ssr=%.2e, x_rank=%d\n", i, ssx, ssr, x_rank);
	      
		  if (ssr/ssx > epsilon) {
			/**
			 * takes the residuals after fitting the regression line (no intercept) 
			 * to the mean value 
			 */
			wresid(resid, N, weights, xbi, resid, beta); 
			x_rank++;
			xbi+=N;
		  }
	      
		}
	    
		/* well, it is question whether we should give printf or just warning here */
		if(x_rank < M){ 
		  if(*trace > 1){
			printf("  glmFit: x_rank=%d, M=%d, X is not full rank\n", x_rank, M);
	        
			printf("  weights:\n");
		   // Rprint_ve(weights, 0, N-1);
	        
			printf("  resid:\n");
		   // Rprint_ve(resid, 0, N-1);
	        
			printf("  y:\n");
		   // Rprint_ve(y, 0, N-1);
		  }
	      
		  //return(0);
		}
	    
		wss = 0.0;
		Nu  = 0;
	    
		for (i=0; i<N; i++) {
	      
		  if (useOffset) {
			mu = invlink(link, z[i] + offset[i] - resid[i]);
		  }else {
			mu = invlink(link, z[i] - resid[i]);
		  }
	      
		  fitted[i] = mu;
		  pi = prior[i];
	      
		  if (!(pi && (weights[i]>0.0))) {
			wi = ri = 0.0;
		  } else {
	        
			Vmu = varfun(family, mu, *phi);
			Nu ++;
	        
			if (link == family) {
			  ri = (y[i] - mu)/Vmu;
			  wi = pi*Vmu;
			}else {
			  D = dlink(link, mu);
			  ri = D*(y[i] - mu);
			  wi = pi/(D*D*Vmu);
			}
			wss += wi*ri*ri;
	        
			weights[i] = wi;
			resid[i]   = ri;
			//if (weights[i] < 0.0001) weights[i] = 0.;
		  }
		}
    
		if(wss > 1.0/epsilon){
		  if(*trace > 1)
			printf("  glmFit: huge wss, indicting failt to fit the model!\n");
		  return(0);
		}
	    
		if (*trace > 5) {
		  printf(" wss=%.3f\n", iter, wss);
		}
    
		convg = (Nu<=0) || (iter && (fabs(wss-wss_last)/(wss_last + 0.1) < conv));
		wss_last = wss;
		iter ++;
	}

	/* assume there is an intercept */
    dfr = Nu - 1 - x_rank;
    
    if (family > 2) {
      *scale = wss_last/(dfr);
    }else{
      *scale = 1.0;
    }

  }
  
  /* ----------------------------------------------------------------*
   * If M>0, include covariates 
   * ----------------------------------------------------------------*/
  
  if (M > 0) {
    convg    = 0;
    iter     = 0;
    wss_last = 0.0;
    
    if (!irls) {
      /* Simple linear Gaussian case */
      xi  = X;
      xbi = Xb;
      x_rank = 0;
      
      for (i=0; i<M; i++, xi+=N) {
        wcenter(xi, N, weights, 1, xbi);
        ssx = wssq(xbi, N, weights);
        
        xbj = Xb;
        
        for (j=0; j<x_rank; j++, xbj+=N)  wresid(xbi, N, weights, xbj, xbi, beta);
        ssr = wssq(xbi, N, weights);
        
        if (ssr/ssx > epsilon) {
          wresid(resid, N, weights, xbi, resid, beta);
          x_rank++;
          xbi+=N;
        }
      }
      
      /* obtain the fitted values */
      for (i=0; i<N; i++) fitted[i] = y[i] - resid[i];
      
      wss_last = wssq(resid, N, weights);
      
    }else{
      
      /* IRLS algorithm */
      while(iter<maxit && !convg) {
        if (*trace > 5) {
          printf("    glmFit: iteration %d: ", iter);
        }

        
        for (i=0; i<N; i++) {
          /**
           * current estimate of eta + (y-mu)/gradient
           *
           * linkfun(link, fitted[i]) = eta
           * resid[i] = (y-mu)/gradien
           */
          z[i] = linkfun(link, fitted[i]) + resid[i];
        }
        
        if (useOffset) {
          for (i=0; i<N; i++) z[i] -= offset[i];
        }
                
        empty = wcenter(z, N, weights, 1, resid);  //removes the mean from z
        
        if (empty == 1) {
          if(*trace > 0)
            printf("  glmFit: empty model, summation of all weights are small!\n");
          return(0);
        }
        
        /**
         * tries to fit the regression line (no intercept) to the residuals
         */
        
        xi  = X;
        xbi = Xb;
        x_rank = 0;
        
        for (i=0; i<M; i++, xi+=N) {
          wcenter(xi, N, weights, 1, xbi);
          ssx = wssq(xbi, N, weights);
          xbj = Xb;
          
          for (j=0; j<x_rank; j++, xbj+=N) wresid(xbi, N, weights, xbj, xbi, beta);
          
          ssr = wssq(xbi, N, weights);
          
          // printf("i=%d, ssx=%.2e, ssr=%.2e, x_rank=%d\n", i, ssx, ssr, x_rank);
          
          if (ssr/ssx > epsilon) {
            /**
             * takes the residuals after fitting the regression line (no intercept) 
             * to the mean value 
             */
            wresid(resid, N, weights, xbi, resid, beta); 
            x_rank++;
            xbi+=N;
          }
          
        }
        
        /* well, it is question whether we should give printf or just warning here */
        if(x_rank < M){ 
          if(*trace > 1){
            printf("  glmFit: x_rank=%d, M=%d, X is not full rank\n", x_rank, M);
            
            printf("  weights:\n");
           // Rprint_ve(weights, 0, N-1);
            
            printf("  resid:\n");
           // Rprint_ve(resid, 0, N-1);
            
            printf("  y:\n");
           // Rprint_ve(y, 0, N-1);
          }
          
          return(0);
        }
        
        wss = 0.0;
        Nu  = 0;
        
        for (i=0; i<N; i++) {
          
          if (useOffset) {
            mu = invlink(link, z[i] + offset[i] - resid[i]);
          }else {
            mu = invlink(link, z[i] - resid[i]);
          }
          
          fitted[i] = mu;
          pi = prior[i];
          
          if (!(pi && (weights[i]>0.0))) {
            wi = ri = 0.0;
          } else {
            
            Vmu = varfun(family, mu, *phi);
            Nu ++;
            
            if (link == family) {
              ri = (y[i] - mu)/Vmu;
              wi = pi*Vmu;
            }else {
              D = dlink(link, mu);
              ri = D*(y[i] - mu);
              wi = pi/(D*D*Vmu);
            }
            wss += wi*ri*ri;
            
            weights[i] = wi;
            resid[i]   = ri;
            /*if (weights[i] < epsilon) weights[i] = 0.;*/
          }
        }
        
        if(wss > 1.0/epsilon){
          if(*trace > 1)
            printf("  glmFit: huge wss, indicting failt to fit the model!\n");
          return(0);
        }
        
        if (*trace > 5) {
          printf(" wss=%.3f\n", iter, wss);
        }
        
        convg = (Nu<=0) || (iter && (fabs(wss-wss_last)/(wss_last + 0.1) < conv));
        wss_last = wss;
        iter ++;
      }
    }
    
    /* assume there is an intercept */
    dfr = Nu - 1 - x_rank;
    
    if (family > 2) {
      *scale = wss_last/(dfr);
    }else{
      *scale = 1.0;
    }
  }
  
  *df_resid = dfr>0? dfr : 0;
  *rank     = x_rank;
  *nIter    = iter;
  
  //free(z);
  delete []z;
  
  return(irls && convg);
}


/**********************************************************************
 
 GLM score test 
 
 Input:
 
 P         Number of new explanatory variables to be added 
 Z         N*P matrix containing covariates
 
 For all other input arguments, see glm_fit, but note that M now coincides 
 with rank -- the number of columns in Xb
 
 Output:
 
 chi2  Score test 
 df    Degrees of freedom for asymptotic chi-squared distribution
 
 **********************************************************************/
            
void MathTools::glm_score_test(int* dims, double *Z, double *resid, 
                    double *weights, double *Xb, double* scaleR,
                    double *chi2, int *df) 
{
  
  int i = 0, j = 0, rank, N, M, P;
  
  double epsilon1 = 1.e-8;   /* First stage singularity test */
  double *Zi = Z;
  double *Zr, *Zri, ssz, ssr, *Zrj, Zrij, ws, wss, wz, test = 0.0;
  double *Xbj, beta;
  double scale = *scaleR;
  
  N = dims[0];
  M = dims[1];
  P = dims[2];
  
  /* Work array */
  Zr  = (double *)calloc(N*P,sizeof(double));
  Zri = Zr;
  
  /* Main algorithm */
  rank = 0;
  
  for (i=0; i<P; i++, Zi+=N) {
    /* Regress each column of Z on X basis */
    wcenter(Zi, N, weights, 1, Zri);
    ssz = wssq(Zri, N, weights);
    Xbj = Xb;
    
    for (j=0; j<M; j++, Xbj+=N){
      wresid(Zri, N, weights, Xbj, Zri, &beta);
    }
    
    ssr = wssq(Zri, N, weights);
    
    if (ssr/ssz > epsilon1) {     /* First singularity test */
      Zrj = Zr;
      for (j=0; j<rank; j++, Zrj+=N){
        wresid(Zri, N, weights, Zrj, Zri, &beta);
      }
      
      /* Sum and sum of squares */
      ws = 0.0, wss = 0.0;
      for (j=0; j<N; j++) {
        Zrij = Zri[j];
        wz   = weights[j]*Zrij;
        ws  += wz*resid[j];
        wss += Zrij*wz;
      }
      
      /* Second singularity test */
      if (wss/ssr > epsilon1) {
        test += ws*ws/(scale*wss);
        rank++;
        Zri  += N;
      }else {
        printf("colinearity in added covaraites Z\n");
      }
      
    }
  }
  
  *chi2 = test;
  *df   = rank;
  
  free(Zr);
}

/**********************************************************************
 *
 * glmNB
 *
 * Originally from MASS/negbin.R 
 * 
 
 Package: MASS
 Priority: recommended
 Version: 7.3-5
 Date: 2010-01-03
 Depends: R (>= 2.10.1), grDevices, graphics, stats, utils
 Suggests: lattice, nlme, survival
 Author: S original by Venables & Ripley. R port by Brian Ripley
 <ripley@stats.ox.ac.uk>, following earlier work by Kurt Hornik
 and Albrecht Gebhardt.
 Maintainer: Brian Ripley <ripley@stats.ox.ac.uk>
 Description: Functions and datasets to support Venables and Ripley,
 'Modern Applied Statistics with S' (4th edition).
 Title: Main Package of Venables and Ripley's MASS
 License: GPL-2 | GPL-3
 URL: http://www.stats.ox.ac.uk/pub/MASS4/
 LazyLoad: yes
 LazyData: yes
 Packaged: 2010-01-03 10:50:27 UTC; ripley
 Repository: CRAN
 Date/Publication: 2010-01-03 14:05:40
 
 **********************************************************************/

/**********************************************************************
 *
 * log likelihood of Poisson
 *
 **********************************************************************/

double MathTools::loglik_Poisson(int N, double* mu, double* y, double* w){
  int i;
  double yi, mui, logL = 0.0;
  
  for (i=0; i<N; i++) {
    yi  = y[i];
    mui = mu[i];
    
    logL += w[i]*(yi*log(mui) - mui - lgammafn(yi + 1.0));
  }
  
  return(logL);
}

/**********************************************************************
 *
 * log likelihood of negative binomial
 *
 **********************************************************************/

double MathTools::loglik_NB(int N, double phi, double* mu, double* y, double* w){
  int i;
  double logL1, logL, yi, mui;
  double th = 1.0/phi;
  double logL0 = th*log(th) - lgammafn(th);
  
  logL = 0.0;
  
  for (i=0; i<N; i++) {
    yi  = y[i];
    mui = mu[i];
    
    if (yi==0) {
      logL1  = th*log(th) - th*log(th + mui);
    }else {
	  double a = lgammafn(th + yi);
	  double b = lgammafn(yi + 1.0);
	  double c = yi*log(mui);
	  double d = (th + yi)*log(th + mui);
      //logL1  = lgammafn(th + yi) - lgammafn(yi + 1.0) + yi*log(mui) - (th + yi)*log(th + mui);
      logL1 = a-b+c-d;
	  logL1 += logL0;
    }

    logL += w[i]*logL1;
  }
  
  return(logL);
}

/**********************************************************************
 *
 * score_info
 *
 * score and Fisher information, i.e., the first and negative second 
 * derivative of likelihood of phi 
 *
 **********************************************************************/

void MathTools::score_info(int N, double theta, double* mu, double* y, double* w, 
                double* score, double* info)
{
  int i;
  double score1=0.0, info1=0.0;
  double wi, mui, yi, scorei, infoi, thMui;
  
  for (i=0; i<N; i++) {
    wi  = w[i];
    yi  = y[i];
    mui = mu[i];
    
    thMui   = theta + mui;
    scorei  = digamma(yi + theta) - digamma(theta) - (theta + yi)/thMui;
    score1 += wi*(scorei - log(thMui) + 1 + log(theta));
    
    infoi   = trigamma(theta) - trigamma(yi + theta) + (mui - yi)/(thMui*thMui);
    info1  += wi*(infoi + 1/thMui - 1/theta);
  }
  
  *score = score1;
  *info  = info1;
}


/**********************************************************************
 *
 * phi_ml
 *
 * MLE of phi (over-dispersion parameter), given mu 
 *
 * Actually we find MLE of 1/phi here and then take inverse
 *
 **********************************************************************/

int MathTools::phi_ml(double* y, double* mu, int N, double* weights, 
            int limit, double eps, double* phi, int initPhi, int trace)
{
  double theta0, del, tmp;
  double score=0.0;
  double info=0.0;
  int i, it=0;
  double minTheta = 1e-5;           
  double maxTheta = 1.0/minTheta;
  int tryPoisson = 0;
  int tryZINB = 0;
  int fail = 0;

  double sumWeights = 0.0;
  
  if(initPhi){
    theta0 = 1.0/(*phi);
  }else{
    theta0 = 0.0;
    for (i=0; i<N; i++) {
      tmp = y[i]/mu[i] - 1.0;
      theta0 += weights[i]*tmp*tmp;
	  sumWeights += weights[i];
    }
    theta0 = sumWeights/theta0;
  }
  
  it  = 0;
  del = 1.0;
  
  if(trace > 5) printf("  phi.ml: initial phi = %.2e\n", 1/theta0);
  
  while(it < limit && fabs(del) > eps) {
    score_info(N, theta0, mu, y, weights, &score, &info);
    del     = score/info;
    theta0 += del;
    it     += 1;
    
    if(trace > 5) printf("  phi.ml: iter %d, phi = %.2e\n", it,  1/theta0);
    
    if (theta0 > maxTheta) {
      theta0 = maxTheta;
      if(trace > 3)
        printf("    Estimate of phi is truncated at %.2e, no overDispersion?\n", maxTheta);
      tryPoisson = 1;
      break;
    }
    
    if(theta0 < minTheta) {
      theta0 = minTheta;
      if(trace > 3)
        printf("    Estimate of phi is truncated at %.2e, too much overDispersion?\n", minTheta);
      tryZINB = 1;
      break;
    }
    
  }
  
  if(it == limit) {
    fail = 1;
    if(trace > 1)
      printf("  phi.ml: iteration limit reached in phi_ml\n");
  }
  
  *phi = 1/theta0;
  
  return(tryPoisson + 2*tryZINB + 4*fail);
}


/**********************************************************************
 *
 * main function of glmNB
 *
 **********************************************************************/

int MathTools::glmNB(int *dims, int *nIter, double *y, double *prior, 
          int *linkR, double *offset, double *X, double *convR, 
          int *rank, double *Xb, double *fitted, double *resid, 
          double *weights, double *phi, double *scale, 
          int *df_resid, int* family, double *twologlik, 
          double *scoreTestP, int *trace, double *beta)
{
  int N, M, maxit, init, useOffset;
  int i, succeed = 0, iter = 0;
  double conv = *convR, initPhi;
  int fam0; 
  int cv = 0;         /* convergence indicator      */
  double nTotal_binom=0.0;  /* only used for binomial link, NOT useful here */
  double del, Lm, Lm0, phi0;
  double score, scoreNum, scoreDen, scorePval, yi, mui;
  
  /* convergence indicator for phi 
   * if phi = 0, NB model is OK. 
   * if phi = 1, suggest we need to use Poisson
   * if phi = 2, suggest we need to use ZINB
   */
  int cvPhi;
  
  N = dims[0];
  M = dims[1];
  maxit = dims[2];
  init  = dims[3];
  useOffset = dims[4];
  
  
  if(*trace > 3) 
    printf("\n  glmNB: N=%d, M=%d, maxit=%d, init=%d, useOffset=%d\n", 
            N, M, maxit, init, useOffset);

  /**
   * if there is no initial values, ignore the parameter *family 
   * and start with a Poission model
   */
  if(!init){
    /* Initial fit */
    
    fam0 = POISSON;
    
    cv = glmFit(&fam0, linkR, dims, nIter, y, prior, offset, 
                X, &nTotal_binom, convR, rank, Xb, fitted, resid, 
                weights, phi, trace, scale, df_resid, beta);
    
    if (cv==0) {
      if(*trace){
        printf("\n  glmNB: fail to converge in initial glmFit using Poission regression\n");
      }
      return(0);
    }

    /* test for overdispersion by Dean's Score test */
    scoreNum = 0.0;
    scoreDen = 0.0;
    
    for (i=0; i<N; i++) {
      yi  = y[i];
      mui = fitted[i];
      scoreNum += (yi - mui)*(yi - mui) - yi;
      scoreDen += mui*mui;
    }
    
    score = scoreNum/sqrt(2.0*scoreDen);
    
    /**
     * double pnorm(double x, double mu, double sigma, int lower_tail, int give_log);
     */
    scorePval = pnorm(score, 0.0, 1.0, 0, 0);
    
    if(*trace > 3) 
      printf("\n  overdispersion score = %.2e, p-value = %.2e\n\n", score, scorePval);

    if(scorePval > *scoreTestP){
      *family    = POISSON;
      Lm         = loglik_Poisson(N, fitted, y, prior);
      *twologlik = 2.0*Lm;
      *phi       = 0.0;
      return(1);
    }else {
      fam0    = NB;
      *family = NB;
    }
    
    /**
     * calculate phi by MLE, without initial values of phi
     */
    cvPhi = phi_ml(y, fitted, N, prior, maxit, conv, phi, 0, *trace);
    
    if(cvPhi==0){
      if(*trace > 3) 
        printf("\n  with mu estimated by Poission model, initial value for phi: %e\n", *phi);
    }else if (cvPhi==1){
      if(*trace > 3) 
        printf("\n  Have to choose Poisson model due to small phi\n");

      *family    = POISSON;
      Lm         = loglik_Poisson(N, fitted, y, prior);
      *twologlik = 2.0*Lm;
      *phi       = 0.0;

      return(1);
    }else if(cvPhi==2){
      if(*trace > 1) 
        printf("\n  The overdispersion parameter is too large: phi=%e\n", *phi);
      
      *family = ZINB;
      //return(0);
    }else { /* estimation of phi fail to converge */
      if(*trace > 1) 
        printf("\n  Estimation of phi fail to converge: phi=%e\n", *phi);

      //return(0);
    }
    
  }else {
    /**
     * if there is initial values, 
     * there must be both initial fitted values and inital phi
     */
    fam0 = *family;
    
    cv = glmFit(&fam0, linkR, dims, nIter, y, prior, offset, 
                X, &nTotal_binom, convR, rank, Xb, fitted, resid, 
                weights, phi, trace, scale, df_resid, beta);
    
    if (cv==0) {
      if(*trace > 1){
        printf("\n  glmNB: fail to converge using initial values, fam0=%d\n", fam0);
      }
      return(0);
    }
    
    /**
     * no need to go further if this is a Poission regression
     */
    
    if(fam0 == POISSON){
      Lm         = loglik_Poisson(N, fitted, y, prior);
      *twologlik = 2.0*Lm;
      *phi       = 0.0;
      return(1);
    }
    
    cvPhi = phi_ml(y, fitted, N, prior, maxit, conv, phi, 0, *trace);
    if(cvPhi==0){
      if(*trace > 3) 
        printf("\n  glmNB: given inital mu, initial value for phi: %e\n", *phi);
    }else {
      if(*trace > 1) 
        printf("\n  glmNB: fail to estimate phi\n");
      //return(0);
    }

  }
  
  del  = 1.0;
  Lm   = loglik_NB(N, *phi, fitted, y, prior);
  Lm0  = Lm + 1.0;

  while (iter < maxit && fabs(Lm0 - Lm) + fabs(del) > conv) {
    
    dims[3] = 1; /* use initial values */
    
    cv = glmFit(&fam0, linkR, dims, nIter, y, prior, offset, 
                X, &nTotal_binom, convR, rank, Xb, fitted, resid, 
                weights, phi, trace, scale, df_resid, beta);
    
    if (cv==0) { 
      if(*trace > 1) 
        printf("\n  glmNB: fail to converge in glmFit\n");
      
      return(0);
    }
    
    phi0  = *phi;
    cvPhi = phi_ml(y, fitted, N, prior, maxit, conv, phi, 1, *trace);
    
    if(cvPhi==0){
      if(*trace > 3) 
        printf("\n  finish phi_ml, cvPhi=%d, phi=%e\n", cvPhi, phi);
    }else if(cvPhi==1){
      if(*trace > 1) 
        printf("  glmNB: trying to fit a NB model while there is no overdispersion\n");
      //return(0);
    }else if(cvPhi==2){
      if(*trace > 1) 
        printf("  glmNB: too large overdispersion\n");
      //return(0);
    }else {
      if(*trace > 1) 
        printf("  glmNB: fail to converge in phi_ml\n");
      //return(0);
    }
    
    del = phi0 - *phi;
    Lm0 = Lm;
    Lm  = loglik_NB(N, *phi, fitted, y, prior);
        
    if (*trace > 3) {
      printf("\n  Phi(%d) = %.2e, logLik = %.2e\n\n", iter, *phi, Lm);
    }
    
    iter++;
  }

  if(iter == maxit) {
    if (*trace) {
      printf("\n  glmNB: Alternation limit reached: iter=%d\n", iter);
    }
  }else {
    succeed = 1;
  }

  *twologlik = 2.0*Lm;
  return(succeed);
}


/**********************************************************************
 *
 * logL_b_TReC: log likelihood of TReC model regarding to b
 *
 **********************************************************************/

//double logL_b_TReC(double b1, int N, int fam, double b0, double phi,  
//                   double* y, double* x, double* mu, double* mu1)
//{
//  int i;
//  double logL, xi, cst0, cst1;
//  
//  cst0 = exp(b1 - b0);
//  cst1 = (1.0 + exp(b1))/(1.0 + exp(b0));
//    
//  for (i=0; i<N; i++) {
//    xi = x[i];
//    
//    if (fabs(xi - 0.0) < 0.01) {
//      mu1[i] = mu[i];
//    }else if (fabs(xi - 2.0) < 0.01) {
//      mu1[i] = mu[i]*cst0;
//    }else if (fabs(xi - 1.0) < 0.01) {
//      mu1[i] = mu[i]*cst1;
//    }else {
//      printf("invalid genotype\n");
//    }
//  }
//  
//  logL = 0.0;
//  if (fam == POISSON) {
//    logL = loglik_Poisson(N, mu1, y);
//  }else if (fam == NB) {
//    logL = loglik_NB(N, phi, mu1, y);
//  }else {
//    printf("wrong family\n");
//  }
//  
//  return(logL);
//}

/**********************************************************************
 *
 * grad_b_TReC: 1st and 2nd gradient of logL_b_TReC
 *
 **********************************************************************/

//*void grad_b_TReC(double b1, int N, int fam, double b0, double phi, 
//                 double* y, double* x, double* mu, double* gr)
//{
//  int i;
//  double ti, mui, xi, grad1, grad2, df_dmu1, df_dmu2, dmu_db;
//  double cst0, cst1, cst2;
//  
//  cst0 = exp(b1 - b0);
//  cst1 = (1.0 + exp(b1))/(1.0 + exp(b0));
//  cst2 = exp(b1)/(1.0 + exp(b1));
//  
//  grad1 = grad2 = 0.0;
//  
//  for (i=0; i<N; i++) {
//    xi = x[i];
//    ti = y[i];
//    
//    if (fabs(xi - 0.0) < 0.01) {
//      continue;
//    }else if (fabs(xi - 2.0) < 0.01) {
//      mui    = mu[i]*cst0;
//      dmu_db = mui;
//    }else if (fabs(xi - 1.0) < 0.01) {
//      mui    = mu[i]*cst1;
//      dmu_db = mui*cst2;
//    }else {
//      printf("invalid genotype\n");
//    }
//    
//    df_dmu1 = ti/mui - (1.0 + phi*ti)/(1.0 + phi*mui);        
//    grad1  += df_dmu1*dmu_db;
//    
//    df_dmu2 = -ti/(mui*mui) + phi*(1.0 + phi*ti)/(1.0 + phi*mui)/(1.0 + phi*mui);        
//    grad2  += (df_dmu2*dmu_db*dmu_db + df_dmu1*dmu_db);
//  }
//  
//  gr[0] = grad1;
//  gr[1] = grad2;
//}*/

/**********************************************************************
 * b_TReC_ml: find the MLE of b in TReC model by 
 * a Newton-Raphson algrorithm
 **********************************************************************/

//*void b_TReC_ml(double* b_xj, int N, int fam, double b0, double phi, 
//               double* y,double* x, double* mu, int limit, 
//               double eps, int trace, int* failR)
//{
//  int it, fail;
//  double logL, b1, del, gr[2];
//  double min_b = -1e5;
//  double max_b =  1e5;
//  
//  it   = 0;
//  del  = 1.0;
//  fail = 0;
//  b1   = b0;
//  gr[0] = 1.0;
//  gr[1] = 1.0;
//    
//  if(trace > 5){
//    printf("b_TReC_ml: limit=%d, eps=%e\n", limit, eps);
//  }
//  
//  while (it < limit && fabs(del) > eps) {
//    
//    grad_b_TReC(b1, N, fam, b0, phi, y, x, mu, gr);
//    
//    del = gr[0]/gr[1];
//    b1 -= del;
//    it += 1;
//    
//    if (it > 10 && fabs(del) > 0.1) {
//      fail = 1;
//      break;
//    }
//    
//    if(trace > 5){
//      printf("      b_TReC_ml: it=%d, b1=%e, del=%e, gr[0]=%e, gr[1]=%e\n", 
//              it, b1, del, gr[0], gr[1]);
//    }
//    
//    if (b1 > max_b) {
//      b1 = max_b;
//      if(trace > 3)
//        printf("    Estimate of b_xj is truncated at %.2e\n", max_b);
//      
//      break;
//    }
//    
//    if(b1 < min_b) {
//      b1 = min_b;
//      if(trace > 3)
//        printf("    Estimate of b_xj is truncated at %.2e\n", min_b);
//      
//      break;
//    }  
//  }
//  
//  if (it == limit) {
//    // fail = 1;
//    
//    if(trace > 1)
//      printf("  b_TReC_ml: iteration limit reached in b_TReC_ml\n");
//  }
//  
//  if (fabs(gr[0]) > 0.01) {
//    fail = 1;
//    
//    if(trace > 1)
//      printf("  b_TReC_ml: 1st derivative = %.2e\n", gr[0]);
//  }
//  
//  if (gr[1] > 1e-16) {
//    fail = 1;
//    
//    if(trace > 1)
//      printf("  b_TReC_ml: 2st derivative = %.2e\n", gr[1]);
//  }
//    
//  *b_xj  = b1;
//  *failR = fail;
//}*/

/**********************************************************************
 *
 * main function of glmNBlog
 * 
 * the last covariate does not follow a log linear model 
 * but a linear model. So iterately estimate all other co
 *
 **********************************************************************/

//int MathTools::glmNBlog(int *dimsNew, int *nIter, double *pY, double *z, 
//             int *linkR, double *offset, double *pX, double *conv, 
//             double *convGLM, int *rank, double *Xb, double *fitted, 
//             double *resid, double *weights, double *phi, double *scale, 
//             int *df_resid, int* family, double *twoLL_trec, 
//             double *scoreTestP, int *trace, double *beta,
//             double *fitted2, double *offsetN)
//{
//  int i, g, k, k1, k2, k3, k4, df1, convBase, convSNPj, succeed, fail;
//  int N     = dimsNew[0];
//  int nX    = dimsNew[1];
//  int maxit = dimsNew[2];
//  int init  = dimsNew[3]; /* whetehr to use initial values */
//  int useOffset = dimsNew[4];
//  double gr[2];
//  
//  double *pZ, pZk, twoLL_trec0, twoLL_trec1, nless5;
//  double bxj_old, bxj_new, btmp, phi_old, paraDiff;
//  
//  succeed = 1;
//  for (k=0; k<N; k++) { fitted2[k] = fitted[k]; }
//  
//  pZ = pX + N*nX;
//  
//  phi_old = *phi;
//  bxj_old = bxj_new = 0.0;
//
//  /*************************************************
//   * must use initial values
//   *************************************************/
//  
//  if(!init){
//    printf("we expect to use initial values\n");  
//  }
//  
//  /*************************************************
//   * phi = 0 for Poisson distribution
//   *************************************************/
//  
//  if(*family==POISSON && phi_old > 0){
//    printf("phi=%e while family is Poisson\n", phi_old);  
//  }
//  
//  /*************************************************
//   * calculate old likelihood       
//   *************************************************/
//  if (*family==NB) {
//    twoLL_trec0 = 2.0*loglik_NB(N, phi_old, fitted2, pY);
//  }else if (*family==POISSON) {
//    twoLL_trec0 = 2.0*loglik_Poisson(N, fitted2, pY);
//  }else {
//    printf("invalid family\n");
//  }
//  
//  for(g=0; g<maxit; g++){
//    
//    paraDiff = 0.0;
//
//    /* --------------------------------------------------------
//     * 1. Given theta, phi, b_0, b_k, and b_u, estimate b_{x_j}
//     * -------------------------------------------------------*/
//        
//    fail = 0;
//
//    b_TReC_ml(&bxj_new, N, *family, bxj_old, phi_old, pY, pZ, 
//              fitted2, maxit, *convGLM, *trace, &fail);
//    
//    if(fail){
//      if(*trace){
//        printf("\n  g=%d, fail b_TReC_ml in glmNBlog, fail=%d\n", g, fail);
//      }
//      
//      succeed=0;
//      break;
//    }
//    
//    if (paraDiff < fabs(bxj_new - bxj_old)) paraDiff = fabs(bxj_new - bxj_old);
//    
//    /**********************************
//     * calculate new likelihood       
//     **********************************/
//    if (*family==NB) {
//      twoLL_trec1 = 2.0*loglik_NB(N, phi_old, fitted2, pY);
//    }else if (*family==POISSON) {
//      twoLL_trec1 = 2.0*loglik_Poisson(N, fitted2, pY);
//    }else {
//      printf("invalid family\n");
//    }
//    
//    // printf("\n===--->>> family=%d, twoLL_trec1=%e\n", *family, twoLL_trec1);
//
//    if (*trace > 1) {
//      printf("\n  g=%d, phi=%e, b_old=%e, b_new=%e\n", g, phi_old, bxj_old, bxj_new);
//      printf("  twoLL_trec(old, new)=(%e, %e)\n", twoLL_trec0, twoLL_trec1);
//    }
//    
//    if((twoLL_trec1 - twoLL_trec0)/(fabs(twoLL_trec0) + 1.0) < -1e-4){
//      
//      nless5 = 0.0;
//      for (k=0; k<N; k++) { 
//        if (pY[k] < 5) nless5 += 1.0;
//      }
//      
//      printf("  g=%d, nless5=%.1f, twoLL(old, new)=(%e, %e)\n", 
//              g, nless5, twoLL_trec0, twoLL_trec1);
//
//      if (nless5 > 0.75*N) {
//        succeed=0;
//        break;
//      }
//      
//      printf("\n  likelihood decreases during update of b_xj in glmNBlog\n");
//    }
//    
//    bxj_old = bxj_new;
//    
//    /* --------------------------------------------------------
//     * 2. Given b_{x_j} and phi, estimate b_0, b_k, b_u, phi
//     * in fact, we do not estimate b_0, b_k, and b_u, 
//     * instead we only need to estimate mu, or fitted2.
//     * -------------------------------------------------------*/
//    
//    twoLL_trec0 = twoLL_trec1;
//    
//    dimsNew[1] = nX;
//    dimsNew[3] = 1; /* use initial values */
//    dimsNew[4] = 1; /* use offset         */
//    
//    for (k=0; k<N; k++){
//      pZk = pZ[k];
//      
//      if (fabs(pZk) < 0.01) {
//        offsetN[k] = offset[k];
//      }else if (fabs(pZk - 1.0) < 0.01) {
//        offsetN[k] = offset[k] + log(0.5*(1 + exp(bxj_old)));
//      }else if (fabs(pZk - 2.0) < 0.01) {
//        offsetN[k] = offset[k] + bxj_old;
//      }else {
//        printf("invalid value of genotype.");
//      }
//    }
//    
//    convSNPj = glmNB(dimsNew, nIter, pY, z, linkR, offsetN, pX, convGLM, 
//                     rank, Xb, fitted2, resid, weights, phi, scale, &df1, 
//                     family, &twoLL_trec1, scoreTestP, trace, &btmp);
//
//    // printf("\n===--->>> family=%d, twoLL_trec1=%e\n", *family, twoLL_trec1);
//
//    if (convSNPj == 0) {
//      if(*trace > 1)
//        printf("\n  g=%d, fail to estimate phi in glmNBlog\n", g);
//      
//      succeed = 0;
//      break;
//    }
//    
//    /* need to substract the extra degree for freedom in the offset */
//    df1 -= 1;
//
//    if (paraDiff < fabs(*phi - phi_old)) paraDiff = fabs(*phi - phi_old);
//    
//    /**********************************
//     * compare likelihood       
//     **********************************/
//    
//    if((twoLL_trec1 - twoLL_trec0)/(fabs(twoLL_trec0) + 1.0) < -1e-4){
//      
//      nless5 = 0.0;
//      for (k=0; k<N; k++) { 
//        if (pY[k] < 5) nless5 += 1.0;
//      }
//      
//      printf("  g=%d, nless5=%.1f, twoLL(old, new)=(%e, %e)\n", 
//              g, nless5, twoLL_trec0, twoLL_trec1);
//
//      if (nless5 > 0.75*N) {
//        succeed=0;
//        break;
//      }
//      
//      warning("\n  likelihood decreases during update of phi in glmNBlog\n");
//      
//      succeed = 0;
//      break;
//    }
//    
//    phi_old = *phi;
//    twoLL_trec0 = twoLL_trec1;
//    
//    if (*trace > 3){      
//      printf("\n  g=%d, parDiff=%e", g, paraDiff);
//      printf("\n  bxj=%e, phi=%e, twoLL_trec1=%e\n", bxj_old, phi_old, twoLL_trec1);
//    }
//    
//    /* --------------------------------------------------------
//     * 3. check convergence
//     * -------------------------------------------------------*/
//    
//    if(paraDiff < *conv){
//      
//      if (*trace > 1){
//        printf("\n  converged in glmNBlog\n");
//        printf("  bxj=%e, phi=%e, df1=%d, twoLL1=%e\n", bxj_old, phi_old, df1, twoLL_trec1);
//      }
//      
//      break;
//    }
//  }
//  
//  if (g >= maxit) {
//    if (*trace){
//      printf("\n  g=%d, paraDiff=%.3e, reach max iteration in glmNBlog", g, paraDiff);
//    }
//    succeed = 0;
//  }
//  
//  if (succeed) {
//    *beta = bxj_old;
//    *twoLL_trec = twoLL_trec1;
//    *df_resid   = df1;
//    for (i=0; i<N; i++) fitted[i] = fitted2[i];    
//  }else {
//    *beta       = 0.0;
//    *twoLL_trec = -1e6;
//    *df_resid   = 0;
//  }
//
//  return(succeed);
//}

double MathTools::lgammafn(double x)
{
    double x0,x2,xp,gl,gl0;
    int n,k;
    static double a[] = {
        8.333333333333333e-02,
       -2.777777777777778e-03,
        7.936507936507937e-04,
       -5.952380952380952e-04,
        8.417508417508418e-04,
       -1.917526917526918e-03,
        6.410256410256410e-03,
       -2.955065359477124e-02,
        1.796443723688307e-01,
       -1.39243221690590};
    
    x0 = x;
    if (x <= 0.0) return 1e308;
    else if ((x == 1.0) || (x == 2.0)) return 0.0;
    else if (x <= 7.0) {
        n = (int)(7-x);
        x0 = x+n;
    }
    x2 = 1.0/(x0*x0);
    xp = 2.0*_PI;
    gl0 = a[9];
    for (k=8;k>=0;k--) {
        gl0 = gl0*x2 + a[k];
    }
    gl = gl0/x0+0.5*log(xp)+(x0-0.5)*log(x0)-x0;
    if (x <= 7.0) {
        for (k=1;k<=n;k++) {
            gl -= log(x0-1.0);
            x0 -= 1.0;
        }
    }
    return gl;
}

void MathTools::dpsifn(double x, int n, int kode, int m, double *ans, int *nz, int *ierr)
{
    const static double bvalues[] = {	/* Bernoulli Numbers */
	 1.00000000000000000e+00,
	-5.00000000000000000e-01,
	 1.66666666666666667e-01,
	-3.33333333333333333e-02,
	 2.38095238095238095e-02,
	-3.33333333333333333e-02,
	 7.57575757575757576e-02,
	-2.53113553113553114e-01,
	 1.16666666666666667e+00,
	-7.09215686274509804e+00,
	 5.49711779448621554e+01,
	-5.29124242424242424e+02,
	 6.19212318840579710e+03,
	-8.65802531135531136e+04,
	 1.42551716666666667e+06,
	-2.72982310678160920e+07,
	 6.01580873900642368e+08,
	-1.51163157670921569e+10,
	 4.29614643061166667e+11,
	-1.37116552050883328e+13,
	 4.88332318973593167e+14,
	-1.92965793419400681e+16
    };

    int i, j, k, mm, mx, nn, np, nx, fn;
    double arg, den, elim, eps, fln, fx, rln, rxsq,
	r1m4, r1m5, s, slope, t, ta, tk, tol, tols, tss, tst,
	tt, t1, t2, wdtol, xdmln, xdmy, xinc, xln = 0.0 /* -Wall */,
	xm, xmin, xq, yint;
    double trm[23], trmr[n_max + 1];

    *ierr = 0;
    if (n < 0 || kode < 1 || kode > 2 || m < 1) {
	*ierr = 1;
	return;
    }
    if (x <= 0.) {
	/* use	Abramowitz & Stegun 6.4.7 "Reflection Formula"
	 *	psi(k, x) = (-1)^k psi(k, 1-x)	-  pi^{n+1} (d/dx)^n cot(x)
	 */
	if (x == (long)x) {
	    /* non-positive integer : +Inf or NaN depends on n */
	    for(j=0; j < m; j++) /* k = j + n : */
		ans[j] = ((j+n)%2==0) ? ML_POSINF : ML_NAN;
	    return;
	}
	/* This could cancel badly */
	dpsifn(1. - x, n, /*kode = */ 1, m, ans, nz, ierr);
	/* ans[j] == (-1)^(k+1) / gamma(k+1) * psi(k, 1 - x)
	 *	     for j = 0:(m-1) ,	k = n + j
	 */

	/* Cheat for now: only work for	 m = 1, n in {0,1,2,3} : */
	if(m > 1 || n > 3) {/* doesn't happen for digamma() .. pentagamma() */
	    /* not yet implemented */
	    *ierr = 4; return;
	}
	x *= _PI; /* pi * x */
	if (n == 0)
	    tt = cos(x)/sin(x);
	else if (n == 1)
	    tt = -1/pow(sin(x),2);
	else if (n == 2)
	    tt = 2*cos(x)/pow(sin(x),3);
	else if (n == 3)
	    tt = -2*(2*pow(cos(x),2) + 1)/pow(sin(x),4);
	else /* can not happen! */
	    tt = ML_NAN;
	/* end cheat */

	s = (n % 2) ? -1. : 1.;/* s = (-1)^n */
	/* t := pi^(n+1) * d_n(x) / gamma(n+1)	, where
	 *		   d_n(x) := (d/dx)^n cot(x)*/
	t1 = t2 = s = 1.;
	for(k=0, j=k-n; j < m; k++, j++, s = -s) {
	    /* k == n+j , s = (-1)^k */
	    t1 *= _PI;/* t1 == pi^(k+1) */
	    if(k >= 2)
		t2 *= k;/* t2 == k! == gamma(k+1) */
	    if(j >= 0) /* by cheat above,  tt === d_k(x) */
		ans[j] = s*(ans[j] + t1/t2 * tt);
	}
	if (n == 0 && kode == 2) /* unused from R, but "wrong": xln === 0 :*/
	    ans[0] += xln;
	return;
    } /* x <= 0 */

    /* else :  x > 0 */
    *nz = 0;
    xln = log(x);
    if(kode == 1 && m == 1) {/* the R case  ---  for very large x: */
	double lrg = 1/(2. * DBL_EPSILON);
	if(n == 0 && x * xln > lrg) {
	    ans[0] = -xln;
	    return;
	}
	else if(n >= 1 && x > n * lrg) {
	    ans[0] = exp(-n * xln)/n; /* == x^-n / n  ==  1/(n * x^n) */
	    return;
	}
    }
    mm = m;
    nx = __min(-Rf_i1mach(15), Rf_i1mach(16));/* = 1021 */
    r1m5 = Rf_d1mach(5);
    r1m4 = Rf_d1mach(4) * 0.5;
    wdtol = __max(r1m4, 0.5e-18); /* 1.11e-16 */

    /* elim = approximate exponential over and underflow limit */
    elim = 2.302 * (nx * r1m5 - 3.0);/* = 700.6174... */
    for(;;) {
	nn = n + mm - 1;
	fn = nn;
	t = (fn + 1) * xln;

	/* overflow and underflow test for small and large x */

	if (fabs(t) > elim) {
	    if (t <= 0.0) {
		*nz = 0;
		*ierr = 2;
		return;
	    }
	}
	else {
	    if (x < wdtol) {
		ans[0] = pow(x, -n-1.0);
		if (mm != 1) {
		    for(k = 1; k < mm ; k++)
			ans[k] = ans[k-1] / x;
		}
		if (n == 0 && kode == 2)
		    ans[0] += xln;
		return;
	    }

	    /* compute xmin and the number of terms of the series,  fln+1 */

	    rln = r1m5 * Rf_i1mach(14);
	    rln = __min(rln, 18.06);
	    fln = __max(rln, 3.0) - 3.0;
	    yint = 3.50 + 0.40 * fln;
	    slope = 0.21 + fln * (0.0006038 * fln + 0.008677);
	    xm = yint + slope * fn;
	    mx = (int)xm + 1;
	    xmin = mx;
	    if (n != 0) {
		xm = -2.302 * rln - __min(0.0, xln);
		arg = xm / n;
		arg = __min(0.0, arg);
		eps = exp(arg);
		xm = 1.0 - eps;
		if (fabs(arg) < 1.0e-3)
		    xm = -arg;
		fln = x * xm / eps;
		xm = xmin - x;
		if (xm > 7.0 && fln < 15.0)
		    break;
	    }
	    xdmy = x;
	    xdmln = xln;
	    xinc = 0.0;
	    if (x < xmin) {
		nx = (int)x;
		xinc = xmin - nx;
		xdmy = x + xinc;
		xdmln = log(xdmy);
	    }

	    /* generate w(n+mm-1, x) by the asymptotic expansion */

	    t = fn * xdmln;
	    t1 = xdmln + xdmln;
	    t2 = t + xdmln;
	    tk = __max(fabs(t), __max(fabs(t1), fabs(t2)));
	    if (tk <= elim) /* for all but large x */
		goto L10;
	}
	nz++; /* underflow */
	mm--;
	ans[mm] = 0.;
	if (mm == 0)
	    return;
    } /* end{for()} */
    nn = (int)fln + 1;
    np = n + 1;
    t1 = (n + 1) * xln;
    t = exp(-t1);
    s = t;
    den = x;
    for(i=1; i <= nn; i++) {
	den += 1.;
	trm[i] = pow(den, (double)-np);
	s += trm[i];
    }
    ans[0] = s;
    if (n == 0 && kode == 2)
	ans[0] = s + xln;

    if (mm != 1) { /* generate higher derivatives, j > n */

	tol = wdtol / 5.0;
	for(j = 1; j < mm; j++) {
	    t /= x;
	    s = t;
	    tols = t * tol;
	    den = x;
	    for(i=1; i <= nn; i++) {
		den += 1.;
		trm[i] /= den;
		s += trm[i];
		if (trm[i] < tols)
		    break;
	    }
	    ans[j] = s;
	}
    }
    return;

  L10:
    tss = exp(-t);
    tt = 0.5 / xdmy;
    t1 = tt;
    tst = wdtol * tt;
    if (nn != 0)
	t1 = tt + 1.0 / fn;
    rxsq = 1.0 / (xdmy * xdmy);
    ta = 0.5 * rxsq;
    t = (fn + 1) * ta;
    s = t * bvalues[2];
    if (fabs(s) >= tst) {
	tk = 2.0;
	for(k = 4; k <= 22; k++) {
	    t = t * ((tk + fn + 1)/(tk + 1.0))*((tk + fn)/(tk + 2.0)) * rxsq;
	    trm[k] = t * bvalues[k-1];
	    if (fabs(trm[k]) < tst)
		break;
	    s += trm[k];
	    tk += 2.;
	}
    }
    s = (s + t1) * tss;
    if (xinc != 0.0) {

	/* backward recur from xdmy to x */

	nx = (int)xinc;
	np = nn + 1;
	if (nx > n_max) {
	    *nz = 0;
	    *ierr = 3;
	    return;
	}
	else {
	    if (nn==0)
		goto L20;
	    xm = xinc - 1.0;
	    fx = x + xm;

	    /* this loop should not be changed. fx is accurate when x is small */
	    for(i = 1; i <= nx; i++) {
		trmr[i] = pow(fx, (double)-np);
		s += trmr[i];
		xm -= 1.;
		fx = x + xm;
	    }
	}
    }
    ans[mm-1] = s;
    if (fn == 0)
	goto L30;

    /* generate lower derivatives,  j < n+mm-1 */

    for(j = 2; j <= mm; j++) {
	fn--;
	tss *= xdmy;
	t1 = tt;
	if (fn!=0)
	    t1 = tt + 1.0 / fn;
	t = (fn + 1) * ta;
	s = t * bvalues[2];
	if (fabs(s) >= tst) {
	    tk = 4 + fn;
	    for(k=4; k <= 22; k++) {
		trm[k] = trm[k] * (fn + 1) / tk;
		if (fabs(trm[k]) < tst)
		    break;
		s += trm[k];
		tk += 2.;
	    }
	}
	s = (s + t1) * tss;
	if (xinc != 0.0) {
	    if (fn == 0)
		goto L20;
	    xm = xinc - 1.0;
	    fx = x + xm;
	    for(i=1 ; i<=nx ; i++) {
		trmr[i] = trmr[i] * fx;
		s += trmr[i];
		xm -= 1.;
		fx = x + xm;
	    }
	}
	ans[mm - j] = s;
	if (fn == 0)
	    goto L30;
    }
    return;

  L20:
    for(i = 1; i <= nx; i++)
	s += 1. / (x + (nx - i)); /* avoid disastrous cancellation, PR#13714 */

  L30:
    if (kode != 2) /* always */
	ans[0] = s - xdmln;
    else if (xdmy != x) {
	xq = xdmy / x;
	ans[0] = s - log(xq);
    }
    return;
} /* dpsifn() */

int MathTools::Rf_i1mach(int i)
{
    switch(i) {

    case  1: return 5;
    case  2: return 6;
    case  3: return 0;
    case  4: return 0;

    case  5: return CHAR_BIT * sizeof(int);
    case  6: return sizeof(int)/sizeof(char);

    case  7: return 2;
    case  8: return CHAR_BIT * sizeof(int) - 1;
    case  9: return INT_MAX;

    case 10: return FLT_RADIX;

    case 11: return FLT_MANT_DIG;
    case 12: return FLT_MIN_EXP;
    case 13: return FLT_MAX_EXP;

    case 14: return DBL_MANT_DIG;
    case 15: return DBL_MIN_EXP;
    case 16: return DBL_MAX_EXP;

    default: return 0;
    }
}

double MathTools::Rf_d1mach(int i)
{
    switch(i) {
    case 1: return DBL_MIN;
    case 2: return DBL_MAX;

    case 3: /* = FLT_RADIX  ^ - DBL_MANT_DIG
	      for IEEE:  = 2^-53 = 1.110223e-16 = .5*DBL_EPSILON */
	return 0.5*DBL_EPSILON;

    case 4: /* = FLT_RADIX  ^ (1- DBL_MANT_DIG) =
	      for IEEE:  = 2^-52 = DBL_EPSILON */
	return DBL_EPSILON;

    case 5: return _LOG10_2;

    default: return 0.0;
    }
}

# define ML_TREAT_psigam(_IERR_)	\
    if(_IERR_ != 0) 			\
	return ML_NAN

#define ISNAN(x) (isnan(x)!=0)


double MathTools::digamma(double x)
{
    double ans;
    int nz, ierr;
    if(isnan(x)) return x;
    dpsifn(x, 0, 1, 1, &ans, &nz, &ierr);
    ML_TREAT_psigam(ierr);
    return -ans;
}

double MathTools::trigamma(double x)
{
    double ans;
    int nz, ierr;
    if(isnan(x)) return x;
    dpsifn(x, 1, 1, 1, &ans, &nz, &ierr);
    ML_TREAT_psigam(ierr);
    return ans;
}

double MathTools::getCumulativeDensityNormal(const double x)
{
  const double c0 = 0.2316419;
  const double c1 = 1.330274429;
  const double c2 = 1.821255978;
  const double c3 = 1.781477937;
  const double c4 = 0.356563782;
  const double c5 = 0.319381530;
  const double c6 = 0.398942280401;
  const double negative = (x < 0 ? 1.0 : 0.0);
  const double xPos = (x < 0.0 ? -x : x);
  const double k = 1.0 / ( 1.0 + (c0 * xPos));
  const double y1 = (((((((c1*k-c2)*k)+c3)*k)-c4)*k)+c5)*k;
  const double y2 = 1.0 - (c6*exp(-0.5*xPos*xPos)*y1);
  return ((1.0-negative)*y2) + (negative*(1.0-y2));
}


double MathTools::pnorm(const double x,const double mean,const double stddev, const bool low_tail, const bool log_p)
{
  double v = getCumulativeDensityNormal((x - mean) / stddev);
  if (!low_tail)
	  v = 1-v;
  if (log_p)
	  v = log(v);
  return v;
}

void MathTools::max(double *v, int nl, int nh, double *val, int *idx)
{
    int i;
    *idx = nl;
    *val = v[nl];

    for(i=nl+1; i<nh; i++){
        if(v[i] > *val){
            *val = v[i];
            *idx = i;
        }
    }
}

/**********************************************************************
 *
 * logsumexp
 *
 * log(sum(exp(v)))
 *
 **********************************************************************/
double MathTools::logsumexp(double* v, const int RN)
{
	double lse = 0; // containing the return value
	int i, idx=0, N=RN;
	double res, val=0.0;

	if(N==0)
	{
		lse = ML_NEGINF;
	}
	else if(N==1)
	{
		lse = v[0];
	}
	else
	{
		max(v, 0, N, &val, &idx);
		if(val==ML_POSINF)
		{
			printf("positive infinite value in v\n");
		}
		res = 0;
		for(i=0; i<N; i++)
		{
			if(i==idx || v[i]==ML_NEGINF)
			{
				continue;
			}
			res = res + exp(v[i] - val);
	    }
		lse = val + log(1+res);
    }
	return lse;
}

double MathTools::fac(int n)
{
    double t=1;
    for (int i=n;i>1;i--)
        t*=i;
    return t;
}

double MathTools::dbinom(int x, int size, double p)
{
	return fac(size)/(fac(size-x)*fac(x))*pow(p,x)*pow(1-p,size-x);
}

double MathTools::dnorm(double x, double mu, double sigma)
{
	return (1.0/(sigma*sqrt(2.0*_PI)))*exp(-0.5*(x-mu)*(x-mu)/(sigma*sigma));
}




///////////////////////// the vesion without a prior

int MathTools::glmFit(bool useLess, int* familyR, int* linkR, int* dims, int* nIter,
           double *y, double *offset, double *z,
           double *X, double *nTotal_binom, double *convR, 
           int *rank, double *Xb, double *fitted, double *resid, 
           double *weights, double *phi, int* trace, 
           double *scale, int *df_resid, double* beta) 
{
  double epsilon = 1e-8;       /* Singularity threshold */
  int N, M, maxit, init, useOffset;
  int i = 0, j=0, Nu, dfr = 0, irls, empty = 0;
  int x_rank = 0, convg = 0, iter = 0;
  
  N = dims[0];
  M = dims[1];
  maxit = dims[2];
  init  = dims[3];
  useOffset = dims[4];
  
  int family  = *familyR;
  int link    = *linkR;
  double conv = *convR;
  
  double mu, ri, wi, D, Vmu, wss, wss_last=0.0, ssx=0.0, ssr=0.0;
  double *xi, *xbi, *xbj, zmu, wsum;
    
  if(family > 6 || family < 1){
    printf("family=%d, ", family);
    printf("Invalid family!\n");
  }
  
  /* Is iteration necessary? */
  irls = !((family==GAUSSIAN) && (link==IDENTITY));
  
  /* ----------------------------------------------------------------*
   * by default, initialize mu (fitted) by y itself, with neccesary 
   * modification, e.g., y + 0.1 to avoid log(0) for Poisson family
   * ----------------------------------------------------------------*/
  if (!init) {
    initialize(family, y, fitted, N, nTotal_binom);
  }
  
  /* ----------------------------------------------------------------*
   * Initialize wi (weights) and (standardized) residual 
   * In IRLS, z_i = eta_i + (y_i - mu_i)*(d_eta/d_mu)_i
   *
   * In the following code:
   * ri = (y_i - mu_i)*(d_eta/d_mu)_i
   * wi = Vmu is the weight,  
   *
   * for those invlaid values, set their weight to be 0.
   * ----------------------------------------------------------------*/
  Nu      = 0;
  wsum    = 0.0;
  
  for (i=0; i<N; i++) {

    mu = fitted[i];
    
    if (!muvalid(family, mu)) {
      wi = ri = 0.0;
    }else {
      Nu ++;
      Vmu = varfun(family, mu, *phi);
      
      if (link == family) {
        ri = (y[i] - mu)/Vmu;
        wi = Vmu;
      }else {
        D  = dlink(link, mu);
        ri = D*(y[i] - mu);
        wi = 1.0/(D*D*Vmu);
      }
    }
    
    weights[i] = wi;
    resid[i]   = ri;
    if (weights[i] < epsilon) weights[i] = 0.;
    
    wsum += weights[i];
  }
  
  /* ----------------------------------------------------------------*
   * If summation of all weights is too small, stop 
   * ----------------------------------------------------------------*/
  
  if (wsum < epsilon) {
    if(*trace)
      printf("  glmFit: empty model, summation of all weights are small!\n");
    
    /* set M = -1, so that no extra computation is done */
    M = -1;
  }
  
  if(*trace > 3){
    printf("\n  glmFit: finish initialization, N=%d, M=%d, family=%d, and irls=%d\n", 
            N, M, family, irls);
  }
  
  /* ----------------------------------------------------------------*
   * If M=0, there is only an intercept 
   * ----------------------------------------------------------------*/
  
  if (M == 0){
    zmu = 0.0;
    
    /* all the observations have the same fitted value  */
    for (i=0; i<N; i++) {
      /**
       * z_i = current estimate of eta + (y-mu)/gradient
       * where
       * linkfun(link, fitted[i]) = eta
       * resid[i] = (y-mu)/gradien
       */
      z[i] = linkfun(link, fitted[i]) + resid[i];
      if (useOffset) { z[i] -= offset[i]; }
      
      zmu += z[i]*weights[i];
    }
    
    zmu = zmu/wsum;
    dfr = Nu;
    
    if (useOffset) {
      for (i=0; i<N; i++) {
        fitted[i] = invlink(link, zmu + offset[i]);
        resid[i]  = y[i] - fitted[i];
      }        
    }else {
      mu = invlink(link, zmu);
      for (i=0; i<N; i++) {
        fitted[i] = mu;
        resid[i]  = y[i] - mu;
      }        
    }
    
    if (family>2){
      *scale = wssq(resid, N, weights)/dfr;
    }else{
      *scale = 1.0;
    }
    
    x_rank = 0;
    
  }else if (M > 0) {
    /* ----------------------------------------------------------------*
     * If M>0, include covariates 
     * ----------------------------------------------------------------*/
    
    convg    = 0;
    iter     = 0;
    wss_last = 0.0;
    
    if (!irls) {
      /* Simple linear Gaussian case */
      xi  = X;
      xbi = Xb;
      x_rank = 0;
      
      for (i=0; i<M; i++, xi+=N) {
        wcenter(xi, N, weights, 1, xbi);
        ssx = wssq(xbi, N, weights);
        
        xbj = Xb;
        
        for (j=0; j<x_rank; j++, xbj+=N)  
          wresid(xbi, N, weights, xbj, xbi, beta);
        
        ssr = wssq(xbi, N, weights);
        
        if (ssr/ssx > epsilon) {
          wresid(resid, N, weights, xbi, resid, beta);
          x_rank++;
          xbi+=N;
        }
      }
      
      /* obtain the fitted values */
      for (i=0; i<N; i++) fitted[i] = y[i] - resid[i];
      
      wss_last = wssq(resid, N, weights);
      
    }else{
      
      /* IRLS algorithm */
      while(iter<maxit && !convg) {
        
        if (*trace > 9) {
          printf("    glmFit: iteration %d: \n", iter);
        }
        
        for (i=0; i<N; i++) {
          /**
           * current estimate of eta + (y-mu)/gradient
           *
           * linkfun(link, fitted[i]) = eta
           * resid[i] = (y-mu)/gradien
           */
          z[i] = linkfun(link, fitted[i]) + resid[i];
        }
        
        if (useOffset) {
          for (i=0; i<N; i++) z[i] -= offset[i];
        }
                
        empty = wcenter(z, N, weights, 1, resid);  //removes the mean from z
        
        if (empty == 1) {
          if(*trace > 0)
            printf("  glmFit: empty model, summation of all weights are small!\n");
          
          break;
        }
        
        if (*trace > 9) {
          printf("    glmFit: iteration %d:, initialized z\n", iter);
        }
        
        
        /**
         * tries to fit the regression line (no intercept) to the residuals
         */
        
        xi  = X;
        xbi = Xb;
        x_rank = 0;
        
        for (i=0; i<M; i++, xi+=N) {
          wcenter(xi, N, weights, 1, xbi);
          ssx = wssq(xbi, N, weights);
          xbj = Xb;
          
          for (j=0; j<x_rank; j++, xbj+=N) wresid(xbi, N, weights, xbj, xbi, beta);
          
          ssr = wssq(xbi, N, weights);
                    
          if (ssr/ssx > epsilon) {
            /**
             * takes the residuals after fitting the regression line (no intercept) 
             * to the mean value 
             */
            wresid(resid, N, weights, xbi, resid, beta); 
            x_rank++;
            xbi+=N;
          }
        }
        
        if (*trace > 9) {
          printf("    glmFit: iteration %d:, got Xb\n", iter);
        }
        
        /* well, it is question whether we should give printf or just warning here */
        if(x_rank < M){ 
          if(*trace > 1){
            printf("  glmFit: x_rank=%d, M=%d, X is not full rank\n", x_rank, M);
          }
          
          break;
        }
        
        wss = 0.0;
        Nu  = 0;
        
        for (i=0; i<N; i++) {
          
          if (useOffset) {
            mu = invlink(link, z[i] + offset[i] - resid[i]);
          }else {
            mu = invlink(link, z[i] - resid[i]);
          }
          
          fitted[i] = mu;
          
          if (weights[i] <= 0.0) {
            wi = ri = 0.0;
          } else {
            
            Vmu = varfun(family, mu, *phi);
            Nu ++;
            
            if (link == family) {
              ri = (y[i] - mu)/Vmu;
              wi = Vmu;
            }else {
              D = dlink(link, mu);
              ri = D*(y[i] - mu);
              wi = 1.0/(D*D*Vmu);
            }
            wss += wi*ri*ri;
            
            weights[i] = wi;
            resid[i]   = ri;
            if (weights[i] < epsilon) weights[i] = 0.;
          }
        }
        
        //if(wss > 1.0/epsilon){
        //  if(*trace > 1)
        //    printf("  glmFit: huge wss, indicting failt to fit the model!\n");
        //  
        //  break;
        //}
        
        if (*trace > 5) {
          printf("    glmFit: iteration %d, wss=%.3f\n", iter, wss);
        }
        
        convg = (Nu<=0) || (iter && (fabs(wss-wss_last)/(wss_last + 0.1) < conv));
        wss_last = wss;
        iter ++;
      }
    }
    
    if (convg) {
      /* assume there is an intercept */
      dfr = Nu - 1 - x_rank;
      
      if (family > 2) {
        *scale = wss_last/(dfr);
      }else{
        *scale = 1.0;
      }
    }else {
      dfr    = 0;
      *scale = 1.0;
    }
  }
  
  *df_resid = dfr>0? dfr : 0;
  *rank     = x_rank;
  *nIter    = iter;
    
  return(irls && convg);
}


double MathTools::loglik_Poisson(int N, double* mu, double* y){
  int i;
  double yi, mui, logL = 0.0;
  
  for (i=0; i<N; i++) {
    yi  = y[i];
    mui = mu[i];
    
    logL += (yi*log(mui) - mui - lgammafn(yi + 1.0));
  }
  
  return(logL);
}

/**********************************************************************
 *
 * log likelihood of negative binomial
 *
 **********************************************************************/

double MathTools::loglik_NB(int N, double phi, double* mu, double* y){
  int i;
  double logL1, logL, yi, mui;
  double th = 1.0/phi;
  double logL0 = th*log(th) - lgammafn(th);
  
  logL = 0.0;
  
  for (i=0; i<N; i++) {
    yi  = y[i];
    mui = mu[i];
    
    if (yi==0) {
      logL1  = th*log(th) - th*log(th + mui);
    }else {
      logL1  = lgammafn(th + yi) - lgammafn(yi + 1.0) + yi*log(mui) - (th + yi)*log(th + mui);
      logL1 += logL0;
    }

    logL += logL1;
  }
  
  return(logL);
}

void MathTools::score_info(int N, double theta, double* mu, double* y, 
                double* score, double* info)
{
  int i;
  double score1=0.0, info1=0.0;
  double mui, yi, scorei, infoi, thMui;
  
  for (i=0; i<N; i++) {
    yi  = y[i];
    mui = mu[i];
    
    thMui   = theta + mui;
    scorei  = digamma(yi + theta) - digamma(theta) - (theta + yi)/thMui;
    score1 += (scorei - log(thMui) + 1 + log(theta));
    
    infoi   = trigamma(theta) - trigamma(yi + theta) + (mui - yi)/(thMui*thMui);
    info1  += (infoi + 1/thMui - 1/theta);
  }
  
  *score = score1;
  *info  = info1;
}

/**********************************************************************
 *
 * phi_ml
 *
 * MLE of phi (over-dispersion parameter), given mu 
 *
 * Actually we find MLE of 1/phi here and then take inverse
 *
 **********************************************************************/

int MathTools::phi_ml(double* y, double* mu, int N, int limit, double eps, 
           double* phi, int initPhi, int trace)
{
  double theta0, del, tmp;
  double score=0.0;
  double info=0.0;
  int i, it=0;
  double minTheta = 1e-5;
  double maxTheta = 1.0/minTheta;
  int tryPoisson = 0;
  int tryZINB = 0;
  int fail = 0;
  
  if(initPhi){
    theta0 = 1.0/(*phi);
  }else{
    theta0 = 0.0;
    for (i=0; i<N; i++) {
      tmp = y[i]/mu[i] - 1.0;
      theta0 += tmp*tmp;
    }
    theta0 = (double)N/theta0;
  }
  
  it  = 0;
  del = 1.0;
  
  if(trace > 5) printf("  phi.ml: initial phi = %.2e\n", 1/theta0);
  
  while(it < limit && fabs(del) > eps) {
    score_info(N, theta0, mu, y, &score, &info);
    del     = score/info;
    theta0 += del;
    it     += 1;
    
    if(trace > 5) printf("  phi.ml: iter %d, phi=%.2e, score=%.2e, info=%.2e\n", 
                          it,  1/theta0, score, info);
    
    if (theta0 > maxTheta) {
      theta0 = maxTheta;
      if(trace > 3)
        printf("    phi is truncated at %.2e, no overDispersion?\n", 1/maxTheta);
      
      tryPoisson = 1;
      break;
    }
    
    if(theta0 < minTheta) {
      theta0 = minTheta;
      if(trace > 3)
        printf("    phi is truncated at %.2e, too much overDispersion?\n", 1/minTheta);
      
      tryZINB = 1;
      break;
    }
    
  }
  
  if(it == limit) {
    fail = 1;
    if(trace > 3)
      printf("  phi.ml: iteration limit reached in phi_ml\n");
  }
  
  *phi = 1/theta0;
  
  return(tryPoisson + 2*tryZINB + 4*fail);
}

/**********************************************************************
 *
 * main function of glmNB
 *
 **********************************************************************/

int MathTools::glmNB(bool useLess,int *dims, int *nIter, double *y, double *z, 
          int *linkR, double *offset, double *X, double *convR, 
          int *rank, double *Xb, double *fitted, double *resid, 
          double *weights, double *phi, double *scale, 
          int *df_resid, int* family, double *twologlik, 
          double *scoreTestP, int *trace, double *beta)
{
  int N, M, maxit, init, useOffset;
  int i, succeed = 0, iter = 0, convged=0;
  double conv = *convR;
  int fam0, cv = 0; /* cv is convergence indicator */
  double nTotal_binom=0.0;  /* only used for binomial link, NOT useful here */
  double del, Lm, Lm0, phi0;
  double score, scoreNum, scoreDen, scorePval, yi, mui;
  
  /* convergence indicator for phi 
   * if cvPhi = 0, NB model is OK. 
   * if cvPhi = 1, suggest we need to use Poisson
   * if cvPhi = 2, suggest we need to use ZINB
   */
  int cvPhi;
  
  N = dims[0];
  M = dims[1];
  maxit = dims[2];
  init  = dims[3];
  useOffset = dims[4];
  
  if(*trace > 3) 
    printf("\n  glmNB: N=%d, M=%d, maxit=%d, init=%d, useOffset=%d\n", 
            N, M, maxit, init, useOffset);

  /**
   * if there is no initial values, ignore the parameter *family 
   * and start with a Poission model
   */
  if(!init){
    /* Initial fit */
    
    fam0 = POISSON;
    
    cv = glmFit(useLess, &fam0, linkR, dims, nIter, y, offset, z, 
                X, &nTotal_binom, convR, rank, Xb, fitted, resid, 
                weights, phi, trace, scale, df_resid, beta);
    
    if (cv==0) {
      if(*trace){
        printf("\n  glmNB: fail to converge in initial glmFit by Poission regression\n");
      }
      
      maxit   = -1;
      succeed = 0;
      return(succeed);
    }else {
      
      /* test for overdispersion by Dean's Score test */
      scoreNum = 0.0;
      scoreDen = 0.0;
      
      for (i=0; i<N; i++) {
        yi  = y[i];
        mui = fitted[i];
        scoreNum += (yi - mui)*(yi - mui) - yi;
        scoreDen += mui*mui;
      }
      
      score = scoreNum/sqrt(2.0*scoreDen);
      
      /**
       * double pnorm(double x, double mu, double sigma, int lower_tail, int give_log);
       */
      scorePval = pnorm(score, 0.0, 1.0, false, false);
      
      if(*trace > 3) 
        printf("\n  overdispersion score = %.2e, p-value = %.2e\n\n", score, scorePval);
      
      if(scorePval > *scoreTestP){
        *family    = POISSON;
        Lm         = loglik_Poisson(N, fitted, y);
        *twologlik = 2.0*Lm;
        *phi       = 0.0;
        
        maxit   = -1;
        succeed = 1;
        return(succeed);

      }else {
        fam0    = NB;
        *family = NB;
      
        /**
         * calculate phi by MLE, without initial values of phi
         */
        cvPhi = phi_ml(y, fitted, N, maxit, conv, phi, 0, *trace);
        
        if(cvPhi==0){
          if(*trace > 3) 
            printf("\n  initial value for phi: %e\n", *phi);
        }else if (cvPhi==1){
          if(*trace > 3) 
            printf("\n  Choose Poisson model due to small phi: %e\n", *phi);
          
          *family    = POISSON;
          Lm         = loglik_Poisson(N, fitted, y);
          *twologlik = 2.0*Lm;
          *phi       = 0.0;
          
          maxit   = -1;
          succeed = 1;
          return(succeed);

        }else if(cvPhi==2){
          if(*trace > 1) 
            printf("\n  The overdispersion parameter is too large: phi=%e\n", *phi);
          
          *family = ZINB;
          
          maxit   = -1;
          succeed = 0;
          return(succeed);

        }else { /* estimation of phi fail to converge */
          if(*trace > 1) 
            printf("\n  Estimation of phi fail to converge: phi=%e\n", *phi);
          
          maxit   = -1;
          succeed = 0;
          return(succeed);
        }
      }
    }
  }else {
    /**
     * if there is initial values, 
     * there must be both initial fitted values and inital phi
     */
    fam0 = *family;
    
    cv = glmFit(useLess, &fam0, linkR, dims, nIter, y, offset, z, 
                X, &nTotal_binom, convR, rank, Xb, fitted, resid, 
                weights, phi, trace, scale, df_resid, beta);
    
    if (cv==0) {
      if(*trace > 1){
        printf("\n  glmNB: fail to converge using initial values, fam0=%d\n", fam0);
      }
      
      maxit   = -1;
      succeed = 0;
      return(succeed);
      
    }else{
      /**
       * no need to go further if this is a Poission regression
       */
      
      if(fam0 == POISSON){
        Lm         = loglik_Poisson(N, fitted, y);
        //printf("\n===--->>> fam0=%d, cv=%d, Lm=%e\n", fam0, cv, Lm);

        *twologlik = 2.0*Lm;
        *phi       = 0.0;

        maxit   = -1;
        succeed = 1;
        return(succeed);

      }else {
        cvPhi = phi_ml(y, fitted, N, maxit, conv, phi, 0, *trace);
        
        if(cvPhi==0){
          if(*trace > 3) 
            printf("\n  glmNB: initial value for phi: %e\n", *phi);
        }else {
          if(*trace > 1) 
            printf("\n  glmNB: fail to estimate phi\n");
          
          maxit   = -1;
          succeed = 0;
          return(succeed);
        }
      }
    
    }
  }
  
  if(maxit > 0){
    Lm   = loglik_NB(N, *phi, fitted, y);
  }else{
    Lm   = 0.0;
  }

  del  = 1.0;
  Lm0  = Lm + 1.0;
  iter = 0;
  convged = 0;
  succeed = 0;
  
  while (iter < maxit && (!convged)) {
    
    dims[3] = 1; /* use initial values */
    
    cv = glmFit(useLess, &fam0, linkR, dims, nIter, y, offset, z, 
                X, &nTotal_binom, convR, rank, Xb, fitted, resid, 
                weights, phi, trace, scale, df_resid, beta);
    
    if (cv==0) { 
      if(*trace > 1) 
        printf("\n  glmNB: fail to converge in glmFit\n");
      
      break;
    }
    
    phi0  = *phi;
    cvPhi = phi_ml(y, fitted, N, maxit, conv, phi, 1, *trace);
    
    if(cvPhi==0){
      if(*trace > 3) 
        printf("\n  finish phi_ml, cvPhi=%d, phi=%e\n", cvPhi, *phi);
    }else {
      if(*trace > 1) 
        printf("  glmNB: fail in phi_ml\n");
      
      break;
    }
    
    del = phi0 - *phi;
    Lm0 = Lm;
    Lm  = loglik_NB(N, *phi, fitted, y);
        
    if (*trace > 3) {
      printf("\n  Phi(%d) = (%e, %e), logLik = (%e, %e)\n\n", 
              iter, phi0, *phi, Lm0, Lm);
    }
    
    convged = fabs(Lm0 - Lm) + fabs(del) < conv;
    
    iter++;
  }

  if(iter == maxit) {
    if (*trace) {
      printf("\n  glmNB: Alternation limit reached: iter=%d\n", iter);
    }
  }else if (convged) {
    succeed = 1;
  }

  *twologlik = 2.0*Lm;
  
  return(succeed);
}

