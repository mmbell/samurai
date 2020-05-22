/*
 *  CostFunction.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "CostFunction.h"
#include <cmath>
#include <iostream>
#include <limits>
#include "timing/gptl.h"
#include "solver.inc"

CostFunction::CostFunction(const Projection& proj, const int& numObs, const int& stateSize) :
	projection(proj)
{
	mObs = numObs;
	nState = stateSize;
}

CostFunction::~CostFunction()
{
	
}

void CostFunction::setNumObservations(const int& numObs)
{
	mObs = numObs;
}

int CostFunction::getNumObservations()
{
	return mObs;
}

void CostFunction::setLengthStateVector(const int& stateSize)
{
	nState = stateSize;
}

int CostFunction::getLengthStateVector()
{
	return nState;
}

bool CostFunction::minimize()
{
  GPTLstart("CostFunction::minimize");

  real ftol = S_CONV_TOL;
  real minimum = 1e34; //used in Samurai CG
  int verbose = S_VERBOSE;

  cout << "Beginning minimize()" << endl;
  cout << "\tSolver conv. tolerance = " << ftol << endl;

  //keep track of linesearch function call count
  ls_cnt = 0; //global var

  //work vector or MT linesearch
  mt_work = new real[nState];

  // choose solver (currState is update by solver)
  if (S_SOLVER == 1) {
    cout << "SOLVER: Samurai Truncated Newton " << endl;
    truncatedNewton(currState, currGradient, ftol);
  } else if (S_SOLVER == 2) {
    cout << "SOLVER: Samurai Conjugate Gradient " << endl;
    conjugateGradient(currState, currGradient, ftol, minimum);
  } else {
    cout << "\tS_SOLVER = " << S_SOLVER << " is not a valid option. Using Samurai TN instead." << endl;
    truncatedNewton(currState, currGradient, ftol);
  }

  if (verbose) {
    if (ls_cnt) cout << "\t\t (Linesearch iterations = " << ls_cnt << " )" << endl;
  }
  
  delete[] mt_work;

  GPTLstop("CostFunction::minimize");
	return true;
}

void CostFunction::truncatedNewton(real* qstate, real* g, const real ftol)
{

  //Newton's method with CG to (inexactly) solve for the search direction
  // qstate = the initial state
  // g = current gradient
  // ftol = desired solving tolerance


  int its, cg_its, total_cg_its;
  int outer_itmax, cg_itmax;
  int j;
  int verbose, neg_curve;
  int ls_ret;

  real cg_tol;
  real f_init, f_val;
  real n_init_grad, n_grad, grad_dot;
  real gg;
  real rr, pAp, rr_m1, r_norm, r0_norm, rel_resid;
  real beta, alpha;
  real initstep;

  real *x = new real[nState];
  real *p = new real[nState];
  real *Ap = new real[nState];
  real *r = new real[nState];


  GPTLstart("CostFunction::TruncNewton");


  cg_tol = S_INNER_CONV_TOL;
  cg_itmax = S_INNER_MAXITER;
  outer_itmax = S_MAXITER;
  verbose = S_VERBOSE;


  //cumulative total of inner cg its  
  total_cg_its = 0;

  //initial step length for newton step
  initstep = 1.0;

  //Newton Step Loop (OUTER LOOP)
  for (its = 0; its < outer_itmax; its++) {

    //get values of the function (into f_val) and the gradient (into g) based on qstate
    //only have to do on the first iteration if using MT line search
    if (its == 0) f_val = funcValueAndGradient(qstate, g);
    
    //calculate the norm of the current gradient
    grad_dot = 0.0;
    for (j = 0; j  < nState; j++) {
      grad_dot +=  g[j]*g[j];
    }
    n_grad = sqrt(grad_dot);

    //collect initial values on first iteration
    if (its == 0) {
      f_init = f_val;
      n_init_grad = n_grad;
      //cout << "Newton: Initial Norm of the Gradient = " << n_init_grad << endl;
    }

    cout << "\tNewton Iteration: " << its << "\tJ = " << f_val << "\tResidual = " << n_grad << endl; //prints cost fcn value       

    //Newton conv. check	 
    gg = n_grad/n_init_grad;

    //cout << "gg = " << gg << endl;

    if ( gg < ftol) {
      cout << "\tMinimum J: " << f_val << endl;
      cout << "\tFound minimum in " << its << " outer Newton iterations." << endl;
      cout << "\t(and a total of " << total_cg_its << " inner CG iterations.)" << endl;
	   
      delete[] x;
      delete[] p;
      delete[] Ap;
      delete[] r;

      GPTLstop("CostFunction::TruncNewton");
      cout << "\tFINAL  ||g(X)||/||g(X0)|| = " << gg << endl;
      return;
    } //end conv. check

    //CG INNER LOOP
    // Use non-preconditioned CG to solve linear system for Newton direction: d_k
    // H_k * d_k = -g_k, where H is the Hessian and g is the gradient
    //init
    for (j = 0; j < nState; j++) {
      x[j] = 0.0;      // zero initial guess (x_0) for CG
      r[j] = -g[j];          //initial residual r0 = b -Ax = -g 
      p[j] = r[j];    //inital direction po = r0 
    }
    r0_norm = n_grad;

    neg_curve = 0; //negative curvature check 

    if (verbose) cout << "\t\tCG iteration 0:   r_norm = " << r0_norm << "     rel_resid = 1.0"  << endl;

    //CG LOOP
    for (cg_its = 0; cg_its < cg_itmax; cg_its ++){

      //update search direction
      if (cg_its == 0) {
				//rr = <r0, r0>
				rr = grad_dot;
				//already have set p0 to r0	
      } else {
				// rr calculated at end of last loop
				//Beta = <r_k, r_k>/<r_k-1, r_k-1>         
				beta = rr/rr_m1;
				// p_k = r_k + beta*p_k-1
				for (j = 0; j < nState; j++) {
	  			p[j] = r[j] + beta*p[j];
				}
     	}

      //Find A*p
      funcHessian(p, Ap);

      // alpha = <r_k, r_k>  / <Ap_k, p_k> 
      pAp = 0.0;
      for (j=0; j < nState; j++) {
				pAp += p[j]*Ap[j];
      }
     
      alpha = rr/pAp;

      //check for negative curvature
      if (pAp < 0) {
				if (verbose) cout << "Negative curvature in CG iteration... pAp = " << pAp << endl;
				neg_curve = 1;
				if (cg_its > 0) {//if first iteration, complete this iteration and then stop 
	  			break;         //else break out now and keep the previous x
				}
      }
      
      //update x and r and compute <r_k+1, r_k+1>
      //x_k+! = x_k + alpha*p_k
      //r_k+1 = r_k - alpha*A*p_k
      rr_m1 = rr;
      rr = 0.0;
      for (j=0; j< nState; j++){
				x[j] = x[j] + alpha*p[j];
				r[j] = r[j] - alpha*Ap[j];  
				rr += r[j]*r[j];
      }
	   
      //CG cconvergence check
      r_norm = sqrt(rr);
      rel_resid = r_norm/r0_norm;
      if (verbose) cout << "\t\tCG iteration " << cg_its + 1 <<  ":  r_norm = " << r_norm << "      rel_resid = " << rel_resid << endl;

      if (rel_resid < cg_tol) {
				break;
      }

      if (neg_curve) break; //this is only at its = 0 => otherwise we already left the loop

    } //end of CG loop 
    total_cg_its += cg_its + 1;

    //Update Newton step 
    //qstate_new = q_state + alpha*d_k
    //newton direction is the final CG iterate; d_k = x
    //Find step length (alpha) via line search

    //Use MT linesearch instead (input state, gradient, search dir, fval, initstep)
    //also returns the gradient
    ls_ret = MTLineSearch(qstate, g, x, &f_val, initstep); 


  } //end of Netwon loop

  // if we got here, then there were too many Newton iterations 
  cout << "Iterations exceeded in Truncated Newton: " << its  << endl;
  GPTLstop("CostFunction::TruncNewton");


  delete[] r;
  delete[] x;
  delete[] p;
  delete[] Ap;

  return;

}
void CostFunction::conjugateGradient(real* q, real* xi, const real ftol, real fret)
{
  GPTLstart("CostFunction::ConjugateGradient");

	int j, its;
	int verbose;

	real gg, gam, fq, dgg, gy, gd, dd, yy, dy, sy, sg, eta;
	real n_init_grad, n_grad;;


	real* g = new real[nState];
	real* h = new real[nState];
	real* y = new real[nState]; //added for alternative beta calculation
	real* s = new real[nState]; //added for alternative beta calculation
	real*q_prev = new  real[nState]; //added for alternative beta calculation

	verbose = S_VERBOSE;

	if (verbose) cout << "\tNCG beta type = " << S_BETA_TYPE << endl;


	fq = funcValue(q); //cost function value
	funcGradient(q, xi); //now xi is the gradient

	//AB - find the norm of the initial gradient - its own loop for now
	n_init_grad = 0.0;
	for (j=0; j< nState; j++) {
	  n_init_grad += xi[j]*xi[j];
	}
	n_init_grad = sqrt(n_init_grad);
	cout << "\tINIT NORM GRADIENT = " << n_init_grad << endl;

	for (j=0; j<nState; j++) {
	 	g[j] = -xi[j]; // g are search directions - start with negative gradient
		xi[j] = h[j] = g[j];  //set xi and h to neg gradient also 
		//(init search direction)
		q_prev[j] = q[j]; //need for other beta calculations
	}

	/* MAIN loop - ITMAX iterations max */
	for (its=0; its< S_MAXITER; its++) {
	  if (verbose) cout << "\t\tIteration: " << its << "\tJ: " << fq << endl; //prints cost fcn value
	  /* call line search to determine step size */
	  /* this is brent's method:
			 fret is updated cost func value (at q) 
			 q is the current guess for the min - updated here to new min state 
	     xi is the search direction - updated here to reflect step taken
		   minimizing f(q + alpha*xi) xi_new = q + alpha*xi
		*/

		dlinmin(q, xi, fret);

		/*convergence check */
		if (S_CG_CONV_TYPE ==1) {
		  if (2.0*fabs(fret-fq) <= ftol*(fabs(fret)+fabs(fq)+S_EPS)) {
			cout << "\tMinimum J: " << fret << endl;
			cout << "\tFound minimum in " << its << " iterations." << endl;

			//AB: extra convergence info
			gg =  (2.0*fabs(fret-fq))/(fabs(fret)+fabs(fq)+S_EPS);
			cout << "FINAL step size convergence value = " << gg << endl;

			//Let's also see what the relative gradient norm is
			n_grad = 0.0;
			for (j=0; j< nState; j++) {
			  n_grad += g[j]*g[j];
			}
			n_grad = sqrt(n_grad);
			dd = n_grad/n_init_grad;
			cout << "(Note: Relative norm of gradient = " << dd << ")" << endl;
			
			delete[] g;
      delete[] h;
      delete[] y;
      delete[] s;
      delete[] q_prev;
			GPTLstop("CostFunction::ConjugateGradient");

			return;
		  }
		}
		else { //check gradient norms 
		  if (its > 0) {
		    gg = n_grad/n_init_grad;
		    if ( gg < ftol) {
		      cout << "\tMinimum J: " << fret << endl;
		      cout << "\tFound minimum in " << its << " iterations." << endl;
		      delete[] g;
		      delete[] h;
		      delete[] y;
		      delete[] s;
		      delete[] q_prev;
		      GPTLstop("CostFunction::ConjugateGradient");

		      cout << "FINAL  ||g(X)||/||g(X0)|| = " << gg << endl;
		      return;
		    }
		  }
		}


		fq = fret;
		/* calculate new gradient (xi) */
		funcGradient(q, xi);

		//AB: for checking rel norm of the gradient
		if (S_CG_CONV_TYPE != 1) {
		  n_grad = 0.0;
		  for (j=0; j< nState; j++) {
		    n_grad += xi[j]*xi[j];
		  }
		  n_grad = sqrt(n_grad);
		  dd = n_grad/n_init_grad;
		  if (verbose) cout << "\t\t\tRelative norm of gradient = " << dd << endl;
		  if (verbose) cout << "\t\t\tNorm of gradient = " << n_grad << endl;
		}

		/****************************/
		/*Calculate Beta (i.e., gam)*/  
		dgg = gg = gy = gd = yy = dy = sy = sg = dd= 0.0;
		switch(S_BETA_TYPE) {
		    case 1: // Polak-Ribiere (PR) - original choice in Samurai
		      for (j=0; j<nState; j++) {
						gg += g[j]*g[j]; // beta denominator 
						dgg += (xi[j]+g[j])*xi[j]; // numerator
		      }
		      if (gg == 0.0) {  //unlikely
		        GPTLstop("CostFunction::ConjugateGradient");
						return;
		      }
		      gam = dgg/gg;
		      //cout << "\tAB Beta: " << gam << endl;
		      break;
		   case 2: // Polak-Ribiere-Polyak (PRP+)  - same as PR +  restart
		      for (j=0; j<nState; j++) {
						gg += g[j]*g[j]; // beta denominator 
						dgg += (xi[j]+g[j])*xi[j]; // numerator
		      }
		      if (gg == 0.0) {  //unlikely
		        GPTLstop("CostFunction::ConjugateGradient");
						return;
		      }
		      gam = dgg/gg;
		      //cout << "\tAB Beta: " << gam << endl;
		      gam = CF_MAX(gam, 0.0); //added to 'restart' if gam is negative 
		      break;
		   case 3: /* Fletcher-Reeves (FR) */
		     for (j=0; j<nState; j++) {
		       gg += g[j]*g[j]; // beta denominator 
		       dgg += xi[j]*xi[j]; // numerator
		     }
		     if (gg == 0.0) {  //unlikely
		       GPTLstop("CostFunction::ConjugateGradient");
		       return;
		     }
		     gam = dgg/gg;
		     break;
		   case 4: /*Dai-Yuan (DY) - this may need a special linesearch...*/
		     /* need to define y */
		     for (j=0; j < nState; j++) {
		       y[j] = xi[j] + g[j];
		       dgg +=  xi[j]*xi[j]; /*DY beta numerator*/ 
		       gg +=  y[j]*h[j]; /*DY denominator*/
		     }
		     if (gg == 0.0) {  //unlikely :)
		       GPTLstop("CostFunction::ConjugateGradient");
		       return;
		     }
		     gam = dgg/gg;
		     break;
    	 case 5: /* Hager-Zhang (HZ) */
		     /* need to define y */
		     for (j=0; j< nState; j++) {
		       y[j] = xi[j] + g[j];
		       dy += y[j]*h[j];  
		       yy += y[j]*y[j];
		       gy += xi[j]*y[j];
		       gd += xi[j]*h[j];
		       gg += g[j]*g[j];
		       dd += h[j]*h[j];
		     }
		     if (dy == 0.0) { 
		       cout << "\tAB - HZ error in dy" << endl;
		       GPTLstop("CostFunction::ConjugateGradient");
		       return;
		     }
		     gam = gy - 2.0*gd*(yy/dy);
		     gam = gam/dy;
		     
		     gg = sqrt(gg);
		     dd = sqrt(dd);
		     eta = dd*CF_MIN(0.01,gg);
		     eta = -1.0/eta;
		     gam = CF_MAX(gam, eta);
		     break;
		   case 6: //Dai-Kou ((DK) 2013
		     for (j=0; j< nState; j++) {
		       y[j] = xi[j] + g[j];
		       s[j] = q[j] - q_prev[j];
		       dy += y[j]*h[j];  
		       yy += y[j]*y[j];
		       gy += xi[j]*y[j];
                       sy += s[j]*y[j];
		       sg += s[j]*xi[j];
		       q_prev[j] = q[j]; //for next iteration
		     }

		     gam = (gy/dy)-((yy/sy)*(sg/dy));
		     break;
		} //end of switch for beta


		/* Update the new search direction (xi) 
		   and h[i] is the prev. search direction */
		for (j=0; j<nState; j++) {
			g[j] = -xi[j];
			xi[j] = h[j] = g[j] + gam*h[j];
		}
	}
	// Got here -> there were too many iterations
	delete[] g;
	delete[] h;
	delete[] y;
	delete[] s;
	delete[] q_prev;

	cout << "Iterations exceeded in inner minimization loop" << endl;
	GPTLstop("CostFunction::ConjugateGradient");
	return;
}

/* line minimization - using derivatives */
void CostFunction::dlinmin(real* &p, real* &xi, real &fret)
{
  GPTLstart("CostFunction::dlinmin");
	const real TOL=2.0e-8;
	int j;
	real xx,xmin,fx,fb,fa,bx,ax;
	int verbose = S_VERBOSE;

	// Fill the temporary state vector
	for (j=0; j<nState; j++) {
		tempState[j] = p[j];
		tempGradient[j] = xi[j];
	}
	ax = 0.0;
	xx = 1.0;
	xmin = 0.0;
	mnbrack(ax,xx,bx,fa,fx,fb);
	fret = dbrent(ax,xx,bx,TOL,xmin);
	for (j=0; j<nState; j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	if (verbose) cout << "\t\t\tLS: alpha = " << xmin << endl;

  GPTLstop("CostFunction::dlinmin");
}

//evaluate function
real CostFunction::f1dim(const real x)
{
  GPTLstart("CostFunction::f1dim");
	int j;
	real f1 = 0.0;
	for (j=0; j<nState; j++) xt[j] = tempState[j] + x*tempGradient[j];
	f1 = funcValue(xt);
	
	
  GPTLstop("CostFunction::f1dim");
	return f1;
}

//evaluate function derivative
real CostFunction::df1dim(const real x)
{
  GPTLstart("CostFunction::df1dim");
	int j;
	real df1 = 0.0;
	for (j=0; j<nState; j++) {
		xt[j] = tempState[j] + x*tempGradient[j];
		df[j] = 0;
	}
	funcGradient(xt, df);
	for (j=0; j<nState; j++) df1 += df[j]*tempGradient[j];


  GPTLstop("CostFunction::df1dim");
	return df1;
}

//eval function value and derivative together
real CostFunction::f1dim_and_df1dim(const real x, real *grad)
{

  GPTLstart("CostFunction::f1dim_and_df1dim");

  int j;
  real df1 = 0.0;
  real f1 = 0.0;
  for (j=0; j<nState; j++) {
    xt[j] = tempState[j] + x*tempGradient[j];
    df[j] = 0;
  }
  f1 = funcValueAndGradient(xt, df);  

  for (j=0; j<nState; j++) {
    df1 += df[j]*tempGradient[j]; 
  }

  *grad = df1;

  GPTLstop("CostFunction::f1dim_and_df1dim");

  return f1;
}


inline void CostFunction::mov3(real &a, real &b, real &c,
							   const real d, const real e, const real f)
{
  GPTLstart("CostFunction::mov3");
	a=d; b=e; c=f;
  GPTLstop("CostFunction::mov3");
}

/* Brent's line seach method - uses derivatives */
real CostFunction::dbrent(const real ax, const real bx, const real cx,
							const real tol, real &xmin)
{
  GPTLstart("CostFunction::dbrent");
	const int DB_ITMAX=100;
	const real DB_ZEPS=numeric_limits<real>::epsilon()*1.0e-3; // Replace with actual machine epsilon
	bool ok1, ok2;
	int iter;
	real a,b,d=0.0,d1,d2,du,dv,dw,dx,e=0.0;
	real fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
	
	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;

	//AB: eval both in one call
	//fw=fv=fx=f1dim(x);
	//dw=dv=dx=df1dim(x);
	fw=fv=fx=f1dim_and_df1dim(x, &dx);
	dw=dv=dx;

	for (iter=0; iter<DB_ITMAX; iter++) {
	  xm = 0.5*(a+b);
	  tol1 = tol*fabs(x)+DB_ZEPS;
	  tol2 = 2.0*tol1;
	  if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
	    xmin = x;
	    GPTLstop("CostFunction::dbrent");
	    ls_cnt += iter;
	    return fx;
	  }
	  if (fabs(e) > tol1) {
	    d1 = 2.0*(b-a);
	    d2 = d1;
	    if (dw != dx) d1 = (w-x)*dx/(dx-dw);
	    if (dv != dx) d2 = (v-x)*dx/(dx-dv);
	    u1 = x+d1;
	    u2 = x+d2;
	    ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
	    ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
	    olde = e;
	    e = d;
	    if (ok1 || ok2) {
	      if (ok1 && ok2)
					d=(fabs(d1) < fabs(d2) ? d1 : d2);
	      else if (ok1)
					d=d1;
	      else
					d=d2;
	      if (fabs(d) <= fabs(0.5*olde)) {
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
		  			d=CF_SIGN(tol1,xm-x);
	      } else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	      }
	    } else {
	      d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	    }
	  } else {
	    d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	  }

	  //AB: combine f1 and df1 calls
	  if (fabs(d) >= tol1) {
	    u=x+d;
	    fu = f1dim_and_df1dim(u, &du);
	  } else {
	    u=x+CF_SIGN(tol1,d);
	    fu = f1dim_and_df1dim(u, &du);
	    if (fu > fx) {
	      xmin = x;
	      GPTLstop("CostFunction::dbrent");
	      ls_cnt += iter;
	      return fx;
	    }
	  }

	  if (fu <= fx) {
	    if (u >= x) a=x; else b=x;
	    	mov3(v,fv,dv,w,fw,dw);
	    	mov3(w,fw,dw,x,fx,dx);
	    	mov3(x,fx,dx,u,fu,du);
	  } else {
	    if (u < x) a=u; else b=u;
	    if (fu <= fw || w == x) {
	      mov3(v,fv,dv,w,fw,dw);
	      mov3(w,fw,dw,u,fu,du);
	    } else if (fu < fv || v == x || v == w) {
	      mov3(v,fv,dv,u,fu,du);
	    }
	  }
	}
	// Too many iterations
	cout << "\tdbrent reached DB_ITMAX: " << DB_ITMAX << endl;
	ls_cnt += iter;
	GPTLstop("CostFunction::dbrent");
	return 0.0;
}

inline void CostFunction::shft3(real &a, real &b, real &c, const real d)
{
  GPTLstart("CostFunction::shft3");
	a=b;
	b=c;
	c=d;
  GPTLstop("CostFunction::shft3");
}

/* Bracketing the minumum of  readl funtion - given 3 points ax, bx, cx */
void CostFunction::mnbrack(real &ax, real &bx, real &cx, 
						   real &fa, real &fb, real &fc)
{
  GPTLstart("CostFunction::mnbrack");
	const real GOLD=1.618034, GLIMIT=100.0, TINY=1.0e-20;
	real ulim,u,r,q,fu;
	
	fa = f1dim(ax);
	fb = f1dim(bx);
	if (fb > fa) {
		CF_SWAP(ax, bx);
		CF_SWAP(fb, fa);
	}
	cx = bx+GOLD*(bx-ax);
	fc = f1dim(cx);
	while (fb > fc) {
		r = (bx - ax)*(fb - fc);
		q = (bx - cx)*(fb - fa);
		u = bx - ((bx-cx)*q-(bx-ax)*r)/
			(2.0*CF_SIGN(CF_MAX(fabs(q-r), TINY), q-r));
		ulim = bx + GLIMIT*(cx-bx);
		if ((bx-u)*(u-cx) > 0.0) {
			fu = f1dim(u);
			if (fu < fc) {
				ax = bx;
				bx = u;
				fa = fb;
				fb = fu;
        GPTLstop("CostFunction::mnbrack");
				return;
			} else if (fu > fb) {
				cx = u;
				fc = fu;
        GPTLstop("CostFunction::mnbrack");
				return;
			}
			u = cx + GOLD*(cx-bx);
			fu = f1dim(u);
		} else if ((cx-u)*(u-ulim) > 0.0) {
			fu = f1dim(u);
			if (fu < fc) {
				shft3(bx,cx,u,u+GOLD*(u-cx));
				shft3(fb,fc,fu,f1dim(u));
			}
		} else if ((u-ulim)*(ulim-cx) >= 0.0) {
			u = ulim;
			fu = f1dim(u);
		} else {
			u = cx + GOLD*(cx-bx);
			fu = f1dim(u);
		}
		shft3(ax,bx,cx,u);
		shft3(fa,fb,fc,fu);
	}
  GPTLstop("CostFunction::mnbrack");
}

//More-and Thuente line search
/* Translation of minpack subroutine cvsrch, whose
   purpose is to find a step which satisfies 
   a sufficient decrease condition and a curvature 
   condition. This is largely taken from O'Leary's matlab version 
   and the PETSc C version*/

int CostFunction::MTLineSearch(real* &x, real* &g, real *s, real *fval, real initstep)
{

  // x = input/output state (x = x + alpha*s)
  // g = input/output gradient
  // s = search direction (input)
  // fval = function values (input/output)
  // initstep = initial step length (should be 1.0 for Newton)

  int j, i;
  int n, bracket, stage1, max_funcs, infoc;
  int reason; // 0 = continue, 1= success, -3 = not search dir, 3 = halted, 4 = halted at max funcs
              // 5 = step at upper bound, 6 = step at lower bound, 7 = halted rtol

  real stx, sty, fx, fy, dgx, dgy, step, dg;
  real ftol, gtol, ftest1, fm, fxm, fym, dgm, dgxm, dgym;
  real rtol;
  real dginit, finit, dgtest;
  real width, width1, xtrapf;

  real stepmin, stepmax;


  int verbose = S_VERBOSE;
  int n_feval = 0;
  
  GPTLstart("CostFunction::MTLineSearch");

  n = nState;

  //compute dginit <g,s> (inital gradient in the search direction)
  dginit = 0.0;
  for (j=0; j< nState; j++){
    dginit += g[j]*s[j];
  }
  
  if (dginit >= 0.0) {
    cout << "Search direction is not a descent direction! dginit = " << dginit << endl;
    reason = -3;
    return reason;
  }

  //Initialize stuff
  reason = 0;
  bracket = 0;
  stage1 = 1;
  finit = *fval;
  infoc = 1;

  stepmin = 1.0e-20; //min allowable step length (MT)
  stepmax = 1.0e20; //max allowable step lengthe (MT)

  ftest1 = 0.0;
  xtrapf = 4.0;
  ftol = 0.0001;  //tol for sufficient decrease
  gtol = 0.9; //tol for curvature condition
  rtol = 1.0e-10; //should reduce if switch to SP
  max_funcs = 30; //max function evals

  dgtest = ftol*dginit;
  width = stepmax - stepmin;
  width1 = 2*width;
 
  //Step, function, and derivative:
  //  stx, fx, dgx - at the best step
  //  sty, fy, dgy - the other endpoint of the interval of uncertainty
  //  step, f, dg - at the current step 
  fx = fy = finit;
  dgx = dgy = dginit;
  stx = sty = 0.0;
  step = initstep;

  //copy the orig x into the work vector
  for (j=0; j<nState; j++) {
    mt_work[j] = x[j];
  }

  //Begin iteration
  for (i=0; i< max_funcs; i++){

    //Set the minimum and maximum steps to correspond to the present interval of uncertainty.
    if (bracket) {
      stepmin = CF_MIN(stx, sty);
      stepmax = CF_MAX(stx, sty);
    } else {
      stepmin = stx;
      stepmax = step + xtrapf*(step - stx);
    }

    // Force the step to be within the bounds
    step = CF_MAX(step, stepmin);
    step = CF_MIN(step, stepmax);

    //in case of unusual termination,go to lowest step size (infoc is set by cstep routine)
    if ((stx!=0) && (((bracket) && (step <= stepmin || step >= stepmax)) || ((bracket) && (stepmax - stepmin <= rtol * stepmax)) || (infoc == 0))) {
      step = stx;
    }

    //  Evaluate the function and gradient at step and compute the directional derivative.
    //the orig x was copied into the work vector
    //x = x-orig + step * s;
    for (j = 0; j<nState; j++) {
      x[j] = mt_work[j] + step*s[j];
    }
    *fval = funcValueAndGradient(x, g); //f is func value at x is g is gradient at x
    ls_cnt ++;
    n_feval++; 

    //compute dg = <g,s>
    dg = 0.0;
    for (j = 0; j< nState; j++){
      dg += g[j]*s[j];
    }
    ftest1 = finit + step*dgtest;

    //to compare with PETSc -tao_ls_monitor
    //cout << "its = " << i+1 << " f = " << *fval << " step = " << step << endl;
    //cout << "stx = " << stx << " fx = " << fx << " dgx = " << dgx << endl;
    //cout << "sty = " << sty << " fy = " << fy << " dgy = " << dgy << endl;


    //test for convergence
    if (((*fval - ftest1 <= 1.0e-10 * fabs(finit)) &&  (fabs(dg) + gtol*dginit <= 0.0))) {
      //converged - yay!
      reason = 1; //success
      break;
    }
    

    // Checks for other convergence problems 
    if (((bracket) && (step <= stepmin || step >= stepmax)) || (!infoc)) {
      cout << "MT: Rounding errors may prevent further progress.  May not be a step satisfying" << endl;
      cout << "sufficient decrease and curvature conditions. Tolerances may be too small." << endl;
      reason = 3;
      break;
    }
    if ((step == stepmax) && (*fval <= ftest1) && (dg <= dgtest)) {
      cout << "MT: Step is at the upper bound, stepmax. " << endl;
      reason = 5;
      break;
    }
    if ((step == stepmin) && (*fval >= ftest1) && (dg >= dgtest)) {
      cout << "MT: Step is at the lower bound, stepmin." << endl;
      reason = 6;
      break;
    }
    if ((bracket) && (stepmax - stepmin <= rtol* stepmax)){
      cout << "MT: Relative width of interval of uncertainty is at most rtol." << endl;
      reason = 7;
      break;
    }

    /*In the first stage, we seek a step for which the modified function
      has a nonpositive value and nonnegative derivative */
    if ((stage1) && (*fval <= ftest1) && (dg >= dginit * CF_MIN(ftol, gtol))) {
         stage1 = 0;
    }

    if ((stage1) && (*fval <= fx) && (*fval > ftest1)) {
      fm   = *fval - step * dgtest;    /* Define modified function */
      fxm  = fx - stx * dgtest;         /* and derivatives */
      fym  = fy - sty * dgtest;
      dgm  = dg - dgtest;
      dgxm = dgx - dgtest;
      dgym = dgy - dgtest;

      // if (dgxm * (ls->step - stx) >= 0.0)
      // Update the interval of uncertainty and compute the new step
      
      infoc = MTcstep(&stx, &fxm, &dgxm, &sty, &fym, &dgym, &step, &fm, &dgm, &bracket, &stepmin, &stepmax);

      fx  = fxm + stx * dgtest; /* Reset the function and */
      fy  = fym + sty * dgtest; /* gradient values */
      dgx = dgxm + dgtest;
      dgy = dgym + dgtest;

    } else {
      // Update the interval of uncertainty and compute the new step
      infoc = MTcstep(&stx, &fx, &dgx, &sty, &fy, &dgy, &step, fval, &dg, &bracket, &stepmin, &stepmax);
    }

    // Force a sufficient decrease in the interval of uncertainty 
    if (bracket) {
      if (fabs(sty - stx) >= 0.66 * width1) {
				step = stx + 0.5*(sty - stx);
      }
      width1 = width;
      width = fabs(sty - stx);
    }

  } //end of i loop for max_funcs
  
  if (i == max_funcs) {
    cout << "\tNumber of line search function evals exceeded the maximum: " << max_funcs << endl;
    reason = 4;
  }

  if (verbose) cout << "\t\tMT LS: Number of function evals = "<< n_feval << endl;

  if (verbose) cout << "\t\tMT LS: step = " << step << endl;


  GPTLstop("CostFunction::MTLineSearch");

  return reason;

}


/*  The purpose of cstep is to compute a safeguarded step for
     a linesearch and to update an interval of uncertainty for
     a minimizer of the function. */

/* The parameter stx contains the step with the least function
    value. The parameter stp contains the current step. It is
    assumed that the derivative at stx is negative in the
    direction of the step. If brackt is set true then a
    minimizer has been bracketed in an interval of uncertainty
    with endpoints stx and sty. */

int CostFunction::MTcstep(real *stx, real *fx, real *dx, real* sty, real *fy, real *dy, real *stp, real *fp, real *dp, int *bracket, real *stepmin, real *stepmax) {

  real gamma1, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
  int  bound;
  int  info = 0; //return this


  //check for errors
  if (*bracket && (*stp <= CF_MIN(*stx,*sty) || (*stp >= CF_MAX(*stx,*sty)))) {
    cout << "MTcstep issue: bad stp in bracket" << endl;
    return info;
  }
  if (*dx * (*stp-*stx) >= 0.0) {
    cout << "MTcstep issue: dx * (stp-stx) >= 0.0" << endl;
    return info;
  }  
  if (*stepmax < *stepmin) {
    cout << "MTcstep issue: stepmax < stepmin" << endl;
    return info;
  }

  //check that derivatives have opposite sign

  sgnd = *dp * (*dx / fabs(*dx));

  /*which case ? */
  if (*fp > *fx) {
    /* Case 1: a higher function value.
      The minimum is bracketed. If the cubic step is closer
     to stx than the quadratic step, the cubic step is taken,
     else the average of the cubic and quadratic steps is taken. */

    info = 1;
    bound = 1;
    theta = 3 * (*fx - *fp) / (*stp - *stx) + *dx + *dp;
    s = CF_MAX(fabs(theta), fabs(*dx));
    s = CF_MAX(s, fabs(*dp));
    gamma1 = s*sqrt(pow(theta/s,2.0) - (*dx/s)*(*dp/s));
    if (*stp < *stx) gamma1 = -gamma1;
    /* Can p be 0?  Check */
    p = (gamma1 - *dx) + theta;
    q = ((gamma1 - *dx) + gamma1) + *dp;
    r = p/q;
    stpc = *stx + r*(*stp - *stx);
    stpq = *stx + ((*dx/((*fx-*fp)/(*stp-*stx)+*dx))*0.5) * (*stp - *stx);

    if (fabs(stpc-*stx) < fabs(stpq-*stx)) {
      stpf = stpc;
    } else {
      stpf = stpc + 0.5*(stpq - stpc);
    }
    *bracket = 1;

  } else if (sgnd < 0.0) {
    /* Case 2: A lower function value and derivatives of
       opposite sign. The minimum is bracketed. If the cubic
       step is closer to stx than the quadratic (secant) step,
       the cubic step is taken, else the quadratic step is taken. */

    info = 2;
    bound = 0;
    theta = 3*(*fx - *fp)/(*stp - *stx) + *dx + *dp;
    s = CF_MAX(fabs(theta), fabs(*dx));
    s = CF_MAX(s, fabs(*dp));
    gamma1 = s*sqrt(pow(theta/s,2.0) - (*dx/s)*(*dp/s));
    if (*stp > *stx) {
      gamma1 = -gamma1;
    }
    p = (gamma1 - *dp) + theta;
    q = ((gamma1 - *dp) + gamma1) + *dx;
    r = p/q;
    stpc = *stp + r*(*stx - *stp);
    stpq = *stp + (*dp/(*dp-*dx))*(*stx - *stp);

    if (fabs(stpc-*stp) > fabs(stpq-*stp)) {
      stpf = stpc;
    } else {
      stpf = stpq;
    }
    *bracket = 1;

  } else if (fabs(*dp) < fabs(*dx)) {
    /* Case 3: A lower function value, derivatives of the
      same sign, and the magnitude of the derivative decreases.
      The cubic step is only used if the cubic tends to infinity
      in the direction of the step or if the minimum of the cubic
      is beyond stp. Otherwise the cubic step is defined to be
      either stepmin or stepmax. The quadratic (secant) step is also
      computed and if the minimum is bracketed then the step
      closest to stx is taken, else the step farthest away is taken. */

    info = 3;
    bound = 1;
    theta = 3*(*fx - *fp)/(*stp - *stx) + *dx + *dp;
    s = CF_MAX(fabs(theta), fabs(*dx));
    s = CF_MAX(s, fabs(*dp));

    /* The case gamma1 = 0 only arises if the cubic does not tend
       to infinity in the direction of the step. */
    gamma1 = s*sqrt(CF_MAX(0.0, pow(theta/s,2.0) - (*dx/s)*(*dp/s)));
    if (*stp > *stx) gamma1 = -gamma1;
    p = (gamma1 - *dp) + theta;
    q = (gamma1 + (*dx - *dp)) + gamma1;
    r = p/q;
    if (r < 0.0 && gamma1 != 0.0) {
      stpc = *stp + r*(*stx - *stp);
    } else if (*stp > *stx) {
      stpc = *stepmax;
    } else     {
      stpc = *stepmin;
    }
    stpq = *stp + (*dp/(*dp-*dx)) * (*stx - *stp);

    if (*bracket) {
      if (fabs(*stp-stpc) < fabs(*stp-stpq)) {
        stpf = stpc;
      } else {
				stpf = stpq;
      }
    } else {
      if (fabs(*stp-stpc) > fabs(*stp-stpq)) {
				stpf = stpc;
      } else {
				stpf = stpq;
      }
    }
  } else {
    /* Case 4: A lower function value, derivatives of the
       same sign, and the magnitude of the derivative does
       not decrease. If the minimum is not bracketed, the step
       is either stpmin or stpmax, else the cubic step is taken. */

    info = 4;
    bound = 0;
    if (*bracket) {
      theta = 3*(*fp - *fy)/(*sty - *stp) + *dy + *dp;
      s = CF_MAX(fabs(theta), fabs(*dy));
      s = CF_MAX(s, fabs(*dp));
      gamma1 = s*sqrt(pow(theta/s,2.0) - (*dy/s)*(*dp/s));
      if (*stp > *sty) {
				gamma1 = -gamma1;
      }      
      p = (gamma1 - *dp) + theta;
      q = ((gamma1 - *dp) + gamma1) + *dy;
      r = p/q;
      stpc = *stp + r*(*sty - *stp);
      stpf = stpc;
    } else if (*stp > *stx) {
      stpf = *stepmax;
    } else {
      stpf = *stepmin;
    }
  }


  /* Update the interval of uncertainty.  This update does not
     depend on the new step or the case analysis above. */

  if (*fp > *fx) {
    *sty = *stp;
    *fy = *fp;
    *dy = *dp;
  } else {
    if (sgnd < 0.0) {
      *sty = *stx;
      *fy = *fx;
      *dy = *dx;
    }
    *stx = *stp;
    *fx = *fp;
    *dx = *dp;
  }

  /* Compute the new step and safeguard it. */
  stpf = CF_MIN(*stepmax,stpf);
  stpf = CF_MAX(*stepmin,stpf);
  *stp = stpf;
  if (*bracket && bound) {
    if (*sty > *stx) {
      *stp = CF_MIN(*stx+0.66*(*sty-*stx),*stp);
    } else {
      *stp = CF_MAX(*stx+0.66*(*sty-*stx),*stp);
    }
  }


  return info;

}
