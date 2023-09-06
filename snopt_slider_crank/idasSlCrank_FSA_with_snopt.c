/* -----------------------------------------------------------------
 * Programmer: Radu Serban and Cosmin Petra @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Simulation of a slider-crank mechanism modelled with 3 generalized
 * coordinates: crank angle, connecting bar angle, and slider location.
 * The mechanism moves under the action of a constant horizontal
 * force applied to the connecting rod and a spring-damper connecting
 * the crank and connecting rod.
 *
 * The equations of motion are formulated as a system of stabilized
 * index-2 DAEs (Gear-Gupta-Leimkuhler formulation).
 *
 * IDAS also computes sensitivities with respect to the problem
 * parameters k (spring constant) and c (damper constant) of the
 * loss that is given by a MSE function with data supplied from another IDAS run
 *   G = \sum_{i=0}^101 0.5 * (y_3^i - \hat{y}_3)^2
 * where \hat{y}_3 is the ground truth data slider position
 * The loss sensitivites are then used to do a parameter identification for the stiffness (K)
 * and the damping (C)
 * For the parameter identification the NLOPT optimiztion package is used with the MMA optimizer
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <idas/idas.h>                 /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

// Includes for the  optmization package
#include <snopt_cwrap.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component i= 1..NEQ */

/* Problem Constants */
#define MAX_ROWS 1000  // Maximum number of rows in the data file
#define NPTS 101 // Number of points in the data

#define NEQ   10
#define NP     2

#define TBEGIN  RCONST(0.0)
#define TEND    RCONST(10.000)

#define RTOLF   RCONST(1.0e-06)
#define ATOLF   RCONST(1.0e-07)

#define RTOLQ   RCONST(1.0e-06)
#define ATOLQ   RCONST(1.0e-08)

#define RTOLFD  RCONST(1.0e-06)
#define ATOLFD  RCONST(1.0e-08)


#define ZERO     RCONST(0.00)
#define QUARTER  RCONST(0.25)
#define HALF     RCONST(0.50)
#define ONE      RCONST(1.00)
#define TWO      RCONST(2.00)
#define FOUR     RCONST(4.00)

typedef struct {
  realtype a;
  realtype J1, J2, m1, m2;
  realtype l0;
  realtype params[2];
  realtype F;
  N_Vector *ground_truth;
} *UserData;



static int fillNVectorsFromFile(const char *fileName, N_Vector *timeVector, N_Vector *sliderPosVector); // Read the ground truth data to add to UserData
static int ressc(realtype tres, N_Vector yy, N_Vector yp,
           N_Vector resval, void *user_data);

static void setIC(N_Vector yy, N_Vector yp, UserData data);
static void force(N_Vector yy, realtype *Q, UserData data);

static int PrintFinalStats(void *mem);
static int check_retval(void *returnvalue, const char *funcname, int opt);
/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */

int main(void)
{
  UserData data;

  void *mem;
  N_Vector yy, yp, id, q, time, sliderPos, *yyS, *ypS, *qS;
  realtype tout, dt, tret;
  realtype pbar[2];
  realtype dp, G, Gm[2], Gp[2];
  int retval, is, nout;
  realtype atolS[NP];
  SUNMatrix A;
  SUNLinearSolver LS;
  SUNContext ctx;

  // Since we don't have a qudrature but a summation we do that computation outside
  // Here is where the variables for those computations is defined
  realtype y3; // Slider position to compute loss
  realtype y3_hat; // Data point
  realtype s3; // Derivative of slider position wrt parameters (k and c)

  A  = NULL;
  LS = NULL;

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(NULL, &ctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return 1;

  id = N_VNew_Serial(NEQ, ctx); // Specifies which equations are differential equations and which are algebraic
  yy = N_VClone(id); // This is the state (y)
  yp = N_VClone(id); // This is state derivative (\dot{y})
  q = N_VNew_Serial(1, ctx); // This *I think* is the quadrature variable associated with the sensitivity equations

  yyS= N_VCloneVectorArray(NP,yy); // these are the sensitivty matrices
  ypS= N_VCloneVectorArray(NP,yp);
  qS = N_VCloneVectorArray(NP, q); // This is the loss sensitivity (\grad{L}{P})

  data = (UserData) malloc(sizeof *data);

  data->a = 0.5;   /* half-length of crank */
  data->J1 = 1.0;  /* crank moment of inertia */
  data->m2 = 1.0;  /* mass of connecting rod */
  data->m1 = 1.0;
  data->J2 = 2.0;  /* moment of inertia of connecting rod */
  data->params[0] = 5.;   /* spring constant */
  data->params[1] = 5.;   /* damper constant */
  data->l0 = 1.0;  /* spring free length */
  data->F = 1.0;   /* external constant force */

  // Create two vectors one for the time and one for the slider position
  time = N_VNew_Serial(NPTS, ctx);
  sliderPos = N_VNew_Serial(NPTS, ctx);
  // Now fill these vectors from the data file that we have
  retval = fillNVectorsFromFile("file3.txt", &time, &sliderPos);

  if (retval != 0) {
      printf("Failed to fill N_Vectors from file.\n");
      return 1;
  }
  // create the vector array for the ground truth data
  data->ground_truth = N_VNewVectorArray(2);
  // assign the vectors that we just filled to these vector arrays
  data->ground_truth[0] = time;
  data->ground_truth[1] = sliderPos;

  // This tells that the 7th to 10th equation is algebraic equation
  N_VConst(ONE, id);
  NV_Ith_S(id, 9) = ZERO;
  NV_Ith_S(id, 8) = ZERO;
  NV_Ith_S(id, 7) = ZERO;
  NV_Ith_S(id, 6) = ZERO;

  // printf("\nSlider-Crank example for IDAS:\n");

  /* Consistent IC*/
  setIC(yy, yp, data);

  // Set the loss to 0
  for (is=0;is<NP;is++) {
    N_VConst(ZERO, yyS[is]);
    N_VConst(ZERO, ypS[is]);
    
  }

  /* IDA initialization */
  mem = IDACreate(ctx);
  retval = IDAInit(mem, ressc, TBEGIN, yy, yp);
  retval = IDASStolerances(mem, RTOLF, ATOLF);
  retval = IDASetUserData(mem, data);
  retval = IDASetId(mem, id);
  retval = IDASetSuppressAlg(mem, SUNTRUE);
  retval = IDASetMaxNumSteps(mem, 20000);

  /* Create dense SUNMatrix for use in linear solves -> This is a "template Jacobian matrix". I think it takes the number of equations */ 
  A = SUNDenseMatrix(NEQ, NEQ, ctx);
  if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(yy, A, ctx);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(mem, LS, A);
  if(check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

  retval = IDASensInit(mem, NP, IDA_SIMULTANEOUS, NULL, yyS, ypS); // Activates forward sensitivity computations
  pbar[0] = data->params[0];pbar[1] = data->params[1]; //DOUBT- What are these pbars? Documentation says they are some sort of scaling factors
  retval = IDASetSensParams(mem, data->params, pbar, NULL);
  retval = IDASensEEtolerances(mem); // Sets same tolerance for sensitivity equations as that for the state equations
  IDASetSensErrCon(mem, SUNTRUE);

  int count_fun_eval = 0;
  int count_fun_with_grad_eval = 0;
  int count_fun_wihtout_grad_eval = 0;


  // Function that does the simulation and provides the necessary jacobian
  void slider_crank_fg(int *Status, int *n, double x[],
                          int *needF, int *nF, double F[],
                          int *needG, int *neG, double G[],
                          char cu[], int *lencu,
                          int iu[], int *leniu,
                          double ru[], int *lenru){
      
      ++count_fun_eval;
      // Run the forward simulation by reinitializing the original problem and the sensitivities
      data->params[0]  = x[0]; // Change the parameters to the new ones
      data->params[1]  = x[1];
      
      // Set the sensitivites to 0
      for (is=0;is<NP;is++) {
        N_VConst(ZERO, yyS[is]);
        N_VConst(ZERO, ypS[is]);
      }
      // Initialize the loss quadratures to 0
      N_VConst(ZERO, q);
      for (is=0; is<NP; is++) N_VConst(ZERO, qS[is]);

      setIC(yy, yp, data); // reset the initial condtions

      IDAReInit(mem, TBEGIN, yy, yp); // Reinitialize the problem
      if(*needG > 0){
        retval = IDASensReInit(mem, IDA_SIMULTANEOUS, yyS, ypS); // reinitialize the sensitivities if the optimizer needs it
      }
      else{
        IDASensToggleOff(mem); // Toggle sensitivity off if the optimizer does not need the sensitivites
      }
      

      /* Perform forward run */
      dt = 0.1; // The states are to be printed every 0.1 seconds
      nout = (TEND - TBEGIN) / dt + 1;
      tout = dt;
      for (int iout=1; iout<nout; iout++) {
        tout = iout*dt;
        retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);
        if (retval < 0) break;


        // Here we get the loss value at this time and add it to the total loss
        y3 = Ith(yy, 3);
        y3_hat = Ith(data->ground_truth[1], iout+1); // +1 is beause inital conditions stored in the data file

        // add to the loss (basically doing "integration" here)
        Ith(q, 1) += HALF * (y3 - y3_hat) * (y3 - y3_hat);
        if(*needG > 0){
          IDAGetSens(mem, &tret, yyS);
          // Get here the sensitivity that comes up in the chain rule for the cost -> Since we are using 
          // sensitivity of state with respect to stiffness (K)
          s3 = Ith(yyS[0],3); 
          
          Ith(qS[0],1) += (y3 - y3_hat) * s3;
          // sensitivity of state with respect to damping (C)
          s3 = Ith(yyS[1],3);
          Ith(qS[1],1) += (y3 - y3_hat) * s3;
        }
      }

      if(*needF > 0){
        F[0] = Ith(q, 1);

      }
   
      if(*needG > 0){
        // Jacobian just two elements
        G[0] = Ith(qS[0],1);
        G[1] = Ith(qS[1],1);
        ++count_fun_with_grad_eval;
      }
      else{
        ++count_fun_wihtout_grad_eval;
      }
  }

  snProblem slider_crank; // Needs to be defined for the optimizer
  int    i, info;
  int    Cold   =  0; // Cold start (don't know what this exactly means)
  int    n      =  2; // Number of parameters to optimzie
  int    nF     =  1; // Length of the cost function vector , in our case its just a scalar
  int    lenA   = 0, neA; // We don't have any A matrix
  int    lenG   = 2, neG; // Number of elements in our jacobian matrix
  int    ObjRow =  0; // Since we only have 1 row
  double ObjAdd =  0; // We don't need to add any constat to the objective row for printing purposes

  int    iGfun[2], jGvar[2];

  int    xstate[2]; // This is the set of initial values for each state -> wierd
  double x[2], xlow[2], xupp[2], xmul[2];

  int    Fstate[1]; // Again set of initial states for the loss function -> needs to be defined for cold start
  double F[1], Flow[1], Fupp[1], Fmul[1]; // Initial value for the loss function -> needs to be defined for cold start

  int    nS, nInf; // "final number of superbasic varaibles"
  double sInf; // "number and sum of the infeasibliities of constarins that lie outside one of their bounds"
  
  // snInit must be called first.
  //   9, 6 are print and summary unit numbers (for Fortran).
  //   6 == standard out
  snInit( &slider_crank, "SlCrank", "SlCrank_snopt.out", 1);
  


  // User workspace allocated and set.
  //   May be accesed in the user-defined function -> Don't think we need this becasue we don't have any thing additional to defince
  slider_crank.leniu = 2;
  slider_crank.iu = (int*)malloc( sizeof(int) * slider_crank.leniu );
  slider_crank.iu[0] = 0;
  slider_crank.iu[1] = 1;
  

  // Set bounds
  xlow[0] = 1e-6;  xupp[0] = 1e20;
  xlow[1] = 1e-6;  xupp[1] = 1e20;

  Flow[0] = -1e20;  Fupp[0] = 1e20;

  for ( i = 0; i < n; i++ ) {
    xstate[i] = 0;
    x[i]      = 0;
    xmul[i]   = 0;
  }

  for ( i = 0; i < nF; i++ ) {
    Fstate[i] = 0;
    F[i]      = 0;
    Fmul[i]   = 0;
  } 

  //initial value of K and C
  x[0] = 10.;
  x[1] = 10.;

  // Set up the Jacobian matrix -> (i,j) elements are non zero
  iGfun[0] = 0;
  jGvar[0] = 0;

  iGfun[1] = 0;
  jGvar[1] = 1;
  neA      = 0;
  neG = 2;
  


  // Read options -> This file isnt there in the tutorials folder so no iddea what this is and where it comes from
  // info = setSpecsfile( &slider_crank, "slCrank.spc" ); 

  // setIntParameter( &slider_crank, "Verify level", 3 );
  setIntParameter( &slider_crank, "Derivative option", 1 );
  info = solveA( &slider_crank, Cold,
		 nF, n, ObjAdd, ObjRow, slider_crank_fg,
		 neA, NULL, NULL, NULL,
		 neG, iGfun, jGvar,
		 xlow, xupp, Flow, Fupp,
		 x, xstate, xmul, F, Fstate, Fmul,
		 &nS, &nInf, &sInf);

  // Deallocate space.
  free( slider_crank.iu );
  deleteSNOPT( &slider_crank );

  

  /* Free memory */
  IDAFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);


  free(data);

  N_VDestroy(id);
  N_VDestroy(yy);
  N_VDestroy(yp);
  N_VDestroy(q);

  N_VDestroyVectorArray(yyS, NP);
  N_VDestroyVectorArray(ypS, NP);
  N_VDestroyVectorArray(qS,  NP);

  SUNContext_Free(&ctx);

  return(0);
}



static int fillNVectorsFromFile(const char *fileName, N_Vector *timeVector, N_Vector *sliderPosVector) {
    FILE *file = fopen(fileName, "r");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    double *time;
    double *sliderPos;
    // Get the data arrays and fill them up
    time = N_VGetArrayPointer(*timeVector);
    sliderPos = N_VGetArrayPointer(*sliderPosVector);
    char line[100];

    // Skip the first line (header)
    fgets(line, sizeof(line), file);

    int row = 0;
    while (fgets(line, sizeof(line), file)) {
        double t, y1, y2, y3;
        if (sscanf(line, "%lf, %lf, %lf, %lf", &t, &y1, &y2, &y3) != 4) {
            printf("Error parsing line %d\n", row + 2);
            fclose(file);
            return 1;
        }

        time[row] = t;
        sliderPos[row] = y3;

        row++;

        if (row > NPTS) {
            printf("Warning: Maximum number of rows reached. rows = %d\n", row);
            break;
        }
    }

    fclose(file);
    

    return 0;
}


static void setIC(N_Vector yy, N_Vector yp, UserData data)
{
  realtype pi;
  realtype a, J1, m2, J2;
  realtype q, p, x;
  realtype Q[3];

  N_VConst(ZERO, yy);
  N_VConst(ZERO, yp);

  pi = FOUR*atan(ONE);

  a = data->a;
  J1 = data->J1;
  m2 = data->m2;
  J2 = data->J2;

  q = pi/TWO;
  p = asin(-a);
  x = cos(p);

  NV_Ith_S(yy,0) = q;
  NV_Ith_S(yy,1) = x;
  NV_Ith_S(yy,2) = p;

  force(yy, Q, data);

  NV_Ith_S(yp,3) = Q[0]/J1;
  NV_Ith_S(yp,4) = Q[1]/m2;
  NV_Ith_S(yp,5) = Q[2]/J2;

}

static void force(N_Vector yy, realtype *Q, UserData data)
{
  realtype a, k, c, l0, F;
  realtype q, x, p;
  realtype qd, xd, pd;
  realtype s1, c1, s2, c2, s21, c21;
  realtype l2, l, ld;
  realtype f, fl;

  a = data->a;
  k = data->params[0];
  c = data->params[1];
  l0 = data->l0;
  F = data->F;

  q = NV_Ith_S(yy,0);
  x = NV_Ith_S(yy,1);
  p = NV_Ith_S(yy,2);

  qd = NV_Ith_S(yy,3);
  xd = NV_Ith_S(yy,4);
  pd = NV_Ith_S(yy,5);

  s1 = sin(q);
  c1 = cos(q);
  s2 = sin(p);
  c2 = cos(p);
  s21 = s2*c1 - c2*s1;
  c21 = c2*c1 + s2*s1;

  l2 = x*x - x*(c2+a*c1) + (ONE + a*a)/FOUR + a*c21/TWO;
  l = sqrt(l2);
  ld = TWO*x*xd - xd*(c2+a*c1) + x*(s2*pd+a*s1*qd) - a*s21*(pd-qd)/TWO;
  ld /= TWO*l;

  f = k*(l-l0) + c*ld;
  fl = f/l;

  Q[0] = - fl * a * (s21/TWO + x*s1) / TWO;
  Q[1] = fl * (c2/TWO - x + a*c1/TWO) + F;
  Q[2] = - fl * (x*s2 - a*s21/TWO) / TWO - F*s2;

}
// This is the residual function
static int ressc(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
{
  UserData data;
  realtype Q[3];
  realtype a, J1, m2, J2;
  realtype *yval, *ypval, *rval;
  realtype q, x, p;
  realtype qd, xd, pd;
  realtype lam1, lam2, mu1, mu2;
  realtype s1, c1, s2, c2;

  data = (UserData) user_data;

  a  = data->a;
  J1 = data->J1;
  m2 = data->m2;
  J2 = data->J2;

  yval = N_VGetArrayPointer(yy);
  ypval = N_VGetArrayPointer(yp);
  rval = N_VGetArrayPointer(rr);

  q = yval[0];
  x = yval[1];
  p = yval[2];

  qd = yval[3];
  xd = yval[4];
  pd = yval[5];

  lam1 = yval[6];
  lam2 = yval[7];

  mu1 = yval[8];
  mu2 = yval[9];

  s1 = sin(q);
  c1 = cos(q);
  s2 = sin(p);
  c2 = cos(p);

  force(yy, Q, data);

  rval[0] = ypval[0] - qd + a*s1*mu1 - a*c1*mu2;
  rval[1] = ypval[1] - xd + mu1;
  rval[2] = ypval[2] - pd + s2*mu1 - c2*mu2;

  rval[3] = J1*ypval[3] - Q[0] + a*s1*lam1 - a*c1*lam2;
  rval[4] = m2*ypval[4] - Q[1] + lam1;
  rval[5] = J2*ypval[5] - Q[2] + s2*lam1 - c2*lam2;

  rval[6] = x - c2 - a*c1;
  rval[7] = -s2 - a*s1;

  rval[8] = a*s1*qd + xd + s2*pd;
  rval[9] = -a*c1*qd - c2*pd;

  return(0);
}


static int PrintFinalStats(void *mem)
{
  int retval;
  long int nst, nni, nnf, nje, nre, nreLS, netf, ncfn;

  retval = IDAGetNumSteps(mem, &nst);
  retval = IDAGetNumResEvals(mem, &nre);
  retval = IDAGetNumJacEvals(mem, &nje);
  retval = IDAGetNumNonlinSolvIters(mem, &nni);
  retval = IDAGetNumErrTestFails(mem, &netf);
  retval = IDAGetNumNonlinSolvConvFails(mem, &nnf);
  retval = IDAGetNumStepSolveFails(mem, &ncfn);
  retval = IDAGetNumLinResEvals(mem, &nreLS);

  printf("\nFinal Run Statistics: \n\n");
  printf("Number of steps                    = %ld\n", nst);
  printf("Number of residual evaluations     = %ld\n", nre+nreLS);
  printf("Number of Jacobian evaluations     = %ld\n", nje);
  printf("Number of nonlinear iterations     = %ld\n", nni);
  printf("Number of error test failures      = %ld\n", netf);
  printf("Number of nonlinear conv. failures = %ld\n", nnf);
  printf("Number of step solver failures     = %ld\n", ncfn);

  return(retval);
}


static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}
