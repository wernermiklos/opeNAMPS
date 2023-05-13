/*==========================================================
 * tensordot.c - contraction routine for multidim (>2) arrays in MatLab.
 *
 * contracts one or more pairs of legs of the two input tensors
 *
 * The calling syntax is:
 *
 *		outtensor = tensordot(A, B, contlegsA, contlegsB, freelegsA, freelegsB)
 *
 *          A:            tensor1   (MatLab multidim array)
 *          B:            tensor2   (MatLab multidim array)
 *          contlegsA:    list of contracted legs in A (row vector)      %int!!!
 *          contlegsB:    list of contracted legs in B (row vector)   the ordering should match the one in contlegsA     %int!!!
 *          freelegsA:    list of free (not contracted) legs in A (row vector). Supposed to be in ascending order.       %int!!!
 *          freelegsB:    list of free (not contracted) legs in B (row vector). Supposed to be in ascending order.       %int!!!
 * This is a MEX-file for MATLAB (R2017a).    
 * Written by M. A. Werner - Budapest 2022
 *
 *========================================================*/

#include "mex.h"
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>

double timediff(struct timeval end, struct timeval start){
    return (double)1000000*(end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
}




void coredotadd(double *Atmp, double *Btmp, double *Ctmp, double fact, mwSize totdim_alpha, mwSize totdim_alpha0, mwSize contdim_alpha, mwSize totdim_beta, mwSize totdim_delta,
             mwSize* alpha_to_A, mwSize* delta_to_A, mwSize* delta_to_B, mwSize* beta_to_B){
    
    mwSize alpha, alpha0, beta, delta, flposB, clposA, posA0, Aind, Cind, j;
    double Bval_tmp;
    
    for (beta = 0; beta < totdim_beta; ++beta){
        flposB = beta_to_B[beta];
        for (delta = 0; delta < totdim_delta; ++delta){
            clposA = delta_to_A[delta];
            if(flposB + delta_to_B[delta] > totdim_delta*totdim_beta){
                mexPrintf("error 1, B: %d\n", flposB + delta_to_B[delta]); return;   
            }
            Bval_tmp = fact*Btmp[flposB + delta_to_B[delta]];
            for(alpha0=0; alpha0 < totdim_alpha; alpha0 += contdim_alpha){    /*could be parfor here? */ 
                posA0 = clposA + alpha_to_A[alpha0];
                for(j=0; j <contdim_alpha; ++j){
                    Aind = posA0 + j;
                    Cind = alpha0 + j + totdim_alpha*beta;
                    /*mexPrintf("%d %d %d %d %d \n",contdim_alpha, j,alpha0 + j, beta, Cind);*/
                    /*if(Aind> totdim_delta*totdim_alpha){
                        mexPrintf("error 2\n"); return;   
                    }
                    if(Cind> totdim_beta*totdim_alpha){
                        mexPrintf("error 3\n"); return;   
                    }*/
                    
                    Ctmp[Cind] += Atmp[Aind]*Bval_tmp;
                }
            }
        }
    }
}

/* The computational routine */
void tensordot(double *Areal, double *Aimag, double *Breal, double *Bimag, double *Creal, double *Cimag, int *contlegsA, int *contlegsB, int *freelegsA, int *freelegsB, 
                             mwSize numcontlegs, mwSize numfreelegsA, mwSize numfreelegsB, const mwSize *truedimsA, const mwSize *truedimsB)
{
    mwSize alpha, alpha0, Aind, Cind;              /* superindex for free legs in A */
    mwSize totdim_alpha = 1;   /* will be set */   
    mwSize beta;               /* superindex for free legs in B */
    mwSize totdim_beta = 1;        
    mwSize gamma;              /* superindex for free legs in C */
    mwSize totdim_gamma = 1;
    mwSize delta;              /* superindex for contracted legs */
    mwSize totdim_delta = 1;
    mwSize i;                  /* loop variable */
    mwSize *dimfactorsA;       /* dimension factors of A, must be allocated and freed */
    mwSize *dimfactorsB;       /* dimension factors of B, must be allocated and freed */
    mwSize *alpha_to_A;        
    mwSize *delta_to_A;
    mwSize *delta_to_B;
    mwSize *beta_to_B;
    mwSize contdim_alpha = 1;
    mwSize totdim_alpha0;
    mwSize tmpsuperind;
    mwSize tmplegind;
    int A_is_complex = 0;
    int B_is_complex = 0;
    int C_is_complex;
    double *Btmp, *Atmp, *Ctmp;   /*temporary pointers (no allocation needed!) */
    double fact;
   
    
    /* first we determine if A/B/C is complex or not */
    if(Aimag != NULL) A_is_complex = 1;
    if(Bimag != NULL) B_is_complex = 1;
    C_is_complex = A_is_complex || B_is_complex;
    
    /* we determine the total dimension of superindices alpha, beta, delta */
    
    for(i=0; i<numfreelegsA; ++i){
        totdim_alpha *= truedimsA[freelegsA[i]-1];                        
        if(freelegsA[i] == i + 1) contdim_alpha = totdim_alpha;   /* until freelegs come as 1 2 3 4 these are continous dimensions... */
    }
    totdim_alpha0 = totdim_alpha / contdim_alpha;
    for(i=0; i<numfreelegsB; ++i){
        totdim_beta *= truedimsB[freelegsB[i]-1];                        
    }
    for(i=0; i<numcontlegs; ++i){
        totdim_delta *= truedimsA[contlegsA[i]-1];                        
    }
    alpha_to_A = (mwSize *)malloc(totdim_alpha*sizeof(mwSize));
    delta_to_A = (mwSize *)malloc(totdim_delta*sizeof(mwSize));
    delta_to_B = (mwSize *)malloc(totdim_delta*sizeof(mwSize));
    beta_to_B = (mwSize *)malloc(totdim_beta*sizeof(mwSize));
           
    
    
    /* we allocate memory for the temporary 'dimfactors' arrays */
    dimfactorsA = (mwSize *)malloc((numfreelegsA + numcontlegs) * sizeof(mwSize));
    dimfactorsB = (mwSize *)malloc((numfreelegsB + numcontlegs) * sizeof(mwSize));
    
    /* We determine dimfactorsA dimfactorsB*/
    dimfactorsA[0] = 1;
    for(i = 1; i < numfreelegsA + numcontlegs; ++i){
        dimfactorsA[i] = dimfactorsA[i-1] * truedimsA[i-1];
    }
    
    dimfactorsB[0] = 1;
    for(i = 1; i < numfreelegsB + numcontlegs; ++i){
        dimfactorsB[i] = dimfactorsB[i-1] * truedimsB[i-1];
    }
    
    for (alpha = 0; alpha < totdim_alpha; ++alpha){
        tmpsuperind = alpha;
        alpha_to_A[alpha] = 0;
        for (i = 0; i < numfreelegsA; ++i){
            tmplegind = tmpsuperind % truedimsA[freelegsA[i]-1];
            alpha_to_A[alpha] += dimfactorsA[freelegsA[i]-1] * tmplegind;
            tmpsuperind = (tmpsuperind - tmplegind) / truedimsA[freelegsA[i]-1];
        }
        /*mexPrintf("alpha: %d, A: %d, Amax: %d\n",alpha, alpha_to_A[alpha], totdim_alpha*totdim_delta);*/
    }
    
    for (delta = 0; delta < totdim_delta; ++delta){
        tmpsuperind = delta;
        delta_to_A[delta] = 0;
        for (i = 0; i < numcontlegs; ++i){
            tmplegind = tmpsuperind % truedimsA[contlegsA[i]-1];
            delta_to_A[delta] += dimfactorsA[contlegsA[i]-1] * tmplegind;
            tmpsuperind = (tmpsuperind - tmplegind) / truedimsA[contlegsA[i]-1];
        }
        /*mexPrintf("delta: %d, A: %d, Amax: %d\n",delta, delta_to_A[delta],totdim_alpha*totdim_delta);*/
    }
    
    for (delta = 0; delta < totdim_delta; ++delta){
        tmpsuperind = delta;
        delta_to_B[delta] = 0;
        for (i = 0; i < numcontlegs; ++i){
            tmplegind = tmpsuperind % truedimsB[contlegsB[i]-1];
            delta_to_B[delta] += dimfactorsB[contlegsB[i]-1] * tmplegind;
            tmpsuperind = (tmpsuperind - tmplegind) / truedimsB[contlegsB[i]-1];
        }
        /*mexPrintf("delta: %d, B: %d, Bmax: %d\n",delta, delta_to_B[delta], totdim_delta*totdim_beta);*/
    }
    
    for (beta = 0; beta < totdim_beta; ++beta){
        tmpsuperind = beta;
        beta_to_B[beta] = 0;
        for (i = 0; i < numfreelegsB; ++i){
            tmplegind = tmpsuperind % truedimsB[freelegsB[i]-1];
            beta_to_B[beta] += dimfactorsB[freelegsB[i]-1] * tmplegind;
            tmpsuperind = (tmpsuperind - tmplegind) / truedimsB[freelegsB[i]-1];
        }
        /*mexPrintf("beta: %d, B: %d, Bmax: %d\n",beta, beta_to_B[beta], totdim_delta*totdim_beta);*/
    }
    
    for (gamma = 0; gamma < totdim_alpha*totdim_beta; ++gamma){
        Creal[gamma] = 0;
        if(C_is_complex) Cimag[gamma] = 0;
    }
    
    /* real-real */
    Btmp = Breal; Atmp = Areal; Ctmp = Creal; fact = 1;
    coredotadd(Atmp, Btmp, Ctmp, fact, totdim_alpha, totdim_alpha0, contdim_alpha, totdim_beta, totdim_delta, alpha_to_A, delta_to_A, delta_to_B, beta_to_B);
    /* real-imag */
    if(B_is_complex){
        Btmp = Bimag; Atmp = Areal; Ctmp = Cimag; fact = 1;
        coredotadd(Atmp, Btmp, Ctmp, fact, totdim_alpha, totdim_alpha0, contdim_alpha, totdim_beta, totdim_delta, alpha_to_A, delta_to_A, delta_to_B, beta_to_B);
    }
    /* imag-real */
    if(A_is_complex){
        Btmp = Breal; Atmp = Aimag; Ctmp = Cimag; fact = 1;
        coredotadd(Atmp, Btmp, Ctmp, fact, totdim_alpha, totdim_alpha0, contdim_alpha, totdim_beta, totdim_delta, alpha_to_A, delta_to_A, delta_to_B, beta_to_B);
    }
    /* imag-imag */
    if(A_is_complex && B_is_complex){
        Btmp = Bimag; Atmp = Aimag; Ctmp = Creal; fact = -1;
        coredotadd(Atmp, Btmp, Ctmp, fact, totdim_alpha, totdim_alpha0, contdim_alpha, totdim_beta, totdim_delta, alpha_to_A, delta_to_A, delta_to_B, beta_to_B);
    }
    
    
    
    free(alpha_to_A);
    free(beta_to_B);
    free(delta_to_A);
    free(delta_to_B);
    free(dimfactorsA);
    free(dimfactorsB);
}




/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *Areal;              /* tensor1 */
    double *Aimag = NULL;
    double *Breal;              /* tensor2 */
    double *Bimag = NULL;
    double *Creal;              /* output tensor */
    double *Cimag = NULL;
    const mwSize *dimsA;           /* array of leg dims in A (matlab screws up 1dim legs) */
    mwSize *truedimsA;       /* true leg dims in A */
    const mwSize *dimsB;           /* array of leg dims in B */
    mwSize *truedimsB;       /* true leg dims in B */
    mwSize *dimsC;          /* array of out dimensions */
    int *contlegsA;          /* array of contracted legIDs in A */
    int *contlegsB;          /* array of contracted legIDs in B */
    int *freelegsA;          /* array of free legIDs in A */
    int *freelegsB;          /* array of free legIDs in B */
    mwSize numcontlegs;          /* number of contracted legs */
    mwSize numfreelegsA;         /* number of free legs in A */
    mwSize numfreelegsB;         /* number of free legs in B */
    mwSize numlegsA;
    mwSize numlegsAnaive;
    mwSize numlegsB;
    mwSize numlegsBnaive;
    mwSize numlegsC;
    mwSize i;
    int complexoutput = 0;
    /*struct timeval starttime, midtime, endtime;*/

    /*gettimeofday(&starttime, NULL);*/
    
    /* TEST NEEDED IF FREELEGS_A AND FREELEGS_B ARE SORTED!!! */   /*STILL!!!*/
    
    /* check for proper number of arguments */
    
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("MyToolbox:tensordot:nrhs","Six inputs required:  A, B, contlegsA, contlegsB, freelegsA, freelegsB");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:tensordot:nlhs","One output required.  (call the function with an explicit rhs: <someOutputVariable> = tensordot(...) )");
    }
    
    /* Check if inputs 3rd-6th inputs are int row vectors */
    if((mxGetM(prhs[2]) > 1) ||
       (mxGetM(prhs[3]) > 1) ||
       (mxGetM(prhs[4]) > 1) ||
       (mxGetM(prhs[5]) > 1) ||
       (!mxIsInt32(prhs[2])) ||
       (!mxIsInt32(prhs[3])) ||
       (!mxIsInt32(prhs[4])) ||
       (!mxIsInt32(prhs[5]))){
        mexErrMsgIdAndTxt("MyToolbox:tensordot:nrhs","contlegsA, contlegsB, freelegsA, freelegsB must be int row vectors");  
    }

    contlegsA = (int *) mxGetData(prhs[2]);
    contlegsB = (int *) mxGetData(prhs[3]);
    freelegsA = (int *) mxGetData(prhs[4]);
    freelegsB = (int *) mxGetData(prhs[5]);
    
    
    numcontlegs = mxGetN(prhs[2]);
    if( numcontlegs != mxGetN(prhs[3])){
         mexErrMsgIdAndTxt("MyToolbox:tensordot:nrhs","contlegsA and contlegsB must have same length.");  
    }
    
    numfreelegsA = mxGetN(prhs[4]);
        
    numfreelegsB = mxGetN(prhs[5]);
    
    numlegsA = numfreelegsA + numcontlegs;
    numlegsAnaive = mxGetNumberOfDimensions(prhs[0]);
    numlegsB = numfreelegsB + numcontlegs;
    numlegsBnaive = mxGetNumberOfDimensions(prhs[1]);
    numlegsC = numfreelegsA + numfreelegsB;
    
    dimsA = mxGetDimensions(prhs[0]);
    dimsB = mxGetDimensions(prhs[1]);
    truedimsA = (mwSize *)malloc(numlegsA * sizeof(mwSize));
    for(i=0; i<numlegsA; ++i){
       if(i < numlegsAnaive) truedimsA[i] = dimsA[i];
       else truedimsA[i] = 1;
       /*mexPrintf("A %d %d\n", i, truedimsA[i]);*/
    }
    truedimsB = (mwSize *)malloc(numlegsB * sizeof(mwSize));
    for(i=0; i<numlegsB; ++i){
       if(i < numlegsBnaive) truedimsB[i] = dimsB[i];
       else truedimsB[i] = 1;
       /*mexPrintf("B %d %d\n", i, truedimsB[i]);*/
    }
    
    if(numlegsC > 0){
        dimsC = (mwSize *)malloc(numlegsC * sizeof(mwSize));
        for(i=0; i<numfreelegsA; ++i) {dimsC[i] = truedimsA[freelegsA[i]-1]; /*mexPrintf("%d %d\n", i, dimsC[i]);*/}
        
        for(i=0; i<numfreelegsB; ++i) {dimsC[numfreelegsA + i] = truedimsB[freelegsB[i]-1]; /*mexPrintf("%d %d\n", i+numfreelegsA, dimsC[i+numfreelegsA]);*/}
    }
    else{
        dimsC = (mwSize *)malloc(sizeof(mwSize));
        dimsC[0] = 1;
    }
    
    
    /* Check if A & B are double or complex matrices, and gather pointers for the data */
    if( mxIsDouble(prhs[0]) &&
        mxIsDouble(prhs[1])) {
        Areal = mxGetPr(prhs[0]);
        Breal = mxGetPr(prhs[1]);
        if(mxIsComplex(prhs[0])){
            complexoutput = 1;
            Aimag = mxGetPi(prhs[0]);
        }
        if(mxIsComplex(prhs[1])){
            complexoutput = 1;
            Bimag = mxGetPi(prhs[1]);
        }
    }
    else{
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "A and B tensors should be double or complex");
    }
    
    
    /* create the output matrix */
    if(numlegsC > 0){
        plhs[0] = mxCreateNumericArray(numlegsC, dimsC, mxDOUBLE_CLASS, complexoutput? mxCOMPLEX : mxREAL);
    }
    else{
        plhs[0] = mxCreateNumericArray(1, dimsC, mxDOUBLE_CLASS, complexoutput? mxCOMPLEX : mxREAL);
    }

    /* get a pointer to the data in the output matrix */
    Creal = mxGetPr(plhs[0]);
    if(complexoutput) Cimag = mxGetPi(plhs[0]);

    /* call the computational routine */
    tensordot(Areal, Aimag, Breal, Bimag, Creal, Cimag, contlegsA, contlegsB, freelegsA, freelegsB, 
              numcontlegs, numfreelegsA, numfreelegsB, truedimsA, truedimsB);
    
    
    free(truedimsA);
    free(truedimsB);
    free(dimsC);
}