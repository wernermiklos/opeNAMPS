/*==========================================================
 * Hashtest.c - example in MATLAB External Interfaces
 *
 * Creates task list for NTdot.
 * 
 *
 * The calling syntax is:
 *
 *		bucketlist = NTdot_createtasklist(keylist)
 *
 * This is a MEX-file for MATLAB.
 * Written by M. A. Werner 2022 *
 *========================================================*/

#define VALUE_AB_INITSIZE 10


#include "mex.h"
#include "matrix.h"
#include "simplehash.h"
/* The computational routine */



/* The content is problem specific */




/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){
    size_t m;
    size_t p = 53;
    hash_bucket *hashtable;
    size_t *values;
    size_t i,j;
    size_t keylength;
    uint16_t *key;
    hash_element *myelement;
    mxArray *buckelements;
    double *tmp;
    
    /* Error check */
    if(nrhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:NTdot_createtasklist:nrhs","One input required:  keylist");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:tensordot:nlhs","One output required ");
    }
    
    m = mxGetN(prhs[0]);
    hashtable = create_hash_table(m,1);
    values = (size_t*)malloc(m*sizeof(size_t));
    keylength = mxGetN(mxGetCell(prhs[0],0));
    
    
    /*for(i=0;i<m;++i) mexPrintf("%d n_elements=%d\n",i,hashtable[i].n_elements);*/
    
    for(i=0; i<m; ++i){
        values[i] = i;
        key = (uint16_t *) mxGetData(mxGetCell(prhs[0],i));
        
        myelement = getcreate_element_ptr(hashtable, key, keylength,m,p);
        myelement->value = (void*) &(values[i]);
    }
    plhs[0] = mxCreateCellMatrix(1,m);
    for(i=0;i<m;++i){
        buckelements = mxCreateDoubleMatrix(1,hashtable[i].n_elements,mxREAL);
        tmp = mxGetPr(buckelements);
        for(j=0; j<hashtable[i].n_elements; ++j) tmp[j] = (double) (*((size_t *)hashtable[i].elements[j].value));
        mxSetCell(plhs[0],i,buckelements);
    }
    free(values);
    clear_hashtable(hashtable,m);
    return;
}
    