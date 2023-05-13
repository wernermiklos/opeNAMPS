/*==========================================================
 * NTdot_createtasklist.c - example in MATLAB External Interfaces
 *
 * Creates task list for NTdot.
 * 
 *
 * The calling syntax is:
 *
 *		[keylist_C, tasklist] = NTdot_createtasklist(keylist_A, keylist_B, matchlayout_A, matchlayout_B, outlayout_A, outlayout_B, outkeylength)
 *
 * This is a MEX-file for MATLAB.
 * Written by M. A. Werner 2022 *
 *========================================================*/

#define VALUE_AB_INITSIZE 10


#include "mex.h"
#include "matrix.h"
#include "simplehash.h"
#include <stdio.h>
/* The computational routine */



/* The content is problem specific */
typedef struct Value_AB {
    size_t *bIDs_A;
    size_t length_A;
    size_t fill_A;
    size_t *bIDs_B;
    size_t length_B;
    size_t fill_B;
} value_AB;

void allocate_ABelement(hash_element* element, size_t initlength){
    element->value = (void *)malloc(sizeof(value_AB));
    ((value_AB*)(element->value))->length_A = initlength;
    ((value_AB*)(element->value))->length_B = initlength;
    ((value_AB*)(element->value))->fill_A = 0;
    ((value_AB*)(element->value))->fill_B = 0;
    ((value_AB*)(element->value))->bIDs_A = (size_t *)malloc(initlength * sizeof(size_t));
    ((value_AB*)(element->value))->bIDs_B = (size_t *)malloc(initlength * sizeof(size_t));
}

void free_ABelement(hash_element* element){
    free(((value_AB*)(element->value))->bIDs_A);
    free(((value_AB*)(element->value))->bIDs_B);
    free(element->value);
}

void add_bID_A_tolist(value_AB* value, size_t bID_A){
    size_t *newbIDs_A;
    size_t i;
    if(value->fill_A == value->length_A){    /* there is no place for the new bID_A -- we allocate more memory */
        newbIDs_A = (size_t *)malloc(2 * value->length_A * sizeof(size_t));
        for(i=0; i<value->fill_A; ++i) newbIDs_A[i] = value->bIDs_A[i];
        free(value->bIDs_A);
        value->bIDs_A = newbIDs_A;
        value->length_A *= 2;
    }
    value->bIDs_A[value->fill_A] = bID_A;
    value->fill_A += 1;
}

void add_bID_B_tolist(value_AB* value, size_t bID_B){
    size_t *newbIDs_B;
    size_t i;
    if(value->fill_B == value->length_B){    /* there is no place for the new bID_B -- we allocate more memory */
        newbIDs_B = (size_t *)malloc(2 * value->length_B * sizeof(size_t));
        for(i=0; i<value->fill_B; ++i) newbIDs_B[i] = value->bIDs_B[i];
        free(value->bIDs_B);
        value->bIDs_B = newbIDs_B;
        value->length_B *= 2;
    }
    value->bIDs_B[value->fill_B] = bID_B;
    value->fill_B += 1;
}

void cleanup_AB(hash_bucket *hash_table, size_t num_of_buckets){
    size_t bucketID,j;
    for(bucketID = 0; bucketID < num_of_buckets; ++bucketID){
        for(j=0;j<hash_table[bucketID].n_elements;++j) free_ABelement(&(hash_table[bucketID].elements[j]));
    }
    return;
}

void create_outkey(uint16_t *outkey_tmp, 
                   uint16_t *keyA, 
                   uint16_t *keyB, 
                   size_t keylength_A, 
                   size_t keylength_B, 
                   int32_t *outlayout_A, 
                   int32_t *outlayout_B){
    size_t i;
    for(i=0; i<keylength_A; ++i) if(outlayout_A[i] > 0) outkey_tmp[outlayout_A[i]-1] = keyA[i];
    for(i=0; i<keylength_B; ++i) if(outlayout_B[i] > 0) outkey_tmp[outlayout_B[i]-1] = keyB[i];
    return;
}

        


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){
    size_t m_matchcomb, m_outkeys, p = 53;                                  /*parameters of the hash function. p=53 is a reasonable prime number */
    size_t *outbIDs;
    size_t numblocks_A, numblocks_B, numblocks_C;
    size_t initbucketsize=1, initABlistlength=1;
    size_t matchkeylength, outkeylength, keylength_A, keylength_B;
    size_t numtasks = 0, numoutblocks=0;
    mxArray *keylist_A;
    mxArray *keylist_B;
    mxArray *tmp_Ckey;
    uint16_t *tmp_Ckey_data;
    int32_t *matchlayout_A;
    int32_t *matchlayout_B;
    int32_t *outlayout_A;
    int32_t *outlayout_B;
    uint16_t *matchkey_tmp;
    uint16_t *outkey_tmp;
    uint16_t *keyA;
    uint16_t *keyB;
    hash_bucket *matchcomb_table;
    hash_bucket *outkey_table;
    hash_element *element_tmp;
    size_t i,j;
    size_t outdims[2];
    size_t bID_A, bID_B, bucketID, elementID, taskID=0;
    value_AB value_tmp;
   
    
    /* Error check */
    if(nrhs!=7) {
        mexErrMsgIdAndTxt("MyToolbox:NTdot_createtasklist:nrhs","Six inputs required:  keylistA, keylistB, matchlayout_A, matchlayoutB, outlayout_A, outlayout_B, outkeylength");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:NTdot_createtasklist:nlhs","Two output required.  (call the function with an explicit rhs: <someOutputVariable> = tensordot(...) )");
    }
    
    if((mxGetM(prhs[2]) > 1) ||
       (mxGetM(prhs[3]) > 1) ||
       (mxGetM(prhs[4]) > 1) ||
       (mxGetM(prhs[5]) > 1) ||
       (mxGetM(prhs[6]) > 1) ||
       (!mxIsInt32(prhs[2])) ||
       (!mxIsInt32(prhs[3])) ||
       (!mxIsInt32(prhs[4])) ||
       (!mxIsInt32(prhs[5])) ||
       (!mxIsInt32(prhs[6]))){
        mexErrMsgIdAndTxt("MyToolbox:NTdot_createtasklist:nrhs","matchlayout_A, matchlayoutB, outlayout_A, outlayout_B, outkeylength must be int32 row vectors (outkeylength is just a number)");  
    }
    
    /* read and digest inputs*/
    numblocks_A = mxGetN(prhs[0]);
    numblocks_B = mxGetN(prhs[1]);
    m_matchcomb = numblocks_A < numblocks_B ? numblocks_A : numblocks_B;
    if(m_matchcomb == 0){  /*one of the tensors is empty, we return empty lists*/
        plhs[0] = mxCreateCellMatrix(1,0);
        plhs[1] = mxCreateDoubleMatrix(0,3,mxREAL);
        return;
    }
    
    matchlayout_A = (int32_t *) mxGetData(prhs[2]);
    matchlayout_B = (int32_t *) mxGetData(prhs[3]);
    matchkeylength = mxGetN(prhs[2]);
    if(mxGetN(prhs[3]) != matchkeylength) mexErrMsgIdAndTxt("MyToolbox:NTdot_createtasklist:nrhs","matchlayoutA and matchlayoutB must have the same length.");
    matchkey_tmp = (uint16_t *)malloc(matchkeylength *sizeof(uint16_t));
    outlayout_A = (int32_t *) mxGetData(prhs[4]);
    outlayout_B = (int32_t *) mxGetData(prhs[5]);
    outkeylength = *((int32_t *) mxGetData(prhs[6]));
    outkey_tmp =  (uint16_t *)malloc(outkeylength *sizeof(uint16_t));
    if(m_matchcomb < 500000) initbucketsize=1000000/m_matchcomb;          /*total init volume of buckets: 10^6 elements*/
    if(m_matchcomb < 5000000) initABlistlength = 10000000/m_matchcomb;    /*init ABlist_length*/
    matchcomb_table = create_hash_table(m_matchcomb,initbucketsize);   
    keylength_A = mxGetN(mxGetCell(prhs[0],0));
    keylength_B = mxGetN(mxGetCell(prhs[1],0));
    
    /* we create the matchcomb_table */
    if(numblocks_A < numblocks_B){
        for(bID_A=0; bID_A<numblocks_A; ++bID_A){
            keyA = (uint16_t *) mxGetData(mxGetCell(prhs[0],bID_A));
            for(j=0; j<matchkeylength; ++j) matchkey_tmp[j] = keyA[matchlayout_A[j]-1];
            element_tmp = getcreate_element_ptr(matchcomb_table, matchkey_tmp, matchkeylength, m_matchcomb, p);
            if(element_tmp->value == NULL) allocate_ABelement(element_tmp, initABlistlength);
            add_bID_A_tolist((value_AB*)(element_tmp->value), bID_A);
        }
        for(bID_B=0; bID_B<numblocks_B; ++bID_B){
            keyB = (uint16_t *) mxGetData(mxGetCell(prhs[1],bID_B));
            for(j=0; j<matchkeylength; ++j) matchkey_tmp[j] = keyB[matchlayout_B[j]-1];
            element_tmp = justget_element_ptr(matchcomb_table, matchkey_tmp, matchkeylength, m_matchcomb, p);
            if(element_tmp != NULL){ 
                numtasks += ((value_AB *)(element_tmp->value))->fill_A;
                add_bID_B_tolist((value_AB*)(element_tmp->value), bID_B);
            }
        }
    }
    else{
        for(bID_B=0; bID_B<numblocks_B; ++bID_B){
            keyB = (uint16_t *) mxGetData(mxGetCell(prhs[1],bID_B));
            for(j=0; j<matchkeylength; ++j) matchkey_tmp[j] = keyB[matchlayout_B[j]-1];
            element_tmp = getcreate_element_ptr(matchcomb_table, matchkey_tmp, matchkeylength, m_matchcomb, p);
            if(element_tmp->value == NULL) allocate_ABelement(element_tmp, initABlistlength);
            add_bID_B_tolist((value_AB*)(element_tmp->value), bID_B);
        }
        for(bID_A=0; bID_A<numblocks_A; ++bID_A){
            keyA = (uint16_t *) mxGetData(mxGetCell(prhs[0],bID_A));
            for(j=0; j<matchkeylength; ++j) matchkey_tmp[j] = keyA[matchlayout_A[j]-1];
            element_tmp = justget_element_ptr(matchcomb_table, matchkey_tmp, matchkeylength, m_matchcomb, p);
            if(element_tmp != NULL){ 
                numtasks += ((value_AB *)(element_tmp->value))->fill_B;
                add_bID_A_tolist((value_AB*)(element_tmp->value), bID_A);
            }
        }
    }
    
    /* we create the task table */
    m_outkeys = numtasks;    /*worst case scenario: all tasks write to different blocks*/
    
    if(numtasks == 0){       /* if no tasks, we finalize and return */
        plhs[0] = mxCreateCellMatrix(1,0);
        plhs[1] = mxCreateDoubleMatrix(0,3,mxREAL);
        free(matchkey_tmp);
        free(outkey_tmp);
        cleanup_AB(matchcomb_table,m_matchcomb);
        clear_hashtable(matchcomb_table,m_matchcomb);
        return;
    }

    if(m_outkeys < 5000000) initbucketsize=10000000/m_outkeys;    /*max 40 MB temporary memory use*/
    outbIDs = (size_t*)malloc(m_outkeys*sizeof(size_t));
    outkey_table = create_hash_table(m_outkeys, initbucketsize);
    plhs[1] = mxCreateDoubleMatrix(numtasks,3,mxREAL);



    for(bucketID=0; bucketID<m_matchcomb; ++bucketID){
        for(elementID=0; elementID<matchcomb_table[bucketID].n_elements; ++elementID){
            value_tmp = *((value_AB *)(matchcomb_table[bucketID].elements[elementID].value));
            for(i=0; i < value_tmp.fill_A; ++i){
                bID_A = value_tmp.bIDs_A[i];
                for(j=0; j<value_tmp.fill_B; ++j){

                    bID_B = value_tmp.bIDs_B[j];
                    create_outkey(outkey_tmp, 
                                  (uint16_t *) mxGetData(mxGetCell(prhs[0],bID_A)),   /*keyA*/
                                  (uint16_t *) mxGetData(mxGetCell(prhs[1],bID_B)),   /*keyB*/
                                  keylength_A,
                                  keylength_B,
                                  outlayout_A,
                                  outlayout_B);
                    mxGetPr(plhs[1])[taskID] = (double)bID_A + 1;
                    mxGetPr(plhs[1])[taskID + numtasks] = (double)bID_B + 1;
                    element_tmp = justget_element_ptr(outkey_table, outkey_tmp, outkeylength, m_outkeys, p);
                    if(element_tmp != NULL) mxGetPr(plhs[1])[taskID + 2*numtasks] = (double)(*((size_t*)(element_tmp->value)));
                    else{
                        element_tmp = getcreate_element_ptr(outkey_table, outkey_tmp, outkeylength, m_outkeys, p);
                        numoutblocks += 1;
                        outbIDs[numoutblocks - 1] = numoutblocks;
                        element_tmp->value = (size_t*) (&(outbIDs[numoutblocks - 1]));
                        mxGetPr(plhs[1])[taskID + 2*numtasks] = (double)(*((size_t*)(element_tmp->value)));
                    }
                    taskID += 1;    
                }
            }
        }
    }
        

    
    /*we create the list of outkeys*/
    plhs[0] = mxCreateCellMatrix(1,numoutblocks);
    outdims[0] = 1; outdims[1] = outkeylength;
    
    

    for (bucketID = 0; bucketID < m_outkeys; ++bucketID){
        for(elementID = 0; elementID < outkey_table[bucketID].n_elements; ++elementID){
            tmp_Ckey = mxCreateCharArray(2, outdims);
            tmp_Ckey_data = (uint16_t *) mxGetData(tmp_Ckey);
            i = *((size_t *) (outkey_table[bucketID].elements[elementID].value)) - 1;
            for(j = 0; j < outkeylength; ++j) tmp_Ckey_data[j] = outkey_table[bucketID].elements[elementID].key[j];
            mxSetCell(plhs[0],i,tmp_Ckey);
        }
    }
    free(matchkey_tmp);
    free(outkey_tmp);
    free(outbIDs);
    cleanup_AB(matchcomb_table,m_matchcomb);
    clear_hashtable(matchcomb_table,m_matchcomb);
    clear_hashtable(outkey_table,m_outkeys);
    
    return;
}
    