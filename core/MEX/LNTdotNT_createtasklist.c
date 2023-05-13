/*==========================================================
 * NTdot_createtasklist.c - example in MATLAB External Interfaces
 *
 * Creates task list for NTdot.
 * 
 *
 * The calling syntax is:
 *
 *		[keylist_C, tasklist] = LNTdotNT_createtasklist(keysA_list, keys_B, matchlayout_A_list, matchlayout_B_list, outlayout_A_list, outlayout_B, outkeylength, no_of_symmetries)
 *
 * This is a MEX-file for MATLAB.
 * Written by M. A. Werner 2023 *
 *========================================================*/

#define VALUE_AB_INITSIZE 10


#include "mex.h"
#include "matrix.h"
#include "simplehash.h"
#include <stdio.h>
/* The computational routine */



/* The content is problem specific */
typedef struct MatchComb_value {
    size_t length;
    size_t fill;
    size_t *bIDs;
} matchcomb_value;


void allocate_matchcomb_element(hash_element* element, size_t initlength){
    element->value = (void *)malloc(sizeof(matchcomb_value));
    ((matchcomb_value*)(element->value))->length = initlength;
    ((matchcomb_value*)(element->value))->fill = 0;
    ((matchcomb_value*)(element->value))->bIDs = (size_t *)malloc(initlength * sizeof(size_t));
}

void free_matchcomb_element(hash_element* element){
    free(((matchcomb_value*)(element->value))->bIDs);
    free(element->value);
}

void add_bID_tolist(matchcomb_value* value, size_t bID){
    size_t *newbIDs;
    size_t i;
    if(value->fill == value->length){    /* there is no place for the new bID -- we allocate more memory */
        newbIDs = (size_t *)malloc(2 * value->length * sizeof(size_t));
        for(i=0; i<value->fill; ++i) newbIDs[i] = value->bIDs[i];
        free(value->bIDs);
        value->bIDs = newbIDs;
        value->length *= 2;
    }
    value->bIDs[value->fill] = bID;
    value->fill += 1;
}

void cleanup_matchcomb(hash_bucket *hash_table, size_t num_of_buckets){
    size_t bucketID,j;
    for(bucketID = 0; bucketID < num_of_buckets; ++bucketID){
        for(j=0;j<hash_table[bucketID].n_elements;++j) free_matchcomb_element(&(hash_table[bucketID].elements[j]));
    }
    return;
}

void create_outkey(uint16_t *outkey_tmp, 
                   uint16_t **keyA_list, 
                   uint16_t *keyB, 
                   size_t keylength_A, 
                   size_t keylength_B, 
                   int32_t **outlayout_A, 
                   int32_t *outlayout_B,
                   size_t no_of_symmetries){
    size_t i, symID;
    for(symID = 0; symID < no_of_symmetries; ++symID) for(i=0; i<keylength_A; ++i) if(outlayout_A[symID][i] > 0) outkey_tmp[outlayout_A[symID][i]-1] = keyA_list[symID][i];
    for(i=0; i<keylength_B; ++i) if(outlayout_B[i] > 0) outkey_tmp[outlayout_B[i]-1] = keyB[i];
    return;
}

/*------------------------------------------------------------------*/        


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){                             
    size_t no_of_symmetries, outkeylength, matchkeylength;
    size_t *numblocks_A, numblocks_B, numblocks_C;
    size_t symID;
    int32_t **matchlayout_A, **matchlayout_B;
    int32_t **outlayout_A;
    int32_t *outlayout_B;
    size_t keylength_A, keylength_B;
    size_t int_tmp;
    hash_bucket **matchcomb_tableA;
    size_t initbucketsizeA, initbucketsize, *init_bIDlistlength;
    size_t bID_A, bID_B, bucketID, elementID, taskID=0, task_start=0, task_end;
    uint16_t *keyA, *keyB;
    size_t i,j;
    hash_element *element_tmp;
    size_t p=53;            /*parameter of the hash function. p=53 is a reasonable prime number */
    size_t numtasks, numtasks_tmp, numoutblocks, aux_ind;
    size_t *blockB_used;
    size_t **bIDs_A;
    size_t *fills_A;
    uint16_t **keyA_list;
    size_t m_outkeys;       /*parameter of the hash function. (outkey_table)*/
    size_t *outbIDs;
    hash_bucket *outkey_table;
    matchcomb_value value_tmp;
    uint16_t *outkey_tmp;
    uint16_t *matchkey_tmp;
    mxArray *tmp_Ckey;
    uint16_t *tmp_Ckey_data;
    size_t outdims[2];
    
    /* Error check */
    if(nrhs!=8) {
        mexErrMsgIdAndTxt("MyToolbox:NTdot_createtasklist:nrhs","Six inputs required:  keysA_list, keysB, matchlayout_A_list, matchlayoutB_list, outlayout_A_list, outlayout_B, outkeylength, no_of_symmetries");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:NTdot_createtasklist:nlhs","Two output required.");
    }
    
    if((mxGetM(prhs[2]) > 1) ||
       (mxGetM(prhs[3]) > 1) ||
       (mxGetM(prhs[4]) > 1) ||
       (mxGetM(prhs[5]) > 1) ||
       (mxGetM(prhs[6]) > 1) ||
       (mxGetM(prhs[7]) > 1)){
        mexErrMsgIdAndTxt("MyToolbox:LNTdotNT_createtasklist:nrhs","matchlayout_A_list, matchlayoutB_list, outlayout_A_list, outlayout_B, outkeylength, no_of_symmetries must be row vectors (outkeylength & no_of_symmetries is just a number)");  
    }

    /* read and digest inputs*/
    no_of_symmetries = *((int32_t *) mxGetData(prhs[7]));
    numblocks_A = (size_t *) malloc(no_of_symmetries*sizeof(size_t));
    int_tmp = 1;
    for(symID = 0; symID < no_of_symmetries; ++symID){
        numblocks_A[symID] = mxGetN(mxGetCell(prhs[0],symID));
        if(numblocks_A[symID] == 0) int_tmp = 0;
    }    
    numblocks_B = mxGetN(prhs[1]);
    blockB_used = (size_t*) malloc(numblocks_B*sizeof(size_t));

    if(int_tmp == 0 || numblocks_B == 0){  /*one of the (sub)tensors is empty, we return empty lists*/
        free(numblocks_A);
        plhs[0] = mxCreateCellMatrix(1,0);
        plhs[1] = mxCreateDoubleMatrix(0,no_of_symmetries + 2,mxREAL);
        return;
    }

    
    matchlayout_A = (int32_t **) malloc(no_of_symmetries*sizeof(int32_t*));
    for(symID = 0; symID < no_of_symmetries; ++symID) matchlayout_A[symID] = (int32_t *) mxGetData(mxGetCell(prhs[2],symID));
    matchlayout_B = (int32_t **) malloc(no_of_symmetries*sizeof(int32_t*));
    for(symID = 0; symID < no_of_symmetries; ++symID) matchlayout_B[symID] = (int32_t *) mxGetData(mxGetCell(prhs[3],symID));
    matchkeylength = mxGetN(mxGetCell(prhs[2],0));
    

    for(symID = 0; symID < no_of_symmetries; ++symID){
        if(mxGetN(mxGetCell(prhs[2],symID) )!= matchkeylength) 
            mexErrMsgIdAndTxt("MyToolbox:LNTdotNT_createtasklist:nrhs","matchlayout_A and matchlayout_B lists element must have the same length.");
        if(mxGetN(mxGetCell(prhs[2],symID)) != matchkeylength) 
            mexErrMsgIdAndTxt("MyToolbox:LNTdotNT_createtasklist:nrhs","matchlayout_A and matchlayout_B lists element must have the same length.");
    }


    matchkey_tmp = (uint16_t *)malloc(matchkeylength *sizeof(uint16_t));

    outlayout_A = (int32_t **) malloc(no_of_symmetries*sizeof(int32_t*));
    for(symID = 0; symID < no_of_symmetries; ++symID) outlayout_A[symID] = (int32_t *) mxGetData(mxGetCell(prhs[4],symID));
    outlayout_B = (int32_t *) mxGetData(prhs[5]);
    outkeylength = *((int32_t *) mxGetData(prhs[6]));
    outkey_tmp =  (uint16_t *)malloc(outkeylength *sizeof(uint16_t));

    init_bIDlistlength = (size_t *) malloc(no_of_symmetries*sizeof(size_t));

    matchcomb_tableA = (hash_bucket **) malloc(no_of_symmetries*sizeof(hash_bucket*));
    for(symID = 0; symID < no_of_symmetries; ++symID){
        initbucketsizeA = 1;
        init_bIDlistlength[symID] = 1;
        if(numblocks_A[symID] < 500000) initbucketsizeA=1000000/numblocks_A[symID];          /*total init volume of buckets: 10^6 elements*/
        if(numblocks_A[symID] < 5000000) init_bIDlistlength[symID] = 10000000/numblocks_A[symID];    /*init bIDlist_length*/
        matchcomb_tableA[symID] = create_hash_table(numblocks_A[symID],initbucketsizeA);
    }   
    keylength_A = mxGetN(mxGetCell(mxGetCell(prhs[0],0),0));
    keylength_B = mxGetN(mxGetCell(prhs[1],0));
    
    /* we create the matchcomb_table */
    for(symID = 0; symID < no_of_symmetries; ++symID){
        for(bID_A=0; bID_A<numblocks_A[symID]; ++bID_A){
            keyA = (uint16_t *) mxGetData(mxGetCell(mxGetCell(prhs[0],symID),bID_A));
            for(j=0; j<matchkeylength; ++j) matchkey_tmp[j] = keyA[matchlayout_A[symID][j]-1];
            element_tmp = getcreate_element_ptr(matchcomb_tableA[symID], matchkey_tmp, matchkeylength, numblocks_A[symID], p);
            if(element_tmp->value == NULL) allocate_matchcomb_element(element_tmp, init_bIDlistlength[symID]);
            add_bID_tolist((matchcomb_value*)(element_tmp->value), bID_A);   
        }    
    }
    numtasks = 0;
    for(bID_B=0; bID_B<numblocks_B; ++bID_B){
        blockB_used[bID_B] = 1;
        numtasks_tmp = 1;
        keyB = (uint16_t *) mxGetData(mxGetCell(prhs[1],bID_B));
        for(symID = 0; symID < no_of_symmetries; ++symID){
            for(j=0; j<matchkeylength; ++j) matchkey_tmp[j] = keyB[matchlayout_B[symID][j]-1];
            element_tmp = justget_element_ptr(matchcomb_tableA[symID], matchkey_tmp, matchkeylength, numblocks_A[symID], p);
            if(element_tmp != NULL){ 
                numtasks_tmp *= ((matchcomb_value *)(element_tmp->value))->fill;
            }
            else{
                blockB_used[bID_B] = 0;
                numtasks_tmp = 0;
                break;
            }
        }
        numtasks += numtasks_tmp;    
    } 
    

    /* we create the task table */
    m_outkeys = numtasks;    /*worst case scenario: all tasks write to different blocks*/
    
    if(numtasks == 0){       /* if no tasks, we finalize and return */
        plhs[0] = mxCreateCellMatrix(1,0);
        plhs[1] = mxCreateDoubleMatrix(0,3,mxREAL);
        for(symID = 0; symID < no_of_symmetries; ++symID){
            cleanup_matchcomb(matchcomb_tableA[symID],numblocks_A[symID]);
            clear_hashtable(matchcomb_tableA[symID],numblocks_A[symID]);
        }
        free(matchcomb_tableA);
        free(matchkey_tmp);
        free(outkey_tmp);
        free(matchlayout_A);
        free(matchlayout_B);
        free(outlayout_A);
        free(numblocks_A);
        free(blockB_used);
        free(init_bIDlistlength);
        return;
    }

    if(m_outkeys < 5000000) initbucketsize=10000000/m_outkeys;    /*max 40 MB temporary memory use*/
    outbIDs = (size_t*)malloc(m_outkeys*sizeof(size_t));
    outkey_table = create_hash_table(m_outkeys, initbucketsize);
    plhs[1] = mxCreateDoubleMatrix(numtasks,2 + no_of_symmetries,mxREAL);


    bIDs_A = (size_t **) malloc(no_of_symmetries*sizeof(size_t*));
    keyA_list = (uint16_t **) malloc(no_of_symmetries*sizeof(uint16_t*));
    fills_A = (size_t *)malloc(no_of_symmetries*sizeof(size_t));
    numoutblocks=0;

    /* filling of the task-list*/
    for(bID_B=0; bID_B<numblocks_B; ++bID_B){
        if(!blockB_used[bID_B]) continue; /* If no matching blocks in A tensor, then next bID_B*/

        keyB = (uint16_t *) mxGetData(mxGetCell(prhs[1],bID_B));
        numtasks_tmp = 1;
        for(symID = 0; symID < no_of_symmetries; ++symID){
            for(j=0; j<matchkeylength; ++j) matchkey_tmp[j] = keyB[matchlayout_B[symID][j]-1];
            element_tmp = justget_element_ptr(matchcomb_tableA[symID], matchkey_tmp, matchkeylength, numblocks_A[symID], p);
            bIDs_A[symID] = ((matchcomb_value *) element_tmp->value)->bIDs;
            fills_A[symID] = ((matchcomb_value *) element_tmp->value)->fill;
            /*keyA_list[symID] = element_tmp->key;*/
            numtasks_tmp *= fills_A[symID];
        }
        task_end = task_start + numtasks_tmp;
        for(taskID = task_start;taskID < task_end; ++taskID){
            aux_ind = taskID-task_start;
            for(symID=0; symID < no_of_symmetries; ++symID){
                bID_A = bIDs_A[symID][aux_ind % fills_A[symID]];
                keyA_list[symID] = (uint16_t *) mxGetData(mxGetCell(mxGetCell(prhs[0],symID),bID_A));
                mxGetPr(plhs[1])[taskID + symID*numtasks] = (double)bID_A + 1;
                aux_ind = (aux_ind - aux_ind % fills_A[symID])/fills_A[symID];
            }
            create_outkey(outkey_tmp, keyA_list,keyB, keylength_A, keylength_B, outlayout_A, outlayout_B, no_of_symmetries);
            element_tmp = justget_element_ptr(outkey_table, outkey_tmp, outkeylength, m_outkeys, p);
            if(element_tmp == NULL){
                element_tmp = getcreate_element_ptr(outkey_table, outkey_tmp, outkeylength, m_outkeys, p);
                numoutblocks += 1;
                outbIDs[numoutblocks - 1] = numoutblocks;
                element_tmp->value = (size_t*) (&(outbIDs[numoutblocks - 1]));
            }
            mxGetPr(plhs[1])[taskID + no_of_symmetries*numtasks] = (double)bID_B + 1;
            mxGetPr(plhs[1])[taskID + (1+no_of_symmetries)*numtasks] = (double)(*((size_t*)(element_tmp->value)));
        }
        task_start = task_end;
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
    free(matchlayout_A);
    free(matchlayout_B);
    free(outlayout_A);
    free(outbIDs);

    free(blockB_used);
    free(init_bIDlistlength);
    free(bIDs_A);
    free(keyA_list);
    free(fills_A);
    for(symID = 0; symID < no_of_symmetries; ++symID){
        cleanup_matchcomb(matchcomb_tableA[symID],numblocks_A[symID]);
        clear_hashtable(matchcomb_tableA[symID],numblocks_A[symID]);
    }
    free(numblocks_A);
    clear_hashtable(outkey_table,m_outkeys);
    
    free(matchcomb_tableA);

    return;
}
    