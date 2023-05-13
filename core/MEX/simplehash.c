#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "simplehash.h"


size_t hash_function(const uint16_t*  key,
                     size_t           keylength,
                     size_t           num_of_buckets,
                     size_t           p)
{                     
    size_t out = 0;
    size_t i,j;
    size_t ptmp = 1;
    for(i=0;i<keylength; ++i){
        out += (ptmp * key[i]) % num_of_buckets;
        ptmp = (ptmp * p) % num_of_buckets;   /* modulo keeps it small */
    }
    return out % num_of_buckets ;
}



int key_cmp(const uint16_t *key1, 
            const uint16_t *key2, 
            size_t keylength)
{
   size_t i;
   int out = 1;
   for(i=0; i<keylength; ++i){
       if(key1[i] != key2[i]){
           out = 0;
           break;
       }
   }
   return out;
}



hash_bucket*  create_hash_table(size_t num_of_buckets, size_t init_bucket_size)
{
    size_t i;
    hash_bucket *out;
    out = (hash_bucket *) malloc (num_of_buckets*sizeof(hash_bucket));
    for(i=0; i<num_of_buckets; ++i){
        out[i].n_elements = 0;
        out[i].bucket_size = init_bucket_size;
        out[i].elements = (hash_element*) malloc(init_bucket_size * sizeof(hash_element));
    }
    return out;
}


hash_element* justget_element_ptr(hash_bucket*     hash_table,
                                  const uint16_t*  key, 
                                  size_t           keylength,
                                  size_t           num_of_buckets,
                                  size_t           p){
    size_t bucketID;
    size_t i;
    bucketID = hash_function(key, keylength, num_of_buckets, p);
    for(i = 0; i < hash_table[bucketID].n_elements; ++i){
           if(key_cmp(hash_table[bucketID].elements[i].key, key, keylength)) return &(hash_table[bucketID].elements[i]);   /* we found the key in the bucket */
    }
    return NULL;
}


hash_element* getcreate_element_ptr(hash_bucket*     hash_table,
                                    const uint16_t*  key, 
                                    size_t           keylength,
                                    size_t           num_of_buckets,
                                    size_t           p){
    size_t bucketID;
    size_t i,j;
    hash_element *newelementlist;
    bucketID = hash_function(key, keylength, num_of_buckets, p);
    for(i = 0; i < hash_table[bucketID].n_elements; ++i){
           if(key_cmp(hash_table[bucketID].elements[i].key, key, keylength)) return &(hash_table[bucketID].elements[i]);   /* we found the key in the bucket */
    }
    /*if we reach this point, we have not found the key in the bucket (the bucket is probably empty or there is a hash collision)*/
    /*If n_sub > 0 :: hash collision, otherwise the bucket was empty (we need to allocate storages) */
    if(hash_table[bucketID].n_elements == hash_table[bucketID].bucket_size){  /* the element list is full: we need to allocate more memory */
        newelementlist = (hash_element*) malloc(2*hash_table[bucketID].bucket_size*sizeof(hash_element));
        for(i=0; i<hash_table[bucketID].n_elements; ++i){
            newelementlist[i] = hash_table[bucketID].elements[i];   /*this is shallow copy (we don't reallocate the value pointers*/
        }
        hash_table[bucketID].bucket_size *= 2;    /* we double the bucket size */
        free(hash_table[bucketID].elements);
        hash_table[bucketID].elements = newelementlist;
    }
    i = hash_table[bucketID].n_elements;
    hash_table[bucketID].n_elements += 1;
    hash_table[bucketID].elements[i].key = (uint16_t *) malloc(keylength * sizeof(uint16_t));
    for(j=0; j<keylength; ++j) hash_table[bucketID].elements[i].key[j] = key[j];
    hash_table[bucketID].elements[i].value = NULL;       /* this is a void pointer, the user must allocate an appropriate variable here. */
    return &(hash_table[bucketID].elements[i]);
}


void clear_hashtable(hash_bucket*     hash_table,
                     size_t           num_of_buckets){                       /*this function deallocates the pointers that were allocated by the general functions. 
                                                                              * hash_table[bucketID].elements[i].value pointers are supposed to be deallocated earlier.*/
    size_t bucketID;
    size_t i;
    for(bucketID = 0; bucketID < num_of_buckets; ++bucketID){
        for(i = 0; i < hash_table[bucketID].n_elements; ++i) free(hash_table[bucketID].elements[i].key);
        free(hash_table[bucketID].elements);
    }
    free(hash_table);
}


