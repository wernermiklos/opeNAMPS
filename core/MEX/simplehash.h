#ifndef SIMPLEHASH_H
#define SIMPLEHASH_H

#include <stdint.h>


typedef struct HashElement {
    uint16_t*  key;
    void*      value;                    /* this will be recasted to the desired type */
} hash_element;


typedef struct HashBucket {
    size_t       n_elements;
    size_t       bucket_size;
    hash_element *elements;
} hash_bucket;


size_t hash_function(const uint16_t*  key,
                     size_t           keylength,
                     size_t           num_of_buckets,
                     size_t           p);

int key_cmp(const uint16_t *key1, const uint16_t *key2, size_t keylength);

hash_bucket*  create_hash_table(size_t    num_of_buckets,
                                size_t    init_bucket_size);

hash_element* getcreate_element_ptr(hash_bucket*     hash_table,
                                    const uint16_t*  key, 
                                    size_t           keylength,
                                    size_t           num_of_buckets,
                                    size_t           p);

hash_element* justget_element_ptr(hash_bucket*     hash_table,
                                  const uint16_t*  key, 
                                  size_t           keylength,
                                  size_t           num_of_buckets,
                                  size_t           p);

void clear_hashtable(hash_bucket*     hash_table,
                     size_t           num_of_buckets);




#endif