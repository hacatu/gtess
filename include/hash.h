#ifndef __HASH_H__
#define __HASH_H__

#include <stddef.h>

typedef struct{
	void *table_a, *table_b;
	uint64_t *flags_a, *flags_b;
	size_t cap, len_a, len_b, i, r, full;
} hashtbl_t;

typedef struct{
	size_t size;
	uint64_t (*hash)(const void*);
	int (*cmp)(const void*, const void*);
	int (*add)(void*, void*);
	void (*del)(void*);
	double load_factor;
} hashtbl_ft;

int hash_init(hashtbl_t*, const hashtbl_ft*, size_t reserve);

void *hash_get(hashtbl_t*, const hashtbl_ft*, const void *key);

void *hash_insert(hashtbl_t*, const hashtbl_ft*, const void *key, int *status);

void *hash_append(hashtbl_t *self, const hashtbl_ft *ft, void *key, int *status);

int hash_remove(hashtbl_t*, const hashtbl_ft*, const void *key);

void hash_delete(hashtbl_t*, const hashtbl_ft*, void *ent);

void hash_clear(hashtbl_t*, const hashtbl_ft*);

void hash_destroy(hashtbl_t*, const hashtbl_ft*);

inline static void *hash_next(hashtbl_t *self, const hashtbl_ft *ft, void *cur){
	uint64_t a = 0, b = 0;
	if(self->table_b){
		if(cur){
			b = (cur - self->table_b)/ft->size;
			if(b++ >= self->len_b){
				goto SEARCH_A;
			}
		}
		for(; b < self->len_b; ++b){
			if((self->flags_b[b >> 5] >> (b&0x1F))&0x100000000ULL){
				return self->table_b + b*ft->size;
			}
		}
	}else
	SEARCH_A: if(cur){
		a = (cur - self->table_a)/ft->size;
		if(a++ >= self->len_a){
			return NULL;
		}
	}
	for(; a < self->len_a; ++a){
		if((self->flags_a[a >> 5] >> (a&0x1F))&0x100000000ULL){
			return self->table_a + a*ft->size;
		}
	}
	return NULL;
}

#endif

