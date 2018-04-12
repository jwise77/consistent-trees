#ifndef CHECK_SYSCALLS_H
#define CHECK_SYSCALLS_H

#include <stdio.h>
FILE *check_fopen(char *filename, char *mode);
void *check_realloc(void *ptr, size_t size, char *reason);

#define stringify(x) #x
#define to_string(x) stringify(x)
#define check_realloc_s(x,y,z) { (x) = check_realloc((x),((int64_t)(y))*((int64_t)(z)), "Reallocating " #x " at " __FILE__ ":" to_string(__LINE__)); }
#define check_malloc_s(x,y,z) check_realloc(NULL,((int64_t)(y))*((int64_t)(z)), "Allocating " #x " at " __FILE__ ":" to_string(__LINE__));

#define check_realloc_var(x,size,cur,new) { if (cur < new) { cur = new; check_realloc_s(x,size,new); } }
#define check_realloc_every(x,size,cur,num) { if (!((cur)%(num))) { check_realloc_s(x,size,(cur)+(num)); } }
#define check_realloc_smart(x,size,cur,new) { if (cur < new) { cur = new*1.05 + 1000; check_realloc_s(x,size,new); } }

#endif /* CHECK_SYSCALLS_H */
