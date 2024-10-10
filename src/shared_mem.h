#include <stdbool.h>          //For bool
#include <semaphore.h>

REAL*
allocate_shm(const long long *n_elems, const char* name);

REAL*
get_shm(const int n_elems, const char* name);

bool
shm_exists(const char* name);

sem_t*
create_sem(const char* name);

sem_t*
get_sem(const char* name);
