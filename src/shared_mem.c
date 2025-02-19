#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>           // For O_* constants
#include <sys/mman.h>        // For shm_open, mmap
#include <sys/stat.h>        // For mode constants
#include <unistd.h>          // For ftruncate
#include "headers_c.h"
#include "shared_mem.h"

#define EXIT_FAILURE 1

const int READ_AND_WRITE_PERMISSIONS_TO_ALL = 0666;
const char* global_shm_name = NULL;
const char* global_sem_name = NULL;

typedef struct
{
	const int fd;
	const off_t size;
} shm;

void cleanup_shared_memory() {
    if (global_shm_name != NULL) {
        if (shm_unlink(global_shm_name) == -1) {
            perror("shm_unlink failed");
        } else {
            //printf("Shared memory segment '%s' removed.\n", global_shm_name);
        }
    }
}
void cleanup_sem() {
    if (global_sem_name != NULL) {
        if (sem_unlink(global_sem_name) == -1) {
            perror("sem_unlink failed");
        } else {
            //printf("Shared memory segment '%s' removed.\n", global_shm_name);
        }
    }
}
void register_exit_cleanup(const char* name)
{
	global_shm_name = name;
	atexit(cleanup_shared_memory);
}
void register_sem_cleanup(const char* name)
{
	global_sem_name = name;
	atexit(cleanup_sem);
}
int open_shm(const char* name, const int create)
{
    const int shm_fd = shm_open(name, create | O_RDWR, READ_AND_WRITE_PERMISSIONS_TO_ALL);
    if (shm_fd == -1) {
        perror("shm_open failed");
        exit(EXIT_FAILURE);
    }
    return shm_fd;
}
shm create_shm(const char* name, const off_t size)
{
    const int shm_fd = open_shm(name,O_CREAT);
    if (ftruncate(shm_fd, size) == -1) {
        perror("ftruncate failed");
        close(shm_fd);
        exit(EXIT_FAILURE);
    }
    return (shm){shm_fd,size};
}
bool shm_exists(const char* name)
{
    const int shm_fd = shm_open(name,  O_RDWR, READ_AND_WRITE_PERMISSIONS_TO_ALL);
    if (shm_fd != -1) close(shm_fd);
    return shm_fd != -1;

}
REAL* map_shm(const shm x)
{
    REAL* addr = (REAL*)mmap(0, x.size, PROT_READ | PROT_WRITE, MAP_SHARED, x.fd, 0);
    if (addr == MAP_FAILED) {
        perror("mmap failed");
        close(x.fd);
        exit(EXIT_FAILURE);
    }
    return addr;
}
REAL* get_shm(const int n_elems, const char* name)
{
    return map_shm(
		    (shm){
		    	open_shm(name,0),
			n_elems*sizeof(REAL)
		    });
}
REAL* FTNIZE(allocate_shm)(const long long *n_elems, const char* name) {

    register_exit_cleanup(name);
    const off_t size = *n_elems*sizeof(REAL);
    const shm shm_res  = create_shm(name, size);
    REAL* res = map_shm(shm_res);
    //TP: create a semaphore with a same name intended for synchronized access to the shared mem file
    create_sem(name);
    close(shm_res.fd);
    return res;
}
sem_t* create_sem(const char* name)
{
    sem_t* sem = sem_open(name, O_CREAT, 0666, 1);
    if (sem == SEM_FAILED) {
        perror("sem_open failed");
        exit(EXIT_FAILURE);
    }
    register_sem_cleanup(name);
    return sem;
}
sem_t* get_sem(const char* name)
{
    sem_t* sem = sem_open(name, 0);
    if (sem == SEM_FAILED) {
        perror("sem_open failed");
        exit(EXIT_FAILURE);
    }
    return sem;
}
