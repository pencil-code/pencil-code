#define _GNU_SOURCE
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <errno.h>
#include <sched.h>
#include <pthread.h>

#include "headers_c.h"

//PTHREAD MUTEX ERROR MESSAGE MEANINGS
//	  [EINVAL]
//The mutex was created with the protocol attribute having the value PTHREAD_PRIO_PROTECT and the calling thread's priority is higher than the mutex's current priority ceiling.
//The pthread_mutex_trylock() function shall fail if:
//
//[EBUSY]
//The mutex could not be acquired because it was already locked.
//The pthread_mutex_lock(), pthread_mutex_trylock(), and pthread_mutex_unlock() functions may fail if:
//
//[EINVAL]
//The value specified by mutex does not refer to an initialized mutex object.
//[EAGAIN]
//[XSI] [Option Start] The mutex could not be acquired because the maximum number of recursive locks for mutex has been exceeded. [Option End]
//The pthread_mutex_lock() function may fail if:
//
//[EDEADLK]
//A deadlock condition was detected or the current thread already owns the mutex.
//The pthread_mutex_unlock() function may fail if:
//
//[EPERM]
//The current thread does not own the mutex.
pthread_mutex_t diag_mutex;
pthread_mutexattr_t mutex_attr;
//pthread_mutex_t diag_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t diag_cond = PTHREAD_COND_INITIALIZER;
#define DIAG_COND_HANDLE (1)
/* ------------------------------------------------------------------------------------ */
//extern "C" 
//{
void
FTNIZE(cond_init)()
{
	pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_ERRORCHECK);
	int res = pthread_mutex_init(&diag_mutex, &mutex_attr);
	if (res != 0)
	{
	  printf("Tried to init mutex: %d\n",res);
  	  printf("EBUSY: %d\n",EBUSY);
  	  printf("EAGAIN: %d\n",EAGAIN);
  	  printf("EPERM: %d\n",EPERM);
  	  printf("ENOMEM: %d\n",ENOMEM);
  	  printf("EINVAL: %d\n",EINVAL);
	  fflush(stdout);
	}
	assert(res == 0);
}
/* ------------------------------------------------------------------------------------ */
void
FTNIZE(cond_wait_single)(const int* cond_handle, volatile bool* flag, volatile bool* val)
{
   switch(*cond_handle)
   {
	case DIAG_COND_HANDLE:
	{
          const int res = pthread_mutex_lock(&diag_mutex);
	  if (res != 0)
	  {
	    printf("lock error message: %d\n",res);
	    printf("EBUSY: %d\n",EBUSY);   // only relevant for trylock
	    printf("EINVAL: %d\n",EINVAL);
	    printf("EAGAIN: %d\n",EAGAIN);
	    printf("EDEADLK: %d\n",EDEADLK);
	    printf("EPERM: %d\n",EPERM);   // only relevant for unlock
	    fflush(stdout);
	    assert(res == 0);
	  }
	  while (*flag != *val) pthread_cond_wait(&diag_cond, &diag_mutex);
	}
	return;
	default:
	  printf("Error - cond_wait: Incorrect cond_handle!!!\n");
	  assert(false); //Incorrect cond var handle
   }
}
/* ------------------------------------------------------------------------------------ */
void
FTNIZE(cond_wait_multi)(const int* cond_handle, volatile bool* flag, volatile bool* val, const int* n)
{
   switch(*cond_handle)
   {
	case DIAG_COND_HANDLE:
		{
          const int res = pthread_mutex_lock(&diag_mutex);
	  if (res != 0)
	  {
	    printf("lock error message: %d\n",res);
	    printf("EBUSY: %d\n",EBUSY);
	    printf("EINVAL: %d\n",EINVAL);
	    printf("EAGAIN: %d\n",EAGAIN);
	    printf("EDEADLK: %d\n",EDEADLK);
	    printf("EPERM: %d\n",EPERM);
	    fflush(stdout);
	    assert(res == 0);
	  }
	  bool condition = true;
	  for (int i = 0; i < *n; ++i)
	    condition &= flag[i] == val[i];

	  while (!condition) 
	  {
		pthread_cond_wait(&diag_cond, &diag_mutex);
	        for (int i = 0; i < *n; ++i)
		  condition &= flag[i] == val[i];
	  }
		}
	  return;
	default:
	  printf("Error - cond_wait: Incorrect cond_handle!!!\n");
	  assert(false); //Incorrect cond var handle
   }
}
/* ---------------------------------------------------------------------------- */
void
FTNIZE(cond_signal)(const int* cond_handle)
{
   switch(*cond_handle)
   {
	case DIAG_COND_HANDLE:
	{
	  pthread_cond_signal(&diag_cond);
	  //pthread_cond_broadcast(&diag_cond);
	  const int res = pthread_mutex_unlock(&diag_mutex);
	  if (res != 0)
	  {
	    printf("unlock error message: %d\n",res);
	    printf("EBUSY: %d\n",EBUSY);
	    printf("EINVAL: %d\n",EINVAL);
	    printf("EAGAIN: %d\n",EAGAIN);
	    printf("EDEADLK: %d\n",EDEADLK);
	    printf("EPERM: %d\n",EPERM);
	    fflush(stdout);
	    assert(res == 0);
	  }
	}
	return;
	default:
	  printf("Error - cond_signal: Incorrect cond_handle!!!\n");
	  assert(false); //Incorrect cond var handle
   }
}
int
FTNIZE(get_cpu_c)()
{
	return sched_getcpu();
}
bool
FTNIZE(set_cpu_c)(int* core_id_in)
{
    const int core_id = *core_id_in;
    // Credit: https://eli.thegreenplace.net/2016/c11-threads-affinity-and-hyperthreading/
    // Create a cpu_set_t object representing a set of CPUs. Clear it and mark
    // only CPU core_id as set.
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);
    const int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
    return (bool)(1-rc);
}

void
FTNIZE(cond_wait)(const int* cond_handle, volatile bool* flag, volatile bool* val)
{
   switch(*cond_handle)
   {
   case DIAG_COND_HANDLE:
           pthread_mutex_lock(&diag_mutex);
	   while(*flag != *val)
		   pthread_cond_wait(&diag_cond, &diag_mutex);
	   return;
   default:
	   printf("Incorrect cond var handle\n");
	   assert(false); //Incorrect cond var handle
   }
}
//}

/* ---------------------------------------------------------------------------- */
