#include "headers_c.h"
#include <assert.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <errno.h>

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
const int DIAG_COND_HANDLE = 1;
/* ---------------------------------------------------------------------------- */
void
FTNIZE(cond_init)()
{
		pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_ERRORCHECK);
		int res = pthread_mutex_init(&diag_mutex, &mutex_attr);
		if(res != 0)
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
void
FTNIZE(cond_wait_single)(const int* cond_handle, volatile bool* flag, volatile bool* val)
{
   switch(*cond_handle)
   {
	case DIAG_COND_HANDLE:
          assert(1==1);
          const int res = pthread_mutex_lock(&diag_mutex);
	  if(res != 0)
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
	  while(*flag != *val) pthread_cond_wait(&diag_cond, &diag_mutex);
	  return;
	default:
	  printf("Error - cond_wait: Incorrect cond_handle!!!\n");
	  assert(false); //Incorrect cond var handle
   }
}
void
FTNIZE(cond_wait_multi)(const int* cond_handle, volatile bool* flag, volatile bool* val, const int* n)
{
   switch(*cond_handle)
   {
	case DIAG_COND_HANDLE:
		//for some stupid syntax reason that I can't figure out now you need this extra assert
	  assert(1==1);
          const int res = pthread_mutex_lock(&diag_mutex);
	  if(res != 0)
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
	  for(int i = 0; i < *n; ++i)
		  condition &= flag[i] == val[i];

	  while(!condition) 
	  {
		pthread_cond_wait(&diag_cond, &diag_mutex);
	        for(int i = 0; i < *n; ++i)
		  condition &= flag[i] == val[i];
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
	  assert(1==1);
	  pthread_cond_signal(&diag_cond);
	  //pthread_cond_broadcast(&diag_cond);
	  const int res = pthread_mutex_unlock(&diag_mutex);
	  if(res != 0)
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
	  return;
	default:
	  printf("Error - cond_signal: Incorrect cond_handle!!!\n");
	  assert(false); //Incorrect cond var handle
   }
}
/* ---------------------------------------------------------------------------- */
