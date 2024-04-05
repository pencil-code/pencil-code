#pragma once
#include "headers_c.h"
#include <assert.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>

pthread_mutex_t diag_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t diag_cond = PTHREAD_COND_INITIALIZER;
const int DIAG_COND_HANDLE = 1;
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

void
FTNIZE(cond_signal)(const int* cond_handle)
{
   switch(*cond_handle)
   {
	   case DIAG_COND_HANDLE:
		   pthread_cond_signal(&diag_cond);
		   pthread_mutex_unlock(&diag_mutex);
		   return;
	   default:
		   printf("Incorrect cond var handle\n");
		   assert(false); //Incorrect cond var handle
   }
}

