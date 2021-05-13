/*                             syscalls_ansi.c
                              -----------------
*/

/* Date:   19-Mar-2010
   Author: Bourdin.KIS (Bourdin@KIS.Uni-Freiburg.de)
   Description:
 ANSI C and standard library callable function wrappers for use in Fortran.
 Written to compensate for inadequatenesses in the Fortran95/2003 standards.
*/
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <dlfcn.h>

#include "headers_c.h"

/* ---------------------------------------------------------------------- */
void FTNIZE(extract_string_c)(char *extract_cmd, char *result)
{
//  Extracts a string (e.g., from a file) by extract_cmd.
//  31-mar-21/MR: coded

    FILE *pipe = popen(extract_cmd, "r");
    if (pipe == NULL) {
        perror("popen");
        exit(EXIT_FAILURE);
    }
    char* ret=fgets(result, sizeof(result), pipe);
    pclose(pipe);
}
/* ---------------------------------------------------------------------- */

int FTNIZE(readlink_c) (const char *linkfile, char *link, const FINT *size)
{
//  Returns the contents of a symbolic link.
//  31-mar-21/MR: coded

  int ret;

  ret=readlink(linkfile,link,*size);
  return ret;
}
/* ---------------------------------------------------------------------- */

int FTNIZE(islink_c) (char *filename)
{
//  Tests whether filename is a symbolic link.
//  21-mar-20/MR: coded

  struct stat fileStat;
  int ret = -1;

  if (lstat(filename, &fileStat) == 0) { 
    ret = (fileStat.st_mode & S_IFMT) == S_IFLNK ? 1 : 0;   
  }
  return ret;
}
/* ---------------------------------------------------------------------- */

void FTNIZE(file_size_c)
     (char *filename, FINT *bytes)
/* Determines the size of a file.
   Returns:
   * positive integer containing the file size of a given file
   * -2 if the file could not be found or opened
   * -1 if retrieving the file size failed
*/
{
  struct stat fileStat;
  int file = -1;

  *bytes = -2;
  file = open (filename, O_RDONLY);
  if(file == -1) return;

  *bytes = -1;
  if(fstat (file, &fileStat) < 0) { close (file); return; }
  close (file);

  *bytes = fileStat.st_size;
}
/* ---------------------------------------------------------------------- */
void FTNIZE(caller2)
     (void (**func)(void*,void*), void* par1, void* par2)
{
  (*func)(par1,par2);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(caller2_str1)
     (void (**func)(char*,int*,void*), char* str, int len, void* par1)
{
  (*func)(str,&len,par1);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(caller3)
     (void (**func)(void*,void*,void*), void* par1, void* par2, void* par3)
{
  (*func)(par1,par2,par3);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(caller1_str)
     (void (**func)(char*,int*), char* str, int len)
{
  (*func)(str,&len);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(caller3_str1)
     (void (**func)(char*,int*,void*,void*), char* str, int len, void* par2, void* par3)
{
  (*func)(str,&len,par2,par3);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(caller4)
     (void (**func)(void*,void*,void*,void*), void* par1, void* par2, void* par3, void* par4)
{
  (*func)(par1,par2,par3,par4);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(caller4_str1)
     (void (**func)(char*,int*,void*,void*,void*), char* str, int len, void* par2, void* par3, void* par4)
{
  (*func)(str,&len,par2,par3,par4);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(caller5)
     (void (**func)(void*,void*,void*,void*,void*), void* par1, void* par2, void* par3, void* par4, void* par5)
{
  (*func)(par1,par2,par3,par4,par5);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(caller5_str5)
     (void (**func)(void*,void*,void*,void*,char*,int*), 
      void* par1, void* par2, void* par3, void* par4, char* str, int len)
{
  (*func)(par1,par2,par3,par4,str,&len);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(caller7_str67)
     (void (**func)(void*,void*,void*,void*,void*,char*,int*,char*,int*), 
      void* par1, void* par2, void* par3, void* par4, void* par5, char* str1, int len1, char* str2, int len2)
{
  (*func)(par1,par2,par3,par4,par5,str1,&len1,str2,&len2);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(caller)
     (void (**func)(void*, ... ), FINT* npar, ... )
{
  va_list ap;
  va_start(ap, npar);

  switch(*npar)
  {
  case 1: (*func)(va_arg(ap,void*)); break;
  case 2: (*func)(va_arg(ap,void*),va_arg(ap,void*)); break;
  //case 2: printf("f,p: %p %p \n",va_arg(ap,void*),va_arg(ap,void*)); break;
  case 3: (*func)(va_arg(ap,void*),va_arg(ap,void*),va_arg(ap,void*)); break;
  case 4: (*func)(va_arg(ap,void*),va_arg(ap,void*),va_arg(ap,void*),va_arg(ap,void*)); break;
  case 5: (*func)(va_arg(ap,void*),va_arg(ap,void*),va_arg(ap,void*),va_arg(ap,void*),va_arg(ap,void*)); break;
  default: return;
  }
  va_end(ap);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(caller0)
     (void (**func)(void))
{
   (*func)(); 
}
/* ---------------------------------------------------------------------- */
void *FTNIZE(dlopen_c)(const char *filename, FINT *flag)
{
 const int ALLNOW=1;
 void *pointer;
 void *p1;
 char *name;

 //dlerror();
 //pointer = dlopen(filename, RTLD_GLOBAL|RTLD_LAZY); //RTLD_LAZY); 

 if (*flag==ALLNOW)
   return dlopen(filename, RTLD_GLOBAL|RTLD_NOW); 
 else
   return dlopen(filename, RTLD_GLOBAL|RTLD_LAZY); 
 return NULL;
}
/* ---------------------------------------------------------------------- */
void *FTNIZE(dlsym_c)(void **handle, const char *symbol)
{
//printf("symbol address= %p\n", dlsym(*handle, symbol));
 return dlsym(*handle, symbol); 
}
/* ---------------------------------------------------------------------- */
void FTNIZE(dlclose_c)(void **handle)
{
 dlclose(*handle);
}
/* ---------------------------------------------------------------------- */
char* FTNIZE(dlerror_c)(void)
{
 char *error=dlerror();
 printf("error = %s\n", error);
 return error;
}       
/* ---------------------------------------------------------------------- */
void FTNIZE(write_binary_file_c)
     (char *filename, FINT *bytes, char *buffer, FINT *result)
/* Writes a given buffer to a binary file.
   Returns:
   * positive integer containing the number of written bytes
   * -2 if the file could not be opened
   * -1 if writing the buffer failed
*/
{
  int file = -1;
  int written = 0;

  *result = -2;
  file = open (filename, O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR);

  if(file == -1) return;

  *result = -1;

  written = (int) write (file, buffer, (size_t) *bytes);
  close (file);
  if (written != *bytes) return;
  *result = written;
}

/* ---------------------------------------------------------------------- */

void FTNIZE(get_pid_c)
     (FINT *pid)
/* Determines the PID of the current process.
   Returns:
   * integer containing the PID of the current process
   * -1 if retrieving the PID failed
*/
{
  pid_t result;

  *pid = -1;
  result = getpid ();
  if (result) *pid = (int) result;
}

/* ---------------------------------------------------------------------- */

void FTNIZE(get_env_var_c)
     (char *name, char *value)
/* Gets the content of an environment variable.
   Returns:
   * string containing the content of the environment variable, if available
   * empty string, if retrieving the environment variable failed
*/
{
  char *env_var;

  env_var = getenv (name);
  if (env_var) strncpy (value, env_var, strlen (env_var));
}

/* ---------------------------------------------------------------------- */

void FTNIZE(directory_exists_c)
     (char *path, FINT *exists)
/* Checks for existence of a directory.
   Returns:
   * 1, if 'path' points to a directory
   * -1, on error
   * 0, otherwise
*/
{
  int status;
  struct stat result;

  *exists = 0;
  status = stat (path, &result);
  if (status == -1) *exists = -1;
  if (S_ISDIR (result.st_mode)) *exists = 1;
}

/* ---------------------------------------------------------------------- */

void FTNIZE(is_nan_c)
     (REAL *value, FINT *result)
/* Determine if value is not a number.
   Returns:
   * 1, if value is not a number
   * 0, if value is a number
   * -1 on failure (value is neither float or double)
*/
{
  *result = -1;

  if (sizeof (*value) == sizeof (double)) *result = isnan ((double) *value);
  /*
    isnanf() is sometimes not available
    if (sizeof (*value) == sizeof (float)) *result = isnanf ((float) *value);
  */
  if (sizeof (*value) == sizeof (float)) *result = !(*value == *value);
}

/* ---------------------------------------------------------------------- */

void FTNIZE(system_c) (char *command)
/* Date:   04-Nov-2011
   Author: MR (matthias.rheinhardt@helsinki.fi)
   Description: function wrapper for ANSI C function system.
*/
{
  int res=system(command);
  if (res == -1) return; // some error handling is missing here [Bourdin.KIS]
}

/* ---------------------------------------------------------------------- */

void FTNIZE(sizeof_real_c)
     (REAL *value, FINT *result)
/* Determine the number of bytes used for value.
   Returns:
   * the number of bytes used for value
*/
{
  *result = sizeof (*value);
}
/* ---------------------------------------------------------------------- */
  void FTNIZE(copy_addr_c)(void *src, void **dest)
  {
    *dest=src;
  }
/* ---------------------------------------------------------------------- */

