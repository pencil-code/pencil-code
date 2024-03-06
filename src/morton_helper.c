#pragma once
#include <stdint.h> // uint64_t
#include <stdio.h>
#include <stdlib.h>
#include "headers_c.h"
typedef struct int3{
	int x, y, z;
} int3;
typedef struct uint3_64 {
  uint64_t x, y, z;
} uint3_64;
/*
 * UINT3_64
 */




static inline uint64_t
mod(const int a, const int b)
{
  const int r = a % b;
  return r < 0 ? (size_t)(r + b) : (size_t)(r);
}



#define MPI_DECOMPOSITION_AXES (3)
static inline uint3_64
morton3D(const uint64_t pid)
{
    uint64_t i, j, k;
    i = j = k = 0;

    if (MPI_DECOMPOSITION_AXES == 3) {
        for (int bit = 0; bit <= 21; ++bit) {
            const uint64_t mask = 0x1l << 3 * bit;
            k |= ((pid & (mask << 0)) >> 2 * bit) >> 0;
            j |= ((pid & (mask << 1)) >> 2 * bit) >> 1;
            i |= ((pid & (mask << 2)) >> 2 * bit) >> 2;
        }
    }
    return (uint3_64){i, j, k};
}

static inline uint64_t
morton1D(const uint3_64 pid)
{
    uint64_t i = 0;

    if (MPI_DECOMPOSITION_AXES == 3) {
        for (int bit = 0; bit <= 21; ++bit) {
            const uint64_t mask = 0x1l << bit;
            i |= ((pid.z & mask) << 0) << 2 * bit;
            i |= ((pid.y & mask) << 1) << 2 * bit;
            i |= ((pid.x & mask) << 2) << 2 * bit;
        }
    }
    else if (MPI_DECOMPOSITION_AXES == 2) {
        for (int bit = 0; bit <= 21; ++bit) {
            const uint64_t mask = 0x1l << bit;
            i |= ((pid.y & mask) << 0) << 1 * bit;
            i |= ((pid.z & mask) << 1) << 1 * bit;
        }
    }
    else if (MPI_DECOMPOSITION_AXES == 1) {
        for (int bit = 0; bit <= 21; ++bit) {
            const uint64_t mask = 0x1l << bit;
            i |= ((pid.z & mask) << 0) << 0 * bit;
        }
    }
    else {
        fprintf(stderr, "Invalid MPI_DECOMPOSITION_AXES\n");
    }

    return i;
}

static inline uint3_64
decompose(const int target_in)
{
    const uint64_t target = (uint64_t) target_in;
    // This is just so beautifully elegant. Complex and efficient decomposition
    // in just one line of code.
    uint3_64 p = morton3D(target - 1);
    return (uint3_64){p.x+1,p.y+1,p.z+1};
}

static inline uint3_64
wrap(const int3 i, const uint3_64 n)
{
    return (uint3_64){
	mod(i.x, n.x),
        mod(i.y, n.y),
        mod(i.z, n.z),
    };
}
static inline int3
getPid3D(const uint64_t pid, const uint3_64 decomp)
{
    const uint3_64 pid3D = morton3D(pid);
    return (int3){(int)pid3D.x, (int)pid3D.y, (int)pid3D.z};
}
int3
FTNIZE(getmortonrank3d)(const int* pid, const int* decomp_x, const int* decomp_y, const int* decomp_z)
{
    return getPid3D((uint64_t)*pid, (uint3_64){(uint64_t)*decomp_x,(uint64_t)*decomp_y,(uint64_t)*decomp_z});
}

static inline int
getPid(const int3 pid_raw, const uint3_64 decomp)
{
    const uint3_64 pid = wrap(pid_raw, decomp);
    return (int)morton1D(pid);
}
int
FTNIZE(getmortonrank)(const int* x, const int* y, const int* z, const int* decomp_x, const int* decomp_y, const int* decomp_z)
{
	return getPid((int3){*x,*y,*z}, (uint3_64){(uint64_t)*decomp_x,(uint64_t)*decomp_y,(uint64_t)*decomp_z});
}

/**
int main(int argc, char* argv[])
{
	int nproc = atoi(argv[1]);

	uint3_64 decomp = decompose(nproc);
	printf("for %d procs ", nproc);
	printf("the morton decomp would be: %ld,%ld,%ld\n", decomp.x,decomp.y,decomp.z);
	return 0;
	
}
**/
