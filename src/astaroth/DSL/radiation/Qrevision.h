input int3 dir, stop

Kernel Qrevision(int3 dir, int3 stop){

  int3 boundary = get_boundary(get_normal())
  int3 sign=sgn(dir)
  stop *= sign

  l,m,n = boundary.x, boundary.y, boundary.z
  Qrad0 = QRAD[l-dir.x][m-dir.y][n-dir.z]    // would not work if the ray is not tied to the thread

  while (l*sign.x < stop.x && m*sign.y < stop.y && n*sign.z < stop.z){

    QRAD[l][m][n] += Qrad0*exp(-TAU[l][m][n])

    l += dir.x
    m += dir.y
    n += dir.z
  }
}


/*calling:

dim3 threadsPerBlock(32, 32);
for xy plane:
dim3 numBlocks(mx/threadsPerBlock.x,my/threadsPerBlock.y);

for xz plane:
dim3 numBlocks(mx/threadsPerBlock.x,mz/threadsPerBlock.z);

for yz plane:
dim3 numBlocks(my/threadsPerBlock.y,mz/threadsPerBlock.z);
*/
