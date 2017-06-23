void initializeCopying();
void finalizeCopying();
void copyOuterHalos(float* grid, float* d_grid);
void copyInnerHalos(float* grid, float* d_grid);
void copyInnerAll(float* grid, float* d_grid);
void copyAll(float* grid, float* d_grid);
