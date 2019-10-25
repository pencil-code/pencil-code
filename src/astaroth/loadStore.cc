//                             loadStore.cc
//                             --------------------

/* Date:   6-Jun-2017
   Author: M. Rheinhardt
   Description: Copier functions for the different "plates" of the halo and the full inner data cube with host-device concurrency.
                Load balance yet to be established.
*/

//C libraries
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
using namespace std;

#include "astaroth.h"
#include "/homeappl/home/mreinhar/git/pencil-code/samples/gputest/src/cparam_c.h"
#include "/homeappl/home/mreinhar/git/pencil-code/samples/gputest/src/cdata_c.h"

const int mxy=mx*my;
extern int halo_yz_size;
static AcReal *halo_yz_buffer; 
extern Node node;

const int BOT=0, TOP=1, TOT=2;
int halo_widths_x[3]={nghost,nghost,2*nghost};    // bottom and top halo width and sum of them
int halo_widths_y[3]={nghost,nghost,2*nghost};
int halo_widths_z[3]={nghost,nghost,2*nghost};

void initLoadStore()
{
//printf("lperi= %d %d %d \n", lperi[0],lperi[1],lperi[2]);
        // halo widths for undivided data cube
        if (!lperi[0]){
                if (lfirst_proc_x) halo_widths_x[BOT]=nghost+1;
                if (llast_proc_x) halo_widths_x[TOP]=nghost+1;
        }
        if (!lyinyang) {
                if (!lperi[1]){
                        if (lfirst_proc_y) halo_widths_y[BOT]=nghost+1;
                        if (llast_proc_y) halo_widths_y[TOP]=nghost+1;
                }
                if (!lperi[2]){
                        if (lfirst_proc_z) halo_widths_z[BOT]=nghost+1;
                        if (llast_proc_z) halo_widths_z[TOP]=nghost+1;
                }
        }
        halo_widths_x[TOT]=halo_widths_x[BOT]+halo_widths_x[TOP];
        halo_widths_y[TOT]=halo_widths_y[BOT]+halo_widths_y[TOP];
        halo_widths_z[TOT]=halo_widths_z[BOT]+halo_widths_z[TOP];

//printf("halo_widths_x= %d %d %d\n",halo_widths_x[BOT],halo_widths_x[TOP],halo_widths_x[TOT]);
//printf("halo_widths_y= %d %d %d\n",halo_widths_y[BOT],halo_widths_y[TOP],halo_widths_x[TOT]);
//printf("halo_widths_z= %d %d %d\n",halo_widths_z[BOT],halo_widths_z[TOP],halo_widths_x[TOT]);
        halo_yz_size = ny*nz*max(halo_widths_x[BOT],halo_widths_x[TOP])*NUM_VTXBUF_HANDLES;
        if (halo_yz_buffer==NULL) halo_yz_buffer=(AcReal*) malloc(halo_yz_size*sizeof(AcReal));    // buffer for yz halos in host

}
/****************************************************************************************************************/
void finalLoadStore()
{
        free(halo_yz_buffer);
}
/****************************************************************************************************************/
void loadOuterFront(AcMesh& mesh, Stream stream)
{
        int3 src={0,0,0};
        int num_vertices=mxy*halo_widths_z[BOT];
        acNodeLoadMeshWithOffset(node, stream, mesh, src, src, num_vertices);

//printf("front:num_vertices= %d \n",num_vertices);
        //!!!cudaHostRegister(mesh, size, cudaHostRegisterDefault);    // time-critical!
}

void loadOuterBack(AcMesh& mesh, Stream stream)
{
        int3 src={0,0,m2-halo_widths_z[TOP]};      // index from m2-halo_widths_z[TOP] to m2-1
        int num_vertices=mxy*halo_widths_z[TOP];
        acNodeLoadMeshWithOffset(node, stream, mesh, src, src, num_vertices);

//printf("back:num_vertices= %d \n",num_vertices);
        //!!!cudaHostRegister(mesh, size, cudaHostRegisterDefault);    // time-critical!
}
 
void loadOuterBot(AcMesh& mesh, Stream stream)
{
        int3 src={0,0,halo_widths_z[BOT]};
        int num_vertices=mx*halo_widths_y[BOT];
//printf("bot:num_vertices= %d %d \n",num_vertices, mz-halo_widths_z[TOT]);
        for (int i=0; i<mz-halo_widths_z[TOT]; i++) {
            acNodeLoadMeshWithOffset(node, stream, mesh, src, src, num_vertices);
            src.z++;

        //!!!cudaHostRegister(mesh, size, cudaHostRegisterDefault);    // time-critical!
        }
}

void loadOuterTop(AcMesh& mesh, Stream stream)
{
        int3 src={0,my-halo_widths_y[TOP],halo_widths_z[BOT]};
        int num_vertices=mx*halo_widths_y[TOP];
//printf("top:num_vertices= %d %d\n",num_vertices,mz-halo_widths_z[TOT]);
        for (int i=0; i<mz-halo_widths_z[TOT]; i++) {
            acNodeLoadMeshWithOffset(node, stream, mesh, src, src, num_vertices);
            src.z++;

        //!!!cudaHostRegister(mesh, size, cudaHostRegisterDefault);    // time-critical!
        }
}

void loadOuterLeft(AcMesh& mesh, Stream stream)
{
    int3 start={0,                      halo_widths_y[BOT],     halo_widths_z[BOT]  };
    int3 end  ={halo_widths_x[BOT]-1,my-halo_widths_y[TOP]-1,mz-halo_widths_z[TOP]-1};

    acLoadYZPlate(start, end, &mesh, halo_yz_buffer);
}

void loadOuterRight(AcMesh& mesh, Stream stream)
{
    int3 start={mx-halo_widths_x[TOP],   halo_widths_y[BOT],     halo_widths_z[BOT]  };
    int3 end  ={mx-1,                 my-halo_widths_y[TOP]-1,mz-halo_widths_z[TOP]-1};

    acLoadYZPlate(start, end, &mesh, halo_yz_buffer);
}

void loadOuterHalos(AcMesh& mesh)
{
    loadOuterFront(mesh,STREAM_DEFAULT);
    loadOuterBack(mesh,STREAM_DEFAULT);
    loadOuterTop(mesh,STREAM_DEFAULT);
    loadOuterBot(mesh,STREAM_DEFAULT);
    loadOuterLeft(mesh,STREAM_DEFAULT);
    loadOuterRight(mesh,STREAM_DEFAULT);
}

