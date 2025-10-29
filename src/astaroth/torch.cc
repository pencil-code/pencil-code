#include "torchfort.h"
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <mpi.h>

//TP: ugly but works
#if AC_DOUBLE_PRECISION
typedef double AcReal;
#define TORCH_PRECISION TORCHFORT_DOUBLE
#else
typedef float  AcReal;
#define TORCH_PRECISION TORCHFORT_FLOAT
#endif

void torch_trainCAPI(int sub_dims[3], AcReal* input, AcReal* label, AcReal* loss_val){

	torchfort_result_t result = torchfort_set_manual_seed(943442);

//	std::ofstream myFile;
//	myFile.open("t.out");
//	myFile << "This is the name of the model: " << std::string(model) << "\n";
	
	//printf("This is the name of the model: %s", model);
	//printf("Calling c API");
	int64_t input_shape[5] = {1, 3, sub_dims[2], sub_dims[1], sub_dims[0]};
	int64_t label_shape[5] = {1, 6, sub_dims[2], sub_dims[1], sub_dims[0]};
  	torchfort_result_t res = torchfort_train("stationary", input, 5, input_shape, label, 5, label_shape, loss_val, TORCH_PRECISION, 0);

 	if (res != TORCHFORT_RESULT_SUCCESS)
 	{
 		fprintf(stderr,"torchfort_train failed!\n");
 	}
}

void torch_inferCAPI(int sub_dims[3], AcReal* input, AcReal* label){

	torchfort_result_t result = torchfort_set_manual_seed(943442);

	int64_t input_shape[5] = {1, 3, sub_dims[2], sub_dims[1], sub_dims[0]};
	int64_t label_shape[5] = {1, 6, sub_dims[2], sub_dims[1], sub_dims[0]};
	torchfort_result_t res = torchfort_inference("stationary", input, 5, input_shape, label, 5, label_shape, TORCH_PRECISION, 0);


 	if (res != TORCHFORT_RESULT_SUCCESS)
 	{
 		fprintf(stderr,"torchfort_train failed!\n");
 	}
}

void torch_createmodel(const char* name, const char* config_fname,MPI_Comm mpi_comm, int device){
	
	torchfort_result_t result = torchfort_set_manual_seed(943442);
	torchfort_result_t res = torchfort_create_distributed_model("cnn", "data/training/unet_torchscript.pt", MPI_COMM_WORLD, device);

 	if (res != TORCHFORT_RESULT_SUCCESS)
 	{
		fprintf(stderr,"torchfort_train failed!\n");
 	}
}
