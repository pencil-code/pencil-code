#include "torchfort.h"
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <mpi.h>

void torch_trainCAPI(int sub_dims[3], float* input, float* label, float* loss_val, bool dble=false){

        torchfort_datatype_t precision = dble ? TORCHFORT_DOUBLE : TORCHFORT_FLOAT;
	torchfort_result_t result = torchfort_set_manual_seed(943442);

//	std::ofstream myFile;
//	myFile.open("t.out");
//	myFile << "This is the name of the model: " << std::string(model) << "\n";
	
	//printf("This is the name of the model: %s", model);
	//printf("Calling c API");
	int64_t input_shape[5] = {5, 3, sub_dims[2], sub_dims[1], sub_dims[0]};
	int64_t label_shape[5] = {5, 6, sub_dims[2], sub_dims[1], sub_dims[0]};
  	torchfort_result_t res = torchfort_train("stationary", input, 5, input_shape, label, 5, label_shape, loss_val, precision, 0);

 	if (res != TORCHFORT_RESULT_SUCCESS)
 	{
 		fprintf(stderr,"torchfort_train failed!\n");
 	}
}

void torch_inferCAPI(int sub_dims[3], float* input, float* label, bool dble=false){

        torchfort_datatype_t precision = dble ? TORCHFORT_DOUBLE : TORCHFORT_FLOAT;
	torchfort_result_t result = torchfort_set_manual_seed(943442);

	int64_t input_shape[5] = {1, 3, sub_dims[2], sub_dims[1], sub_dims[0]};
	int64_t label_shape[5] = {1, 6, sub_dims[2], sub_dims[1], sub_dims[0]};
	torchfort_result_t res = torchfort_inference("stationary", input, 5, input_shape, label, 5, label_shape, precision, 0);


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
