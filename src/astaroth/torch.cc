
#include "torchfort.h"
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <mpi.h>

extern "C" void torch_trainCAPI(float* input, float* label, float* loss_val){

	torchfort_result_t result = torchfort_set_manual_seed(943442);

//	std::ofstream myFile;
//	myFile.open("t.out");


//	myFile << "This is the name of the model: " << std::string(model) << "\n";
	
	//printf("This is the name of the model: %s", model);
	//printf("Calling c API");
	int64_t input_arr[5] = {1, 3, 38, 38, 38};
	int64_t label_arr[5] = {1, 6, 38, 38, 38};
 torchfort_result_t res = torchfort_train("stationary", input, 5, input_arr, label, 5, label_arr, loss_val, TORCHFORT_FLOAT, 0);

 if(res != TORCHFORT_RESULT_SUCCESS)
 {
 		fprintf(stderr,"torchfort_train failed!\n");
 }
}


extern "C" void torch_inferCAPI(float* input, float* label){

	torchfort_result_t result = torchfort_set_manual_seed(943442);
	int64_t input_shape[5] = {1, 3, 38, 38, 38};
	int64_t label_shape[5] = {1, 6, 38, 38, 38};

	torchfort_result_t res = torchfort_inference("stationary", input, 5, input_shape, label, 5, label_shape, TORCHFORT_FLOAT, 0);


 if(res != TORCHFORT_RESULT_SUCCESS)
 {
 		fprintf(stderr,"torchfort_train failed!\n");
 }
}




extern "C" void torch_createmodel(const char* name, const char* config_fname,MPI_Comm mpi_comm, int device){
	
	torchfort_result_t result = torchfort_set_manual_seed(943442);
	torchfort_result_t res = torchfort_create_distributed_model("cnn", "data/training/unet_torchscript.pt", MPI_COMM_WORLD, 0);

 	if(res != TORCHFORT_RESULT_SUCCESS)
 	{
		fprintf(stderr,"torchfort_train failed!\n");
 	}
}
