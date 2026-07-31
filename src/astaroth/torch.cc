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

/***********************************************************************************************/
bool
initialize_torch()
{
	const torchfort_result_t res = torchfort_set_manual_seed(943442);
	return res == TORCHFORT_RESULT_SUCCESS;
}
/***********************************************************************************************/
bool torch_train_CAPI(int sub_dims[3], AcReal* input, AcReal* label, AcReal* loss_val,
		     const int input_fields, const int output_fields, const char* model_name){


	int64_t input_shape[5] = {1, input_fields,  sub_dims[2], sub_dims[1], sub_dims[0]};
	int64_t label_shape[5] = {1, output_fields, sub_dims[2], sub_dims[1], sub_dims[0]};
  	const torchfort_result_t res = torchfort_train(model_name, input, 5, input_shape, label, 5, label_shape, loss_val, TORCH_PRECISION, 0);
	return res == TORCHFORT_RESULT_SUCCESS;
}

/***********************************************************************************************/
bool torch_infer_CAPI(int sub_dims[3], AcReal* input, AcReal* label, 
		     const int input_fields, const int output_fields, const char* model_name, bool subsample){

	torchfort_result_t res = torchfort_set_manual_seed(943442);

	int64_t input_shape[5] = {1, input_fields, sub_dims[2], sub_dims[1], sub_dims[0]};
	
	int64_t label_shape[5] = {1, output_fields, sub_dims[2], sub_dims[1], sub_dims[0]};

/*
	if (subsample) {
    int64_t new_vals[] = {1, output_fields, sub_dims[2]/7, sub_dims[1]/7, sub_dims[0]/7};
    std::copy(std::begin(new_vals), std::end(new_vals), std::begin(label_shape));
	}
*/
	res = torchfort_inference(model_name, input, 5, input_shape, label, 5, label_shape, TORCH_PRECISION, 0);
	return res == TORCHFORT_RESULT_SUCCESS;
}

/***********************************************************************************************/
bool torch_create_model_CAPI(const char* name, const char* config_fname, int device){
	
	const torchfort_result_t res = torchfort_create_model(name, config_fname, device);
	return res == TORCHFORT_RESULT_SUCCESS;
}

/***********************************************************************************************/
bool torch_create_distributed_model_CAPI(const char* name, const char* config_fname, MPI_Comm mpi_comm, int device){
	
	const torchfort_result_t res = torchfort_create_distributed_model(name, config_fname, mpi_comm, device);
	return res == TORCHFORT_RESULT_SUCCESS;
}
/***********************************************************************************************/
bool torch_load_CAPI(const char* name, const char* fname){
	const torchfort_result_t res = torchfort_load_model(name, fname);
	return res == TORCHFORT_RESULT_SUCCESS;
}
/***********************************************************************************************/
bool torch_load_checkpoint_CAPI(const char* name, const char* checkpoint_dir, int64_t* step_train, int64_t* step_inference){
	const torchfort_result_t res = torchfort_load_checkpoint(name, checkpoint_dir, step_train, step_inference);
	return res == TORCHFORT_RESULT_SUCCESS;
}
/***********************************************************************************************/
bool torch_save_model_CAPI(const char* name, const char* fname){
	const torchfort_result_t res = torchfort_save_model(name, fname);
	return res == TORCHFORT_RESULT_SUCCESS;
}
/***********************************************************************************************/
bool torch_save_checkpoint_CAPI(const char* name, const char* checkpoint_dir){
	const torchfort_result_t res = torchfort_save_checkpoint(name, checkpoint_dir);
	return res == TORCHFORT_RESULT_SUCCESS;
}
/***********************************************************************************************/
