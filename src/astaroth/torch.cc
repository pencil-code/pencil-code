#include "torchfort.h"
#include <iostream>
#include <vector>
#include <utility>
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
	return res != TORCHFORT_RESULT_SUCCESS;
}
/***********************************************************************************************/
bool torch_wandb_log_double(const char* name, const char* metric_name, int64_t step, double value){
  const torchfort_result_t res = torchfort_wandb_log_double(name, metric_name, step, value);
}
/***********************************************************************************************/
bool torch_train_CAPI(int sub_dims[3], AcReal* input, AcReal* label, AcReal* loss_val,
		     const int input_fields, const int output_fields, const char* model_name){


	int64_t input_shape[5] = {1, input_fields,  sub_dims[2], sub_dims[1], sub_dims[0]};
	int64_t label_shape[5] = {1, output_fields, sub_dims[2], sub_dims[1], sub_dims[0]};
  	const torchfort_result_t res = torchfort_train(model_name, input, 5, input_shape, label, 5, label_shape, loss_val, TORCH_PRECISION, 0);
	return res != TORCHFORT_RESULT_SUCCESS;
}
/***********************************************************************************************/
bool torch_train_multiarg_CAPI(int sub_dims[3], const std::vector<std::pair<AcReal*, int>>& inputs, 
        const std::vector<std::pair<AcReal*, int>>& outputs, AcReal* loss_val,  const char* model_name){
    

    torchfort_tensor_list_t inputs_tensor;
    torchfort_tensor_list_create(&inputs_tensor);

    torchfort_tensor_list_t outputs_tensor;
    torchfort_tensor_list_create(&outputs_tensor);

    for(std::size_t i = 0; i < inputs.size(); ++i) {
      if (inputs[i].first == NULL) {
        fprintf(stderr, "ERROR: ASTAROTH device pointers are NULL! (input: %p), %d\n", (void*)inputs[i].first, i);
        fflush(stderr);
      }
	    int64_t input_shape[5] = {1, inputs[i].second,  sub_dims[2], sub_dims[1], sub_dims[0]};
        torchfort_tensor_list_add_tensor(inputs_tensor, inputs[i].first, 5, input_shape, TORCH_PRECISION);
    }
    
    for(std::size_t i = 0; i < outputs.size(); ++i) {
      if (outputs[i].first == NULL) {
        fprintf(stderr, "ERROR: ASTAROTH device pointers are NULL! (outputs: %p), %d\n", (void*)outputs[i].first, i);
        fflush(stderr);
      }
	    int64_t output_shape[5] = {1, outputs[i].second,  sub_dims[2], sub_dims[1], sub_dims[0]};
        torchfort_tensor_list_add_tensor(outputs_tensor, outputs[i].first, 5, output_shape, TORCH_PRECISION);
    }

    const torchfort_result_t res = torchfort_train_multiarg(model_name, inputs_tensor, outputs_tensor, loss_val, nullptr, 0);
    torchfort_result_t res_destroy = torchfort_tensor_list_destroy(inputs_tensor);
    res_destroy = torchfort_tensor_list_destroy(outputs_tensor);
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
	return res != TORCHFORT_RESULT_SUCCESS;
}

/***********************************************************************************************/
bool torch_create_model_CAPI(const char* name, const char* config_fname, int device){
	
	const torchfort_result_t res = torchfort_create_model(name, config_fname, device);
	return res != TORCHFORT_RESULT_SUCCESS;
}

/***********************************************************************************************/
bool torch_create_distributed_model_CAPI(const char* name, const char* config_fname, MPI_Comm mpi_comm, int device){
	
	const torchfort_result_t res = torchfort_create_distributed_model(name, config_fname, mpi_comm, device);
	return res != TORCHFORT_RESULT_SUCCESS;
}
/***********************************************************************************************/
bool torch_load_CAPI(const char* name, const char* fname){
	const torchfort_result_t res = torchfort_load_model(name, fname);
	return res != TORCHFORT_RESULT_SUCCESS;
}
/***********************************************************************************************/
bool torch_load_checkpoint_CAPI(const char* name, const char* checkpoint_dir, int64_t* step_train, int64_t* step_inference){
	const torchfort_result_t res = torchfort_load_checkpoint(name, checkpoint_dir, step_train, step_inference);
	return res != TORCHFORT_RESULT_SUCCESS;
}
/***********************************************************************************************/
bool torch_save_model_CAPI(const char* name, const char* fname){
	const torchfort_result_t res = torchfort_save_model(name, fname);
	return res != TORCHFORT_RESULT_SUCCESS;
}
/***********************************************************************************************/
bool torch_save_checkpoint_CAPI(const char* name, const char* checkpoint_dir){
	const torchfort_result_t res = torchfort_save_checkpoint(name, checkpoint_dir);
	return res != TORCHFORT_RESULT_SUCCESS;
}
/***********************************************************************************************/
