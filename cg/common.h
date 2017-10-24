#ifndef TIRAMISU_CG_COMMON_H
#define TIRAMISU_CG_COMMON_H
#include "HalideRuntime.h"

#include <cstdio>
#include <cstdlib>
#include <string>

// wrap a vector in a halide_buffer_t
inline void wrap_vector(halide_buffer_t* ret, float* data, int size) {

  ret->device = 0;
  ret->device_interface = 0;
  ret->host = (uint8_t*)data;
  ret->flags = 0;
  ret->type = halide_type_of<float>();
  ret->dimensions = 1;

  ret->dim = (halide_dimension_t*)malloc(sizeof(halide_dimension_t));
  ret->dim[0].min = 0;
  ret->dim[0].extent = size;
  ret->dim[0].stride = 1;

  ret->set_host_dirty();
}

// wrap a matrix in a halide_buffer_t
inline void wrap_matrix(halide_buffer_t* ret, float* data, int size) {

  ret->device = 0;
  ret->device_interface = 0;
  ret->host = (uint8_t*)data;
  ret->flags = 0;
  ret->type = halide_type_of<float>();
  ret->dimensions = 2;

  ret->dim = (halide_dimension_t*)malloc(2*sizeof(halide_dimension_t));
  ret->dim[0].min = 0;
  ret->dim[0].extent = size;
  ret->dim[0].stride = 1;

  ret->dim[1].min = 0;
  ret->dim[1].extent = size;
  ret->dim[1].stride = size;


  ret->set_host_dirty();
}

void print_vec(std::string name, float* vec, int size) {
  printf("%s: ", name.c_str());
  for (int i=0; i<size; i++)
    printf("%3.3f ", vec[i]);
  printf("\n");
}

void print_vec(std::string name, halide_buffer_t* vec, int size) {
  print_vec(name, (float*)(vec->host), size);
}

void generate_test_data(float** mat, float** b) {
  float *mat_ret = (float*)malloc(sizeof(float) * 9);
  float *b_ret = (float*)malloc(sizeof(float) * 3);

  mat_ret[0] = 5;
  mat_ret[3] = 1;
  mat_ret[6] = 0.5;
  mat_ret[1] = 0.4;
  mat_ret[4] = 5;
  mat_ret[7] = 0.3;
  mat_ret[2] = 0.6;
  mat_ret[5] = 0.8;
  mat_ret[8] = 5;

  b_ret[0] = 1.f;
  b_ret[1] = 2.f;
  b_ret[2] = 3.f;

  *mat = mat_ret;
  *b = b_ret;

  // expected result: {0.074, 0.362, 0.533}
} 

void generate_random_test_data(float** mat, float** b, int size) {
  float *mat_ret = (float*)malloc(sizeof(float) * size * size);
  float *b_ret = (float*)malloc(sizeof(float) * size);

  srand(size * (uint64_t)mat_ret);

  for (int i=0; i<size; i++)
    b_ret[i] = 10 * i + 10*(float)rand()/RAND_MAX;

  for (int i=0; i<size; i++)
    for (int j=0; j<size; j++)
      mat_ret[i+size*j] = (float)rand()/RAND_MAX;

  for (int i=0; i<size; i++)
    mat_ret[i+size*i] = size*10;

  *mat = mat_ret;
  *b = b_ret;
}
#endif
