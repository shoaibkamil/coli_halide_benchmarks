/**

  var tol = 1e-12;
  var maxiters=100;
  var r = b - (A*xguess);
  var p = r;
  var iter = 0;
  var x = xguess;
  var normr2 = r' * r;
  while (normr2 > tol) and (iter < maxiters)
    Ap = A * p;
    denom = p' * Ap;
    alpha = normr2 / denom;
    x = x + alpha * p;
    normr2old = normr2;
    r = r - alpha * Ap;
    normr2 = r' * r;
    beta = normr2 / normr2old;
    p = r + beta * p;
    iter = iter + 1;
  end
**/
#include <cstdio>
#include "Halide.h"

#include "saxpy.h"
#include "sdot.h"
#include "sgemv.h"
#include "sveccopy.h"
#include "svecsub.h"

using namespace Halide;

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

  ret->set_device_dirty();
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


  ret->set_device_dirty();
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
  mat_ret[1] = 1;
  mat_ret[2] = 0.5;
  mat_ret[3] = 0.4;
  mat_ret[4] = 5;
  mat_ret[5] = 0.3;
  mat_ret[6] = 0.6;
  mat_ret[7] = 0.8;
  mat_ret[8] = 5;

  b_ret[0] = 1.f;
  b_ret[1] = 2.f;
  b_ret[2] = 3.f;

  *mat = mat_ret;
  *b = b_ret;

  // expected result: {0.074, 0.362, 0.533}
} 

// CG
// Inputs: Matrix A, RHS vector b, tolerance, maxiters
// Temporaries: p, r
// Outputs: x, where A*x = b
void cg_1(halide_buffer_t* A, halide_buffer_t* b, float tol, int maxiters,
          halide_buffer_t* x,
          halide_buffer_t* p, halide_buffer_t* r, halide_buffer_t* tmp) {
  // allocate temporaries
  int iter = 0;
  float normr2;
  float alpha;
  float normr2old;
  float denom;
  float beta;

  halide_buffer_t normr2_buf;
  halide_buffer_t denom_buf;

  // wrap 0-dim temporaries
  wrap_vector(&normr2_buf, &normr2, 1);
  wrap_vector(&denom_buf, &denom, 1);

  // calculate residual
  sgemv(1.0f, A, x, 0.0f, tmp, tmp);
  print_vec("tmp", tmp, 3);
  svecsub(tmp, b, r);

  // calculate norm squared
  sdot(r, r, &normr2_buf);

  // copy r to p
  sveccopy(r, p);

  // main loop
  while (normr2 > tol && iter < maxiters) {
    // calculate denom
    sgemv(1.0f, A, p, 0.0f, tmp, tmp);
    sdot(r, tmp, &denom_buf);

    alpha = normr2 / denom;

    // update x
    // TODO: make sure the order here is correct
    saxpy(alpha, p, x, x);

    // update normr2
    normr2old = normr2;
    saxpy(-alpha, tmp, r, r);
    sdot(r, r, &normr2_buf);

    // update p
    beta = normr2 / normr2old;
    saxpy(beta, r, p, p);
    iter++;
  }
}



int main(int argc, char* argv[]) {

  float tol = 1e-4;
  int maxiters = 100;

  // generate test data
  float* mat;
  float* b;
  halide_buffer_t mat_buf;
  halide_buffer_t b_buf;

  generate_test_data(&mat, &b);
  wrap_vector(&b_buf, b, 3);
  wrap_matrix(&mat_buf, mat, 3);

  // create x and temporaries
  float x[3] = {0};
  float p[3] = {0};
  float r[3] = {0};
  float tmp[3] = {0};
  halide_buffer_t x_buf;
  halide_buffer_t p_buf;
  halide_buffer_t r_buf;
  halide_buffer_t tmp_buf;

  wrap_vector(&x_buf, x, 3);
  wrap_vector(&p_buf, p, 3);
  wrap_vector(&r_buf, r, 3);
  wrap_vector(&tmp_buf, tmp, 3);

  // call CG
  cg_1(&mat_buf, &b_buf, tol, maxiters,
       &x_buf,
       &p_buf, &r_buf, &tmp_buf);

  print_vec("x", x, 3);

  return 0;
}
