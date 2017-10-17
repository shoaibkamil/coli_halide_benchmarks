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
#include <cstdlib>
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
  svecsub(b, tmp, r);

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
    saxpy(alpha, p, x, x);

    // update normr2
    normr2old = normr2;
    saxpy(-alpha, tmp, r, r);

    sdot(r, r, &normr2_buf);

    // update p
    beta = normr2 / normr2old;
    saxpy(beta, p, r, p);
    iter++;
  }
}

int test_small() {

  float tol = 1e-5;
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

  // check answer
  float normr2;
  halide_buffer_t normr2_buf;
  wrap_vector(&normr2_buf, &normr2, 1);
  sgemv(1.0f, &mat_buf, &x_buf, 0.0f, &tmp_buf, &tmp_buf);
  svecsub(&tmp_buf, &b_buf, &r_buf);
  sdot(&r_buf, &r_buf, &normr2_buf);
  print_vec("answer", &tmp_buf, 3);
  print_vec("expected", &b_buf, 3);

  if (normr2 >= tol) {
    printf("|r|^2 = %f\n", normr2);
    printf("FAILED\n");
    return -1;
  }

  printf("SUCCESS\n");
  return 0;
}

int test_rand(int size) {

  float tol = 1e-5;
  int maxiters = 100;

  // generate test data
  float* mat;
  float* b;
  halide_buffer_t mat_buf;
  halide_buffer_t b_buf;

  generate_random_test_data(&mat, &b, size);
  wrap_vector(&b_buf, b, size);
  wrap_matrix(&mat_buf, mat, size);

  // create x and temporaries
  float x[size];
  float p[size];
  float r[size];
  float tmp[size];
  for (int i=0; i<size; i++)
    x[i] = p[i] = r[i] = tmp[i] = 0.0f;

  halide_buffer_t x_buf;
  halide_buffer_t p_buf;
  halide_buffer_t r_buf;
  halide_buffer_t tmp_buf;

  wrap_vector(&x_buf, x, size);
  wrap_vector(&p_buf, p, size);
  wrap_vector(&r_buf, r, size);
  wrap_vector(&tmp_buf, tmp, size);

  // call CG
  cg_1(&mat_buf, &b_buf, tol, maxiters,
       &x_buf,
       &p_buf, &r_buf, &tmp_buf);

  // check answer
  float normr2;
  halide_buffer_t normr2_buf;
  wrap_vector(&normr2_buf, &normr2, 1);
  sgemv(1.0f, &mat_buf, &x_buf, 0.0f, &tmp_buf, &tmp_buf);
  svecsub(&tmp_buf, &b_buf, &r_buf);
  sdot(&r_buf, &r_buf, &normr2_buf);
  print_vec("answer", &tmp_buf, size);
  print_vec("expected", &b_buf, size);

  if (normr2 >= tol) {
    printf("|r|^2 = %f\n", normr2);
    printf("FAILED\n");
    return -1;
  }

  printf("SUCCESS\n");
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    printf("Usage: %s <size>\n", argv[0]);
    return 1;
  }
  return test_rand(atoi(argv[1]));
}
