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

#include "Halide.h"

#include "saxpy.h"
#include "sdot.h"
#include "sgemv.h"
#include "sveccopy.h"
#include "svecsub.h"

using namespace Halide;

// wrap a vector in a halide_buffer_t
void wrap_vector(halide_buffer_t* ret, float* data, int size) {

  ret->device = 0;
  ret->device_interface = 0;
  ret->host = (uint8_t*)data;
  ret->flags = 0;
  ret->type = halide_type_of<float>();
  ret->dimensions = 2;

  ret->dim = (halide_dimension_t*)malloc(sizeof(halide_dimension_t));
  ret->dim[0].min = 0;
  ret->dim[0].extent = size;
  ret->dim[0].stride = 1;

  ret->set_device_dirty();
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



}
