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

using namespace Halide;

// CG
// Inputs: Matrix A, RHS vector b, tolerance, maxiters
// Outputs: x, where A*x = b
void cg_1(halide_buffer_t* A, halide_buffer_t* b, float tol, int maxiters,
          halide_buffer_t* x) {
  // allocate temporaries
  int iter = 0;

  // calculate residual
  sgemv(A, x, tmp);
  svecsub(tmp, b, r);

  // calculate norm squared
  sdot(r, r, normr2);

  // copy r to p
  sveccopy(r, p);

  // main loop
  while (normr2 > tol && iter < maxiters) {
    // calculate denom
    sgemv(A, p, tmp);
    sdot(r, tmp, denom);

    alpha = normr2 / denom;

    // update x
    // TODO: make sure the order here is correct
    saxpy(x, p, alpha, x);

    // update normr2
    float normr2old = normr2;
    saxpy(r, tmp, -alpha, r);
    sdot(r, r, normr2);

    // update p
    float beta = normr2 / normr2old;
    saxpy(p, r, beta, p);
    iter++;
  }
}




int main(int argc, char* argv[]) {

  float tol = 1e-4;
  int maxiters = 100;



}
