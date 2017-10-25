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

#include "common.h"
#include "mkl.h"

// CG
// Inputs: Matrix A, RHS vector b, tolerance, maxiters
// Temporaries: p, r
// Outputs: x, where A*x = b
void cg_mkl(int rows, float* A, float* b, float tol, int maxiters,
          float* x,
          float* p, float* r, float* tmp) {
  // allocate temporaries
  int iter = 0;
  float normr2;
  float alpha;
  float normr2old;
  float denom;
  float beta;


  // calculate residual
  cblas_sgemv(CblasColMajor, CblasNoTrans, rows, rows,
    1.0f, A, rows, x, 1, 0.0f, tmp, 1);

  // first copy b into r
  // then compute r = -tmp + r
  cblas_scopy(rows, b, 1, r, 1);
  cblas_saxpy(rows, -1.0f, tmp, 1, r, 1);

  // calculate norm squared
  normr2 = cblas_sdot(rows, r, 1, r, 1);

  // copy r to p
  cblas_scopy(rows, r, 1, p, 1);

  // main loop
  while (normr2 > tol && iter < maxiters) {
    // calculate denom
    cblas_sgemv(CblasColMajor, CblasNoTrans, rows, rows,
      1.0f, A, rows, p, 1, 0.0f, tmp, 1);

    denom = cblas_sdot(rows, r, 1, tmp, 1);

    alpha = normr2 / denom;

    // update x
    cblas_saxpy(rows, alpha, p, 1, x, 1);

    // update normr2
    normr2old = normr2;
    cblas_saxpy(rows, -alpha, tmp, 1, r, 1);

    normr2 = cblas_sdot(rows, r, 1, r, 1);

    // update p
    beta = normr2 / normr2old;
    cblas_sscal(rows, beta, p, 1);
    cblas_saxpy(rows, 1.0f, r, 1, p, 1);

    iter++;
  }
}

int test_small() {

  float tol = 1e-5;
  int maxiters = 100;

  // generate test data
  float* mat;
  float* b;

  generate_test_data(&mat, &b);

  // create x and temporaries
  float x[3] = {0};
  float p[3] = {0};
  float r[3] = {0};
  float tmp[3] = {0};

  // call CG
  cg_mkl(3, mat, b, tol, maxiters,
       x,
       p, r, tmp);

  print_vec("x", x, 3);

  // check answer
  float normr2;
  cblas_sgemv(CblasColMajor, CblasNoTrans, 3, 3,
    1.0f, mat, 3, x, 1, 0.0f, tmp, 1);
  cblas_scopy(3, b, 1, r, 1);
  cblas_saxpy(3, -1.0f, tmp, 1, r, 1);
  normr2 = cblas_sdot(3, r, 1, r, 1);

  print_vec("answer", tmp, 3);
  print_vec("expected", b, 3);

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

  generate_random_test_data(&mat, &b, size);

  // create x and temporaries
  float x[size];
  float p[size];
  float r[size];
  float tmp[size];
  for (int i=0; i<size; i++)
    x[i] = p[i] = r[i] = tmp[i] = 0.0f;

  // call CG
  cg_mkl(size, mat, b, tol, maxiters,
       x,
       p, r, tmp);

  // check answer
  float normr2;
  cblas_sgemv(CblasColMajor, CblasNoTrans, size, size,
    1.0f, mat, size, x, 1, 0.0f, tmp, 1);
  cblas_scopy(size, b, 1, r, 1);
  cblas_saxpy(size, -1.0f, tmp, 1, r, 1);
  normr2 = cblas_sdot(size, r, 1, r, 1);

  print_vec("answer", tmp, size);
  print_vec("expected", b, size);

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
