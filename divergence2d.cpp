#include "Halide.h"


using namespace Halide;

int main(int argc, char* argv[]) {
  ImageParam in{Float(32), 2, "input"};
  Param<float> alpha;
  Param<float> beta;

  Func divergence2d{"divergence2d"};
  Var x, y;

  divergence2d(x, y) = alpha * (in(x+1, y) + in(x-1, y)) +
                       beta  * (in(x, y+1) + in(x, y-1));

  divergence2d.parallel(y).vectorize(x, 8);

  divergence2d.compile_to_lowered_stmt("divergence2d.html", {in, alpha, beta}, HTML);


  return 0;
}
