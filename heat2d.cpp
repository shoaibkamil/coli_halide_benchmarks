#include "Halide.h"

using namespace Halide;

int main(int argc, char* argv[]) {
  ImageParam in{Float(32), 2, "input"};
  Param<float> alpha;
  Param<float> beta;

  Func heat2d{"heat2d"};
  Var x, y;

  heat2d(x, y) = alpha * in(x, y) +
                 beta * (in(x+1, y) + in(x-1, y)+
                         in(x, y+1) + in(x, y-1));

  heat2d.parallel(y).vectorize(x, 8);

  heat2d.compile_to_lowered_stmt("heat2d.html", {in, alpha, beta}, HTML);


  return 0;
}
