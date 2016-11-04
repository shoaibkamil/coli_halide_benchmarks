#include "Halide.h"

#define RADIUS 3

using namespace Halide;

int main(int argc, char* argv[]) {
  ImageParam in{Float(32), 2, "input"};
  // kernel is 2*radius x 2*radius
  ImageParam kernel{Float(32), 2, "kernel"};

  Func filter2D_nordom{"filter2D_nordom"};
  Var x, y;

  for (int i=-RADIUS; i<RADIUS; i++) {
    for (int j=-RADIUS; j<RADIUS; j++)  {
      e += in(x+i, y+j) * kernel(RADIUS+i, RADIUS+j);
    }
  }

  filter2D_nordom(x, y) = e;

  filter2D_nordom.parallel(y).vectorize(x, 8);

  filter2D_nordom.compile_to_lowered_stmt("filter2D_nordom.html", {in, kernel}, HTML);


  return 0;
}
