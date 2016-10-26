#include "Halide.h"


using namespace Halide;


int main(int argc, char* argv[]) {
  ImageParam in{UInt(8), 2, "input"};

  Func histogram{"histogram"};
  Var x,y;

  RDom r(in);

  histogram(x) = 0;
  histogram(clamp(in(r.x, r.y), 0, 255)) += 1;

  // TODO: use an rfactor schedule

  histogram.compile_to_lowered_stmt("histogram.html", {in}, HTML);


  return 0;
}

