#include "Halide.h"

using namespace Halide;

int main(int argc, char* argv[]) {
  ImageParam in{Float(32), 1, "input"};

  Func slide1d{"slide1d"};
  Var x;
  RDom r{2, in.width()-1};

  slide1d(x) = undef<float>();

  // slide1d(x) = in(x) + in(x-1) + in(x+1)
  // -->
  // slide1d(x) = ( in(x-1) + in(x-2) + in(x)) - in(x-2) + in(x+1)

  slide1d(1) = in(1) + in(0) + in(2);
  slide1d(r.x) = slide1d(r.x-1) - in(r.x-2) + in(r.x+1);



  slide1d.compile_to_lowered_stmt("slide1d.html", {in}, HTML);


  return 0;
}
