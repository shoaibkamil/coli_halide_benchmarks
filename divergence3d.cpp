#include "Halide.h"


using namespace Halide;

int main(int argc, char* argv[]) {
  ImageParam in1{Float(32), 3, "input1"};
  ImageParam in2{Float(32), 3, "input2"};
  ImageParam in3{Float(32), 3, "input3"};
  Param<float> alpha;
  Param<float> beta;
  Param<float> gamma;

  Func divergence3d{"divergence3d"};
  Var x, y, z;

  divergence3d(x, y, z) = alpha * (in1(x+1, y, z) + in1(x-1, y, z)) +
                       beta  * (in2(x, y+1, z) + in2(x, y-1, z)) +
                       gamma  * (in3(x, y, z+1) + in3(x, y, z-1));

  divergence3d.parallel(y).vectorize(x, 8);

  divergence3d.compile_to_lowered_stmt("divergence3d.html", {in1, in2, in3, alpha, beta, gamma}, HTML);


  return 0;
}
