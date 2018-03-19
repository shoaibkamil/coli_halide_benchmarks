#include "Halide.h"

using namespace Halide;

int main(int argc, char* argv[]) {
  ImageParam in{Float(32), 3, "input"};
  Param<float> alpha;
  Param<float> beta;
  Param<int> t;

  Func heat3d{"heat3d"};
  Var x, y, z;
  Var k;

  RDom st{1, in.dim(0).extent()-1, 1, in.dim(1).extent()-1, 1, in.dim(2).extent()-1, 1, t, "spacetime"};

  heat3d(x, y, z, k) = undef<float>();
  heat3d(st.x, st.y, st.z, 0) = alpha * in(st.x, st.y, st.z) +
                 beta * (in(st.x+1, st.y, st.z) + in(st.x-1, st.y, st.z)+
                         in(st.x, st.y+1, st.z) + in(st.x, st.y-1, st.z) +
                         in(st.x, st.y, st.z+1) + in(st.x, st.y, st.z-1));


  heat3d(st.x, st.y, st.z, (st.w+1)%2) = alpha * heat3d(st.x, st.y, st.z, st.w%2) +
                  beta * (heat3d(st.x+1, st.y, st.z, st.w%2) + heat3d(st.x-1, st.y, st.z, st.w%2) +
                          heat3d(st.x, st.y+1, st.z, st.w%2) + heat3d(st.x, st.y-1, st.z, st.w%2) +
                          heat3d(st.x, st.y, st.z+1, st.w%2) + heat3d(st.x, st.y, st.z-1, st.w%2));


  //heat3d.parallel(z).vectorize(x, 8);
  heat3d.update(1).allow_race_conditions().parallel(st.z).vectorize(st.x, 8);

  heat3d.compile_to_lowered_stmt("heat3d.html", {in, alpha, beta, t}, HTML);


  return 0;
}
