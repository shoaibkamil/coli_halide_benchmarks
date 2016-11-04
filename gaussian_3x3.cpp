#include "Halide.h"

using namespace Halide;

#define kernelX_length 7
#define kernelY_length 7

int main(int argc, char* argv[]) {
  ImageParam in{Float(32), 2, "input"};
  ImageParam kernelX{Float(32), 1, "kernelx"};
  ImageParam kernelY{Float(32), 1, "kernely"};

  

  Func gaussian{"gaussian"};
  Func gaussian_x{"gaussian_x"};
  Var x,y;

//  gaussian_x(x, y) = sum(in(clamp(x+krx.x-kernelX_length/2, 0, in.width()-1), y) * kernelX(krx.x));
//  gaussian(x, y) = sum(gaussian_x(x, clamp(y+kry.x-kernelY_length/2, 0, in.height()-1)) * kernelY(kry.x));

  Expr e,f;
  e = 0.0f;
  for (int i=0; i<kernelX_length; i++) {
    e += in(x+i-kernelX_length/2, y) * kernelX(i);
  }
  gaussian_x(x, y) = e;

  f = 0.0f;
  for (int i=0; i<kernelX_length; i++) {
    f += gaussian_x(x, y+i-kernelX_length/2) * kernelY(i);
  }

  gaussian(x, y) = f;

  Var x_inner, y_inner, x_outer, y_outer, tile_index;
  gaussian.tile(x, y, x_outer, y_outer, x_inner, y_inner, 4, 4)
          .fuse(x_outer, y_outer, tile_index).compute_root().parallel(tile_index);
  gaussian_x.compute_at(gaussian, y_inner);

  gaussian.compile_to_lowered_stmt("gaussian_3x3.html",
    {in, kernelX, kernelY}, HTML);


  return 0;
}

