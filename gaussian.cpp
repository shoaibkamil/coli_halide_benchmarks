#include "Halide.h"

using namespace Halide;

int main(int argc, char* argv[]) {
  ImageParam in{Float(32), 2, "input"};
  ImageParam kernelX{Float(32), 1, "kernelx"};
  ImageParam kernelY{Float(32), 1, "kernely"};
  Param<uint8_t> kernelX_length{"kernelX_length"};
  Param<uint8_t> kernelY_length{"kernelY_length"};

  Func gaussian{"gaussian"};
  Func gaussian_x{"gaussian_x"};
  Var x,y;
  RDom krx{0, kernelX_length};
  RDom kry{0, kernelY_length};

  gaussian_x(x, y) = sum(in(clamp(x+krx.x-kernelX_length/2, 0, in.width()-1), y) * kernelX(krx.x));
  gaussian(x, y) = sum(gaussian_x(x, clamp(y+kry.x-kernelY_length/2, 0, in.height()-1)) * kernelY(kry.x));

  Var x_inner, y_inner, x_outer, y_outer, tile_index;
  gaussian.tile(x, y, x_outer, y_outer, x_inner, y_inner, 4, 4)
          .fuse(x_outer, y_outer, tile_index).compute_root().parallel(tile_index);
  gaussian_x.compute_at(gaussian, y_inner);

  gaussian.compile_to_lowered_stmt("gaussian.html",
    {in, kernelX, kernelY, kernelX_length, kernelY_length}, HTML);


  return 0;
}

