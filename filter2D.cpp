#include "Halide.h"

using namespace Halide;

int main(int argc, char* argv[]) {
  ImageParam in{Float(32), 2, "input"};
  ImageParam kernel{Float(32), 2, "kernel"};

  Param<uint8_t> kernel_rows{"kernel_rows"};
  Param<uint8_t> kernel_cols{"kernel_cols"};

  Func filter2D{"filter2D"};
  Var x, y;
  RDom r(0, kernel_rows, 0, kernel_cols);

  Expr row = clamp(x + r.x - kernel_rows/2, 0, in.width());
  Expr col = clamp(y + r.y - kernel_cols/2, 0, in.height());

  filter2D(x, y) = sum(in(col, row) * kernel(r.x, r.y));

  filter2D.parallel(y).vectorize(x, 8);

  filter2D.compile_to_lowered_stmt("filter2D.html", {in, kernel, kernel_rows, kernel_cols}, HTML);


  return 0;
}
