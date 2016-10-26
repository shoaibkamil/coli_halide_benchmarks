HALIDE_PATH = ../Halide
HALIDE_LIB = -L${HALIDE_PATH}/bin -lHalide
HALIDE_INC = -I${HALIDE_PATH}/include

.cpp:
	c++ -std=c++11 ${HALIDE_INC} -o $@ $< ${HALIDE_LIB}
