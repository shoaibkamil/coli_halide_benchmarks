HALIDE_PATH = ../Halide
HALIDE_LIB = -L${HALIDE_PATH}/bin -lHalide
HALIDE_INC = -I${HALIDE_PATH}/include

all: cvt_color filter2D histogram gaussian

cvt_color: cvt_color.cpp
	c++ -std=c++11 ${HALIDE_INC} -o $@ $< ${HALIDE_LIB}

filter2D: filter2D.cpp
	c++ -std=c++11 ${HALIDE_INC} -o $@ $< ${HALIDE_LIB}

histogram: histogram.cpp
	c++ -std=c++11 ${HALIDE_INC} -o $@ $< ${HALIDE_LIB}

gaussian: gaussian.cpp
	c++ -std=c++11 ${HALIDE_INC} -o $@ $< ${HALIDE_LIB}

clean:
	rm -f cvt_color filter2D histogram gaussian *.html
