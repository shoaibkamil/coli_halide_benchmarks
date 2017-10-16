
HALIDE_ROOT=/Users/kamil/documents/work/adobe/halide
HALIDE_INCLUDE_DIR=${HALIDE_ROOT}/include
HALIDE_LIBRARY_DIR=${HALIDE_ROOT}/lib

all: generators

generators: blas_l1_generators blas_l2_generators vecutil_generators

blas_l1_generators: blas_l1_generators.cpp
	c++ -fno-rtti -std=c++11 ./blas_l1_generators.cpp ${HALIDE_ROOT}/tools/GenGen.cpp \
		-o blas_l1_generators \
		-I ${HALIDE_INCLUDE_DIR} -L ${HALIDE_LIBRARY_DIR} -lHalide -lcurses -lz

blas_l2_generators: blas_l2_generators.cpp
	c++ -fno-rtti -std=c++11 ./blas_l2_generators.cpp ${HALIDE_ROOT}/tools/GenGen.cpp \
		-o blas_l2_generators \
		-I ${HALIDE_INCLUDE_DIR} -L ${HALIDE_LIBRARY_DIR} -lHalide -lcurses -lz

vecutil_generators: vecutil_generators.cpp
	c++ -fno-rtti -std=c++11 ./vecutil_generators.cpp ${HALIDE_ROOT}/tools/GenGen.cpp \
		-o vecutil_generators \
		-I ${HALIDE_INCLUDE_DIR} -L ${HALIDE_LIBRARY_DIR} -lHalide -lcurses -lz

clean:
	rm -f blas_l1_generators blas_l2_generators vecutil_generators
