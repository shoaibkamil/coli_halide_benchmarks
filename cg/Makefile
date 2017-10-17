
HALIDE_ROOT=/Users/kamil/documents/work/adobe/halide
HALIDE_INCLUDE_DIR=${HALIDE_ROOT}/include
HALIDE_LIBRARY_DIR=${HALIDE_ROOT}/lib

OBJS = saxpy.o sdot.o sgemv.o sveccopy.o svecsub.o runtime.a
HEADERS = saxpy.h sdot.h sgemv.h sveccopy.h svecsub.h

cg: objs
	c++ -fno-rtti -std=c++11 cg.cpp ${OBJS} -I ${HALIDE_INCLUDE_DIR} -o cg \
		-L ${HALIDE_LIBRARY_DIR} -lHalide -lcurses -lz


objs: ${OBJS}

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

saxpy.o: blas_l1_generators
	./blas_l1_generators -g saxpy -e h,o,stmt -o . target=host-no_asserts-no_bounds_query vectorize=true add_to_y=true scale_x=true

sdot.o: blas_l1_generators
	./blas_l1_generators -g sdot -e h,o,stmt -o . target=host-no_asserts-no_bounds_query vectorize=true

sgemv.o: blas_l2_generators
	./blas_l2_generators -g sgemv -e h,o,stmt -o . target=host-no_asserts-no_bounds_query vectorize=true parallel=true

sveccopy.o: vecutil_generators
	./vecutil_generators -g sveccopy -e h,o,stmt -o . target=host-no_asserts-no_bounds_query vectorize=true

svecsub.o: vecutil_generators
	./vecutil_generators -g svecsub -e h,o,stmt -o . target=host-no_asserts-no_bounds_query vectorize=true

runtime.a: vecutil_generators
	./vecutil_generators -o . -r runtime target=host-no_asserts-no_bounds_query

clean:
	rm -f blas_l1_generators blas_l2_generators vecutil_generators ${OBJS} *.stmt cg ${HEADERS}
