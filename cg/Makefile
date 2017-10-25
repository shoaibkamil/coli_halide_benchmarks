
HALIDE_ROOT=/Users/skamil/documents/work/adobe/halide
HALIDE_INCLUDE_DIR=${HALIDE_ROOT}/include
HALIDE_LIBRARY_DIR=${HALIDE_ROOT}/lib

MKL_COMPILE_OPTS=-DMKL_ILP64 -m64 -I${MKLROOT}/include -I${HALIDE_INCLUDE_DIR}
MKL_LINK_OPTS=-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

OBJS = saxpy.o sdot.o sgemv.o sveccopy.o svecsub.o runtime.a
HEADERS = saxpy.h sdot.h sgemv.h sveccopy.h svecsub.h

cg: objs common.h cg.cpp
	c++ -g -fno-rtti -std=c++11 cg.cpp ${OBJS} -I ${HALIDE_INCLUDE_DIR} -o cg \
		-L ${HALIDE_LIBRARY_DIR} -lHalide -lcurses -lz

cg_mkl: cg_mkl.cpp common.h
	c++ -g -std=c++11 cg_mkl.cpp ${MKL_COMPILE_OPTS} -o cg_mkl ${MKL_LINK_OPTS}

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
	rm -rf blas_l1_generators blas_l2_generators vecutil_generators ${OBJS} *.stmt cg ${HEADERS} *dSYM
