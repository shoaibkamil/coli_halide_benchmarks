#include <vector>
#include "Halide.h"

using namespace Halide;

// Generator class for vectory copy
template<class T>
class CopyGenerator :
        public Generator<CopyGenerator<T>> {
  public:
    typedef Generator<CopyGenerator<T>> Base;
    using Base::target;
    using Base::get_target;
    using Base::natural_vector_size;
    template<typename T2> using Input = typename Base::template Input<T2>;
    template<typename T2> using Output = typename Base::template Output<T2>;

    GeneratorParam<bool> vectorize_ = {"vectorize", true};
    GeneratorParam<int>  block_size_ = {"block_size", 1024};

    // Standard ordering of parameters in AXPY functions.
    Input<Buffer<T>> x_ = {"x", 1};

    Output<Buffer<T>> result_ = {"result", 1};

    void Schedule(Func result, Expr width) {
        Var i("i"), o("o");
    }

    template<class Arg>
    Expr calc(Arg i) {
        return x_(i);
    }

    void generate() {
        assert(get_target().has_feature(Target::NoAsserts));
        assert(get_target().has_feature(Target::NoBoundsQuery));

        const int vec_size = vectorize_? natural_vector_size(type_of<T>()): 1;
        Expr size = x_.width();
        Expr size_vecs = (size / vec_size) * vec_size;
        Expr size_tail = size - size_vecs;

        Var  i("i");
        RDom vecs(0, size_vecs, "vec");
        RDom tail(size_vecs, size_tail, "tail");
        result_(i) = undef(type_of<T>());
        result_(vecs) = calc(vecs);
        result_(tail) = calc(tail);

        if (vectorize_) {
            Var ii("ii");
            result_.update().vectorize(vecs, vec_size);
        }

        result_.bound(i, 0, x_.width());
        result_.dim(0).set_bounds(0, x_.width());

        x_.dim(0).set_min(0);
    }
};

// Generator class for vector subtract
template<class T>
class SubtractGenerator :
        public Generator<SubtractGenerator<T>> {
  public:
    typedef Generator<SubtractGenerator<T>> Base;
    using Base::target;
    using Base::get_target;
    using Base::natural_vector_size;
    template<typename T2> using Input = typename Base::template Input<T2>;
    template<typename T2> using Output = typename Base::template Output<T2>;

    GeneratorParam<bool> vectorize_ = {"vectorize", true};
    GeneratorParam<int>  block_size_ = {"block_size", 1024};
    GeneratorParam<bool> add_to_y_ = {"add_to_y", true};

    // Standard ordering of parameters in AXPY functions.
    Input<Buffer<T>> x_ = {"x", 1};
    Input<Buffer<T>> y_ = {"y", 1};

    Output<Buffer<T>> result_ = {"result", 1};

    void Schedule(Func result, Expr width) {
        Var i("i"), o("o");
    }

    template<class Arg>
    Expr calc(Arg i) {
       return  x_(i) - y_(i);
    }

    void generate() {
        assert(get_target().has_feature(Target::NoAsserts));
        assert(get_target().has_feature(Target::NoBoundsQuery));

        const int vec_size = vectorize_? natural_vector_size(type_of<T>()): 1;
        Expr size = x_.width();
        Expr size_vecs = (size / vec_size) * vec_size;
        Expr size_tail = size - size_vecs;

        Var  i("i");
        RDom vecs(0, size_vecs, "vec");
        RDom tail(size_vecs, size_tail, "tail");
        result_(i) = undef(type_of<T>());
        result_(vecs) = calc(vecs);
        result_(tail) = calc(tail);

        if (vectorize_) {
            Var ii("ii");
            result_.update().vectorize(vecs, vec_size);
        }

        result_.bound(i, 0, x_.width());
        result_.dim(0).set_bounds(0, x_.width());

        x_.dim(0).set_min(0);
        y_.dim(0).set_bounds(0, x_.width());
    }
};

HALIDE_REGISTER_GENERATOR(CopyGenerator<float>, sveccopy);
HALIDE_REGISTER_GENERATOR(SubtractGenerator<float>, svecsub);
