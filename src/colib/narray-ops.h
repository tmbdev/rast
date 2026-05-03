// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

/// \file narray-ops.h
/// \brief Mathematical operators on arrays

#ifndef h_narray_ops_
#define h_narray_ops_

#include <math.h>
#include <stdlib.h>
#include "colib/narray.h"

namespace narray_ops {

    namespace {
        template <class T>
        inline T abs_(T x) { return x>0?x:-x; }
        template <class T,class S>
        inline T max_(T x,S y) { return x>y?x:y; }
        template <class T,class S>
        inline T min_(T x,S y) { return x<y?x:y; }
    }

    /// min/max

    template <class T,class S>
    T amin(colib::narray<T> &out) {
        T result = out.at1d(0);
        for(int i=1;i<out.length1d();i++)
            result = min_(result,out.at1d(i));
        return result;
    }

    template <class T,class S>
    T amax(colib::narray<T> &out) {
        T result = out.at1d(0);
        for(int i=1;i<out.length1d();i++)
            result = max_(result,out.at1d(i));
        return result;
    }

    template <class T,class S>
    void max(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = max_(out.at1d(i),in);
    }

    template <class T,class S>
    void min(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = min_(out.at1d(i),in);
    }

    template <class T,class S>
    void max(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = max_(out.at1d(i),in.at1d(i));
    }

    template <class T,class S>
    void min(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = min_(out.at1d(i),in.at1d(i));
    }

    template <class T,class S,class R>
    void greater(colib::narray<T> &out,S in,R no,R yes) {
        for(int i=0;i<out.length1d();i++)
            if(out.at1d(i)>in)
                out.at1d(i) = yes;
            else
                out.at1d(i) = no;
    }

    template <class T,class S,class R>
    void less(colib::narray<T> &out,S in,R no,R yes) {
        for(int i=0;i<out.length1d();i++)
            if(out.at1d(i)<in)
                out.at1d(i) = yes;
            else
                out.at1d(i) = no;
    }

    // arithmetic operations

    // add

    template <class T,class S>
    void add(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) += in;
    }

    template <class T,class S>
    void operator+=(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) += in;
    }

    template <class T,class S>
    void add(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) += in.at1d(i);
    }

    template <class T,class S>
    void operator+=(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) += in.at1d(i);
    }

    template <class T,class S>
    void add(colib::narray<T> &out,colib::narray<S> &in1,colib::narray<S> &in2) {
        CHECK_ARG(samedims(in1,in2));
        makelike(out,in1);
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = in1.at1d(i) + in2.at1d(i);
    }

    // sub

    template <class T,class S>
    void sub(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) -= in;
    }

    template <class T,class S>
    void operator-=(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) -= in;
    }

    template <class T,class S>
    void sub(S in,colib::narray<T> &out) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = in - out.at1d(i);
    }

    template <class T,class S>
    void sub(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) -= in.at1d(i);
    }

    template <class T,class S>
    void operator-=(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) -= in.at1d(i);
    }

    template <class T,class S>
    void sub(colib::narray<T> &out,colib::narray<S> &in1,colib::narray<S> &in2) {
        CHECK_ARG(samedims(in1,in2));
        makelike(out,in1);
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = in1.at1d(i) - in2.at1d(i);
    }

    // mul

    template <class T,class S>
    void mul(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) *= in;
    }

    template <class T,class S>
    void operator*=(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) *= in;
    }

    template <class T,class S>
    void mul(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) *= in.at1d(i);
    }

    template <class T,class S>
    void operator*=(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) *= in.at1d(i);
    }

    template <class T,class S>
    void mul(colib::narray<T> &out,colib::narray<S> &in1,colib::narray<S> &in2) {
        CHECK_ARG(samedims(in1,in2));
        makelike(out,in1);
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = in1.at1d(i) * in2.at1d(i);
    }

    // div

    template <class T,class S>
    void div(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) /= in;
    }

    template <class T,class S>
    void div(S value,colib::narray<T> &out) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = value / out.at1d(i);
    }

    template <class T,class S>
    void operator/=(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) /= in;
    }

    template <class T,class S>
    void div(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) /= in.at1d(i);
    }

    template <class T,class S>
    void operator/=(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) /= in.at1d(i);
    }

    template <class T,class S>
    void div(colib::narray<T> &out,colib::narray<S> &in1,colib::narray<S> &in2) {
        CHECK_ARG(samedims(in1,in2));
        makelike(out,in1);
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = in1.at1d(i) / in2.at1d(i);
    }

    // mod

    template <class T,class S>
    void mod(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) %= in;
    }

    template <class T,class S>
    void mod(S value,colib::narray<T> &out) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = value % out.at1d(i);
    }

    template <class T,class S>
    void operator%=(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) %= in;
    }

    template <class T,class S>
    void mod(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) %= in.at1d(i);
    }

    template <class T,class S>
    void mod(colib::narray<T> &out,colib::narray<S> &in1,colib::narray<S> &in2) {
        CHECK_ARG(samedims(in1,in2));
        makelike(out,in1);
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = in1.at1d(i) % in2.at1d(i);
    }

    template <class T,class S>
    void operator%=(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) %= in.at1d(i);
    }

    // pow

    template <class T,class S>
    void pow(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = ::pow(out.at1d(i),in);
    }

    template <class T,class S>
    void pow(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = ::pow(out.at1d(i),in.at1d(i));
    }

    template <class T,class S>
    void pow(colib::narray<T> &out,colib::narray<S> &in1,colib::narray<S> &in2) {
        CHECK_ARG(samedims(in1,in2));
        makelike(out,in1);
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = ::pow(in1.at1d(i),in2.at1d(i));
    }

    ////////////////////////////////////////////////////////////////
    /// bit operations
    ////////////////////////////////////////////////////////////////

    // and

    template <class T,class S>
    void  bits_and(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) &= in.at1d(i);
    }

    template <class T,class S>
    void operator&=(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) &= in.at1d(i);
    }

    template <class T,class S>
    void  bits_and(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) &= in;
    }

    template <class T,class S>
    void operator&=(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) &= in;
    }

    // or

    template <class T,class S>
    void  bits_or(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) |= in.at1d(i);
    }

    template <class T,class S>
    void operator|=(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) |= in.at1d(i);
    }

    template <class T,class S>
    void  bits_or(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) |= in;
    }

    template <class T,class S>
    void operator|=(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) |= in;
    }

    // xor

    template <class T,class S>
    void  bits_xor(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) ^= in.at1d(i);
    }

    template <class T,class S>
    void operator^=(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) &= in.at1d(i);
    }

    template <class T,class S>
    void  bits_xor(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) ^= in;
    }

    template <class T,class S>
    void operator^=(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) &= in;
    }

    // lshift

    template <class T,class S>
    void  bits_lshift(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) <<= in.at1d(i);
    }

    template <class T,class S>
    void operator<<=(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) <<= in.at1d(i);
    }

    template <class T,class S>
    void  bits_lshift(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) <<= in;
    }

    template <class T,class S>
    void operator<<=(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) <<= in;
    }

    // rshift

    template <class T,class S>
    void  bits_rshift(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) >>= in.at1d(i);
    }

    template <class T,class S>
    void operator>>=(colib::narray<T> &out,colib::narray<S> &in) {
        CHECK_ARG(samedims(out,in));
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) >>= in.at1d(i);
    }

    template <class T,class S>
    void  bits_rshift(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) >>= in;
    }

    template <class T,class S>
    void operator>>=(colib::narray<T> &out,S in) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) >>= in;
    }

    /// misc arithmetic operations

    template <class T>
    void neg(colib::narray<T> &out) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = -out.at1d(i);
    }

    template <class T>
    void invert(colib::narray<T> &out) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = ~out.at1d(i);
    }

    template <class T>
    void abs(colib::narray<T> &out) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = abs_(out.at1d(i));
    }

    template <class T>
    void log(colib::narray<T> &out) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = ::log(out.at1d(i));
    }

    template <class T>
    void exp(colib::narray<T> &out) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = ::exp(out.at1d(i));
    }

    template <class T>
    void sin(colib::narray<T> &out) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = ::sin(out.at1d(i));
    }

    template <class T>
    void cos(colib::narray<T> &out) {
        for(int i=0;i<out.length1d();i++)
            out.at1d(i) = ::cos(out.at1d(i));
    }
}

#endif
