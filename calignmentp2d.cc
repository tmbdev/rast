/* Copyright (c) 1990-1995 by Thomas M. Breuel */

////////////////////////////////////////////////////////////////
// TODO:
// - increase default allocation size for narray
////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fpu_control.h>

#include "misc.h"
#include "narray.h"
#include "vec2.h"
#include "smartptr.h"
#include "narray-util.h"
#include "trie.h"
using namespace colib;

#include "rast.h"
#include "util.h"

namespace lumo_calignmentp2d {

bool verbose = igetenv("verbose", 1);

typedef vec2 Ipoint;
typedef vec2 Mpoint;
typedef vec2 Msource;

inline int urand48() { return abs(int(lrand48())); }

template <class T>
void shuffle(narray<T> &narray) {
  int n = narray.length();
  for (int i = 0; i < n - 1; i++) {
    int j = urand48() % (n - i) + i;
    if (i != j) swap(narray[i], narray[j]);
  }
}

static bool use_trie = true;

struct CAlignmentP2D : AlignmentP2D {
  float eps;
  narray<vec2> image;
  narray<vec2> model;
  Trie2<int> itrie;
  enum { offset = 1000 };
  float minscale;
  float maxscale;
  void computeImagePointTable() {
    if (!use_trie) return;
    itrie.init(eps, 2000, 2000, -1000, -1000);
    for (int i = 0; i < image.length(); i++) {
      vec2 p = image.at(i);
      itrie.add(p[0], p[1], i);
    }
  }

  vec2 transl;
  vec2 rotation;

  float evaluateAlignment(int i0, int i1, int m0, int m1) {
    vec2 pi0 = image.at(i0);
    vec2 pi1 = image.at(i1);
    vec2 pm0 = model.at(m0);
    vec2 pm1 = model.at(m1);
    rotation = cdiv(pi1 - pi0, pm1 - pm0);
    float scale = rotation.magnitude();
    if (scale < minscale || scale > maxscale) return 0.0;
    transl = pi0 - cmul(rotation, pm0);
    ASSERT(distance(cmul(rotation, pm0) + transl, pi0) < 0.01);
    return evaluateAlignment0();
  }

  float evaluateAlignment0() {
    float total = 0;
    for (int m = 0; m < model.length(); m++) {
      vec2 p = model.at(m);
      vec2 tp = cmul(rotation, p) + transl;
      int mi = -1;
      if (use_trie) {
        narray<int> candidates;
        itrie.query(candidates, tp[0] - eps, tp[1] - eps, tp[0] + eps,
                    tp[1] + eps);
        for (int ii = 0; ii < candidates.length(); ii++) {
          int i = candidates[ii];
          if (distance(image.at(i), tp) < eps) {
            mi = i;
            total++;
            break;
          }
        }
      } else {
        for (int i = 0; i < image.length(); i++) {
          if (distance(image.at(i), tp) < eps) {
            mi = i;
            total++;
            break;
          }
        }
      }
    }
    return total;
  }

  struct Solution {
    int i0, i1, m0, m1;
    float quality;
    vec2 translation;
    vec2 rotation;
  } solution;

  void searchForBestAlignment() {
    solution.quality = 0;
    for (int i0 = 0; i0 < image.length(); i0++) {
      for (int i1 = i0 + 1; i1 < image.length(); i1++) {
        for (int m0 = 0; m0 < model.length(); m0++) {
          for (int m1 = 0; m1 < model.length(); m1++) {
            if (m0 == m1) continue;
            float quality = evaluateAlignment(i0, i1, m0, m1);
            if (quality > solution.quality) {
              solution.i0 = i0;
              solution.i1 = i1;
              solution.m0 = m0;
              solution.m1 = m1;
              solution.quality = quality;
              solution.translation = transl;
              solution.rotation = rotation;
            }
          }
        }
      }
    }
  }

  void compute() {
    computeImagePointTable();
    searchForBestAlignment();
  }
  void clear_ipoints() { image.clear(); }
  void add_ipoint(float x, float y) { image.push(vec2(x, y)); }
  void clear_mpoints() { model.clear(); }
  void add_mpoint(float x, float y) { model.push(vec2(x, y)); }
  void set_srange(float min, float max) {
    minscale = min;
    maxscale = max;
  }
  void set_epsilon(float e) { eps = e; }

  float quality() { return solution.quality; }
  float translation(int dim) { return solution.translation[dim]; }
  float angle() { return normangleOf(angleOf(solution.rotation)); }
  float scale() { return norm(solution.rotation); }
};
}

AlignmentP2D *makeAlignmentP2D() {
  return new lumo_calignmentp2d::CAlignmentP2D();
}
