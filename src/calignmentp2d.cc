// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "rast.h"
#include "trie.h"
#include "util.h"
#include "vec2.h"

namespace lumo_calignmentp2d {

// Diagnostic flag: parsed once at startup from the `verbose` env-var; read-only
// thereafter. Falls under the "Debug and diagnostic settings" carve-out from
// the project's globals rule.
bool verbose = igetenv("verbose", 1);

inline int urand48() { return std::abs(int(lrand48())); }

static bool use_trie = true;  // build-time switch; read-only at runtime

struct CAlignmentP2D : AlignmentP2D {
  float eps{0.0f};
  std::vector<vec2> image;
  std::vector<vec2> model;
  Trie2<int> itrie;
  float minscale{1.0f};
  float maxscale{1.0f};

  vec2 transl;
  vec2 rotation;

  void computeImagePointTable() {
    if (!use_trie) return;
    itrie.init(eps, 2000, 2000, -1000, -1000);
    for (std::size_t i = 0; i < image.size(); i++) {
      const vec2 &p = image[i];
      itrie.add(p[0], p[1], int(i));
    }
  }

  float evaluateAlignment(int i0, int i1, int m0, int m1) {
    vec2 pi0 = image[i0];
    vec2 pi1 = image[i1];
    vec2 pm0 = model[m0];
    vec2 pm1 = model[m1];
    rotation = cdiv(pi1 - pi0, pm1 - pm0);
    float scale = rotation.magnitude();
    if (scale < minscale || scale > maxscale) return 0.0;
    transl = pi0 - cmul(rotation, pm0);
    assert(distance(cmul(rotation, pm0) + transl, pi0) < 0.01);
    return evaluateAlignment0();
  }

  float evaluateAlignment0() {
    float total = 0;
    for (std::size_t m = 0; m < model.size(); m++) {
      vec2 p = model[m];
      vec2 tp = cmul(rotation, p) + transl;
      if (use_trie) {
        std::vector<int> candidates;
        itrie.query(candidates, tp[0] - eps, tp[1] - eps, tp[0] + eps, tp[1] + eps);
        for (int i : candidates) {
          if (distance(image[i], tp) < eps) {
            total++;
            break;
          }
        }
      } else {
        for (std::size_t i = 0; i < image.size(); i++) {
          if (distance(image[i], tp) < eps) {
            total++;
            break;
          }
        }
      }
    }
    return total;
  }

  struct Solution {
    int i0{0}, i1{0}, m0{0}, m1{0};
    float quality{0.0f};
    vec2 translation;
    vec2 rotation;
  } solution;

  void searchForBestAlignment() {
    solution.quality = 0;
    const int ni = int(image.size()), nm = int(model.size());
    for (int i0 = 0; i0 < ni; i0++) {
      for (int i1 = i0 + 1; i1 < ni; i1++) {
        for (int m0 = 0; m0 < nm; m0++) {
          for (int m1 = 0; m1 < nm; m1++) {
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

  void compute() override {
    computeImagePointTable();
    searchForBestAlignment();
  }
  void clear_ipoints() override { image.clear(); }
  void add_ipoint(float x, float y) override { image.push_back(vec2(x, y)); }
  void clear_mpoints() override { model.clear(); }
  void add_mpoint(float x, float y) override { model.push_back(vec2(x, y)); }
  void set_srange(float min, float max) override {
    minscale = min;
    maxscale = max;
  }
  void set_epsilon(float e) override { eps = e; }

  float quality() override { return solution.quality; }
  float translation(int dim) override { return solution.translation[dim]; }
  float angle() override { return normangleOf(angleOf(solution.rotation)); }
  float scale() override { return norm(solution.rotation); }
};

}  // namespace lumo_calignmentp2d

std::unique_ptr<AlignmentP2D> makeAlignmentP2D() {
  return std::make_unique<lumo_calignmentp2d::CAlignmentP2D>();
}
