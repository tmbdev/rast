// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include "heap.h"
#include "rast.h"
#include "util.h"
#include "vec2.h"

namespace lumo_clinesp2d {

inline float min3(float a, float b, float c) { return std::min(a, std::min(b, c)); }

struct Point {
  vec2 p;
  float a;
  float w;
};

struct LineRegion {
  float th0, ux0, uy0;
  float th1, ux1, uy1;
  float thm, uxm, uym;
  float r0, r1, rm;
  float rerr, therr;

  void print(FILE *stream = stdout) {
    std::fprintf(stream, "<LineRegion %g %g   %g %g>", th0, th1, r0, r1);
  }

  void set(float th0_, float th1_, float r0_, float r1_) {
    if (th1_ <= th0_) throw "parameters (th)";
    if (r1_ <= r0_) throw "parameters (r)";
    this->th0 = th0_;
    this->th1 = th1_;
    this->r0 = r0_;
    this->r1 = r1_;
    ux0 = std::cos(th0_);
    uy0 = std::sin(th0_);
    ux1 = std::cos(th1_);
    uy1 = std::sin(th1_);
    thm = (th0_ + th1_) / 2;
    uxm = std::cos(thm);
    uym = std::sin(thm);
    double factor = std::cos((th1_ - th0_) / 2);
    this->rm = std::max(0.0, r0_ * factor);
    rerr = (r1_ - r0_) / 2.0f;
    therr = (th1_ - th0_) / 2.0f;
  }

  float dist(float x, float y) const {
    float dot0 = ux0 * x + uy0 * y;
    float dot1 = ux1 * x + uy1 * y;
    float d0 = dot0 - r1;
    float d1 = dot1 - r1;
    float dotm = uxm * x + uym * y;
    float d2 = -(dot0 - r0);
    float d3 = -(dotm - rm);
    float d4 = -(dot1 - r0);
    float upper = 0.0;
    if (d0 >= 0 && d1 >= 0) upper = std::min(d0, d1);
    float lower = 0.0;
    if (d2 >= 0 && d3 >= 0 && d4 >= 0) lower = min3(d2, d3, d4);
    return std::max(upper, lower);
  }

  float adist(float a, bool unoriented) const {
    float diff = a - thm;
    if (unoriented) {
      while (diff < -M_PI / 4) diff += M_PI / 2;
      while (diff > M_PI / 4) diff -= M_PI / 2;
      diff = std::fabs(diff);
      diff -= therr;
      if (diff < 0) return 0.0;
      return diff;
    } else {
      while (diff < -M_PI / 2) diff += M_PI;
      while (diff > M_PI / 2) diff -= M_PI;
      diff = std::fabs(diff);
      diff -= therr;
      if (diff < 0) return 0.0;
      return diff;
    }
  }

  void split(LineRegion *result, float thscale) const {
    float dr = r1 - r0;
    float da = thscale * (th1 - th0);
    if (da > dr) {
      float m = (th0 + th1) / 2;
      result[0].set(th0, m, r0, r1);
      result[1].set(m, th1, r0, r1);
    } else {
      float m = (r0 + r1) / 2;
      result[0].set(th0, th1, r0, m);
      result[1].set(th0, th1, m, r1);
    }
  }
};

using Matches = std::vector<int>;
using CMatches = std::shared_ptr<Matches>;

struct State {
  int generation{0};
  float weight{0.0f};
  LineRegion region;
  CMatches matches{std::make_shared<Matches>()};
};

using CState = std::shared_ptr<State>;

struct CLinesP2D : LinesP2D {
  float eps{2.0f}, aeps{0.05f};
  float tol{0.1f}, atol{0.001f};
  float minweight{0.0f};
  float maxoffset{3000.0f};
  int maxresults{1};
  int generation{0};
  int verbose{0};
  bool unoriented{true};
  bool lsq{true};

  std::vector<Point> points;
  std::vector<bool> used;

  Heap<CState> queue;
  std::vector<CState> results;

  void filter(CState &state) {
    LineRegion &region = state->region;
    Matches &input = *state->matches;
    auto cresult = std::make_shared<Matches>();
    Matches &result = *cresult;
    float weight = 0.0;
    float eps2 = eps * eps;
    float aeps2 = aeps * aeps;
    for (int index : input) {
      if (used[index]) continue;
      const Point &p = points[index];
      float q = 1.0;
      float da = region.adist(p.a, unoriented);
      if (lsq) {
        q *= std::max(0.0, 1.0 - da * da / aeps2);
      } else {
        if (da > aeps) q = 0.0;
      }
      if (q == 0.0) continue;
      float d = region.dist(p.p[0], p.p[1]);
      if (lsq) {
        q *= std::max(0.0, 1.0 - d * d / eps2);
      } else {
        if (d > eps) q = 0.0;
      }
      if (q == 0.0) continue;
      q *= p.w;
      weight += q;
      result.push_back(index);
    }
    state->weight = weight;
    state->matches = cresult;
  }

  void seedAll(std::vector<int> &indices) {
    indices.clear();
    indices.reserve(points.size());
    for (int j = 0; j < int(points.size()); j++) indices.push_back(j);
  }

  void compute(float a0, float a1, float d0, float d1) override {
    if (verbose) std::fprintf(stderr, "[#segments %zu]\n", points.size());
    generation = 0;
    used.assign(points.size(), false);

    auto initial = std::make_shared<State>();
    seedAll(*initial->matches);
    initial->region.set(a0, a1, d0, d1);
    filter(initial);
    queue.insert(initial, initial->weight);
    compute_();
  }

  void compute() override {
    if (verbose) std::fprintf(stderr, "[#segments %zu]\n", points.size());
    generation = 0;
    used.assign(points.size(), false);

    for (int i = 0; i < 8; i++) {
      auto initial = std::make_shared<State>();
      seedAll(*initial->matches);
      if (unoriented)
        initial->region.set(i * M_PI / 4, (i + 1) * M_PI / 4, 0.0, maxoffset);
      else
        initial->region.set(i * M_PI / 4, (i + 1) * M_PI / 4, -maxoffset, maxoffset);
      filter(initial);
      queue.insert(initial, initial->weight);
    }
    compute_();
  }

  void compute_() {
    for (int iter = 0;; iter++) {
      if (queue.empty()) break;
      CState state = queue.extractMax();
      if (state->generation != generation) {
        filter(state);
        state->generation = generation;
        queue.insert(state, state->weight);
        continue;
      }
      LineRegion &pregion = state->region;
      if (verbose > 1) {
        pregion.print(stderr);
        std::fprintf(stderr, "(%g)\n", state->weight);
      } else if (verbose > 0) {
        if (iter % 1000 == 0 && !queue.empty()) {
          std::fprintf(stderr, "[%d %zu %g %zu]\n", iter, queue.length(), queue.topPriority(),
                       results.size());
        }
      }

      if (state->weight < minweight) continue;

      if (pregion.rerr < tol && pregion.therr < atol) {
        results.push_back(state);
        Matches &matches = *state->matches;
        for (int idx : matches) {
          used[idx] = true;
          const Point &p = points[idx];
          if (verbose) std::printf("P %g %g %g\n", p.p[0], p.p[1], p.a);
        }
        generation++;
        if (int(results.size()) >= maxresults) break;
        continue;
      }

      LineRegion regions[2];
      state->region.split(regions, tol / atol);
      for (int i = 0; i < 2; i++) {
        auto child = std::make_shared<State>();
        child->matches = state->matches;
        child->weight = state->weight;
        child->region = regions[i];
        filter(child);
        queue.insert(child, child->weight);
      }
    }
  }

  void clear_ipoints() override { points.clear(); }
  void add_ipoint(float x, float y, float a, float w = 1.0) override {
    points.push_back({vec2(x, y), a, w});
  }

  void set_maxresults(int n) override { maxresults = n; }
  void set_minweight(float value) override { minweight = value; }
  void set_maxoffset(float value) override { maxoffset = value; }
  void set_lsq(bool value) override { lsq = value; }
  void set_unoriented(bool value) override { unoriented = value; }
  void set_error(float eps_, float aeps_) override {
    this->eps = eps_;
    this->aeps = aeps_;
  }
  void set_tolerance(float tol_, float atol_) override {
    this->tol = tol_;
    this->atol = atol_;
  }
  void set_breakpenalty(float, float) override {}
  void set_verbose(int value) override { verbose = value; }

  int nresults() override { return int(results.size()); }
  int nmatches(int i) override { return int(results[i]->matches->size()); }
  float weight(int i) override { return results[i]->weight; }
  float angle(int i) override { return results[i]->region.thm; }
  float offset(int i) override { return results[i]->region.rm; }
};

}  // namespace lumo_clinesp2d

std::unique_ptr<LinesP2D> makeLinesP2D() {
  return std::make_unique<lumo_clinesp2d::CLinesP2D>();
}
