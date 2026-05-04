// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <vector>

#include "heap.h"
#include "rast.h"
#include "util.h"
#include "vec2.h"

namespace lumo_crastp2d {

struct Ipoint {
  vec2 p;
  float a;
};

struct Msource {
  vec2 p;
  float a;
  float eps;
  float aeps;
};

static float angle_diff(float a1, float a2) {
  float d = a1 - a2;
  while (d < -kPi) d += kTwoPi;
  while (d > kPi) d -= kTwoPi;
  return std::fabs(d);
}

static float unoriented_angle_diff(float a1, float a2) {
  float d = a1 - a2;
  while (d < -kHalfPi) d += kPi;
  while (d > kHalfPi) d -= kPi;
  return std::fabs(d);
}

struct Region {
  std::vector<float> low;
  std::vector<float> high;
  vec2 translation() const { return vec2((high[0] + low[0]) / 2.0f, (high[1] + low[1]) / 2.0f); }
  float angle() const { return (high[2] + low[2]) / 2.0f; }
  float scale() const { return (high[3] + low[3]) / 2.0f; }
  vec2 rotation() const {
    float a = angle();
    float s = scale();
    return vec2(s * std::cos(a), s * std::sin(a));
  }
  float tdelta() const {
    return 1.5f * std::max((high[0] - low[0]) / 2.0f, (high[1] - low[1]) / 2.0f);
  }
  float adelta() const { return (high[2] - low[2]) / 2.0f; }
  float smax() const { return high[3]; }
  float sdelta() const { return (high[3] - low[3]) / 2.0f; }
};

struct IMPair {
  short msource;
  short ipoint;
  IMPair() : msource(0), ipoint(0) {}
  IMPair(int ms, int ip) : msource(short(ms)), ipoint(short(ip)) {}
};

using Pairs = std::vector<IMPair>;
using CPairs = std::shared_ptr<Pairs>;

class CRastP2D;

struct State {
  int depth{0};
  int generation{0};
  Region region;
  CPairs parent_matches{std::make_shared<Pairs>()};
  CPairs matches{std::make_shared<Pairs>()};
  float lbound{0.0f};
  float ubound{0.0f};

  void print(FILE *stream) {
    std::fprintf(stream, "<%d u %g l %g low %g %g %g %g high %g %g %g %g>", depth, ubound, lbound,
                 region.low[0], region.low[1], region.low[2], region.low[3], region.high[0],
                 region.high[1], region.high[2], region.high[3]);
  }

  void set(int depth_, const Region &oregion, CPairs omatches) {
    this->depth = depth_;
    region = oregion;
    parent_matches = std::move(omatches);
  }

  void init(const std::vector<Msource> &msources, const std::vector<Ipoint> &ipoints) {
    depth = 0;
    region.low = {0.0f, 0.0f, 0.0f};
    region.high = {0.0f, 0.0f, 0.0f};
    lbound = 0.0f;
    ubound = float(msources.size());
    Pairs &omatches = *parent_matches;
    omatches.clear();
    for (std::size_t i = 0; i < msources.size(); i++) {
      for (std::size_t j = 0; j < ipoints.size(); j++) {
        omatches.emplace_back(int(i), int(j));
      }
    }
  }
  void eval(CRastP2D &env);
};

using CState = std::shared_ptr<State>;

struct CRastP2D : RastP2D {
  std::vector<float> splitscale{1.0f, 1.0f, 500.0f, 500.0f};

  bool final(const Region &r, float delta) {
    for (std::size_t i = 0; i < r.low.size(); i++) {
      float v = (r.high[i] - r.low[i]) * splitscale[i];
      if (v > delta) return false;
    }
    return true;
  }

  void split(Region &left, Region &right, const Region &r) {
    int dim = int(r.low.size());
    int mi = -1;
    float mv = 0.0;
    for (int i = 0; i < dim; i++) {
      float v = (r.high[i] - r.low[i]) * splitscale[i];
      if (v < mv) continue;
      mv = v;
      mi = i;
    }
    float meanv = (r.high[mi] + r.low[mi]) / 2.0f;
    left.low = r.low;
    left.high = r.high;
    left.high[mi] = meanv;
    right.low = r.low;
    right.low[mi] = meanv;
    right.high = r.high;
  }

  std::vector<Ipoint> ipoints;
  std::vector<Msource> msources;
  std::vector<bool> used;

  Heap<CState> queue;
  std::vector<CState> results;

  bool verbose{false};
  float tolerance{1e-3f};
  float min_q{3.0f};
  int maxresults{1};
  std::vector<float> tlow{-1000.0f, -1000.0f, 0.0f, 0.9f};
  std::vector<float> thigh{1000.0f, 1000.0f, float(kTwoPi), 1.1f};
  int generation{1};
  bool use_lsq{false};
  bool unoriented{true};

  double priority(CState state) {
    double priority = state->ubound + 1e-4 * state->lbound;
    assert(priority < state->ubound + 1);
    return priority;
  }

  int n_nodes{0};
  int n_transforms{0};
  int n_distances{0};

  void start_match() {
    n_nodes = 0;
    n_transforms = 0;
    n_distances = 0;
    results.clear();
    queue.clear();
    used.assign(ipoints.size(), false);
    auto initial_state = std::make_shared<State>();
    initial_state->init(msources, ipoints);
    initial_state->region.low = tlow;
    initial_state->region.high = thigh;
    initial_state->eval(*this);
    initial_state->generation = generation;
    queue.insert(initial_state, initial_state->ubound);
    for (int iter = 0;; iter++) {
      if (queue.empty()) break;
      CState top = queue.extractMax();
      if (top->generation != generation) {
        top->eval(*this);
        top->generation = generation;
        queue.insert(top, priority(top));
        continue;
      }
      if (verbose && iter % 10000 == 0) {
        float q = !results.empty() ? results[0]->ubound : 0.0f;
        std::fprintf(stderr, "# %10d result %6g queue %7zu", iter, q, 1 + queue.length());
        std::fprintf(stderr, "   ");
        top->print(stderr);
        std::fprintf(stderr, "\n");
      }
      if (top->ubound == top->lbound || final(top->region, tolerance)) {
        results.push_back(top);
        const Pairs &matches = *top->matches;
        for (const auto &pr : matches) used[pr.ipoint] = true;
        generation++;
        if (int(results.size()) >= maxresults) return;
        continue;
      }
      Region subregions[2];
      CState substates[2];
      split(subregions[0], subregions[1], top->region);
      for (int i = 0; i < 2; i++) {
        substates[i] = std::make_shared<State>();
        substates[i]->set(top->depth + 1, subregions[i], top->matches);
        substates[i]->eval(*this);
        substates[i]->generation = generation;
        if (substates[i]->ubound < min_q) continue;
        queue.insert(substates[i], priority(substates[i]));
      }
    }
  }

  void clear_msources() override { msources.clear(); }
  void add_msource(float x, float y, float a, float eps, float aeps) override {
    msources.push_back({vec2(x, y), a, eps, aeps});
  }

  void clear_ipoints() override { ipoints.clear(); }
  void add_ipoint(float x, float y, float a) override { ipoints.push_back({vec2(x, y), a}); }

  void set_maxresults(int n) override { maxresults = n; }
  void set_verbose(bool value) override { verbose = value; }
  void set_tolerance(float value) override {
    if (value < 1e-3)
      throw std::invalid_argument("tolerance too small; would fail to converge occasionally");
    tolerance = value;
  }
  void set_min_q(float min_q_) override { this->min_q = min_q_; }
  void set_xrange(float x0, float x1) override {
    tlow[0] = x0;
    thigh[0] = x1;
  }
  void set_yrange(float y0, float y1) override {
    tlow[1] = y0;
    thigh[1] = y1;
  }
  void set_arange(float a0, float a1) override {
    tlow[2] = a0;
    thigh[2] = a1;
  }
  void set_srange(float s0, float s1) override {
    tlow[3] = s0;
    thigh[3] = s1;
  }
  void set_lsq(bool value) override { use_lsq = value; }
  void set_unoriented(bool value) override { unoriented = value; }
  void match() override { start_match(); }

  int nresults() const override { return int(results.size()); }
  float ubound(int rank) const override { return results[rank]->ubound; }
  float lbound(int rank) const override { return results[rank]->lbound; }
  float translation(int rank, int dim) const override { return results[rank]->region.translation()[dim]; }
  float angle(int rank) const override { return results[rank]->region.angle(); }
  float scale(int rank) const override { return results[rank]->region.scale(); }
};

void State::eval(CRastP2D &env) {
  env.n_nodes++;
  const std::vector<Msource> &msources = env.msources;
  const std::vector<Ipoint> &ipoints = env.ipoints;
  const std::vector<bool> &used = env.used;

  Pairs &nmatches = *matches;
  nmatches.clear();
  lbound = 0.0;
  ubound = 0.0;

  vec2 translation = region.translation();
  vec2 rotation = region.rotation();
  float tdelta = region.tdelta();
  float angle = region.angle();
  float adelta = region.adelta();
  float sdelta = region.sdelta();
  float smax = region.smax();

  const Pairs &omatches = *parent_matches;
  int n = int(omatches.size());
  for (int i = 0; i < n;) {
    env.n_transforms++;
    int msource_index = omatches[i].msource;
    const Msource &msource = msources[msource_index];
    vec2 tmpoint = cmul(rotation, msource.p) + translation;
    float eps2 = sqr(msource.eps);
    float nmsource = norm(msource.p);
    float strict = msource.eps;
    float delta = tdelta + nmsource * smax * adelta + nmsource * sdelta;
    float loose = strict + delta;

    float tangle = angle + msource.a;
    float astrict = msource.aeps;
    float aloose = astrict + adelta;

    float llbound = 0.0;
    float lubound = 0.0;

    for (; i < n && omatches[i].msource == msource_index; i++) {
      env.n_distances++;
      int ipoint_index = omatches[i].ipoint;
      if (used[ipoint_index]) continue;
      const Ipoint &ipoint = ipoints[ipoint_index];
      float adiff;
      if (env.unoriented)
        adiff = unoriented_angle_diff(ipoint.a, tangle);
      else
        adiff = angle_diff(ipoint.a, tangle);
      if (adiff > aloose) continue;
      float err = distance(tmpoint, ipoint.p);
      if (err > loose) continue;
      if (env.use_lsq) {
        float ud = std::max(0.0f, err - delta);
        float uq = std::max(0.0f, 1.0f - ud * ud / eps2);
        lubound = std::max(lubound, uq);
        float ld = err;
        float lq = std::max(0.0f, 1.0f - ld * ld / eps2);
        llbound = std::max(llbound, lq);
      } else {
        if (err < strict && adiff < astrict) llbound = 1.0;
        lubound = 1.0;
      }
      nmatches.emplace_back(msource_index, ipoint_index);
    }
    lbound += llbound;
    ubound += lubound;
  }
}

}  // namespace lumo_crastp2d

std::unique_ptr<RastP2D> makeRastP2D() {
  return std::make_unique<lumo_crastp2d::CRastP2D>();
}
