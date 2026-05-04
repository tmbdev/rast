// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <utility>
#include <vector>

#include "heap.h"
#include "rast.h"
#include "util.h"
#include "vec2.h"

namespace lumo_crastrs2d {

struct Blob {
  vec2 p;
  float r{0};
  Blob() = default;
  Blob(vec2 p_, float r_) { set(p_, r_); }
  void set(vec2 p_, float r_) {
    this->p = p_;
    this->r = r_;
  }
  float within(float eps, float delta, vec2 a) const {
    float d = p.distance(a);
    return (d >= eps - delta && d <= eps + delta);
  }
  float lsq(float eps, float delta, vec2 a) const {
    float d = p.distance(a);
    return std::max(0.0, 1.0 - sqr(std::max(0.0f, (d - delta) / (eps * eps))));
  }
};

struct Segment {
  vec2 p, q;
  float a{0};
  vec2 dir;
  float l0{0}, l1{0};
  vec2 normal;
  float d{0};
  Segment() = default;
  Segment(vec2 p_, vec2 q_) { set(p_, q_); }
  void set(vec2 p_, vec2 q_) {
    this->p = p_;
    this->q = q_;
    dir = (q_ - p_).normalized();
    l0 = dir * p_;
    l1 = dir * q_;
    if (l1 < l0) std::swap(l0, l1);
    normal = ::normal(dir);
    d = normal * p_;
    a = angleOf(dir);
  }
  float length() const { return (q - p).magnitude(); }
  vec2 sub(float l) const { return p + dir * l; }
  float within(float eps, float delta, vec2 a_) const {
    float la = dir * a_;
    if (la < l0 - eps - delta || la > l1 + eps + delta) return 0;
    float err = std::max(0.0f, std::fabs(normal * a_ - d) - delta);
    return err < eps;
  }
  float lsq(float eps, float delta, vec2 a_) const {
    float la = dir * a_;
    float q_ = 1.0;
    float dl = la - (l0 - delta);
    if (dl < -eps) return 0;
    if (dl < 0) q_ = 1.0 - sqr(dl * dl / (eps * eps));
    float dr = la - (l1 + delta);
    if (dr > eps) return 0;
    if (dr > 0) q_ = 1.0 - sqr(dr * dr / (eps * eps));
    q_ *= std::max(0.0, 1.0 - sqr(std::max(0.0f, std::fabs(normal * a_ - d) - delta)) /
                            (eps * eps));
    return q_;
  }
};

using Ipoint = Segment;
using Msource = Segment;

static float angle_diff(float a1, float a2) {
  float d = a1 - a2;
  while (d < -M_PI) d += 2 * M_PI;
  while (d > M_PI) d -= 2 * M_PI;
  return std::fabs(d);
}

static float unoriented_angle_diff(float a1, float a2) {
  float d = a1 - a2;
  while (d < -M_PI / 2) d += M_PI;
  while (d > M_PI / 2) d -= M_PI;
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
  unsigned char msource;
  unsigned char l;
  short ipoint;
  IMPair() : msource(0), l(0), ipoint(0) {}
  IMPair(int ms, int ip, float l_) : msource((unsigned char)ms), l((unsigned char)l_),
                                     ipoint((short)ip) {}
};

using Pairs = std::vector<IMPair>;
using CPairs = std::shared_ptr<Pairs>;

class CRastRS2D;

struct State {
  int depth{0};
  int generation{0};
  Region region;
  CPairs parent_matches{std::make_shared<Pairs>()};
  CPairs matches{std::make_shared<Pairs>()};
  float lbound{0.0f};
  float ubound{0.0f};

  void print(FILE *stream = stdout) {
    std::fprintf(stream, "<%d [%g:%g] (%g:%g %g:%g) %g:%g %g:%g>", depth, lbound, ubound,
                 region.low[0], region.high[0], region.low[1], region.high[1], region.low[2],
                 region.high[2], region.low[3], region.high[3]);
  }

  void set(int depth_, const Region &oregion, CPairs omatches) {
    this->depth = depth_;
    region = oregion;
    parent_matches = std::move(omatches);
  }

  void init(std::vector<Msource> &msources, std::vector<Ipoint> &ipoints, float sdist) {
    if (msources.size() > 255) throw "too many line segments in model";
    if (ipoints.size() > 32000) throw "too many line segments in image";
    depth = 0;
    region.low = {0.0f, 0.0f, 0.0f};
    region.high = {0.0f, 0.0f, 0.0f};
    lbound = 0.0f;
    ubound = float(msources.size());
    Pairs &omatches = *parent_matches;
    omatches.clear();
    for (std::size_t i = 0; i < msources.size(); i++) {
      for (float l = 0.0, ml = msources[i].length(); l < ml / 2; l += sdist) {
        for (std::size_t j = 0; j < ipoints.size(); j++) {
          omatches.emplace_back(int(i), int(j), l);
          omatches.emplace_back(int(i), int(j), ml - l);
        }
      }
    }
  }
  void eval(CRastRS2D &env);
};

using CState = std::shared_ptr<State>;

struct CRastRS2D : RastRS2D {
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
  std::vector<float> thigh{1000.0f, 1000.0f, float(2 * M_PI), 1.1f};
  int generation{1};
  bool use_lsq{false};
  bool unoriented{true};
  float eps{3.0f};
  float aeps{0.05f};
  float sdist{3.0f};
  float qtolerance{1e-4f};

  double priority(CState state) {
    double priority = state->ubound + 1e-4 * state->lbound;
    if (priority >= state->ubound + 1) throw "error";
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
    initial_state->init(msources, ipoints, sdist);
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
      if ((1.0 - qtolerance) * top->ubound <= top->lbound || final(top->region, tolerance)) {
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
  void add_mseg(float x, float y, float x1, float y1) override {
    Msource ms;
    ms.set(vec2(x, y), vec2(x1, y1));
    msources.push_back(ms);
  }

  void clear_ipoints() override { ipoints.clear(); }
  void add_iseg(float x, float y, float x1, float y1) override {
    Ipoint ip;
    ip.set(vec2(x, y), vec2(x1, y1));
    ipoints.push_back(ip);
  }

  void set_maxresults(int n) override { maxresults = n; }
  void set_verbose(bool value) override { verbose = value; }
  void set_eps(float eps_, float aeps_, float sdist_) override {
    this->eps = eps_;
    this->aeps = aeps_;
    this->sdist = sdist_;
  }
  void set_tolerance(float value) override {
    if (value < 1e-3) throw "tolerance too small; would fail to converge occasionally";
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
  void set_qtolerance(float value) override { qtolerance = value; }
  void match() override { start_match(); }

  int nresults() override { return int(results.size()); }
  float ubound(int rank) override { return results[rank]->ubound; }
  float lbound(int rank) override { return results[rank]->lbound; }
  float translation(int rank, int dim) override { return results[rank]->region.translation()[dim]; }
  float angle(int rank) override { return results[rank]->region.angle(); }
  float scale(int rank) override { return results[rank]->region.scale(); }
};

void State::eval(CRastRS2D &env) {
  env.n_nodes++;
  std::vector<Msource> &msources = env.msources;
  std::vector<Ipoint> &ipoints = env.ipoints;
  std::vector<bool> &used = env.used;

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
  float eps = env.eps;
  float aeps = env.aeps;
  for (int i = 0; i < n;) {
    env.n_transforms++;
    int msource_index = omatches[i].msource;
    Msource &msource = msources[msource_index];
    vec2 tmpoint0 = cmul(rotation, msource.p) + translation;
    vec2 tmpoint1 = cmul(rotation, msource.q) + translation;
    Segment tmseg(tmpoint0, tmpoint1);
    float nmsource = std::max(norm(msource.p), norm(msource.q));
    float delta = tdelta + nmsource * smax * adelta + nmsource * sdelta;

    float tangle = angle + msource.a;
    float aloose = aeps + adelta;

    for (; i < n;) {
      if (omatches[i].msource != msource_index) break;
      float start_l = omatches[i].l;
      float llbound = 0.0;
      float lubound = 0.0;
      for (; i < n; i++) {
        float l = omatches[i].l;
        if (l != start_l) break;
        vec2 tmpoint = tmseg.sub(l);
        env.n_distances++;
        int ipoint_index = omatches[i].ipoint;
        if (used[ipoint_index]) continue;
        Ipoint &ipoint = ipoints[ipoint_index];
        float adiff;
        if (env.unoriented)
          adiff = unoriented_angle_diff(ipoint.a, tangle);
        else
          adiff = angle_diff(ipoint.a, tangle);
        if (adiff > aloose) continue;
        if (!env.use_lsq) {
          float uq = ipoint.within(eps, delta, tmpoint);
          if (uq <= 0.0) continue;
          float lq = ipoint.within(eps, 0.0, tmpoint);
          lubound = std::max(lubound, uq);
          llbound = std::max(llbound, lq);
          nmatches.emplace_back(msource_index, ipoint_index, l);
        } else {
          float uq = ipoint.lsq(eps, delta, tmpoint);
          if (uq <= 0.0) continue;
          float lq = ipoint.lsq(eps, 0.0, tmpoint);
          lubound = std::max(lubound, uq);
          llbound = std::max(llbound, lq);
          nmatches.emplace_back(msource_index, ipoint_index, l);
        }
      }
      lbound += llbound;
      ubound += lubound;
    }
  }
}

}  // namespace lumo_crastrs2d

std::unique_ptr<RastRS2D> makeRastRS2D() {
  return std::make_unique<lumo_crastrs2d::CRastRS2D>();
}
