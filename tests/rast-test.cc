// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

// Simple regression tests for the RAST library.
//
// These tests exercise only a small part of the library code. They are
// intended to catch major bloopers and crashes only ("regression to the
// non-working state"), not subtle numerical problems.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <memory>
#include <vector>

#include "rast.h"
#include "util.h"
#include "vec2.h"

#define begin_trials(NAME, N, FREQ)                                                                \
  if (1) {                                                                                         \
    const char *NAME__ = (NAME);                                                                   \
    std::fprintf(stderr, "================ STARTING %s (%d)\n", NAME__, (N));                      \
    for (int TRIAL = 0; TRIAL < (N); TRIAL++) {                                                    \
      try {                                                                                        \
        if (FREQ > 0 && TRIAL > 0 && TRIAL % FREQ == 0)                                            \
          std::fprintf(stderr, " %d", TRIAL);

#define end_trials                                                                                 \
  }                                                                                                \
  catch (const char *s) {                                                                          \
    std::fprintf(stderr, "\n%s: trial %d failed: %s\n", NAME__, TRIAL, s);                         \
  }                                                                                                \
  catch (const std::exception &e) {                                                                \
    std::fprintf(stderr, "\n%s: trial %d failed: %s\n", NAME__, TRIAL, e.what());                  \
  }                                                                                                \
  catch (...) {                                                                                    \
    std::fprintf(stderr, "\n%s: trial %d failed with unknown exception\n", NAME__, TRIAL);         \
  }                                                                                                \
  }                                                                                                \
  std::fprintf(stderr, "\n================ FINISHED %s\n", NAME__);                                \
  }

#undef assert
#define assert(X)                                                                                  \
  do {                                                                                             \
    if (!(X)) throw "ASSERTION FAILED: " #X;                                                       \
  } while (0)

struct Segment {
  vec2 u, v;
  Segment() = default;
  Segment(const vec2 &u_, const vec2 &v_) : u(u_), v(v_) {}
  vec2 sample(float lambda) const { return u + (v - u) * lambda; }
  float ldist(vec2 p) const {
    vec2 n = normal((v - u).normalized());
    float offset = n * u;
    return std::fabs(n * p - offset);
  }
  float offset() const {
    vec2 n = normal((v - u).normalized());
    return n * u;
  }
  float angle() const {
    vec2 n = normal(v - u);
    return std::atan2(n[1], n[0]);
  }
};

bool within(float x, float y, float eps) { return std::fabs(x - y) <= eps; }

float uniform(float low, float high) { return float(drand48()) * (high - low) + low; }
int iuniform(int low, int high) { return static_cast<int>(lrand48() % (high - low)) + low; }
vec2 vuniform(float x0, float y0, float x1, float y1) {
  return vec2(uniform(x0, x1), uniform(y0, y1));
}

void test_linesp2d_1() {
  begin_trials("linesp2d_1", 1000, 10) {
    std::unique_ptr<LinesP2D> lines(makeLinesP2D());
    vec2 u, v;
    do {
      u = vuniform(0, 0, 512, 512);
      v = vuniform(0, 0, 512, 512);
    } while (distance(u, v) < 1.0f);
    Segment seg(u, v);
    float tol = 0.1f, atol = 1e-4f;
    lines->set_tolerance(tol, atol);
    float offset = seg.offset();
    float angle = seg.angle();
    int npoints = iuniform(2, 200);
    for (int i = 0; i < npoints; i++) {
      vec2 p = seg.sample(static_cast<float>(i) / static_cast<float>(npoints));
      lines->add_ipoint(p[0], p[1], angle, 1.0f);
    }
    lines->compute();
    if (lines->nresults() < 1) throw "didn't get any results";
    float loffset = lines->offset(0);
    float langle = lines->angle(0);

    offset = std::fabs(offset);
    loffset = std::fabs(loffset);
    angle = normalize_orientation(angle);
    langle = normalize_orientation(langle);
    float s = 4;
    assert(within(loffset, offset, s * tol));
    assert(within(langle, angle, s * tol));
  }
  end_trials;
}

void test_liness2d_1() {
  begin_trials("liness2d_1", 1000, 10) {
    std::unique_ptr<LinesS2D> lines(makeLinesS2D());
    vec2 u, v;
    do {
      u = vuniform(0, 0, 512, 512);
      v = vuniform(0, 0, 512, 512);
    } while (distance(u, v) < 1.0f);
    Segment seg(u, v);
    float tol = 0.1f, atol = 1e-4f;
    lines->set_tolerance(tol, atol);
    float offset = seg.offset();
    float angle = seg.angle();
    int npoints = iuniform(2, 200);
    for (int i = 0; i < npoints; i++) {
      vec2 p = seg.sample(static_cast<float>(i) / static_cast<float>(npoints));
      vec2 q = seg.sample(0.5f + static_cast<float>(i) / static_cast<float>(npoints));
      lines->add_iseg(p[0], p[1], q[0], q[1], angle, 1.0f);
    }
    lines->compute();
    if (lines->nresults() < 1) throw "didn't get any results";
    float loffset = lines->offset(0);
    float langle = lines->angle(0);

    offset = std::fabs(offset);
    loffset = std::fabs(loffset);
    angle = normalize_orientation(angle);
    langle = normalize_orientation(langle);
    float s = 4;
    assert(within(loffset, offset, s * tol));
    assert(within(langle, angle, s * tol));
  }
  end_trials;
}

void test_rastp2d_1() {
  // bounded error
  begin_trials("rastp2d_1", 1000, 10) {
    std::unique_ptr<InstanceP2D> instance(makeInstanceP2D());
    std::unique_ptr<RastP2D> rast(makeRastP2D());
    instance->set_nclutter(0);
    instance->set_nmodel_total(20);
    instance->set_nmodel_unoccluded(20);
    instance->set_error(0.0);
    instance->set_aerror(0.0);
    instance->set_srange(0.5, 2.0f);
    instance->generate();
    float tol = 1e-2f;
    float eps = 1.0f;
    float aeps = 0.1f;
    rast->set_min_q(0);
    rast->set_tolerance(tol);
    rast->set_srange(0.5, 2.0f);
    assert(instance->nimage() == instance->nmodel());
    for (int i = 0; i < instance->nimage(); i++) {
      float x, y, a;
      instance->get_image(x, y, a, i);
      rast->add_ipoint(x, y, a);
    }
    for (int i = 0; i < instance->nmodel(); i++) {
      float x, y, a;
      instance->get_model(x, y, a, i);
      rast->add_msource(x, y, a, eps, aeps);
    }
    rast->match();
    assert(rast->nresults() > 0);
    assert(rast->ubound(0) == instance->nmodel());
    assert(within(rast->translation(0, 0), instance->get_param(0), 2.0f));
    assert(within(rast->translation(0, 1), instance->get_param(1), 2.0f));
    assert(within(rast->angle(0), instance->get_param(2), 0.07f));
    assert(within(rast->scale(0), instance->get_param(3), 0.05f));
  }
  end_trials;
}

void test_rastp2d_2() {
  // robust least square
  begin_trials("rastp2d_2", 1000, 10) {
    std::unique_ptr<InstanceP2D> instance(makeInstanceP2D());
    std::unique_ptr<RastP2D> rast(makeRastP2D());
    instance->set_nclutter(0);
    instance->set_nmodel_total(20);
    instance->set_nmodel_unoccluded(20);
    instance->set_error(0.0);
    instance->set_aerror(0.0);
    instance->set_srange(0.5, 2.0f);
    instance->generate();
    float tol = 1e-2f;
    float eps = 1.0f;
    float aeps = 0.1f;
    rast->set_min_q(0);
    rast->set_tolerance(tol);
    rast->set_srange(0.5, 2.0f);
    rast->set_lsq(true);
    assert(instance->nimage() == instance->nmodel());
    for (int i = 0; i < instance->nimage(); i++) {
      float x, y, a;
      instance->get_image(x, y, a, i);
      rast->add_ipoint(x, y, a);
    }
    for (int i = 0; i < instance->nmodel(); i++) {
      float x, y, a;
      instance->get_model(x, y, a, i);
      rast->add_msource(x, y, a, eps, aeps);
    }
    rast->match();
    assert(rast->nresults() > 0);
    assert(within(rast->translation(0, 0), instance->get_param(0), 0.1f));
    assert(within(rast->translation(0, 1), instance->get_param(1), 0.1f));
    assert(within(rast->angle(0), instance->get_param(2), 0.05f));
    assert(within(rast->scale(0), instance->get_param(3), 0.05f));
  }
  end_trials;
}

struct seg2 {
  vec2 u, v;
  seg2() = default;
  seg2(vec2 u_, vec2 v_) : u(u_), v(v_) {}
  seg2 transform(vec2 rot, vec2 tr) { return seg2(cmul(rot, u) + tr, cmul(rot, v) + tr); }
  float angle() const { return angleOf(v - u); }
  float length() const { return (v - u).magnitude(); }
};

void test_rasts2d_1() {
  begin_trials("rasts2d_1", 100, 1) {
    vec2 tr(urand(0.0, 512.0), urand(0.0, 512.0));
    float alpha = urand(0.0, 2 * M_PI);
    float scale = urand(0.95, 1.05);
    vec2 rot = vec2(scale * std::cos(alpha), scale * std::sin(alpha));
    int nmodel = 20;
    std::vector<seg2> model;
    for (int i = 0; i < nmodel; i++) {
      vec2 u(urand(-100.0, 100.0), urand(-100.0, 100.0));
      vec2 v(urand(-100.0, 100.0), urand(-100.0, 100.0));
      model.push_back(seg2(u, v));
    }
    std::vector<seg2> image;
    for (int i = int(model.size()) - 1; i >= 0; i--) {
      image.push_back(model[i].transform(rot, tr));
    }
    std::unique_ptr<RastS2D> rast(makeRastS2D());
    for (auto &s : image) rast->add_iseg(s.u[0], s.u[1], s.v[0], s.v[1]);
    for (auto &s : model) rast->add_mseg(s.u[0], s.u[1], s.v[0], s.v[1]);
    rast->set_srange(0.95, 1.05);
    rast->set_arange(0.0, 2 * M_PI);
    float eps = 1.0, aeps = 0.01;
    rast->set_eps(eps, aeps);
    rast->set_verbose(false);
    rast->set_lsq(true);
    rast->set_tolerance(0.1);
    rast->match();
    assert(rast->nresults() > 0);
    assert(within(rast->translation(0, 0), tr[0], 2.0f * eps));
    assert(within(rast->translation(0, 1), tr[1], 2.0f * eps));
    assert(within(rast->angle(0), alpha, 2.0f * aeps));
  }
  end_trials;
}

int main() {
  srand48(0);
  test_rasts2d_1();
  test_rastp2d_2();
  test_rastp2d_1();
  test_linesp2d_1();
  test_liness2d_1();
}
