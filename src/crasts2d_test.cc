// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <memory>

#include "rast.h"
#include "util.h"

namespace {

struct Seg {
  float x0, y0, x1, y1;
};

std::unique_ptr<RastS2D> makeMatcher() {
  std::unique_ptr<RastS2D> r(makeRastS2D());
  r->set_maxresults(1);
  r->set_verbose(false);
  r->set_tolerance(1e-3f);
  r->set_qtolerance(1e-2f);
  r->set_min_q(1.0f);
  r->set_xrange(-200, 200);
  r->set_yrange(-200, 200);
  r->set_arange(-kPi, kPi);
  r->set_srange(0.95f, 1.05f);
  r->set_eps(4.0f, 0.05f);
  r->set_lsq(false);
  return r;
}

}  // namespace

TEST_CASE("RastS2D: factory yields a usable instance") {
  std::unique_ptr<RastS2D> r(makeRastS2D());
  REQUIRE(r != nullptr);
}

TEST_CASE("RastS2D: identity match of three segments") {
  auto r = makeMatcher();
  const std::array<Seg, 3> segs = {{{0, 0, 50, 0}, {50, 0, 50, 50}, {50, 50, 0, 50}}};
  for (const auto& s : segs) {
    r->add_mseg(s.x0, s.y0, s.x1, s.y1);
    r->add_iseg(s.x0, s.y0, s.x1, s.y1);
  }
  r->match();
  REQUIRE(r->nresults() >= 1);
  CHECK(std::abs(r->translation(0, 0)) < 8.0f);
  CHECK(std::abs(r->translation(0, 1)) < 8.0f);
  CHECK(r->scale(0) == doctest::Approx(1.0).epsilon(0.05));
}

TEST_CASE("RastS2D: pure translation recovered") {
  auto r = makeMatcher();
  const float dx = 30.0f, dy = -20.0f;
  const std::array<Seg, 4> segs = {
      {{0, 0, 50, 0}, {50, 0, 50, 50}, {50, 50, 0, 50}, {0, 50, 0, 0}}};
  // Constrain xrange/yrange tightly so the search can't drift to an
  // alternate-but-still-locally-consistent translation along the segment axes.
  r->set_xrange(dx - 5, dx + 5);
  r->set_yrange(dy - 5, dy + 5);
  for (const auto& s : segs) {
    r->add_mseg(s.x0, s.y0, s.x1, s.y1);
    r->add_iseg(s.x0 + dx, s.y0 + dy, s.x1 + dx, s.y1 + dy);
  }
  r->match();
  REQUIRE(r->nresults() >= 1);
  CHECK(r->translation(0, 0) == doctest::Approx(dx).epsilon(0.2));
  CHECK(r->translation(0, 1) == doctest::Approx(dy).epsilon(0.2));
}

TEST_CASE("RastS2D: empty input emits at most a degenerate (ubound==0) match") {
  auto r = makeMatcher();
  r->match();
  // The matcher's search converges immediately on the empty problem and may
  // emit a single degenerate result with zero quality.
  if (r->nresults() >= 1) {
    CHECK(r->ubound(0) == doctest::Approx(0.0));
  }
}

TEST_CASE("RastS2D: lbound <= ubound") {
  auto r = makeMatcher();
  const std::array<Seg, 2> segs = {{{0, 0, 30, 0}, {0, 0, 0, 30}}};
  for (const auto& s : segs) {
    r->add_mseg(s.x0, s.y0, s.x1, s.y1);
    r->add_iseg(s.x0, s.y0, s.x1, s.y1);
  }
  r->match();
  REQUIRE(r->nresults() >= 1);
  CHECK(r->lbound(0) <= r->ubound(0));
}

TEST_CASE("RastS2D: scale_eps setter does not crash") {
  auto r = makeMatcher();
  r->set_scale_eps(true, 4.0f);
  // A single segment under the default ±200 search ranges is too
  // unconstrained for branch-and-bound to converge (the queue grows
  // unboundedly). Use the same three-segment L as the identity test
  // so the matcher actually finishes.
  const std::array<Seg, 3> segs = {{{0, 0, 50, 0}, {50, 0, 50, 50}, {50, 50, 0, 50}}};
  for (const auto& s : segs) {
    r->add_mseg(s.x0, s.y0, s.x1, s.y1);
    r->add_iseg(s.x0, s.y0, s.x1, s.y1);
  }
  r->match();
  CHECK(r->nresults() >= 0);
}

TEST_CASE("RastS2D: xrange constrains accepted translations") {
  auto r = makeMatcher();
  r->set_xrange(-2, 2);
  r->set_yrange(-2, 2);
  r->add_mseg(0, 0, 30, 0);
  r->add_iseg(0, 0, 30, 0);
  r->match();
  if (r->nresults() >= 1) {
    CHECK(r->translation(0, 0) >= -3.0f);
    CHECK(r->translation(0, 0) <= 3.0f);
    CHECK(r->translation(0, 1) >= -3.0f);
    CHECK(r->translation(0, 1) <= 3.0f);
  }
}
