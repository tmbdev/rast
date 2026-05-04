// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <utility>

#include "rast.h"

namespace {

using FloatPair = std::pair<float, float>;

std::unique_ptr<RastP2D> makeMatcher() {
  std::unique_ptr<RastP2D> r(makeRastP2D());
  r->set_maxresults(1);
  r->set_verbose(false);
  r->set_tolerance(1e-3f);
  r->set_min_q(3.0f);
  r->set_xrange(-1000, 1000);
  r->set_yrange(-1000, 1000);
  r->set_arange(0.0f, 2 * M_PI);
  r->set_srange(0.8f, 1.2f);
  r->set_lsq(false);
  r->set_unoriented(true);
  return r;
}

}  // namespace

TEST_CASE("RastP2D: factory yields a usable instance") {
  std::unique_ptr<RastP2D> r(makeRastP2D());
  REQUIRE(r != nullptr);
}

TEST_CASE("RastP2D: identity match recovers near-zero translation and unit scale") {
  auto r = makeMatcher();
  const float eps = 5.0f, aeps = 0.1f;
  const std::array<FloatPair, 4> points = {{{0, 0}, {30, 0}, {0, 30}, {15, 20}}};
  for (const auto& p : points) {
    r->add_msource(p.first, p.second, 0.0f, eps, aeps);
    r->add_ipoint(p.first, p.second, 0.0f);
  }
  r->match();
  REQUIRE(r->nresults() >= 1);
  CHECK(std::abs(r->translation(0, 0)) < 5.0f);
  CHECK(std::abs(r->translation(0, 1)) < 5.0f);
  CHECK(r->scale(0) == doctest::Approx(1.0).epsilon(0.05));
}

TEST_CASE("RastP2D: pure translation is recovered") {
  auto r = makeMatcher();
  const float eps = 5.0f, aeps = 0.1f;
  const float dx = 50.0f, dy = -30.0f;
  const std::array<FloatPair, 4> points = {{{0, 0}, {30, 0}, {0, 30}, {15, 20}}};
  for (const auto& p : points) {
    r->add_msource(p.first, p.second, 0.0f, eps, aeps);
    r->add_ipoint(p.first + dx, p.second + dy, 0.0f);
  }
  r->match();
  REQUIRE(r->nresults() >= 1);
  CHECK(r->translation(0, 0) == doctest::Approx(dx).epsilon(0.2));
  CHECK(r->translation(0, 1) == doctest::Approx(dy).epsilon(0.2));
}

TEST_CASE("RastP2D: ubound is at most msources count") {
  auto r = makeMatcher();
  const std::array<FloatPair, 3> points = {{{0, 0}, {10, 0}, {0, 10}}};
  for (const auto& p : points) {
    r->add_msource(p.first, p.second, 0.0f, 5.0f, 0.1f);
    r->add_ipoint(p.first, p.second, 0.0f);
  }
  r->match();
  REQUIRE(r->nresults() >= 1);
  CHECK(r->ubound(0) <= 3);
  CHECK(r->lbound(0) <= r->ubound(0));
}

TEST_CASE("RastP2D: empty input emits at most a degenerate (ubound==0) match") {
  auto r = makeMatcher();
  r->match();
  if (r->nresults() >= 1) {
    CHECK(r->ubound(0) == doctest::Approx(0.0));
  }
}

TEST_CASE("RastP2D: synthetic InstanceP2D match (deterministic)") {
  srand48(7);
  std::unique_ptr<InstanceP2D> inst(makeInstanceP2D());
  inst->set_image_size(512);
  inst->set_model_size(100);
  inst->set_nclutter(20);
  inst->set_nmodel_total(15);
  inst->set_nmodel_unoccluded(12);
  inst->set_error(2.0f);
  inst->set_aerror(0.05f);
  inst->set_srange(0.95f, 1.05f);
  inst->generate();

  auto r = makeMatcher();
  r->set_xrange(0, 512);
  r->set_yrange(0, 512);
  r->set_srange(0.95f, 1.05f);
  r->set_min_q(3.0f);
  for (int i = 0; i < inst->nmodel(); i++) {
    float x = 0, y = 0, a = 0;
    inst->get_model(x, y, a, i);
    r->add_msource(x, y, a, 4.0f, 0.1f);
  }
  for (int i = 0; i < inst->nimage(); i++) {
    float x = 0, y = 0, a = 0;
    inst->get_image(x, y, a, i);
    r->add_ipoint(x, y, a);
  }
  r->match();
  REQUIRE(r->nresults() >= 1);
  CHECK(r->translation(0, 0) == doctest::Approx(inst->get_param(0)).epsilon(0.05));
  CHECK(r->translation(0, 1) == doctest::Approx(inst->get_param(1)).epsilon(0.05));
}

TEST_CASE("RastP2D: maxresults limits results") {
  auto r = makeMatcher();
  r->set_maxresults(3);
  const std::array<FloatPair, 3> points = {{{0, 0}, {10, 0}, {0, 10}}};
  for (const auto& p : points) {
    r->add_msource(p.first, p.second, 0.0f, 5.0f, 0.1f);
    r->add_ipoint(p.first, p.second, 0.0f);
  }
  r->match();
  CHECK(r->nresults() <= 3);
}

TEST_CASE("RastP2D: clear_msources / clear_ipoints reset state") {
  auto r = makeMatcher();
  r->add_msource(0, 0, 0, 5.0f, 0.1f);
  r->add_ipoint(0, 0, 0);
  r->clear_msources();
  r->clear_ipoints();
  r->match();
  // After clearing all inputs, the matcher converges immediately on the empty
  // problem; any emitted result is a degenerate ubound==0 match.
  if (r->nresults() >= 1) {
    CHECK(r->ubound(0) == doctest::Approx(0.0));
  }
}

TEST_CASE("RastP2D: tolerance below 1e-3 is rejected") {
  auto r = makeMatcher();
  CHECK_THROWS_AS(r->set_tolerance(1e-6f), std::invalid_argument);
}
