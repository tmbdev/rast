// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <memory>
#include <utility>

#include "rast.h"

namespace {

using FloatPair = std::pair<float, float>;

constexpr float kEps = 0.5f;

std::unique_ptr<AlignmentP2D> makeAligner(float epsilon = kEps) {
  std::unique_ptr<AlignmentP2D> a(makeAlignmentP2D());
  a->set_epsilon(epsilon);
  a->set_srange(0.5f, 2.0f);
  return a;
}

}  // namespace

TEST_CASE("AlignmentP2D: factory yields a usable instance") {
  std::unique_ptr<AlignmentP2D> a(makeAlignmentP2D());
  REQUIRE(a != nullptr);
}

TEST_CASE("AlignmentP2D: identity match recovers translation 0 and scale 1") {
  auto a = makeAligner();
  const std::array<FloatPair, 4> points = {{{0, 0}, {10, 0}, {0, 10}, {7, 5}}};
  for (const auto& p : points) {
    a->add_mpoint(p.first, p.second);
    a->add_ipoint(p.first, p.second);
  }
  a->compute();
  CHECK(a->quality() >= 4);
  CHECK(a->translation(0) == doctest::Approx(0.0).epsilon(0.05));
  CHECK(a->translation(1) == doctest::Approx(0.0).epsilon(0.05));
  CHECK(a->scale() == doctest::Approx(1.0).epsilon(0.05));
}

TEST_CASE("AlignmentP2D: pure translation is recovered") {
  auto a = makeAligner();
  const float dx = 17.0f, dy = -23.0f;
  const std::array<FloatPair, 4> points = {{{0, 0}, {10, 0}, {0, 10}, {7, 5}}};
  for (const auto& p : points) {
    a->add_mpoint(p.first, p.second);
    a->add_ipoint(p.first + dx, p.second + dy);
  }
  a->compute();
  CHECK(a->quality() >= 4);
  CHECK(a->translation(0) == doctest::Approx(dx).epsilon(0.02));
  CHECK(a->translation(1) == doctest::Approx(dy).epsilon(0.02));
  CHECK(a->scale() == doctest::Approx(1.0).epsilon(0.05));
}

TEST_CASE("AlignmentP2D: pure rotation is recovered") {
  auto a = makeAligner();
  const float theta = static_cast<float>(M_PI / 6);  // 30 degrees
  const float c = std::cos(theta), s = std::sin(theta);
  const std::array<FloatPair, 4> model = {{{0, 0}, {10, 0}, {0, 10}, {7, 5}}};
  for (const auto& p : model) {
    a->add_mpoint(p.first, p.second);
    a->add_ipoint(c * p.first - s * p.second, s * p.first + c * p.second);
  }
  a->compute();
  CHECK(a->quality() >= 4);
  CHECK(a->scale() == doctest::Approx(1.0).epsilon(0.05));
  CHECK(std::cos(a->angle() - theta) == doctest::Approx(1.0).epsilon(0.01));
}

TEST_CASE("AlignmentP2D: scaled match recovers scale") {
  auto a = makeAligner(/*epsilon=*/1.0f);
  const float scale = 1.5f;
  const std::array<FloatPair, 4> points = {{{0, 0}, {10, 0}, {0, 10}, {7, 5}}};
  for (const auto& p : points) {
    a->add_mpoint(p.first, p.second);
    a->add_ipoint(p.first * scale, p.second * scale);
  }
  a->compute();
  CHECK(a->quality() >= 4);
  CHECK(a->scale() == doctest::Approx(scale).epsilon(0.05));
}

TEST_CASE("AlignmentP2D: scale outside srange is rejected") {
  auto a = makeAligner();
  a->set_srange(0.95f, 1.05f);  // tight band around 1.0
  const std::array<FloatPair, 3> points = {{{0, 0}, {10, 0}, {0, 10}}};
  for (const auto& p : points) {
    a->add_mpoint(p.first, p.second);
    a->add_ipoint(p.first * 2.0f, p.second * 2.0f);
  }
  a->compute();
  CHECK(a->quality() < 3);
}

TEST_CASE("AlignmentP2D: clear_mpoints / clear_ipoints reset state") {
  auto a = makeAligner();
  a->add_mpoint(0, 0);
  a->add_mpoint(10, 0);
  a->add_ipoint(0, 0);
  a->add_ipoint(10, 0);
  a->clear_mpoints();
  a->clear_ipoints();
  a->add_mpoint(0, 0);
  a->add_mpoint(5, 0);
  a->add_ipoint(0, 0);
  a->add_ipoint(5, 0);
  a->compute();
  CHECK(a->translation(0) == doctest::Approx(0.0).epsilon(0.1));
  CHECK(a->translation(1) == doctest::Approx(0.0).epsilon(0.1));
}
