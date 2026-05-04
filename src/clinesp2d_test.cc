// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#include <doctest/doctest.h>

#include <cmath>
#include <memory>

#include "rast.h"
#include "util.h"

namespace {

std::unique_ptr<LinesP2D> makeFinder() {
  std::unique_ptr<LinesP2D> finder(makeLinesP2D());
  finder->set_maxresults(1);
  finder->set_error(2.0f, 0.1f);
  finder->set_tolerance(0.1f, 0.001f);
  finder->set_minweight(3.0f);
  finder->set_maxoffset(1000.0f);
  finder->set_lsq(false);
  finder->set_unoriented(true);
  return finder;
}

}  // namespace

TEST_CASE("LinesP2D: factory yields a usable instance") {
  std::unique_ptr<LinesP2D> finder(makeLinesP2D());
  REQUIRE(finder != nullptr);
}

TEST_CASE("LinesP2D: horizontal line at y=cy is recovered") {
  auto finder = makeFinder();
  const float cy = 250.0f;
  for (float x = -200; x <= 200; x += 5.0f) {
    finder->add_ipoint(x, cy, /*a=*/0.0f, /*w=*/1.0f);
  }
  finder->compute();
  REQUIRE(finder->nresults() >= 1);
  // For a horizontal line in normal-form representation, theta is the angle
  // of the normal vector and offset = r in r = x*cos(theta) + y*sin(theta).
  // Here y = cy means theta = pi/2 and offset = cy.
  const float theta = finder->angle(0);
  const float offset = finder->offset(0);
  CHECK(std::cos(theta - static_cast<float>(M_PI_2)) == doctest::Approx(1.0).epsilon(0.05));
  CHECK(offset == doctest::Approx(cy).epsilon(0.05));
}

TEST_CASE("LinesP2D: angled line slope ~0.3") {
  auto finder = makeFinder();
  const float cx = 256.0f, cy = 256.0f;
  const float dx = 1.0f, dy = 0.3f;
  const float len = std::sqrt(dx * dx + dy * dy);
  const float lineAngle = std::atan2(dy, dx);
  for (float t = -100.0f; t <= 100.0f; t += 2.0f) {
    finder->add_ipoint(cx + t * dx / len, cy + t * dy / len, lineAngle, 1.0f);
  }
  finder->compute();
  REQUIRE(finder->nresults() >= 1);
  CHECK(finder->weight(0) > 10);
}

TEST_CASE("LinesP2D: empty input gives no results") {
  auto finder = makeFinder();
  finder->compute();
  CHECK(finder->nresults() == 0);
}

TEST_CASE("LinesP2D: weight tracks number of inliers") {
  auto finder = makeFinder();
  finder->set_minweight(0.0f);
  for (float x = 0; x < 50; x += 5.0f) {
    finder->add_ipoint(x, 100.0f, 0.0f, 1.0f);
  }
  finder->compute();
  REQUIRE(finder->nresults() >= 1);
  CHECK(finder->weight(0) >= 5.0f);
}

TEST_CASE("LinesP2D: maxresults limits the result count") {
  auto finder = makeFinder();
  finder->set_maxresults(3);
  for (float x = 0; x < 50; x += 5.0f) {
    finder->add_ipoint(x, 100.0f, 0.0f, 1.0f);
    finder->add_ipoint(100.0f, x, kHalfPi, 1.0f);
  }
  finder->compute();
  CHECK(finder->nresults() <= 3);
}

TEST_CASE("LinesP2D: nmatches reports inlier count") {
  auto finder = makeFinder();
  for (float x = 0; x < 50; x += 5.0f) {
    finder->add_ipoint(x, 100.0f, 0.0f, 1.0f);
  }
  finder->compute();
  REQUIRE(finder->nresults() >= 1);
  CHECK(finder->nmatches(0) >= 5);
}

TEST_CASE("LinesP2D: clear_ipoints resets") {
  auto finder = makeFinder();
  for (float x = 0; x < 50; x += 5.0f) {
    finder->add_ipoint(x, 100.0f, 0.0f, 1.0f);
  }
  finder->clear_ipoints();
  finder->compute();
  CHECK(finder->nresults() == 0);
}
