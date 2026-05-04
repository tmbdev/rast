// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#include <doctest/doctest.h>

#include <cmath>
#include <memory>

#include "rast.h"

namespace {

std::unique_ptr<LinesS2D> makeFinder() {
  std::unique_ptr<LinesS2D> finder(makeLinesS2D());
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

TEST_CASE("LinesS2D: factory yields a usable instance") {
  std::unique_ptr<LinesS2D> finder(makeLinesS2D());
  REQUIRE(finder != nullptr);
}

TEST_CASE("LinesS2D: horizontal segments along y=cy are recovered") {
  auto finder = makeFinder();
  const float cy = 200.0f;
  for (float x = 0; x < 400; x += 20.0f) {
    finder->add_iseg(x, cy, x + 18.0f, cy, /*a=*/0.0f, /*w=*/1.0f);
  }
  finder->compute();
  REQUIRE(finder->nresults() >= 1);
  const float theta = finder->angle(0);
  const float offset = finder->offset(0);
  CHECK(std::cos(theta - static_cast<float>(M_PI_2)) == doctest::Approx(1.0).epsilon(0.05));
  CHECK(offset == doctest::Approx(cy).epsilon(0.5));
}

TEST_CASE("LinesS2D: empty input gives no results") {
  auto finder = makeFinder();
  finder->compute();
  CHECK(finder->nresults() == 0);
}

TEST_CASE("LinesS2D: maxresults limits the result count") {
  auto finder = makeFinder();
  finder->set_maxresults(2);
  for (float x = 0; x < 100; x += 20.0f) {
    finder->add_iseg(x, 100.0f, x + 18.0f, 100.0f, 0.0f, 1.0f);
    finder->add_iseg(50.0f, x, 50.0f, x + 18.0f, M_PI_2, 1.0f);
  }
  finder->compute();
  CHECK(finder->nresults() <= 2);
}

TEST_CASE("LinesS2D: vertical line is recovered") {
  auto finder = makeFinder();
  const float cx = 100.0f;
  for (float y = 0; y < 300; y += 20.0f) {
    finder->add_iseg(cx, y, cx, y + 18.0f, M_PI_2, 1.0f);
  }
  finder->compute();
  REQUIRE(finder->nresults() >= 1);
  const float theta = finder->angle(0);
  const float offset = finder->offset(0);
  CHECK(std::cos(theta - 0.0f) == doctest::Approx(1.0).epsilon(0.05));
  CHECK(offset == doctest::Approx(cx).epsilon(0.5));
}

TEST_CASE("LinesS2D: clear_ipoints resets") {
  auto finder = makeFinder();
  for (float x = 0; x < 100; x += 20.0f) {
    finder->add_iseg(x, 100.0f, x + 18.0f, 100.0f, 0.0f, 1.0f);
  }
  finder->clear_ipoints();
  finder->compute();
  CHECK(finder->nresults() == 0);
}

TEST_CASE("LinesS2D: weight is positive when a line is found") {
  auto finder = makeFinder();
  finder->set_minweight(0.0f);
  for (float x = 0; x < 200; x += 20.0f) {
    finder->add_iseg(x, 50.0f, x + 18.0f, 50.0f, 0.0f, 1.0f);
  }
  finder->compute();
  REQUIRE(finder->nresults() >= 1);
  CHECK(finder->weight(0) > 0.0f);
}
