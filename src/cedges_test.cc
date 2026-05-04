// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#include <doctest/doctest.h>

#include <cmath>
#include <cstdlib>
#include <memory>
#include <vector>

#include "cedges.h"

using iupr_cedges::EdgeDetector;
using iupr_cedges::makeEdgeDetector;

namespace {

// Build an in-memory grayscale image with a black square on a white background.
// Layout matches set_image: row-major, w*h bytes, origin top-left.
std::vector<unsigned char> makeSquareImage(int w, int h, int box_x0, int box_y0, int box_x1,
                                           int box_y1) {
  std::vector<unsigned char> data(static_cast<std::size_t>(w * h), 255);
  for (int y = box_y0; y < box_y1; y++) {
    for (int x = box_x0; x < box_x1; x++) {
      data[static_cast<std::size_t>(y * w + x)] = 0;
    }
  }
  return data;
}

}  // namespace

TEST_CASE("EdgeDetector: factory yields a usable instance") {
  std::unique_ptr<EdgeDetector> det = makeEdgeDetector();
  REQUIRE(det != nullptr);
}

TEST_CASE("EdgeDetector: dim() reflects set_image() dimensions") {
  auto det = makeEdgeDetector();
  const int w = 64, h = 48;
  auto img = makeSquareImage(w, h, 16, 12, 48, 36);
  det->set_image(img.data(), w, h);
  CHECK(det->dim(0) == w);
  CHECK(det->dim(1) == h);
}

TEST_CASE("EdgeDetector: compute() finds edges around a black square") {
  auto det = makeEdgeDetector();
  det->set_gauss(1.0f, 1.0f);
  det->set_noise(0.3f, 1.5f, 3.0f);
  det->set_poly(2.0f, 1.0f);
  const int w = 64, h = 64;
  auto img = makeSquareImage(w, h, 16, 16, 48, 48);
  det->set_image(img.data(), w, h);
  det->compute();

  std::vector<unsigned char> out(static_cast<std::size_t>(w * h), 0);
  det->get_eimage(out.data(), w, h);

  int edge_count = 0;
  for (unsigned char p : out) {
    if (p) edge_count++;
  }
  // Some edge pixels must show up around the square's perimeter.
  CHECK(edge_count > 0);
}

TEST_CASE("EdgeDetector: get_eimage rejects mismatched dimensions") {
  auto det = makeEdgeDetector();
  const int w = 32, h = 32;
  auto img = makeSquareImage(w, h, 8, 8, 24, 24);
  det->set_image(img.data(), w, h);
  det->compute();

  std::vector<unsigned char> wrong(16 * 16, 0);
  CHECK_THROWS_AS(det->get_eimage(wrong.data(), 16, 16), std::length_error);
}

TEST_CASE("EdgeDetector: nextchain iterates exhaustively") {
  auto det = makeEdgeDetector();
  det->set_gauss(1.0f, 1.0f);
  det->set_noise(0.3f, 1.5f, 3.0f);
  det->set_poly(2.0f, 1.0f);
  const int w = 64, h = 64;
  auto img = makeSquareImage(w, h, 16, 16, 48, 48);
  det->set_image(img.data(), w, h);
  det->compute();

  int chains = 0;
  while (det->nextchain()) {
    chains++;
    REQUIRE(det->npoints() >= 0);
    if (det->npoints() > 0) {
      float x = -1.0f, y = -1.0f;
      det->point(0, x, y);
      CHECK(x >= 0.0f);
      CHECK(y >= 0.0f);
      CHECK(x < static_cast<float>(w));
      CHECK(y < static_cast<float>(h));
    }
    REQUIRE(det->nsegments() >= 0);
  }
  // The square's outline produces at least one chain.
  CHECK(chains >= 1);
}

TEST_CASE("EdgeDetector: gradient_magnitude is non-negative") {
  auto det = makeEdgeDetector();
  const int w = 32, h = 32;
  auto img = makeSquareImage(w, h, 8, 8, 24, 24);
  det->set_image(img.data(), w, h);
  det->compute();
  for (int x = 0; x < w; x++) {
    for (int y = 0; y < h; y++) {
      CHECK(det->gradient_magnitude(x, y) >= 0.0f);
    }
  }
}

TEST_CASE("EdgeDetector: clear() resets the detector for reuse") {
  auto det = makeEdgeDetector();
  const int w = 32, h = 32;
  auto img = makeSquareImage(w, h, 8, 8, 24, 24);
  det->set_image(img.data(), w, h);
  det->compute();
  det->clear();
  CHECK(det->dim(0) == 0);
  CHECK(det->dim(1) == 0);
}
