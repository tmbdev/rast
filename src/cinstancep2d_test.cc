// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#include <doctest/doctest.h>

#include <cmath>
#include <cstdlib>
#include <memory>
#include <stdexcept>

#include "rast.h"

namespace {

std::unique_ptr<InstanceP2D> makeFixedInstance(int seed = 42) {
  srand48(seed);
  std::unique_ptr<InstanceP2D> inst(makeInstanceP2D());
  inst->set_image_size(512);
  inst->set_model_size(100);
  inst->set_nclutter(20);
  inst->set_nmodel_total(10);
  inst->set_nmodel_unoccluded(8);
  inst->set_error(2.0f);
  inst->set_aerror(0.05f);
  inst->set_srange(0.9f, 1.1f);
  return inst;
}

}  // namespace

TEST_CASE("InstanceP2D: factory yields a usable instance") {
  std::unique_ptr<InstanceP2D> inst(makeInstanceP2D());
  REQUIRE(inst != nullptr);
}

TEST_CASE("InstanceP2D: nimage = nmodel_unoccluded + nclutter") {
  auto inst = makeFixedInstance();
  inst->generate();
  CHECK(inst->nimage() == 8 + 20);
  CHECK(inst->nmodel() == 10);
}

TEST_CASE("InstanceP2D: parameters are within configured ranges") {
  auto inst = makeFixedInstance();
  inst->generate();
  const float tx = inst->get_param(0);
  const float ty = inst->get_param(1);
  const float angle = inst->get_param(2);
  const float scale = inst->get_param(3);
  CHECK(tx >= 0.0f);
  CHECK(tx <= 512.0f);
  CHECK(ty >= 0.0f);
  CHECK(ty <= 512.0f);
  CHECK(angle >= 0.0f);
  CHECK(angle <= 2 * M_PI);
  CHECK(scale >= 0.9f);
  CHECK(scale <= 1.1f);
}

TEST_CASE("InstanceP2D: same seed produces same parameters") {
  auto a = makeFixedInstance(123);
  a->generate();
  const float a_tx = a->get_param(0);
  const float a_ty = a->get_param(1);

  auto b = makeFixedInstance(123);
  b->generate();
  CHECK(b->get_param(0) == doctest::Approx(a_tx));
  CHECK(b->get_param(1) == doctest::Approx(a_ty));
}

TEST_CASE("InstanceP2D: different seeds usually produce different parameters") {
  auto a = makeFixedInstance(1);
  a->generate();
  auto b = makeFixedInstance(2);
  b->generate();
  // At least one of (tx, ty, angle, scale) should differ.
  const bool same = (a->get_param(0) == b->get_param(0)) && (a->get_param(1) == b->get_param(1)) &&
                    (a->get_param(2) == b->get_param(2)) && (a->get_param(3) == b->get_param(3));
  CHECK_FALSE(same);
}

TEST_CASE("InstanceP2D: get_param out-of-range throws") {
  auto inst = makeFixedInstance();
  inst->generate();
  CHECK_THROWS_AS(inst->get_param(4), std::out_of_range);
  CHECK_THROWS_AS(inst->get_param(-1), std::out_of_range);
}

TEST_CASE("InstanceP2D: model points fit within the configured model_size") {
  auto inst = makeFixedInstance();
  inst->set_model_size(50);
  inst->generate();
  for (int i = 0; i < inst->nmodel(); i++) {
    float x = 0, y = 0, a = 0;
    inst->get_model(x, y, a, i);
    CHECK(std::abs(x) <= 50.0f);
    CHECK(std::abs(y) <= 50.0f);
  }
}
