// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>
#include <stdexcept>
#include <vector>

#include "rast.h"
#include "util.h"
#include "vec2.h"

namespace lumo_cinstancep2d {

struct avec {
  vec2 p;
  float a;
};

using Msource = avec;
using Ipoint = avec;

inline int urand48() { return std::abs(int(lrand48())); }

template <class T>
static void shuffle(std::vector<T> &v) {
  for (std::size_t i = 0; i + 1 < v.size(); i++) {
    std::size_t j = urand48() % (v.size() - i) + i;
    if (i != j) std::swap(v[i], v[j]);
  }
}

struct CInstanceP2D : InstanceP2D {
  int image_size{512};
  int model_size{100};

  int nclutter{50};      // configurable via set_nclutter() / CLI nclutter env-var
  int nmodel_total{20};
  int nmodel_unoccluded{10};
  float error{5.0f};
  float aerror{0.1f};
  float minscale{1.0f};
  float maxscale{1.0f};

  vec2 translation;
  float angle{0.0f};
  float scale{1.0f};

  std::vector<Msource> msources;
  std::vector<Ipoint> ipoints;

  CInstanceP2D() = default;

  float get_param(int i) const override {
    switch (i) {
      case 0: return translation[0];
      case 1: return translation[1];
      case 2: return angle;
      case 3: return scale;
      default: throw std::out_of_range("InstanceP2D::get_param: index out of range");
    }
  }

  // urand returns double; cast to float at the boundary.
  static float frand(double low, double high) { return static_cast<float>(urand(low, high)); }

  void generate() override {
    float dx = frand(0.0, image_size);
    float dy = frand(0.0, image_size);
    translation = vec2(dx, dy);
    angle = frand(0.0, 2.0 * M_PI);
    scale = frand(minscale, maxscale);
    vec2 rotation = vec2(scale * std::cos(angle), scale * std::sin(angle));
    msources.clear();
    ipoints.clear();
    for (int i = 0; i < nmodel_total; i++) {
      Msource m;
      m.p = vec2(frand(-model_size, model_size), frand(-model_size, model_size));
      m.a = frand(0.0, 2 * M_PI);
      msources.push_back(m);
    }
    for (int i = 0; i < nmodel_unoccluded; i++) {
      Ipoint p;
      p.p = cmul(rotation, msources[i].p) + translation + randomUniformVectorFromCircle(error);
      p.a = msources[i].a + angle + frand(-aerror, aerror);
      ipoints.push_back(p);
    }
    shuffle(msources);
    for (int i = 0; i < nclutter; i++) {
      Ipoint p;
      p.p = vec2(frand(0, image_size), frand(0, image_size));
      p.a = frand(0.0, 2 * M_PI);
      ipoints.push_back(p);
    }
    shuffle(ipoints);
  }

  void set_image_size(int r) override { image_size = r; }
  void set_model_size(int r) override { model_size = r; }
  void set_nclutter(int v) override { nclutter = v; }
  void set_nmodel_total(int v) override { nmodel_total = v; }
  void set_nmodel_unoccluded(int v) override { nmodel_unoccluded = v; }
  void set_error(float v) override { error = v; }
  void set_aerror(float v) override { aerror = v; }
  void set_srange(float min, float max) override {
    minscale = min;
    maxscale = max;
  }

  int nimage() const override { return int(ipoints.size()); }
  void get_image(float &x, float &y, float &a, int i) const override {
    x = ipoints[i].p[0];
    y = ipoints[i].p[1];
    a = ipoints[i].a;
  }

  int nmodel() const override { return int(msources.size()); }
  void get_model(float &x, float &y, float &a, int i) const override {
    x = msources[i].p[0];
    y = msources[i].p[1];
    a = msources[i].a;
  }
};

}  // namespace lumo_cinstancep2d

std::unique_ptr<InstanceP2D> makeInstanceP2D() {
  return std::make_unique<lumo_cinstancep2d::CInstanceP2D>();
}
