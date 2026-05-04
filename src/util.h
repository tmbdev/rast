// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#ifndef RAST_SRC_UTIL_H_
#define RAST_SRC_UTIL_H_

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <sys/time.h>

#include "vec2.h"

template <class T>
inline T sqr(T x) {
  return x * x;
}

constexpr float kPi = static_cast<float>(M_PI);
constexpr float kTwoPi = 2.0f * kPi;

inline float normangleOf(float a) {
  while (a < 0) a += kTwoPi;
  while (a >= kTwoPi) a -= kTwoPi;
  return a;
}

inline float normalize_orientation(float a) {
  while (a < 0) a += kPi;
  while (a >= kPi) a -= kPi;
  return a;
}

inline float normalize_angle_centered(float a) {
  while (a < -kPi) a += kTwoPi;
  while (a >= kPi) a -= kTwoPi;
  return a;
}

inline double urand(double low, double high) { return drand48() * (high - low) + low; }

inline vec2 randomUniformVectorFromCircle(float epsilon) {
  for (;;) {
    vec2 v(static_cast<float>(urand(-1.0, 1.0)), static_cast<float>(urand(-1.0, 1.0)));
    if (v.magnitude() < 1.0f) return v * epsilon;
  }
}

// Environment-variable parameter readers. Per the project's coding guidelines,
// reading lowercase env-vars at CLI invocation time (e.g. `verbose=1 ./rast ...`)
// is the chosen mechanism for simple programs. These helpers wrap `std::getenv`
// with a default value and optional logging via `verbose_params`.

inline int igetenv(const char *name, int dflt) {
  int result = std::getenv(name) ? std::atoi(std::getenv(name)) : dflt;
  int where = 0;
  if (std::strcmp(name, "verbose_params") != 0) where = igetenv("verbose_params", 0);
  switch (where) {
    case 1: std::fprintf(stdout, "__param__ %s = %d\n", name, result); break;
    case 2: std::fprintf(stderr, "__param__ %s = %d\n", name, result); break;
    default:;
  }
  return result;
}

inline float fgetenv(const char *name, float dflt) {
  float result = std::getenv(name) ? float(std::atof(std::getenv(name))) : dflt;
  int where = igetenv("verbose_params", 0);
  switch (where) {
    case 1: std::fprintf(stdout, "__param__ %s = %g\n", name, double(result)); break;
    case 2: std::fprintf(stderr, "__param__ %s = %g\n", name, double(result)); break;
    default:;
  }
  return result;
}

inline double dgetenv(const char *name, double dflt) {
  double result = std::getenv(name) ? std::atof(std::getenv(name)) : dflt;
  int where = igetenv("verbose_params", 0);
  switch (where) {
    case 1: std::fprintf(stdout, "__param__ %s = %g\n", name, result); break;
    case 2: std::fprintf(stderr, "__param__ %s = %g\n", name, result); break;
    default:;
  }
  return result;
}

inline const char *sgetenv(const char *name, const char *dflt) {
  const char *result = std::getenv(name) ? std::getenv(name) : dflt;
  int where = igetenv("verbose_params", 0);
  switch (where) {
    case 1: std::fprintf(stdout, "__param__ %s = %s\n", name, result); break;
    case 2: std::fprintf(stderr, "__param__ %s = %s\n", name, result); break;
    default:;
  }
  return result;
}

inline int clamp(int i, int n) {
  if (i < 0) return 0;
  if (i >= n) return n - 1;
  return i;
}

inline int mkseed() {
  int seed = igetenv("seed", -1);
  if (seed != -1) return seed;
  return 0;
}

inline double now() {
  struct timeval tv;
  gettimeofday(&tv, nullptr);
  return tv.tv_sec + 1e-6 * tv.tv_usec;
}

#endif  // RAST_SRC_UTIL_H_
