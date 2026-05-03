// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#ifndef RAST_SRC_VEC2_H_
#define RAST_SRC_VEC2_H_

#include <cassert>
#include <cmath>

// 2D float vector used by the RAST algorithms. There is no standard library
// equivalent, so this stays as a local type.

struct vec2 {
  float data[2]{0.0f, 0.0f};

  vec2() = default;
  explicit vec2(float v0, float v1) {
    data[0] = v0;
    data[1] = v1;
  }

  float operator[](int i) const {
    assert(unsigned(i) < 2);
    return data[i];
  }
  float &operator[](int i) {
    assert(unsigned(i) < 2);
    return data[i];
  }
  float operator()(int i) const { return (*this)[i]; }
  float &operator()(int i) { return (*this)[i]; }

  vec2 operator+(const vec2 &other) const { return vec2{data[0] + other[0], data[1] + other[1]}; }
  vec2 operator-(const vec2 &other) const { return vec2{data[0] - other[0], data[1] - other[1]}; }
  float operator*(const vec2 &other) const { return data[0] * other[0] + data[1] * other[1]; }
  vec2 operator*(float scale) const { return vec2{data[0] * scale, data[1] * scale}; }
  vec2 operator/(float scale) const { return vec2{data[0] / scale, data[1] / scale}; }

  float magnitude_squared() const { return data[0] * data[0] + data[1] * data[1]; }
  float magnitude() const { return std::sqrt(magnitude_squared()); }
  float angle() const { return std::atan2(data[1], data[0]); }
  vec2 normalized() const { return *this * (1.0f / magnitude()); }
  float distance(const vec2 &b) const { return (*this - b).magnitude(); }
};

inline vec2 normal(const vec2 &v) { return vec2{-v[1], v[0]}; }

inline vec2 cmul(const vec2 &a, const vec2 &b) {
  return vec2{a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]};
}

inline vec2 cdiv(const vec2 &a, const vec2 &b) {
  const float n = b[0] * b[0] + b[1] * b[1];
  return vec2{(a[0] * b[0] + a[1] * b[1]) / n, (a[1] * b[0] - a[0] * b[1]) / n};
}

inline float angleOf(const vec2 &v) { return std::atan2(v[1], v[0]); }
inline float distance(const vec2 &a, const vec2 &b) { return (a - b).magnitude(); }
inline float norm(const vec2 &v) { return v.magnitude(); }

#endif  // RAST_SRC_VEC2_H_
