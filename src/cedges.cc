// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

/* Compile with: g++ -std=c++23 -DMAIN cedges.cc -o cedges */

/*
    A self-contained implementation of something like the Canny edge detector
    in C++23. Reads PGM and PPM files and outputs edge maps, edge chains,
    sampled edge chains, and polygonal approximations.

    The implementation uses a FIR filter and difference-of-gaussians instead
    of Canny's IIR filter; I have code for the latter as well but haven't
    found it to be noticeably better, and its behavior around image borders
    is less predictable.

    Compile with -DMAIN to get the main program. Without that, you get a
    class iupr_cedges::CEdges (factory: makeEdgeDetector()) that you can
    call from within your own programs.

    If you find this code useful and use it in some research leading to a
    publication, I would appreciate being mentioned in the Acknowledgements
    section.
*/

/*
Usage:
   cedges < input > output

Input image must be in .pgm or .ppm format

Parameters (via environment):
   esx,esy    3.0       Gaussian sigma
   efrac      0.3       Canny fraction
   elow       2.0       Canny low threshold
   ehigh      4.0       Canny high threshold
   eminlength 5.0       min length of segment in polygonal approximation
   emaxdist   1.5       max deviation in polygonal approximation
   espacing   4         interval at which edges are sampled
   eformat    segments  segments: output line segments
                        sampled: output edge samples
                        poly: output polygons
                        chain: output chains of edge pixels
                        map: output edge map in PBM format
*/

#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <vector>

// Keep the legacy paren-call syntax (arr(i,j)) on top of mdarray. The Kokkos
// reference implementation defaults to bracket-only when C++23's
// multidimensional subscript is detected; opt back into both.
#define MDSPAN_USE_PAREN_OPERATOR 1
#include <experimental/mdarray>
#include <experimental/mdspan>

#include "cedges.h"

namespace iupr_cedges {

namespace stdex = std::experimental;

template <class T>
inline T sqr(T x) {
  return x * x;
}

template <class T>
using Image = stdex::mdarray<T, stdex::dextents<int, 2>>;

using ByteImage = Image<unsigned char>;
using FloatImage = Image<float>;

template <class T>
static void fill_image(Image<T>& a, T value) {
  std::fill(a.data(), a.data() + a.size(), value);
}

template <class S, class D>
static void copy_image(Image<D>& dst, const Image<S>& src) {
  dst = Image<D>(src.extent(0), src.extent(1));
  const int n = static_cast<int>(src.size());
  for (int i = 0; i < n; i++) {
    dst.data()[i] = static_cast<D>(src.data()[i]);
  }
}

struct vec2 {
  enum { N = 2 };
  float data[2];
  explicit vec2(float v0 = 0.0f, float v1 = 0.0f) {
    data[0] = v0;
    data[1] = v1;
  }
  int length() const { return N; }
  float at(int i) const {
    assert(static_cast<unsigned>(i) < static_cast<unsigned>(N));
    return data[i];
  }
  float operator()(int i) const {
    assert(static_cast<unsigned>(i) < static_cast<unsigned>(N));
    return data[i];
  }
  float& operator()(int i) {
    assert(static_cast<unsigned>(i) < static_cast<unsigned>(N));
    return data[i];
  }
  float operator[](int i) const {
    assert(static_cast<unsigned>(i) < static_cast<unsigned>(N));
    return data[i];
  }
  float& operator[](int i) {
    assert(static_cast<unsigned>(i) < static_cast<unsigned>(N));
    return data[i];
  }
  vec2 operator+(const vec2& other) const { return vec2(at(0) + other(0), at(1) + other(1)); }
  vec2 operator-(const vec2& other) const { return vec2(at(0) - other(0), at(1) - other(1)); }
  float operator*(const vec2& other) const { return at(0) * other(0) + at(1) * other(1); }
  vec2 operator*(float scale) const { return vec2(at(0) * scale, at(1) * scale); }
  vec2 operator/(float scale) const { return vec2(at(0) / scale, at(1) / scale); }
  float magnitude() const { return std::sqrt(sqr(data[0]) + sqr(data[1])); }
  float angle() const { return std::atan2(data[1], data[0]); }
  float magnitude_squared() const { return sqr(data[0]) + sqr(data[1]); }
  vec2 normalized() const { return operator*(1.0f / magnitude()); }
  vec2 normal() const { return vec2(-data[1], data[0]); }
  float distance(const vec2& b) const { return (*this - b).magnitude(); }
  vec2 cmul(const vec2& b) const {
    const vec2& a = *this;
    return vec2(a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]);
  }
  vec2 cdiv(const vec2& b) const {
    const vec2& a = *this;
    float n = sqr(b[0]) + sqr(b[1]);
    return vec2((a[0] * b[0] + a[1] * b[1]) / n, (a[1] * b[0] - a[0] * b[1]) / n);
  }
};

static int scanint(std::istream& in) {
  int value;
  if (!(in >> value)) throw std::runtime_error("read_pnm: number format error in image");
  return value;
}

static int getbyte(std::istream& in) {
  int value = in.get();
  if (value < 0) throw std::runtime_error("read_pnm: image short");
  return value;
}

static void read_pnm(ByteImage& image, std::istream& in) {
  int cp = in.get();
  if (cp != 'P') throw std::runtime_error("read_pnm: unknown format");
  int ctype = in.get() - '0';
  if (ctype != 2 && ctype != 3 && ctype != 5 && ctype != 6)
    throw std::runtime_error("read_pnm: cannot handle this type");
  int params[3] = {0, 0, 0};
  int nparams = 0;
  while (nparams < 3) {
    int c = in.get();
    if (c < 0) throw std::runtime_error("read_pnm: unexpected eof");
    if (std::isdigit(static_cast<unsigned char>(c))) {
      in.unget();
      if (!(in >> params[nparams])) throw std::runtime_error("read_pnm: bad number format");
      nparams++;
      c = in.get();
      if (c < 0 || !std::isspace(static_cast<unsigned char>(c)))
        throw std::runtime_error("read_pnm: bad header format");
      continue;
    }
    if (c == '#') {
      for (;;) {
        c = in.get();
        if (c == '\n') break;
        if (c < 0) throw std::runtime_error("read_pnm: unexpected eof");
      }
    }
  }
  image = ByteImage(params[0], params[1]);
  for (int j = image.extent(1) - 1; j >= 0; j--) {
    for (int i = 0; i < image.extent(0); i++) {
      int value = 0;
      switch (ctype) {
      case 2:
        value = scanint(in);
        break;
      case 3:
        value = (scanint(in) + scanint(in) + scanint(in)) / 3;
        break;
      case 5:
        value = getbyte(in);
        break;
      case 6:
        value = (getbyte(in) + getbyte(in) + getbyte(in)) / 3;
        break;
      default:
        throw std::runtime_error("bad type");
      }
      image(i, j) = static_cast<unsigned char>(value);
    }
  }
}

static void write_pnm(std::ostream& out, const ByteImage& image) {
  out << "P5\n" << image.extent(0) << " " << image.extent(1) << " 255\n";
  for (int j = image.extent(1) - 1; j >= 0; j--) {
    for (int i = 0; i < image.extent(0); i++) {
      out.put(static_cast<char>(image(i, j)));
    }
  }
}

static void read_pnm(ByteImage& image, const char* file) {
  std::ifstream stream(file, std::ios::binary);
  if (!stream) throw std::runtime_error("read_pnm: open failed");
  read_pnm(image, stream);
}

static void write_pnm(const char* file, const ByteImage& image) {
  std::ofstream stream(file, std::ios::binary);
  if (!stream) throw std::runtime_error("write_pnm: open failed");
  write_pnm(stream, image);
}

static void gauss1d(std::vector<float>& out, const std::vector<float>& in, float sigma) {
  out.assign(in.size(), 0.0f);
  const int range = 1 + static_cast<int>(3.0f * sigma);
  std::vector<float> mask(static_cast<std::size_t>(2 * range + 1));
  for (int i = 0; i <= range; i++) {
    const float y = std::exp(-static_cast<float>(i * i) / (2.0f * sigma * sigma));
    mask[static_cast<std::size_t>(range + i)] = y;
    mask[static_cast<std::size_t>(range - i)] = y;
  }
  float total = 0.0f;
  for (float v : mask) total += v;
  for (float& v : mask) v /= total;

  const int n = static_cast<int>(in.size());
  const int m = static_cast<int>(mask.size());
  for (int i = 0; i < n; i++) {
    float sum = 0.0f;
    for (int j = 0; j < m; j++) {
      int index = i + j - range;
      if (index < 0) index = 0;
      if (index >= n) index = n - 1;
      sum += in[static_cast<std::size_t>(index)] * mask[static_cast<std::size_t>(j)];
    }
    out[static_cast<std::size_t>(i)] = sum;
  }
}

static void gauss2d(FloatImage& a, float sx, float sy) {
  const int w = a.extent(0);
  const int h = a.extent(1);
  std::vector<float> r, s;
  // smooth along axis 1 for each i (the dim-1 stride is 1)
  r.resize(static_cast<std::size_t>(h));
  for (int i = 0; i < w; i++) {
    for (int j = 0; j < h; j++) r[static_cast<std::size_t>(j)] = a(i, j);
    gauss1d(s, r, sy);
    for (int j = 0; j < h; j++) a(i, j) = s[static_cast<std::size_t>(j)];
  }
  // smooth along axis 0 for each j
  r.resize(static_cast<std::size_t>(w));
  for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) r[static_cast<std::size_t>(i)] = a(i, j);
    gauss1d(s, r, sx);
    for (int i = 0; i < w; i++) a(i, j) = s[static_cast<std::size_t>(i)];
  }
}

static int isign(float x) { return x >= 0.0f ? 1 : -1; }

static void nonmaxsup(ByteImage& out, const FloatImage& gradm, const FloatImage& gradx,
                      const FloatImage& grady) {
  const int w = gradm.extent(0);
  const int h = gradm.extent(1);
  out = ByteImage(w, h);
  fill_image(out, static_cast<unsigned char>(0));
  for (int i = 1; i < w - 1; i++) {
    for (int j = 1; j < h - 1; j++) {
      float dx = gradx(i, j);
      float ux = std::fabs(dx);
      float dy = grady(i, j);
      float uy = std::fabs(dy);
      int bx = isign(dx);
      int by = isign(dy);
      int ax = bx * (ux > uy);
      int ay = by * (ux <= uy);
      float vx, vy;
      if (ax) {
        vy = ux;
        vx = uy;
      } else {
        vx = ux;
        vy = uy;
      }
      float c = gradm(i, j);
      float u = gradm(i - ax, j - ay);
      float d = gradm(i - bx, j - by);
      if (vy * c <= (vx * d + (vy - vx) * u)) continue;
      u = gradm(i + ax, j + ay);
      d = gradm(i + bx, j + by);
      if (vy * c < (vx * d + (vy - vx) * u)) continue;
      out(i, j) = 255;
    }
  }
}

static float masked_fractile(const FloatImage& gradm, const ByteImage& mask, float frac,
                             int bins = 1000) {
  const bool use_mask = mask.size() > 0;
  std::vector<int> hist(static_cast<std::size_t>(bins), 0);
  float minv = 1e30f, maxv = -1e30f;
  int count = 0;
  const int total = static_cast<int>(gradm.size());
  const float* g = gradm.data();
  const unsigned char* m = mask.data();
  for (int i = 0; i < total; i++) {
    if (use_mask && !m[i]) continue;
    count++;
    if (maxv < g[i]) maxv = g[i];
    if (minv > g[i]) minv = g[i];
  }
  if (count < 2) return minv;
  if (maxv == minv) return minv;
  const int limit = static_cast<int>(frac * static_cast<float>(count));
  const float scale = static_cast<float>(bins) / (maxv - minv);
  for (int i = 0; i < total; i++) {
    if (use_mask && !m[i]) continue;
    int bin = static_cast<int>(scale * (g[i] - minv));
    hist[static_cast<std::size_t>(std::min(bins - 1, bin))]++;
  }
  int i = 0, sum = 0;
  for (; i < bins && sum < limit; i++) {
    sum += hist[static_cast<std::size_t>(i)];
  }
  return (maxv - minv) * static_cast<float>(i) / static_cast<float>(bins) + minv;
}

static void masked_fill(const ByteImage& in, ByteImage& out, int x, int y) {
  const int w = in.extent(0);
  const int h = in.extent(1);
  if (x < 0 || x >= w) return;
  if (y < 0 || y >= h) return;
  if (!in(x, y)) return;
  if (out(x, y)) return;
  out(x, y) = 1;
  masked_fill(in, out, x + 1, y);
  masked_fill(in, out, x, y + 1);
  masked_fill(in, out, x - 1, y);
  masked_fill(in, out, x, y - 1);
  masked_fill(in, out, x + 1, y + 1);
  masked_fill(in, out, x - 1, y + 1);
  masked_fill(in, out, x - 1, y + 1);
  masked_fill(in, out, x + 1, y - 1);
}

static void hysteresis_thresholding(ByteImage& out, const FloatImage& gradm, const ByteImage& mask,
                                    float tlow, float thigh) {
  const int w = gradm.extent(0);
  const int h = gradm.extent(1);
  ByteImage low(w, h);
  fill_image(low, static_cast<unsigned char>(0));
  out = ByteImage(w, h);
  fill_image(out, static_cast<unsigned char>(0));
  const bool use_mask = mask.size() > 0;
  for (int i = 0; i < w; i++)
    for (int j = 0; j < h; j++) {
      if (use_mask && !mask(i, j)) continue;
      if (gradm(i, j) > tlow) low(i, j) = 1;
    }
  for (int i = 0; i < w; i++)
    for (int j = 0; j < h; j++) {
      if (gradm(i, j) > thigh) masked_fill(low, out, i, j);
    }
}

static void thin(ByteImage& uci) {
  enum { OFF = 0, ON = 1, SKEL = 2, DEL = 3 };

  static const unsigned char ttable[256] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, /* 00 */
      0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, /* 10 */
      0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, /* 20 */
      0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, /* 30 */
      0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, /* 40 */
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, /* 50 */
      0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, /* 60 */
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, /* 70 */
      0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, /* 80 */
      1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, /* 90 */
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, /* a0 */
      1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, /* b0 */
      0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, /* c0 */
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, /* d0 */
      0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, /* e0 */
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0  /* f0 */
  };

  static const int nx[] = {1, 1, 0, -1, -1, -1, 0, 1};
  static const int ny[] = {0, 1, 1, 1, 0, -1, -1, -1};

  const int w = uci.extent(0) - 1;
  const int h = uci.extent(1) - 1;

  const int total = static_cast<int>(uci.size());
  unsigned char* d = uci.data();
  for (int i = 0; i < total; i++) {
    d[i] = d[i] ? ON : OFF;
  }

  int flag;
  do {
    flag = 0;
    for (int j = 0; j < 8; j += 2) {
      for (int x = 1; x < w; x++)
        for (int y = 1; y < h; y++) {
          if (uci(x, y) != ON) continue;
          if (uci(x + nx[j], y + ny[j]) != OFF) continue;
          int b = 0;
          for (int i = 7; i >= 0; i--) {
            b <<= 1;
            b |= (uci(x + nx[i], y + ny[i]) != OFF);
          }
          if (ttable[b])
            uci(x, y) = SKEL;
          else {
            uci(x, y) = DEL;
            flag = 1;
          }
        }
      if (!flag) continue;
      for (int x = 1; x < w; x++)
        for (int y = 1; y < h; y++)
          if (uci(x, y) == DEL) uci(x, y) = OFF;
    }
  } while (flag);

  for (int i = 0; i < total; i++) {
    d[i] = (d[i] == SKEL) ? 255 : 0;
  }
}

static float point_line_dist(vec2 p, vec2 a, vec2 b) {
  vec2 delta = b - a;
  float mag = delta.magnitude();
  // if the distance is small, just return the point distance;
  // that's the right thing for approx_chain
  if (mag < 1e-4f) return a.distance(p);
  vec2 normal = delta.normal() / mag;
  float offset = normal * a;
  return std::fabs(normal * p - offset);
}

static void approx_chain(std::vector<int>& breaks, const std::vector<vec2>& chain, int i0, int i1,
                         float maxdist) {
  float md = 0.0f;
  int mi = -1;
  vec2 a = chain[static_cast<std::size_t>(i0)];
  vec2 b = chain[static_cast<std::size_t>(i1)];
  for (int i = i0; i <= i1; i++) {
    float d = point_line_dist(chain[static_cast<std::size_t>(i)], a, b);
    if (d <= md) continue;
    md = d;
    mi = i;
  }
  if (md < maxdist) return;
  assert(mi != i0 && mi != i1);
  approx_chain(breaks, chain, i0, mi, maxdist);
  breaks.push_back(mi);
  approx_chain(breaks, chain, mi, i1, maxdist);
}

struct ChainTracer {
  enum { OFF = 0, ON = 1, DONE = 2 };

  ByteImage bi;
  int sx = 0;
  int sy = 0;
  int x = 0;
  int y = 0;
  int w = 0;
  int h = 0;

  int count_neighbors() const {
    int nn = 0;
    if (bi(x + 1, y)) nn++;
    if (bi(x + 1, y + 1)) nn++;
    if (bi(x, y + 1)) nn++;
    if (bi(x - 1, y + 1)) nn++;
    if (bi(x - 1, y)) nn++;
    if (bi(x - 1, y - 1)) nn++;
    if (bi(x, y - 1)) nn++;
    if (bi(x + 1, y - 1)) nn++;
    return nn;
  }

  bool nextpixel() {
    if (bi(x + 1, y) == ON) {
      x = x + 1;
      return true;
    }
    if (bi(x + 1, y + 1) == ON) {
      x = x + 1;
      y = y + 1;
      return true;
    }
    if (bi(x, y + 1) == ON) {
      y = y + 1;
      return true;
    }
    if (bi(x - 1, y + 1) == ON) {
      x = x - 1;
      y = y + 1;
      return true;
    }
    if (bi(x - 1, y) == ON) {
      x = x - 1;
      return true;
    }
    if (bi(x - 1, y - 1) == ON) {
      x = x - 1;
      y = y - 1;
      return true;
    }
    if (bi(x, y - 1) == ON) {
      y = y - 1;
      return true;
    }
    if (bi(x + 1, y - 1) == ON) {
      x = x + 1;
      y = y - 1;
      return true;
    }
    return false;
  }

  bool nextstart() {
    for (; sx < w; sx++)
      for (sy = 0; sy < h; sy++) {
        if (!bi(sx, sy)) continue;
        if (bi(sx, sy) == DONE) continue;
        x = sx;
        y = sy;
        return true;
      }
    return false;
  }

  void set_image(const ByteImage& image) {
    const int iw = image.extent(0);
    const int ih = image.extent(1);
    bi = ByteImage(iw, ih);
    const int total = static_cast<int>(image.size());
    for (int i = 0; i < total; i++) bi.data()[i] = image.data()[i] ? ON : OFF;
    for (int i = 0; i < iw; i++) {
      bi(i, 0) = OFF;
      bi(i, ih - 1) = OFF;
    }
    for (int j = 0; j < ih; j++) {
      bi(0, j) = OFF;
      bi(iw - 1, j) = OFF;
    }
    w = iw - 1;
    h = ih - 1;
    x = 0;
    y = 0;
    sx = 0;
    sy = 0;
  }

  bool get_chain(std::vector<vec2>& points, bool close = false) {
    points.clear();
    if (!nextstart()) return false;
    do {
      points.push_back(vec2(static_cast<float>(x), static_cast<float>(y)));
      bi(x, y) = DONE;
    } while (nextpixel());
    x = sx;
    y = sy;
    if (nextpixel()) {
      // sweep up the other direction (if any)
      std::reverse(points.begin(), points.end());
      do {
        points.push_back(vec2(static_cast<float>(x), static_cast<float>(y)));
        bi(x, y) = DONE;
      } while (nextpixel());
    }
    if (close && !points.empty() && points.front().distance(points.back()) < 2.0f) {
      // close circular chains
      points.push_back(points.front());
    }
    return true;
  }

  std::vector<int> breaks;
  std::vector<vec2> chain;

  bool get_poly(std::vector<vec2>& poly, float maxdist = 1.0f, bool close = false) {
    poly.clear();
    breaks.clear();
    if (!get_chain(chain, close)) return false;
    if (chain.size() < 3) {
      for (const vec2& v : chain) poly.push_back(v);
    } else {
      breaks.push_back(0);
      approx_chain(breaks, chain, 0, static_cast<int>(chain.size()) - 1, maxdist);
      breaks.push_back(static_cast<int>(chain.size()) - 1);
      for (int i : breaks) poly.push_back(chain[static_cast<std::size_t>(i)]);
    }
    return true;
  }

  bool started() const { return bi.size() > 0; }

  void clear() {
    bi = ByteImage();
    breaks.clear();
    chain.clear();
  }
};

struct CEdges : EdgeDetector {
  float p_sx = 3.0f;
  float p_sy = 3.0f;
  float p_frac = 0.3f;
  float p_tlow = 2.0f;
  float p_thigh = 4.0f;
  float p_minlength = 5.0f;
  float p_maxdist = 0.5f;
  float noise = 0.0f;

  FloatImage image;
  ByteImage uedges;
  FloatImage smoothed;
  FloatImage gradm;
  FloatImage gradx;
  FloatImage grady;
  ByteImage edges;
  ChainTracer tracer;

  std::vector<vec2> chain;
  std::vector<int> breaks;

  struct Sample {
    vec2 c;
    vec2 g;
    int n;
  };
  std::vector<Sample> samples;

  void set_gauss(float sx, float sy) override {
    p_sx = sx;
    p_sy = sy;
  }
  void set_noise(float frac, float low, float high) override {
    p_frac = frac;
    p_tlow = low;
    p_thigh = high;
  }
  void set_poly(float minlength, float maxdist) override {
    p_minlength = minlength;
    p_maxdist = maxdist;
  }

  void clear() override {
    image = FloatImage();
    uedges = ByteImage();
    smoothed = FloatImage();
    gradm = FloatImage();
    gradx = FloatImage();
    grady = FloatImage();
    edges = ByteImage();
    tracer.clear();
    chain.clear();
    breaks.clear();
    samples.clear();
    noise = 0.0f;
  }

  void load_pnm(std::istream& in) {
    ByteImage bimage;
    read_pnm(bimage, in);
    copy_image(image, bimage);
  }
  void load_pnm(const char* file) override {
    ByteImage bimage;
    read_pnm(bimage, file);
    copy_image(image, bimage);
  }
  void save_pnm(std::ostream& out) const { write_pnm(out, edges); }
  void save_pnm(const char* file) override { write_pnm(file, edges); }

  int dim(int i) const override { return image.extent(i); }

  void set_image(unsigned char* p, int w, int h) override {
    image = FloatImage(w, h);
    const int n = static_cast<int>(image.size());
    for (int i = 0; i < n; i++) image.data()[i] = static_cast<float>(p[i]);
  }
  void set_pixmap(unsigned char* p, int w, int h) override {
    image = FloatImage(w, h);
    for (int yy = h - 1; yy >= 0; yy--)
      for (int xx = 0; xx < w; xx++) image(xx, yy) = static_cast<float>(*p++);
  }

  void compute() override {
    const int w = image.extent(0);
    const int h = image.extent(1);

    copy_image(smoothed, image);
    gauss2d(smoothed, p_sx, p_sy);

    gradm = FloatImage(w, h);
    gradx = FloatImage(w, h);
    grady = FloatImage(w, h);
    fill_image(gradm, 0.0f);
    fill_image(gradx, 0.0f);
    fill_image(grady, 0.0f);
    for (int i = w - 2; i >= 0; i--)
      for (int j = h - 2; j >= 0; j--) {
        float v = smoothed(i, j);
        float dx = smoothed(i + 1, j) - v;
        float dy = smoothed(i, j + 1) - v;
        gradx(i, j) = dx;
        grady(i, j) = dy;
        gradm(i, j) = std::sqrt(sqr(dx) + sqr(dy));
      }

    nonmaxsup(uedges, gradm, gradx, grady);

    {
      const int n = static_cast<int>(uedges.size());
      unsigned char* ud = uedges.data();
      for (int i = 0; i < n; i++)
        if (ud[i]) ud[i] = 255;
    }
    thin(uedges);

    ByteImage temp;
    noise = masked_fractile(gradm, temp, p_frac);
    const float low = p_tlow * noise;
    const float high = p_thigh * noise;
    hysteresis_thresholding(edges, gradm, uedges, low, high);
    {
      const int n = static_cast<int>(edges.size());
      unsigned char* ed = edges.data();
      for (int i = 0; i < n; i++) ed[i] = ed[i] ? 255 : 0;
    }

    tracer.set_image(edges);
  }

  void get_eimage(unsigned char* p, int w, int h) const override {
    if (edges.extent(0) != w || edges.extent(1) != h)
      throw std::length_error("output image has the wrong size");
    const int n = static_cast<int>(edges.size());
    for (int i = 0; i < n; i++) p[i] = edges.data()[i];
  }
  void get_epixmap(unsigned char* out, int w, int h) const override {
    if (edges.extent(0) != w || edges.extent(1) != h)
      throw std::length_error("output image has the wrong size");
    for (int yy = h - 1; yy >= 0; yy--)
      for (int xx = 0; xx < w; xx++) *out++ = edges(xx, yy);
  }
  float gradient_magnitude(int x, int y) const override { return gradm(x, y); }
  float gradient_angle(int x, int y) const override {
    return std::atan2(grady(x, y), gradx(x, y));
  }

  bool nextchain() override {
    if (!tracer.get_chain(chain)) return false;
    breaks.clear();
    breaks.push_back(0);
    approx_chain(breaks, chain, 0, static_cast<int>(chain.size()) - 1, p_maxdist);
    breaks.push_back(static_cast<int>(chain.size()) - 1);
    return true;
  }
  int npoints() const override { return static_cast<int>(chain.size()); }
  void point(int index, float& x, float& y) const override {
    x = chain[static_cast<std::size_t>(index)][0];
    y = chain[static_cast<std::size_t>(index)][1];
  }
  int nsegments() const override { return static_cast<int>(breaks.size()) - 1; }
  void segment(int i, float& x0, float& y0, float& x1, float& y1, float& angle, float& magnitude,
               int& n) const override {
    i++;
    int i0 = breaks[static_cast<std::size_t>(i - 1)];
    int i1 = breaks[static_cast<std::size_t>(i)];
    vec2 a = chain[static_cast<std::size_t>(i0)];
    vec2 b = chain[static_cast<std::size_t>(i1)];
    vec2 g(0.0f, 0.0f);
    for (int j = i0; j <= i1; j++) {
      int xx = static_cast<int>(chain[static_cast<std::size_t>(j)][0]);
      int yy = static_cast<int>(chain[static_cast<std::size_t>(j)][1]);
      g = g + vec2(gradx(xx, yy), grady(xx, yy));
    }
    x0 = a[0];
    y0 = a[1];
    x1 = b[0];
    y1 = b[1];
    magnitude = g.magnitude();
    angle = g.angle();
    while (angle < 0.0f) angle += 2.0f * static_cast<float>(M_PI);
    n = i1 - i0 + 1;
  }

  void sampleat(int spacing) {
    samples.clear();
    for (std::size_t k = 1; k < breaks.size(); k++) {
      const int start = breaks[k - 1];
      const int end = breaks[k];
      for (int j = start + spacing / 2; j < end; j += spacing) {
        const int j0 = std::max(0, j - spacing / 2);
        const int j1 = std::min(end, j + spacing / 2);
        vec2 c(0.0f, 0.0f);
        vec2 g(0.0f, 0.0f);
        for (int kk = j0; kk <= j1; kk++) {
          int xx = static_cast<int>(chain[static_cast<std::size_t>(kk)][0]);
          int yy = static_cast<int>(chain[static_cast<std::size_t>(kk)][1]);
          c = c + vec2(static_cast<float>(xx), static_cast<float>(yy));
          g = g + vec2(gradx(xx, yy), grady(xx, yy));
        }
        c = c / static_cast<float>(j1 - j0 + 1);
        Sample s;
        s.c = c;
        s.g = g;
        s.n = j1 - j0 + 1;
        samples.push_back(s);
      }
    }
  }
  int nsamples() const { return static_cast<int>(samples.size()); }
  void sample(int i, float& x, float& y, float& angle, float& mag, int& n_pts) const {
    const Sample& s = samples[static_cast<std::size_t>(i)];
    x = s.c(0);
    y = s.c(1);
    angle = std::atan2(s.g(1), s.g(0));
    mag = std::hypot(s.g(0), s.g(1));
    n_pts = s.n;
  }
};

std::unique_ptr<EdgeDetector> makeEdgeDetector() { return std::make_unique<CEdges>(); }

}  // namespace iupr_cedges

#ifdef MAIN

namespace iupr_cedges {

inline int igetenv(const char* name, int dflt) {
  const char* env = std::getenv(name);
  int result = env ? std::atoi(env) : dflt;
  int where = 0;
  if (std::strcmp(name, "verbose_params") != 0) where = igetenv("verbose_params", 0);
  switch (where) {
  case 1:
    std::cout << "__param__ " << name << " = " << result << "\n";
    break;
  case 2:
    std::cerr << "__param__ " << name << " = " << result << "\n";
    break;
  default:;
  }
  return result;
}

inline float fgetenv(const char* name, float dflt) {
  const char* env = std::getenv(name);
  float result = env ? static_cast<float>(std::atof(env)) : dflt;
  int where = igetenv("verbose_params", 0);
  switch (where) {
  case 1:
    std::cout << "__param__ " << name << " = " << result << "\n";
    break;
  case 2:
    std::cerr << "__param__ " << name << " = " << result << "\n";
    break;
  default:;
  }
  return result;
}

inline const char* sgetenv(const char* name, const char* dflt) {
  const char* env = std::getenv(name);
  const char* result = env ? env : dflt;
  int where = igetenv("verbose_params", 0);
  switch (where) {
  case 1:
    std::cout << "__param__ " << name << " = " << result << "\n";
    break;
  case 2:
    std::cerr << "__param__ " << name << " = " << result << "\n";
    break;
  default:;
  }
  return result;
}

const char* usage =
    "Usage:\n"
    "   cedges < input > output\n"
    "\n"
    "Input image must be in .pgm or .ppm format\n"
    "\n"
    "Parameters (via environment):\n"
    "   esx,esy    3.0       Gaussian sigma\n"
    "   efrac      0.3       Canny fraction\n"
    "   elow       2.0       Canny low threshold\n"
    "   ehigh      4.0       Canny high threshold\n"
    "   eminlength 5.0       min length of segment in polygonal approximation\n"
    "   emaxdist   1.5       max deviation in polygonal approximation\n"
    "   espacing   4         interval at which edges are sampled\n"
    "   eformat    segments  segments: output line segments\n"
    "                        sampled: output edge samples\n"
    "                        poly: output polygons\n"
    "                        chain: output chains of edge pixels\n"
    "                        map: output edge map in PBM format\n";

}  // namespace iupr_cedges

int main(int argc, char** argv) {
  using namespace iupr_cedges;
  auto cedges = std::make_unique<CEdges>();
  cedges->set_gauss(fgetenv("esx", 3.0f), fgetenv("esy", 3.0f));
  cedges->set_noise(fgetenv("efrac", 0.3f), fgetenv("elow", 2.0f), fgetenv("ehigh", 4.0f));
  const float eminlength = fgetenv("eminlength", 5.0f);
  cedges->set_poly(eminlength, fgetenv("emaxdist", 1.5f));
  const char* p_format = sgetenv("eformat", "segments");
  const int p_spacing = igetenv("espacing", 4);

  if (argc > 2 || (argc == 2 && (!std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "-?"))) ||
      (argc == 1 && isatty(0))) {
    std::cerr << usage;
    return 1;
  }

  try {
    if (argc == 2)
      cedges->load_pnm(argv[1]);
    else
      cedges->load_pnm(std::cin);
  } catch (const std::exception& e) {
    std::cerr << "problem reading image: " << e.what() << "\n";
    return 2;
  } catch (...) {
    std::cerr << "problem reading image\n";
    return 2;
  }

  cedges->compute();

  if (!std::strcmp(p_format, "map")) {
    ByteImage out_image(cedges->dim(0), cedges->dim(1));
    cedges->get_eimage(out_image.data(), out_image.extent(0), out_image.extent(1));
    write_pnm(std::cout, out_image);
  } else if (!std::strcmp(p_format, "chain")) {
    while (cedges->nextchain()) {
      for (int i = 0; i < cedges->npoints(); i++) {
        float x, y;
        cedges->point(i, x, y);
        std::cout << x << " " << y << "\n";
      }
      std::cout << "\n";
    }
  } else if (!std::strcmp(p_format, "poly")) {
    while (cedges->nextchain()) {
      for (int i = 0; i < cedges->nsegments(); i++) {
        float x, y, x1, y1, a, w;
        int n;
        cedges->segment(i, x, y, x1, y1, a, w, n);
        if (i == 0) std::cout << x << " " << y << "\n";
        std::cout << x1 << " " << y1 << "\n";
      }
      std::cout << "\n";
    }
  } else if (!std::strcmp(p_format, "segments")) {
    while (cedges->nextchain()) {
      int count = 0;
      for (int i = 0; i < cedges->nsegments(); i++) {
        float x, y, x1, y1, a, w;
        int n;
        cedges->segment(i, x, y, x1, y1, a, w, n);
        if (std::sqrt(sqr(x1 - x) + sqr(y1 - y)) < eminlength) continue;
        std::cout << x << " " << y << " " << x1 << " " << y1 << "  " << a << " " << w << " " << n
                  << "\n";
        count++;
      }
      if (count > 0) std::cout << "\n";
    }
  } else if (!std::strcmp(p_format, "sampled")) {
    while (cedges->nextchain()) {
      cedges->sampleat(p_spacing);
      for (int i = 0; i < cedges->nsamples(); i++) {
        float x, y, a, w;
        int n;
        cedges->sample(i, x, y, a, w, n);
        std::cout << x << " " << y << "  " << a << " " << w << " " << n << "\n";
      }
      std::cout << "\n";
    }
  } else {
    std::cerr << p_format << ": unknown format\n";
    return 1;
  }
  return 0;
}

#endif
