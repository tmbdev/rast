/* Copyright 1992-1995 by Thomas M. Breuel */
/* This code comes without any warranty whatsoever.  Use at your own risk. */
/* Licensed under the terms of the GNU General Public License, Version 2 */
/* Contact info: http://www.lumo.com/ */

/* Compile with: g++ -g -O2 -DMAIN cedges.cc */

/* this is a test */

/*
    A self-contained implementation of something like the Canny edge detector
    in C++. Reads PGM and PPM files and outputs edge maps, edge chains,
    sampled edge chains, and polygonal approximations.

    The implementation uses a FIR filter and difference-of-gaussians instead
    of Canny's IIR filter; I have code for the latter as well but haven't
    found it to be noticeably better, and its behavior around image borders
    is less predictable.

    Compile with -DMAIN in order to get the main program.

    Without that, you get a class lumo_cedges::CEdges that you can call from
    within your own programs.

    If you find this code useful and use it in some research leading
    to a publication, I would appreciate being mentioned in the
    Acknowledgements section.

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

#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>

#include "cedges.h"

namespace iupr_cedges {
// template <class T> inline T signbit(T a) { return a<0?1:0; }
template <class T>
inline T sqr(T x) {
  return x * x;
}
template <class T>
inline T min(T a, T b) {
  return a < b ? a : b;
}
template <class T>
inline T max(T a, T b) {
  return a > b ? a : b;
}
template <class T>
inline void swap(T &a, T &b) {
  T temp = a;
  a = b;
  b = temp;
}

template <class T>
struct autodel {
  T *data;
  autodel() { data = 0; }
  ~autodel() {
    if (data) delete data;
  }
  autodel(T *other) { data = other; }
  void operator=(T *other) {
    if (data) delete data;
    data = other;
  }
  T &operator*() { return *data; }
  T *operator->() { return data; }
  operator bool() { return !!data; }
};

struct Stdio {
  FILE *stream;
  operator FILE *() { return stream; }
  Stdio() { stream = 0; }
  Stdio(const char *file, const char *mode) {
    stream = fopen(file, mode);
    if (!stream) throw "open failed";
  }
  ~Stdio() {
    if (stream) {
      fclose(stream);
      stream = 0;
    }
  }
  void open(const char *file, const char *mode) {
    if (stream) {
      fclose(stream);
      stream = 0;
    }
    stream = fopen(file, mode);
    if (!stream) throw "open failed";
  }
  void close() {
    if (stream) {
      fclose(stream);
      stream = 0;
    }
  }
  void popen(const char *file, const char *mode) {
    if (stream) {
      fclose(stream);
      stream = 0;
    }
    stream = ::popen(file, mode);
    if (!stream) throw "open failed";
  }
  void pclose() {
    if (stream) {
      ::pclose(stream);
      stream = 0;
    }
  }
};

#define ASSERT(X)                            \
  do {                                       \
    if (!(X)) throw "ASSERTION FAILED: " #X; \
  } while (0)
#define RANGE(X, N) ASSERT(unsigned(X) < unsigned(N))

template <class T>
struct Array {
  T *data;
  int total;
  int size;
  int dims[2];

  explicit Array() {
    total = 4;
    size = 0;
    dims[0] = dims[1] = 0;
    data = new T[total];
  }

  explicit Array(int d0, int d1 = 1) {
    total = size = d0 * d1;
    dims[0] = d0;
    dims[1] = d1;
    data = new T[total];
  }

  ~Array() {
    delete[] data;
    total = size = 0;
    dims[0] = dims[1] = 0;
    data = 0;
  }

  void operator=(T value) {
    for (int i = 0; i < total; i++) data[i] = value;
  }

  void resize(int d0, int d1 = 1) {
    if (d0 * d1 > size) {
      delete[] data;
      total = size = d0 * d1;
      data = new T[total];
    }
    dims[0] = d0;
    dims[1] = d1;
    size = d0 * d1;
  }

  template <class S>
  void copy(Array<S> &other) {
    resize(other.dim(0), other.dim(1));
    for (int i = 0, n = length(); i < n; i++) operator[](i) = other[i];
  }

  void reshape(int d0, int d1 = 1) {
    if (d0 * d1 != size) throw "incompatible reshape";
    dims[0] = d0;
    dims[1] = d1;
  }

  int dim(int i) { return dims[i]; }

  T &operator()(int i0) {
#ifndef UNSAFE
    if (unsigned(i0) >= unsigned(dims[0])) throw "index error";
#endif
    return data[i0];
  }

  T &operator()(int i0, int i1) {
#ifndef UNSAFE
    if (unsigned(i0) >= unsigned(dims[0])) throw "index error";
    if (unsigned(i1) >= unsigned(dims[1])) throw "index error";
#endif
    return data[i1 + dims[1] * i0];
  }

  int length() const { return size; }
  T &operator[](int i0) const {
    if (unsigned(i0) >= unsigned(size)) throw "index error";
    return data[i0];
  }
  T &push() {
    if (size == total) {
      int ntotal = total * 2 + 1;
      T *ndata = new T[ntotal];
      for (int i = 0; i < total; i++) ndata[i] = data[i];
      delete[] data;
      total = ntotal;
      data = ndata;
    }
    return data[size++];
  }
  void push(const T &obj) { push() = obj; }
  void push(T &obj) { push() = obj; }
  T &pop() {
    if (size < 1) throw "index error (pop)";
    return data[--size];
  }
  void clear() { size = 0; }

  void reset() {
    resize(4, 0);
    size = 0;
    dims[0] = dims[1] = 0;
  }

  void getd0(Array<T> &out, int i) {
    out.resize(dims[1]);
    for (int j = 0; j < dims[1]; j++) {
      out(j) = operator()(i, j);
    }
  }
  void getd1(Array<T> &out, int j) {
    out.resize(dims[0]);
    for (int i = 0; i < dims[0]; i++) {
      out(i) = operator()(i, j);
    }
  }
  void putd0(Array<T> &in, int i) {
    for (int j = 0; j < dims[1]; j++) {
      operator()(i, j) = in(j);
    }
  }
  void putd1(Array<T> &in, int j) {
    for (int i = 0; i < dims[0]; i++) {
      operator()(i, j) = in(i);
    }
  }

 private:
  Array(Array<T> &);
  void operator=(Array<T> &);
};

typedef Array<unsigned char> ByteArray;
typedef Array<float> FloatArray;
typedef Array<int> IntArray;

template <class T>
T min(Array<T> &data) {
  T result = data[0];
  for (int i = 1, n = data.length(); i < n; i++) {
    T value = data[i];
    if (value < result) result = value;
  }
  return result;
}

template <class T>
T max(Array<T> &data) {
  T result = data[0];
  for (int i = 1, n = data.length(); i < n; i++) {
    T value = data[i];
    if (value > result) result = value;
  }
  return result;
}

template <class T>
int count(Array<T> &data) {
  int total = 0;
  for (int i = 0, n = data.length(); i < n; i++) {
    if (data[i]) total++;
  }
  return total;
}

template <class T>
void reverse(Array<T> &data) {
  for (int i = 0, n = data.length(), n2 = n / 2; i < n2; i++) {
    swap(data[i], data[n - i - 1]);
  }
}

struct vec2 {
  enum { N = 2 };
  float data[2];
  explicit vec2(float v0 = 0, float v1 = 0) {
    data[0] = v0;
    data[1] = v1;
  }
  int length() const { return N; }
  float at(int i) const {
    RANGE(i, N);
    return data[i];
  }
  float operator()(int i) const {
    RANGE(i, N);
    return data[i];
  }
  float &operator()(int i) {
    RANGE(i, N);
    return data[i];
  }
  float operator[](int i) const {
    RANGE(i, N);
    return data[i];
  }
  float &operator[](int i) {
    RANGE(i, N);
    return data[i];
  }
  vec2 operator+(const vec2 &other) const {
    return vec2(at(0) + other(0), at(1) + other(1));
  }
  vec2 operator-(const vec2 &other) const {
    return vec2(at(0) - other(0), at(1) - other(1));
  }
  float operator*(const vec2 &other) const {
    return at(0) * other(0) + at(1) * other(1);
  }
  vec2 operator*(float scale) const {
    return vec2(at(0) * scale, at(1) * scale);
  }
  vec2 operator/(float scale) const {
    return vec2(at(0) / scale, at(1) / scale);
  }
  float magnitude() const { return sqrt(sqr(data[0]) + sqr(data[1])); }
  float angle() const { return atan2(data[1], data[0]); }
  float magnitude_squared() const { return sqr(data[0]) + sqr(data[1]); }
  vec2 normalized() const { return operator*(1.0 / magnitude()); }
  vec2 normal() const { return vec2(-data[1], data[0]); }
  float distance(const vec2 &b) {
    const vec2 &a = *this;
    return (a - b).magnitude();
  }
  inline vec2 cmul(const vec2 &b) const {
    const vec2 &a = *this;
    return vec2(a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]);
  }
  inline vec2 cdiv(const vec2 &b) const {
    const vec2 &a = *this;
    double n = sqr(b[0]) + sqr(b[1]);
    return vec2((a[0] * b[0] + a[1] * b[1]) / n,
                (a[1] * b[0] - a[0] * b[1]) / n);
  }
};

typedef Array<vec2> Vec2Array;

#if 0
    static float chainlength(Vec2Array &ps) {
        float total = 0.0;
        for(int i=1;i<ps.length();i++) {
            total += ps[i].distance(ps[i-1]);
        }
        return total;
    }
#endif

static int scanint(FILE *file) {
  int value;
  if (fscanf(file, "%d", &value) != 1)
    throw "read_pnm: number format error in image";
  return value;
}

static int getbyte(FILE *file) {
  int value = getc(file);
  if (value == -1) throw "read_pnm: image short";
  return value;
}

static void read_pnm(ByteArray &image, FILE *file) {
  char cp = getc(file);
  if (cp != 'P') throw "read_pnm: unknown foramt";
  int ctype = getc(file) - '0';
  if (ctype != 2 && ctype != 3 && ctype != 5 && ctype != 6)
    throw "read_pnm: cannot handle this type";
  int params[3];
  int nparams = 0;
  while (nparams < 3) {
    int c = getc(file);
    if (c < 0) throw "read_pnm: unexpected eof";
    if (isdigit(c)) {
      ungetc(c, file);
      if (fscanf(file, "%d", &params[nparams]) != 1)
        throw "read_pnm: bad number format";
      nparams++;
      c = getc(file);
      if (!isspace(c)) throw "read_pnm: bad header format";
      continue;
    }
    if (c == '#') {
      for (;;) {
        c = getc(file);
        if (c == '\n') break;
        if (c == -1) throw "read_pnm: unexpected eof";
      }
    }
  }
  image.resize(params[0], params[1]);
  for (int j = image.dim(1) - 1; j >= 0; j--) {
    for (int i = 0; i < image.dim(0); i++) {
      int value;
      switch (ctype) {
        case 2:
          value = scanint(file);
          break;
        case 3:
          value = (scanint(file) + scanint(file) + scanint(file)) / 3;
          break;
        case 5:
          value = getbyte(file);
          break;
        case 6:
          value = (getbyte(file) + getbyte(file) + getbyte(file)) / 3;
          break;
        default:
          throw "bad type";
      }
      image(i, j) = value;
    }
  }
}

static void write_pnm(FILE *file, ByteArray &image) {
  fprintf(file, "P5\n%d %d %d\n", image.dim(0), image.dim(1), 255);
  for (int j = image.dim(1) - 1; j >= 0; j--) {
    for (int i = 0; i < image.dim(0); i++) {
      fputc(image(i, j), file);
    }
  }
}

static void read_pnm(ByteArray &image, const char *file) {
  Stdio stream(file, "r");
  read_pnm(image, stream);
}

static void write_pnm(const char *file, ByteArray &image) {
  Stdio stream(file, "w");
  write_pnm(stream, image);
}

#if 0
    static void display(ByteArray &image) {
        Stdio stream;
        stream.popen("display","w");
        write_pnm(stream,image);
        stream.pclose();
    }
#endif

struct scaled_array {
  ByteArray data;
  scaled_array(FloatArray &image, float minv = 1e30, float maxv = -1e30) {
    int w = image.dim(0), h = image.dim(1), n = image.length();
    data.resize(w, h);
    if (minv > maxv) {
      minv = min(image);
      maxv = max(image);
    }
    if (maxv == minv) maxv = minv + 1.0;
    float scale = 256.0 * maxv - minv;
    for (int i = 0; i < n; i++) {
      int value = int(scale * (image[i] - minv));
      if (value < 0) value = 0;
      if (value > 255) value = 255;
      data[i] = value;
    }
  }
  operator ByteArray &() { return data; }
};

static void gauss1d(FloatArray &out, FloatArray &in, float sigma) {
  out.resize(in.length());
  // make a normalized mask
  int range = 1 + int(3.0 * sigma);
  FloatArray mask(2 * range + 1);
  for (int i = 0; i <= range; i++) {
    double y = exp(-i * i / 2.0 / sigma / sigma);
    mask(range + i) = mask(range - i) = y;
  }
  float total = 0.0;
  for (int i = 0; i < mask.length(); i++) total += mask[i];
  for (int i = 0; i < mask.length(); i++) mask[i] /= total;

  // apply it
  int n = in.length();
  for (int i = 0; i < n; i++) {
    double total = 0.0;
    for (int j = 0; j < mask.length(); j++) {
      int index = i + j - range;
      if (index < 0) index = 0;
      if (index >= n) index = n - 1;
      total += in[index] * mask[j];  // it's symmetric
    }
    out[i] = total;
  }
}

static void gauss2d(FloatArray &a, float sx, float sy) {
  FloatArray r, s;
  for (int i = 0; i < a.dim(0); i++) {
    a.getd0(r, i);
    gauss1d(s, r, sy);
    a.putd0(s, i);
  }
  for (int j = 0; j < a.dim(1); j++) {
    a.getd1(r, j);
    gauss1d(s, r, sx);
    a.putd1(s, j);
  }
}

template <class T>
inline T isign(T x) {
  return ((x) >= 0 ? 1 : -1);
}

static void nonmaxsup(ByteArray &out, FloatArray &gradm, FloatArray &gradx,
                      FloatArray &grady) {
  int w = gradm.dim(0), h = gradm.dim(1);
  out.resize(w, h);
  out = 0;
  for (int i = 1; i < w - 1; i++) {
    for (int j = 1; j < h - 1; j++) {
      float dx = gradx(i, j);
      float ux = fabs(dx);
      float dy = grady(i, j);
      float uy = fabs(dy);
      int bx = int(isign(dx));
      int by = int(isign(dy));
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

static float masked_fractile(FloatArray &gradm, const ByteArray &mask,
                             float frac, int bins = 1000) {
  bool use_mask = mask.length() > 0;
  IntArray hist(bins);
  hist = 0;
  float minv = 1e30, maxv = -1e30;
  int count = 0;
  for (int i = 0, n = gradm.length(); i < n; i++) {
    if (use_mask && !mask[i]) continue;
    count++;
    if (maxv < gradm[i]) maxv = gradm[i];
    if (minv > gradm[i]) minv = gradm[i];
  }
  if (count < 2) return minv;
  if (maxv == minv) return minv;
  int limit = int(frac * count);
  float scale = bins / (maxv - minv);
  for (int i = 0, n = gradm.length(); i < n; i++) {
    if (use_mask && !mask[i]) continue;
    int bin = int(scale * (gradm[i] - minv));
    hist[min(bins - 1, bin)]++;
  }
  int i = 0, sum = 0;
  for (; i < bins && sum < limit; i++) {
    sum += hist[i];
  }
  return (maxv - minv) * i / bins + minv;
}

static void masked_fill(ByteArray &in, ByteArray &out, int x, int y) {
  int w = in.dim(0), h = in.dim(1);
  if (x < 0) return;
  if (x >= w) return;
  if (y < 0) return;
  if (y >= h) return;
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

static void hysteresis_thresholding(ByteArray &out, FloatArray &gradm,
                                    ByteArray &mask, float tlow, float thigh) {
  int w = gradm.dim(0), h = gradm.dim(1);
  ByteArray low(w, h);
  low = 0;
  out.resize(w, h);
  out = 0;
  for (int i = 0; i < w; i++)
    for (int j = 0; j < h; j++) {
      if (mask.dim(0) > 0 && !mask(i, j)) continue;
      if (gradm(i, j) > tlow) low(i, j) = 1;
    }
  for (int i = 0; i < w; i++)
    for (int j = 0; j < h; j++) {
      if (gradm(i, j) > thigh) masked_fill(low, out, i, j);
    }
}

static void thin(ByteArray &uci) {
  enum { OFF = 0, ON = 1, SKEL = 2, DEL = 3 };

  static char ttable[256] = {
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

  static int nx[] = {1, 1, 0, -1, -1, -1, 0, 1};
  static int ny[] = {0, 1, 1, 1, 0, -1, -1, -1};

  int w = uci.dim(0) - 1;
  int h = uci.dim(1) - 1;

  for (int i = 0, n = uci.length(); i < n; i++) {
    if (uci[i])
      uci[i] = ON;
    else
      uci[i] = OFF;
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

  for (int i = 0, n = uci.length(); i < n; i++) {
    if (uci[i] == SKEL)
      uci[i] = 255;
    else
      uci[i] = 0;
  }
};

inline float point_line_dist(vec2 p, vec2 a, vec2 b) {
  vec2 delta = b - a;
  float mag = delta.magnitude();
  // if the distance is small, just return the point distance;
  // that's the right thing for approx_chain
  if (mag < 1e-4) return a.distance(p);
  vec2 normal = delta.normal() / mag;
  float offset = normal * a;
  return fabs(normal * p - offset);
}

static void approx_chain(IntArray &breaks, Vec2Array &chain, int i0, int i1,
                         float maxdist) {
  float md = 0.0;
  int mi = -1;
  vec2 a = chain[i0];
  vec2 b = chain[i1];
  for (int i = i0; i <= i1; i++) {
    float d = point_line_dist(chain[i], a, b);
    // printf("%d %g %g %g\n",i,d,a.distance(chain[i]),b.distance(chain[i]));
    if (d <= md) continue;
    md = d;
    mi = i;
  }
  // printf("i0 %d i1 %d mi %d md %g maxdist %g\n",i0,i1,mi,md,maxdist);
  if (md < maxdist) return;
  ASSERT(mi != i0 && mi != i1);
  approx_chain(breaks, chain, i0, mi, maxdist);
  breaks.push(mi);
  approx_chain(breaks, chain, mi, i1, maxdist);
}

struct ChainTracer {
  enum { OFF = 0, ON = 1, DONE = 2 };

  ByteArray bi;
  int sx, sy;
  int x, y;
  int w, h;

  int count_neighbors() {
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
      y = y;
      return true;
    }
    if (bi(x + 1, y + 1) == ON) {
      x = x + 1;
      y = y + 1;
      return true;
    }
    if (bi(x, y + 1) == ON) {
      x = x;
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
      y = y;
      return true;
    }
    if (bi(x - 1, y - 1) == ON) {
      x = x - 1;
      y = y - 1;
      return true;
    }
    if (bi(x, y - 1) == ON) {
      x = x;
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

  void set_image(ByteArray &image) {
    w = image.dim(0);
    h = image.dim(1);
    x = 0;
    y = 0;
    sx = 0;
    sy = 0;
    bi.resize(w, h);
    for (int i = 0, n = image.length(); i < n; i++) bi[i] = image[i] ? ON : OFF;
    for (int i = 0; i < w; i++) bi(i, 0) = bi(i, h - 1) = OFF;
    for (int j = 0; j < h; j++) bi(0, j) = bi(w - 1, j) = OFF;
    w--;
    h--;
  }

  bool get_chain(Vec2Array &points, bool close = false) {
    points.clear();
    if (!nextstart()) return false;
    do {
      points.push() = vec2(x, y);
      bi(x, y) = DONE;
    } while (nextpixel());
    x = sx;
    y = sy;
    if (nextpixel()) {
      // sweep up the other direction (if any)
      reverse(points);
      do {
        points.push() = vec2(x, y);
        bi(x, y) = DONE;
      } while (nextpixel());
    }
    if (close && points[0].distance(points[points.length() - 1]) < 2.0) {
      // close circular chains
      points.push(points[0]);
    }
    return true;
  }

  IntArray breaks;
  Vec2Array chain;

  bool get_poly(Vec2Array &poly, float maxdist = 1.0, bool close = false) {
    poly.clear();
    breaks.clear();
    if (!get_chain(chain, close)) return false;
    if (chain.length() < 3) {
      for (int i = 0; i < chain.length(); i++) poly.push(chain[i]);
    } else {
      breaks.push(0);
      approx_chain(breaks, chain, 0, chain.length() - 1, maxdist);
      breaks.push(chain.length() - 1);
      for (int i = 0; i < breaks.length(); i++) poly.push(chain[breaks[i]]);
    }
    return true;
  }

  bool started() { return bi.dim(0) > 0; }

  void clear() {
    bi.clear();
    breaks.clear();
    chain.clear();
  }
};

struct CEdges : EdgeDetector {
  float p_sx;
  float p_sy;
  float p_frac;
  float p_tlow;
  float p_thigh;
  float p_minlength;
  float p_maxdist;

  CEdges() {
    p_sx = 3.0;
    p_sy = p_sx;
    p_frac = 0.3;
    p_tlow = 2.0;
    p_thigh = 4.0;
    p_minlength = 5.0;
    p_maxdist = 0.5;
  }

  ~CEdges() {}

  void set_gauss(float sx, float sy) {
    p_sx = sx;
    p_sy = sy;
  }

  void set_noise(float frac, float low, float high) {
    p_frac = frac;
    p_tlow = low;
    p_thigh = high;
  }

  void set_poly(float minlength, float maxdist) {
    p_minlength = minlength;
    p_maxdist = maxdist;
  }

  FloatArray image;
  ByteArray uedges;
  FloatArray smoothed;
  FloatArray gradm, gradx, grady;
  ByteArray edges;
  ChainTracer tracer;
  float noise;

  void clear() {
    image.clear();
    uedges.clear();
    smoothed.clear();
    gradm.clear();
    gradx.clear();
    grady.clear();
    edges.clear();
    tracer.clear();
    noise = 0.0;
  }
  void load_pnm(FILE *file) {
    ByteArray bimage;
    read_pnm(bimage, file);
    image.copy(bimage);
  }
  void load_pnm(char *file) {
    ByteArray bimage;
    read_pnm(bimage, file);
    image.copy(bimage);
  }
  void save_pnm(FILE *file) { write_pnm(file, edges); }
  void save_pnm(char *file) { write_pnm(file, edges); }
  int dim(int i) { return image.dim(i); }
  void set_image(unsigned char *p, int w, int h) {
    image.resize(w, h);
    for (int i = 0; i < image.length(); i++) image[i] = p[i];
  }
  void set_pixmap(unsigned char *p, int w, int h) {
    image.resize(w, h);
    for (int y = h - 1; y >= 0; y--)
      for (int x = 0; x < w; x++) image(x, y) = *p++;
  }
  void compute() {
    int w = image.dim(0);
    int h = image.dim(1);

    smoothed.copy(image);
    gauss2d(smoothed, p_sx, p_sy);

    gradm.resize(w, h);
    gradx.resize(w, h);
    grady.resize(w, h);
    gradm = 0.0;
    gradx = 0.0;
    grady = 0.0;
    for (int i = w - 2; i >= 0; i--)
      for (int j = h - 2; j >= 0; j--) {
        float v = smoothed(i, j);
        float dx = smoothed(i + 1, j) - v;
        float dy = smoothed(i, j + 1) - v;
        gradx(i, j) = dx;
        grady(i, j) = dy;
        gradm(i, j) = sqrt(sqr(dx) + sqr(dy));
      }

    nonmaxsup(uedges, gradm, gradx, grady);

    for (int i = 0, n = uedges.length(); i < n; i++)
      if (uedges[i]) uedges[i] = 255;
    thin(uedges);

    ByteArray temp(0, 0);
    noise = masked_fractile(gradm, temp, p_frac);
    float low = p_tlow * noise;
    float high = p_thigh * noise;
    hysteresis_thresholding(edges, gradm, uedges, low, high);
    for (int i = 0; i < edges.length(); i++) edges[i] = 255 * !!edges[i];

    tracer.set_image(edges);
  }
  void get_eimage(unsigned char *p, int w, int h) {
    if (edges.dim(0) != w || edges.dim(1) != h)
      throw "output image has the wrong size";
    for (int i = 0; i < edges.length(); i++) p[i] = edges[i];
  }
  void get_epixmap(unsigned char *image, int w, int h) {
    if (edges.dim(0) != w || edges.dim(1) != h)
      throw "output image has the wrong size";
    for (int y = h - 1; y >= 0; y--)
      for (int x = 0; x < w; x++) *image++ = edges(x, y);
  }
  float gradient_magnitude(int x, int y) { return gradm(x, y); }
  float gradient_angle(int x, int y) { return atan2(grady(x, y), gradx(x, y)); }

  Vec2Array chain;
  IntArray breaks;

  bool nextchain() {
    if (!tracer.get_chain(chain)) return false;
    breaks.clear();
    breaks.push(0);
    approx_chain(breaks, chain, 0, chain.length() - 1, p_maxdist);
    breaks.push(chain.length() - 1);
    return true;
  }
  int npoints() { return chain.length(); }
  void point(int index, float &x, float &y) {
    x = chain[index][0];
    y = chain[index][1];
  }
  int nsegments() { return breaks.length() - 1; }
  void segment(int i, float &x0, float &y0, float &x1, float &y1, float &angle,
               float &magnitude, int &n) {
    i++;
    int i0 = breaks[i - 1];
    int i1 = breaks[i];
    vec2 a = chain[i0];
    vec2 b = chain[i1];
    vec2 g = vec2(0.0, 0.0);
    for (int j = i0; j <= i1; j++) {
      int x = int(chain[j][0]);
      int y = int(chain[j][1]);
      g = g + vec2(gradx(x, y), grady(x, y));
    }
    // g = g / (i1-i0+1);
    x0 = a[0];
    y0 = a[1];
    x1 = b[0];
    y1 = b[1];
    magnitude = g.magnitude();
    angle = g.angle();
    while (angle < 0) angle += 2 * M_PI;
    n = i1 - i0 + 1;
  }
  struct Sample {
    vec2 c;
    vec2 g;
    int n;
  };
  Array<Sample> samples;
  void sampleat(int spacing) {
    samples.clear();
    for (int i = 1; i < breaks.length(); i++) {
      int start = breaks[i - 1];
      int end = breaks[i];
      for (int j = start + spacing / 2; j < end; j += spacing) {
        int j0 = max(0, j - spacing / 2);
        int j1 = min(end, j + spacing / 2);
        vec2 c = vec2(0, 0);
        vec2 g = vec2(0, 0);
        for (int k = j0; k <= j1; k++) {
          int x = int(chain[k][0]);
          int y = int(chain[k][1]);
          c = c + vec2(x, y);
          g = g + vec2(gradx(x, y), grady(x, y));
        }
        c = c / (j1 - j0 + 1);
        // g = g / (j1-j0+1);
        Sample &sample = samples.push();
        sample.c = c;
        sample.g = g;
        sample.n = j1 - j0 + 1;
      }
    }
  }
  int nsamples() { return samples.length(); }
  void sample(int i, float &x, float &y, float &angle, float &mag,
              int &npoints) {
    Sample &s = samples(i);
    x = s.c(0);
    y = s.c(1);
    angle = atan2(s.g(1), s.g(0));
    mag = hypot(s.g(0), s.g(1));
    npoints = s.n;
  }
};

EdgeDetector *makeEdgeDetector() { return new CEdges(); }
}

#ifdef MAIN

namespace iupr_cedges {
inline int igetenv(const char *name, int dflt) {
  int result = getenv(name) ? atoi(getenv(name)) : dflt;
  int where = 0;
  if (strcmp(name, "verbose_params")) where = igetenv("verbose_params", 0);
  switch (where) {
    case 1:
      fprintf(stdout, "__param__ %s = %d\n", name, result);
      break;
    case 2:
      fprintf(stderr, "__param__ %s = %d\n", name, result);
      break;
    default:;
  }
  return result;
}
inline float fgetenv(const char *name, float dflt) {
  float result = getenv(name) ? atof(getenv(name)) : dflt;
  int where = igetenv("verbose_params", 0);
  switch (where) {
    case 1:
      fprintf(stdout, "__param__ %s = %g\n", name, result);
      break;
    case 2:
      fprintf(stderr, "__param__ %s = %g\n", name, result);
      break;
    default:;
  }
  return result;
}
inline double dgetenv(const char *name, double dflt) {
  double result = getenv(name) ? atof(getenv(name)) : dflt;
  int where = igetenv("verbose_params", 0);
  switch (where) {
    case 1:
      fprintf(stdout, "__param__ %s = %g\n", name, result);
      break;
    case 2:
      fprintf(stderr, "__param__ %s = %g\n", name, result);
      break;
    default:;
  }
  return result;
}
inline const char *sgetenv(const char *name, const char *dflt) {
  const char *result = getenv(name) ? getenv(name) : dflt;
  int where = igetenv("verbose_params", 0);
  switch (where) {
    case 1:
      fprintf(stdout, "__param__ %s = %s\n", name, result);
      break;
    case 2:
      fprintf(stderr, "__param__ %s = %s\n", name, result);
      break;
    default:;
  }
  return result;
}

const char *usage =
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
}

int main(int argc, char **argv) {
  using namespace iupr_cedges;
  autodel<CEdges> cedges = new CEdges();
  cedges->set_gauss(fgetenv("esx", 3.0), fgetenv("esy", 3.0));
  cedges->set_noise(fgetenv("efrac", 0.3), fgetenv("elow", 2.0),
                    fgetenv("ehigh", 4.0));
  float eminlength = fgetenv("eminlength", 5.0);
  cedges->set_poly(eminlength, fgetenv("emaxdist", 1.5));
  const char *p_format = sgetenv("eformat", "segments");
  int p_spacing = igetenv("espacing", 4);

  if (argc > 2 ||
      (argc == 2 && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "-?"))) ||
      (argc == 1 && isatty(0))) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }

  try {
    if (argc == 2)
      cedges->load_pnm(argv[1]);
    else
      cedges->load_pnm(stdin);
  } catch (const char *s) {
    fprintf(stderr, "problem reading image: %s\n", s);
    exit(2);
  } catch (...) {
    fprintf(stderr, "problem reading image\n");
    exit(2);
  }

  cedges->compute();

  if (!strcmp(p_format, "map")) {
    ByteArray edges(cedges->dim(0), cedges->dim(1));
    cedges->get_eimage(&edges(0, 0), edges.dim(0), edges.dim(1));
    write_pnm(stdout, edges);
  } else if (!strcmp(p_format, "chain")) {
    while (cedges->nextchain()) {
      for (int i = 0; i < cedges->npoints(); i++) {
        float x, y;
        cedges->point(i, x, y);
        printf("%g %g\n", x, y);
      }
      printf("\n");
    }
  } else if (!strcmp(p_format, "poly")) {
    while (cedges->nextchain()) {
      for (int i = 0; i < cedges->nsegments(); i++) {
        float x, y, x1, y1, a, w;
        int n;
        cedges->segment(i, x, y, x1, y1, a, w, n);
        if (i == 0) printf("%g %g\n", x, y);
        printf("%g %g\n", x1, y1);
      }
      printf("\n");
    }
  } else if (!strcmp(p_format, "segments")) {
    while (cedges->nextchain()) {
      int count = 0;
      for (int i = 0; i < cedges->nsegments(); i++) {
        float x, y, x1, y1, a, w;
        int n;
        cedges->segment(i, x, y, x1, y1, a, w, n);
        if (sqrt(sqr(x1 - x) + sqr(y1 - y)) < eminlength) continue;
        printf("%g %g %g %g  %g %g %d\n", x, y, x1, y1, a, w, n);
        count++;
      }
      if (count > 0) printf("\n");
    }
  } else if (!strcmp(p_format, "sampled")) {
    while (cedges->nextchain()) {
      cedges->sampleat(p_spacing);
      for (int i = 0; i < cedges->nsamples(); i++) {
        float x, y, a, w;
        int n;
        cedges->sample(i, x, y, a, w, n);
        printf("%g %g  %g %g %d\n", x, y, a, w, n);
      }
      printf("\n");
    }
  } else {
    fprintf(stderr, "%s: unknown format\n", p_format);
    exit(1);
  }
}

#endif
