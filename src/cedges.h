// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#ifndef RAST_SRC_CEDGES_H_
#define RAST_SRC_CEDGES_H_

/// \file cedges.h
/// \brief Public interface for the self-contained Canny-style edge detector.
///
/// One implementation, accessed through `iupr_cedges::makeEdgeDetector()`.
/// The detector reads a grayscale or RGB image (PGM/PPM or a buffer),
/// runs gradient + non-maximum suppression + hysteresis thresholding,
/// and exposes the resulting edges as either a chain-by-chain stream or
/// a polygonal approximation.
///
/// Single-threaded: a given EdgeDetector instance must be touched by at
/// most one thread at a time.

#include <memory>

namespace iupr_cedges {

/// \brief Canny-style edge detector with chain and polygon outputs.
///
/// Typical usage:
///   1. configure parameters via `set_*()`
///   2. supply input (`load_pnm` or `set_image` / `set_pixmap`)
///   3. call `compute()`
///   4. iterate chains with `nextchain()` and read out via `point()` /
///      `segment()` / similar.
struct EdgeDetector {
  EdgeDetector() = default;
  virtual ~EdgeDetector() = default;
  EdgeDetector(const EdgeDetector &) = delete;
  EdgeDetector &operator=(const EdgeDetector &) = delete;
  EdgeDetector(EdgeDetector &&) = delete;
  EdgeDetector &operator=(EdgeDetector &&) = delete;

  /// Set the Gaussian sigma for the smoothing pre-filter.
  virtual void set_gauss(float sx, float sy) = 0;
  /// Set the Canny noise model: `frac` = noise fraction estimate; `low`
  /// and `high` are the hysteresis thresholds.
  virtual void set_noise(float frac, float low, float high) = 0;
  /// Set the polygonal-approximation parameters: `minlength` is the
  /// shortest segment kept; `maxdist` is the max chord-deviation tolerated.
  virtual void set_poly(float minlength, float maxdist) = 0;
  /// Reset internal state; required between consecutive image runs.
  virtual void clear() = 0;

  /// Read input from a PGM or PPM file. Throws `std::runtime_error` on
  /// I/O or format errors.
  virtual void load_pnm(const char *file) = 0;
  /// Write the edge map as a PNM file.
  virtual void save_pnm(const char *file) = 0;

  /// Image dimension along axis `i` (0 = width, 1 = height).
  virtual int dim(int i) const = 0;
  /// Set the input image from a caller-owned grayscale buffer (size w*h).
  /// Buffer must remain valid until `compute()` returns.
  virtual void set_image(unsigned char *p, int w, int h) = 0;
  /// Set the input image from a caller-owned RGB buffer (size w*h*3).
  virtual void set_pixmap(unsigned char *p, int w, int h) = 0;

  /// Run the detection pipeline (smoothing, gradients, NMS, hysteresis).
  virtual void compute() = 0;

  /// Copy out the binary edge map into the caller-owned buffer.
  /// Buffer dimensions must match `dim(0) x dim(1)`; throws otherwise.
  virtual void get_eimage(unsigned char *p, int w, int h) const = 0;
  /// Copy out a colorized visualization of the edges.
  virtual void get_epixmap(unsigned char *image, int w, int h) const = 0;

  virtual float gradient_magnitude(int x, int y) const = 0;
  virtual float gradient_angle(int x, int y) const = 0;

  /// Advance to the next edge chain. Returns false when the chains are
  /// exhausted. After this call, `npoints()` and `nsegments()` describe
  /// the current chain.
  virtual bool nextchain() = 0;
  /// Number of points in the current chain.
  virtual int npoints() const = 0;
  /// Read out the `index`-th point of the current chain.
  virtual void point(int index, float &x, float &y) const = 0;
  /// Number of segments in the polygonal approximation of the current chain.
  virtual int nsegments() const = 0;
  /// Read out the `i`-th segment of the polygonal approximation.
  virtual void segment(int i, float &x0, float &y0, float &x1, float &y1, float &angle,
                       float &magnitude, int &n) const = 0;
};

/// Allocate a default-configured edge detector.
std::unique_ptr<EdgeDetector> makeEdgeDetector();

}  // namespace iupr_cedges

#endif  // RAST_SRC_CEDGES_H_
