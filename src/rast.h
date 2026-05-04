// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#ifndef RAST_SRC_RAST_H_
#define RAST_SRC_RAST_H_

/// \file rast.h
/// \brief Public umbrella header for the RAST family of branch-and-bound
///        feature-matching algorithms (Recognition by Adaptive Subdivisions
///        of Transformation space).
///
/// Each algorithm is exposed as an abstract interface plus a `make*()`
/// factory returning `std::unique_ptr` to a freshly-allocated implementation.
/// Implementations are private to their respective `.cc` files; clients see
/// only the interfaces here.
///
/// All interfaces are single-threaded: a given instance must be touched by
/// at most one thread at a time. Distinct instances are independent.
///
/// All result accessors are valid only after the corresponding `compute()`
/// or `match()` call has returned. Calling them before, or with `rank` /
/// `i` outside `[0, nresults())`, is undefined behavior.

#include <memory>

/// \internal Helper macro: declare an abstract interface with a virtual
/// destructor and deleted copy/move (preventing slicing through base
/// references). Used to keep each interface body terse and consistent.
#define RAST_INTERFACE(NAME)                          \
  NAME() = default;                                   \
  virtual ~NAME() = default;                          \
  NAME(const NAME &) = delete;                        \
  NAME &operator=(const NAME &) = delete;             \
  NAME(NAME &&) = delete;                             \
  NAME &operator=(NAME &&) = delete

/// \brief Find lines through a set of oriented points using the RAST
///        branch-and-bound search over (theta, r) space.
///
/// Each input point carries a 2D location, an orientation angle, and an
/// optional weight. The algorithm partitions the (theta, r) parameter
/// space and returns the top-`maxresults` lines whose accumulated point
/// weight is maximal.
struct LinesP2D {
  RAST_INTERFACE(LinesP2D);

  /// Cap the number of results emitted by `compute()`. Default 1.
  virtual void set_maxresults(int n) = 0;
  virtual void set_breakpenalty(float eps, float cost) = 0;
  /// Set spatial and angular error tolerances used during scoring.
  virtual void set_error(float eps, float aeps) = 0;
  /// Set termination tolerance on the parameter-space rectangle.
  virtual void set_tolerance(float tol, float atol) = 0;
  virtual void set_verbose(int value) = 0;
  /// Drop solutions whose accumulated weight is below this threshold.
  virtual void set_minweight(float value) = 0;
  /// Maximum line offset r searched.
  virtual void set_maxoffset(float value) = 0;
  /// Switch between hard-threshold scoring (false) and least-squares (true).
  virtual void set_lsq(bool value) = 0;
  /// Treat input angles modulo pi (true) instead of 2pi (false).
  virtual void set_unoriented(bool value) = 0;

  /// Discard all previously-added image points.
  virtual void clear_ipoints() = 0;
  /// Append an image point at (x, y) with orientation `a` (radians) and
  /// weight `w`.
  virtual void add_ipoint(float x, float y, float a, float w = 1.0) = 0;

  /// Run the line search over the default (theta, r) range.
  virtual void compute() = 0;
  /// Run the search restricted to theta in [a0, a1] and r in [d0, d1].
  virtual void compute(float a0, float a1, float d0, float d1) = 0;

  /// Number of lines found by the most recent `compute()`.
  virtual int nresults() const = 0;
  /// Accumulated weight of the line at `rank` (higher = better).
  virtual float weight(int rank) const = 0;
  /// Theta of the line at `rank` (radians).
  virtual float angle(int rank) const = 0;
  /// Offset r of the line at `rank`.
  virtual float offset(int rank) const = 0;
  /// Number of points contributing to the line at `rank`.
  virtual int nmatches(int rank) const = 0;
};

/// \brief Allocate a default-configured LinesP2D instance.
std::unique_ptr<LinesP2D> makeLinesP2D();

/// \brief Find lines through a set of line segments using the RAST search.
///
/// Like LinesP2D but with segments as inputs: each segment contributes its
/// two endpoints to the line score, and the algorithm rewards segments
/// whose entire length lies within the (theta, r) tolerance.
struct LinesS2D {
  RAST_INTERFACE(LinesS2D);

  virtual void set_maxresults(int n) = 0;
  virtual void set_breakpenalty(float eps, float cost) = 0;
  virtual void set_error(float eps, float aeps) = 0;
  virtual void set_tolerance(float tol, float atol) = 0;
  virtual void set_verbose(int value) = 0;
  virtual void set_minweight(float value) = 0;
  virtual void set_maxoffset(float value) = 0;
  virtual void set_lsq(bool value) = 0;
  virtual void set_unoriented(bool value) = 0;

  virtual void clear_ipoints() = 0;
  /// Append an image segment from (x, y) to (x1, y1) with orientation `a`.
  virtual void add_iseg(float x, float y, float x1, float y1, float a, float w = 1.0) = 0;

  virtual void compute() = 0;

  virtual int nresults() const = 0;
  virtual float weight(int rank) const = 0;
  virtual float angle(int rank) const = 0;
  virtual float offset(int rank) const = 0;
};

std::unique_ptr<LinesS2D> makeLinesS2D();

/// \brief Generate synthetic instances of the 2D point recognition problem.
///
/// Produces a model point set, a transformed-and-cluttered image point set,
/// and the ground-truth (translation, angle, scale) used to relate them.
/// Intended for benchmarking / testing of the matchers below.
struct InstanceP2D {
  RAST_INTERFACE(InstanceP2D);

  virtual void set_image_size(int r) = 0;
  virtual void set_model_size(int r) = 0;
  /// Number of clutter points added to the image.
  virtual void set_nclutter(int v) = 0;
  /// Total model-point count.
  virtual void set_nmodel_total(int v) = 0;
  /// Of those, the count present (un-occluded) in the image.
  virtual void set_nmodel_unoccluded(int v) = 0;
  virtual void set_error(float v) = 0;
  virtual void set_aerror(float v) = 0;
  virtual void set_srange(float min, float max) = 0;

  /// Generate a fresh instance using the configured parameters and the
  /// current `drand48()` state. Call `srand48(seed)` before this for
  /// reproducibility.
  virtual void generate() = 0;

  /// Number of image points, including clutter.
  virtual int nimage() const = 0;
  /// Read out the i-th image point (0-based).
  virtual void get_image(float &x, float &y, float &a, int i) const = 0;

  virtual int nmodel() const = 0;
  virtual void get_model(float &x, float &y, float &a, int i) const = 0;

  /// Ground-truth parameter `i`: 0=tx, 1=ty, 2=angle, 3=scale.
  /// Throws `std::out_of_range` if `i` is outside [0, 4).
  virtual float get_param(int i) const = 0;
};

std::unique_ptr<InstanceP2D> makeInstanceP2D();

/// \brief Match oriented model points against oriented image points by
///        branch-and-bound search over (tx, ty, theta, scale).
///
/// Each model point may carry its own spatial and angular tolerance,
/// allowing stronger features to enforce tighter matches.
struct RastP2D {
  RAST_INTERFACE(RastP2D);

  virtual void set_maxresults(int n) = 0;
  virtual void set_verbose(bool value) = 0;
  /// Set termination tolerance. Throws `std::invalid_argument` if value
  /// is below 1e-3 (numerically unstable below that point).
  virtual void set_tolerance(float value) = 0;
  /// Drop subregions whose upper bound on quality falls below `min_q`.
  virtual void set_min_q(float min_q) = 0;
  virtual void set_xrange(float x0, float x1) = 0;
  virtual void set_yrange(float y0, float y1) = 0;
  virtual void set_arange(float a0, float a1) = 0;
  virtual void set_srange(float s0, float s1) = 0;
  virtual void set_lsq(bool value) = 0;
  virtual void set_unoriented(bool value) = 0;

  virtual void clear_msources() = 0;
  /// Append a model point with per-source spatial (`eps`) and angular
  /// (`aeps`) tolerances.
  virtual void add_msource(float x, float y, float a, float eps, float aeps) = 0;

  virtual void clear_ipoints() = 0;
  virtual void add_ipoint(float x, float y, float a) = 0;

  /// Run the branch-and-bound match. The set of matches is committed to
  /// `nresults()` etc. on return.
  virtual void match() = 0;

  virtual int nresults() const = 0;
  /// Upper bound on the match quality at `rank`.
  virtual float ubound(int rank) const = 0;
  /// Lower bound on the match quality at `rank`.
  virtual float lbound(int rank) const = 0;
  /// Component `dim` (0 or 1) of the recovered translation at `rank`.
  virtual float translation(int rank, int dim) const = 0;
  virtual float angle(int rank) const = 0;
  virtual float scale(int rank) const = 0;
};

std::unique_ptr<RastP2D> makeRastP2D();

/// \brief Match model line segments to image line segments by RAST search.
struct RastS2D {
  RAST_INTERFACE(RastS2D);

  virtual void set_maxresults(int n) = 0;
  virtual void set_verbose(bool value) = 0;
  /// Throws `std::invalid_argument` if value is below 1e-3.
  virtual void set_tolerance(float value) = 0;
  virtual void set_eps(float eps, float aeps) = 0;
  virtual void set_scale_eps(bool value, float ieps) = 0;
  virtual void set_lsq(bool value) = 0;
  virtual void set_qtolerance(float value) = 0;
  virtual void set_min_q(float min_q) = 0;
  virtual void set_xrange(float x0, float x1) = 0;
  virtual void set_yrange(float y0, float y1) = 0;
  virtual void set_arange(float a0, float a1) = 0;
  virtual void set_srange(float s0, float s1) = 0;

  virtual void clear_msources() = 0;
  virtual void add_mseg(float x0, float y0, float x1, float y1) = 0;

  virtual void clear_ipoints() = 0;
  virtual void add_iseg(float x0, float y0, float x1, float y1) = 0;

  virtual void match() = 0;

  virtual int nresults() const = 0;
  virtual float ubound(int rank) const = 0;
  virtual float lbound(int rank) const = 0;
  virtual float translation(int rank, int dim) const = 0;
  virtual float angle(int rank) const = 0;
  virtual float scale(int rank) const = 0;
};

std::unique_ptr<RastS2D> makeRastS2D();

/// \brief Match model segments to image segments using RAST + segment
///        sampling. Each model segment is sampled at intervals of
///        `sdist` and treated as multiple constraint points.
///
/// `set_eps()` takes a third parameter `sdist` controlling the sampling
/// density. Throws `std::length_error` if the model has more than 255
/// segments or the image more than 32000.
struct RastSS2D {
  RAST_INTERFACE(RastSS2D);

  virtual void set_maxresults(int n) = 0;
  virtual void set_verbose(bool value) = 0;
  virtual void set_tolerance(float value) = 0;
  virtual void set_eps(float eps, float aeps, float sdist) = 0;
  virtual void set_lsq(bool value) = 0;
  virtual void set_qtolerance(float value) = 0;
  virtual void set_min_q(float min_q) = 0;
  virtual void set_xrange(float x0, float x1) = 0;
  virtual void set_yrange(float y0, float y1) = 0;
  virtual void set_arange(float a0, float a1) = 0;
  virtual void set_srange(float s0, float s1) = 0;

  virtual void clear_msources() = 0;
  virtual void add_mseg(float x0, float y0, float x1, float y1) = 0;

  virtual void clear_ipoints() = 0;
  virtual void add_iseg(float x0, float y0, float x1, float y1) = 0;

  virtual void match() = 0;

  virtual int nresults() const = 0;
  virtual float ubound(int rank) const = 0;
  virtual float lbound(int rank) const = 0;
  virtual float translation(int rank, int dim) const = 0;
  virtual float angle(int rank) const = 0;
  virtual float scale(int rank) const = 0;
};

std::unique_ptr<RastSS2D> makeRastSS2D();

/// \brief RAST + sampling segment matcher with a Blob-aware scoring path.
///
/// Same input shape as `RastSS2D`. The internal scoring includes a
/// blob-distance term in addition to the segment-distance terms; the
/// public interface does not yet expose blob inputs separately.
struct RastRS2D {
  RAST_INTERFACE(RastRS2D);

  virtual void set_maxresults(int n) = 0;
  virtual void set_verbose(bool value) = 0;
  virtual void set_tolerance(float value) = 0;
  virtual void set_eps(float eps, float aeps, float sdist) = 0;
  virtual void set_lsq(bool value) = 0;
  virtual void set_qtolerance(float value) = 0;
  virtual void set_min_q(float min_q) = 0;
  virtual void set_xrange(float x0, float x1) = 0;
  virtual void set_yrange(float y0, float y1) = 0;
  virtual void set_arange(float a0, float a1) = 0;
  virtual void set_srange(float s0, float s1) = 0;

  virtual void clear_msources() = 0;
  virtual void add_mseg(float x0, float y0, float x1, float y1) = 0;

  virtual void clear_ipoints() = 0;
  virtual void add_iseg(float x0, float y0, float x1, float y1) = 0;

  virtual void match() = 0;

  virtual int nresults() const = 0;
  virtual float ubound(int rank) const = 0;
  virtual float lbound(int rank) const = 0;
  virtual float translation(int rank, int dim) const = 0;
  virtual float angle(int rank) const = 0;
  virtual float scale(int rank) const = 0;
};

std::unique_ptr<RastRS2D> makeRastRS2D();

/// \brief Brute-force exact alignment of two 2D point sets.
///
/// Picks the (model, image) point pair that anchors a similarity transform
/// and counts inliers. O(n^2 * m^2) over input sizes; intended for small
/// problems and for ground-truth comparison against the RAST matchers.
struct AlignmentP2D {
  RAST_INTERFACE(AlignmentP2D);

  /// Inlier distance threshold.
  virtual void set_epsilon(float e) = 0;
  virtual void set_srange(float min, float max) = 0;

  virtual void clear_mpoints() = 0;
  virtual void add_mpoint(float x, float y) = 0;

  virtual void clear_ipoints() = 0;
  virtual void add_ipoint(float x, float y) = 0;

  virtual void compute() = 0;

  /// Number of model points that aligned with image points within `epsilon`.
  virtual float quality() const = 0;
  /// Component `dim` (0 or 1) of the recovered translation.
  virtual float translation(int dim) const = 0;
  virtual float angle() const = 0;
  virtual float scale() const = 0;
};

std::unique_ptr<AlignmentP2D> makeAlignmentP2D();

#endif  // RAST_SRC_RAST_H_
