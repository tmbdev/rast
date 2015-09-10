/* Copyright (c) 1990-1995 by Thomas M. Breuel */

// line finding using the RAST algorithm

struct LinesP2D {
  // parameters
  virtual void set_maxresults(int n) = 0;
  virtual void set_breakpenalty(float eps, float cost) = 0;
  virtual void set_error(float eps, float aeps) = 0;
  virtual void set_tolerance(float tol, float atol) = 0;
  virtual void set_verbose(int value) = 0;
  virtual void set_minweight(float value) = 0;
  virtual void set_maxoffset(float value) = 0;
  virtual void set_lsq(bool value) = 0;
  virtual void set_unoriented(bool value) = 0;

  // image points
  virtual void clear_ipoints() = 0;
  virtual void add_ipoint(float x, float y, float a, float w = 1.0) = 0;

  // compute matches
  virtual void compute() = 0;
  virtual void compute(float a0, float a1, float d0, float d1) = 0;

  // get results
  virtual int nresults() = 0;
  virtual float weight(int rank) = 0;
  virtual float angle(int rank) = 0;
  virtual float offset(int rank) = 0;
  virtual int nmatches(int rank) = 0;

  virtual ~LinesP2D() {}
};

LinesP2D *makeLinesP2D();

// line finding using the RAST algorithm

struct LinesS2D {
  // parameters
  virtual void set_maxresults(int n) = 0;
  virtual void set_breakpenalty(float eps, float cost) = 0;
  virtual void set_error(float eps, float aeps) = 0;
  virtual void set_tolerance(float tol, float atol) = 0;
  virtual void set_verbose(int value) = 0;
  virtual void set_minweight(float value) = 0;
  virtual void set_maxoffset(float value) = 0;
  virtual void set_lsq(bool value) = 0;
  virtual void set_unoriented(bool value) = 0;

  // image points
  virtual void clear_ipoints() = 0;
  virtual void add_iseg(float x, float y, float x1, float y1, float a,
                        float w = 1.0) = 0;

  // compute matches
  virtual void compute() = 0;

  // get results
  virtual int nresults() = 0;
  virtual float weight(int rank) = 0;
  virtual float angle(int rank) = 0;
  virtual float offset(int rank) = 0;

  virtual ~LinesS2D() {}
};

LinesS2D *makeLinesS2D();

// generate instances of the 2D recognition problem

struct InstanceP2D {
  // parameters
  virtual void set_image_size(int r) = 0;
  virtual void set_model_size(int r) = 0;
  virtual void set_nclutter(int v) = 0;
  virtual void set_nmodel_total(int v) = 0;
  virtual void set_nmodel_unoccluded(int v) = 0;
  virtual void set_error(float v) = 0;
  virtual void set_aerror(float v) = 0;
  virtual void set_srange(float min, float max) = 0;  // range of scales

  // generate instance
  virtual void generate() = 0;

  // get image data
  virtual int nimage() = 0;
  virtual void get_image(float &x, float &y, float &a, int i) = 0;

  // get model data
  virtual int nmodel() = 0;
  virtual void get_model(float &x, float &y, float &a, int i) = 0;

  virtual ~InstanceP2D() {}
  virtual float get_param(int i) = 0;
};

InstanceP2D *makeInstanceP2D();

// point matching using the RAST algorithm

struct RastP2D {
  // parameters
  virtual void set_maxresults(int n) = 0;
  virtual void set_verbose(bool value) = 0;
  virtual void set_tolerance(float value) = 0;
  virtual void set_min_q(float min_q) = 0;
  virtual void set_xrange(float x0, float x1) = 0;
  virtual void set_yrange(float y0, float y1) = 0;
  virtual void set_arange(float a0, float a1) = 0;  // range of angles
  virtual void set_srange(float s0, float s1) = 0;  // range of scales
  virtual void set_lsq(bool value) = 0;
  virtual void set_unoriented(bool value) = 0;

  // set model points
  virtual void clear_msources() = 0;
  virtual void add_msource(float x, float y, float a, float eps,
                           float aeps) = 0;

  // set image points
  virtual void clear_ipoints() = 0;
  virtual void add_ipoint(float x, float y, float a) = 0;

  // perform match
  virtual void match() = 0;

  // get results
  virtual int nresults() = 0;
  virtual float ubound(int rank) = 0;
  virtual float lbound(int rank) = 0;
  virtual float translation(int rank, int dim) = 0;
  virtual float angle(int rank) = 0;
  virtual float scale(int rank) = 0;

  virtual ~RastP2D() {}
};

RastP2D *makeRastP2D();

// line segment matching using the RAST algorithm

struct RastS2D {
  // parameters
  virtual void set_maxresults(int n) = 0;
  virtual void set_verbose(bool value) = 0;
  virtual void set_tolerance(float value) = 0;
  virtual void set_eps(float eps, float aeps) = 0;
  virtual void set_scale_eps(bool value, float ieps) = 0;
  virtual void set_lsq(bool value) = 0;
  virtual void set_qtolerance(float value) = 0;
  virtual void set_min_q(float min_q) = 0;
  virtual void set_xrange(float x0, float x1) = 0;
  virtual void set_yrange(float y0, float y1) = 0;
  virtual void set_arange(float a0, float a1) = 0;  // range of angles
  virtual void set_srange(float s0, float s1) = 0;  // range of scales

  // set model segments
  virtual void clear_msources() = 0;
  virtual void add_mseg(float x0, float y0, float x1, float y1) = 0;

  // set image segments
  virtual void clear_ipoints() = 0;
  virtual void add_iseg(float x0, float y0, float x1, float y1) = 0;

  // perform match
  virtual void match() = 0;

  // get results
  virtual int nresults() = 0;
  virtual float ubound(int rank) = 0;
  virtual float lbound(int rank) = 0;
  virtual float translation(int rank, int dim) = 0;
  virtual float angle(int rank) = 0;
  virtual float scale(int rank) = 0;

  virtual ~RastS2D() {}
};

RastS2D *makeRastS2D();

// line segment matching using the RAST algorithm and sampling

struct RastSS2D {
  // parameters
  virtual void set_maxresults(int n) = 0;
  virtual void set_verbose(bool value) = 0;
  virtual void set_tolerance(float value) = 0;
  virtual void set_eps(float eps, float aeps, float sdist) = 0;
  virtual void set_lsq(bool value) = 0;
  virtual void set_qtolerance(float value) = 0;
  virtual void set_min_q(float min_q) = 0;
  virtual void set_xrange(float x0, float x1) = 0;
  virtual void set_yrange(float y0, float y1) = 0;
  virtual void set_arange(float a0, float a1) = 0;  // range of angles
  virtual void set_srange(float s0, float s1) = 0;  // range of scales

  // set model segments
  virtual void clear_msources() = 0;
  virtual void add_mseg(float x0, float y0, float x1, float y1) = 0;

  // set image segments
  virtual void clear_ipoints() = 0;
  virtual void add_iseg(float x0, float y0, float x1, float y1) = 0;

  // perform match
  virtual void match() = 0;

  // get results
  virtual int nresults() = 0;
  virtual float ubound(int rank) = 0;
  virtual float lbound(int rank) = 0;
  virtual float translation(int rank, int dim) = 0;
  virtual float angle(int rank) = 0;
  virtual float scale(int rank) = 0;

  virtual ~RastSS2D() {}
};

RastSS2D *makeRastSS2D();

// line+blob segment matching using the RAST algorithm and sampling

struct RastRS2D {
  // parameters
  virtual void set_maxresults(int n) = 0;
  virtual void set_verbose(bool value) = 0;
  virtual void set_tolerance(float value) = 0;
  virtual void set_eps(float eps, float aeps, float sdist) = 0;
  virtual void set_lsq(bool value) = 0;
  virtual void set_qtolerance(float value) = 0;
  virtual void set_min_q(float min_q) = 0;
  virtual void set_xrange(float x0, float x1) = 0;
  virtual void set_yrange(float y0, float y1) = 0;
  virtual void set_arange(float a0, float a1) = 0;  // range of angles
  virtual void set_srange(float s0, float s1) = 0;  // range of scales

  // set model segments
  virtual void clear_msources() = 0;
  virtual void add_mseg(float x0, float y0, float x1, float y1) = 0;

  // set image segments
  virtual void clear_ipoints() = 0;
  virtual void add_iseg(float x0, float y0, float x1, float y1) = 0;

  // perform match
  virtual void match() = 0;

  // get results
  virtual int nresults() = 0;
  virtual float ubound(int rank) = 0;
  virtual float lbound(int rank) = 0;
  virtual float translation(int rank, int dim) = 0;
  virtual float angle(int rank) = 0;
  virtual float scale(int rank) = 0;

  virtual ~RastRS2D() {}
};

// point matching using alignment

struct AlignmentP2D {
  // set parameters
  virtual void set_epsilon(float e) = 0;
  virtual void set_srange(float min, float max) = 0;  // range of scales

  // set model points
  virtual void clear_mpoints() = 0;
  virtual void add_mpoint(float x, float y) = 0;

  // set image points
  virtual void clear_ipoints() = 0;
  virtual void add_ipoint(float x, float y) = 0;

  // compute optimal match
  virtual void compute() = 0;

  // get optimal match
  virtual float quality() = 0;
  virtual float translation(int dim) = 0;
  virtual float angle() = 0;
  virtual float scale() = 0;

  virtual ~AlignmentP2D() {}
};

AlignmentP2D *makeAlignmentP2D();
