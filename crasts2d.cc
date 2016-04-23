/* Copyright (c) 1990-1995 by Thomas M. Breuel */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "colib/misc.h"
#include "colib/narray.h"
#include "colib/vec2.h"
using namespace colib;

#include "util.h"
#include "rast.h"

namespace lumo_crasts2d {

struct Segment {
  /// two vertices of segment
  vec2 p, q;
  /// angle of dir 
  float a;
  /// directional vector
  vec2 dir;
  /// l0 = dir * p; l1 = dir * q;
  float l0, l1;
  /// normal vector
  vec2 normal;
  /// d = normal * p;
  float d;
  float weight;
  Segment() {}
  Segment(vec2 p, vec2 q) { set(p, q); }
  /// set the segment's vertices with p & q inputs and calculate all segment's parameters
  void set(vec2 p, vec2 q) {
    this->p = p;
    this->q = q;
    dir = (q - p).normalized();
    l0 = dir * p;
    l1 = dir * q;
    if (l1 < l0) swap(l0, l1);
    normal = ::normal(dir);
    d = normal * p;
    a = angleOf(dir);
  }
  float length() { return (q - p).magnitude(); }
  /// sub = p + dir * l;
  vec2 sub(float l) { return p + dir * l; }
  /// ### Check if a point is within the segment i.e. within rectangular region defined by the two vertices  
  /// along segment's dir: \f$ l_0-\epsilon-\Delta < l_a < l_1+\epsilon+\Delta \f$ 
  /// along segment's normal: \f$ err=max(0., |normal*a -d| - \Delta) < \epsilon \f$ 
  float within(float eps, float delta, vec2 a) {
    float la = dir * a;
    if (la < l0 - eps - delta || la > l1 + eps + delta) return 0;
    float err = max(0.0, fabs(normal * a - d) - delta);
    return err < eps;
    ;
  }
#if 0
	float lsq(float eps,float delta,vec2 a) {
	    float la = dir * a;
	    if(la<l0-eps-delta || la>l1+eps+delta) return 0;
	    float eps2 = sqr(eps);
	    float q = max(0.0,1.0 - sqr(max(0.0,fabs(normal * a - d) - delta))/eps2);
	    return q;
	}
#else
  ///   
  float lsq(float eps, float delta, vec2 a) {
    float la = dir * a;
    float q = 1.0;
    float dl = la - (l0 - delta);
    if (dl < -eps) return 0;
    if (dl < 0) q = 1.0 - sqr(dl * dl / (eps * eps));
    float dr = la - (l1 + delta);
    if (dr > eps) return 0;
    if (dr > 0) q = 1.0 - sqr(dr * dr / (eps * eps));
    q *= max(0.0,
             1.0 - sqr(max(0.0, fabs(normal * a - d) - delta)) / (eps * eps));
    return q; 
  }
#endif

  /**
   * @brief check whether point a & b are within two boundaries of the segment
   * @image within.png
   * @return (high - low) + 2 * delta
   * @n      where:
   * @n             low = min(la,lb);
   * @n             high = max(la,lb);
   */
  float within(float eps, float delta, vec2 a, vec2 b) {
    float la = dir * a;
    float lb = dir * b;
    float low = min(la, lb);
    float high = max(la, lb);
    if (low < l0 - eps - delta || high > l1 + eps + delta) return 0;
    float erra = max(0.0, fabs(normal * a - d) - delta);
    if (erra > eps) return 0;
    float errb = max(0.0, fabs(normal * b - d) - delta);
    if (errb > eps) return 0;
    return (high - low) + 2 * delta;
  }
  /**
   * @brief check whether point a & b are within two boundaries of the segment
   * @return uq or lq = qa * qb * ((high - low) + 2 * delta)
   * @n      where:
   * @n             qa = max(0.0, 1.0 - sqr(erra) / sqr(eps))
   * @n             qb = max(0.0, 1.0 - sqr(errb) / sqr(eps))
   */
  float lsq(float eps, float delta, vec2 a, vec2 b) {
    float la = dir * a;
    float lb = dir * b;
    float low = min(la, lb);
    float high = max(la, lb);
    if (low < l0 - eps - delta || high > l1 + eps + delta) return 0;
    float erra = max(0.0, fabs(normal * a - d) - delta);
    if (erra > eps) return 0;
    float errb = max(0.0, fabs(normal * b - d) - delta);
    if (errb > eps) return 0;
    float qa = max(0.0, 1.0 - sqr(erra) / sqr(eps));
    float qb = max(0.0, 1.0 - sqr(errb) / sqr(eps));
    return qa * qb * ((high - low) + 2 * delta);
  }
};

typedef Segment Ipoint;

typedef Segment Msource;

static float angle_diff(float a1, float a2) {
  float d = a1 - a2;
  while (d < -M_PI) d += 2 * M_PI;
  while (d > M_PI) d -= 2 * M_PI;
  return fabs(d);
}

static float unoriented_angle_diff(float a1, float a2) {
  float d = a1 - a2;
  while (d < -M_PI / 2) d += M_PI;
  while (d > M_PI / 2) d -= M_PI;
  return fabs(d);
}

struct Region {
  //vector<float> low; -> deprecated! use namespace annotation instead
  rastUtils::vector<float> low;
  //vector<float> high; -> deprecated! use namespace annotation instead
  rastUtils::vector<float> high;
  
  /// @brief represents a region of parameters. Translation value of this region i.e. the center of region 
  vec2 translation() {
    return vec2((high(0) + low(0)) / 2.0, (high(1) + low(1)) / 2.0);
  }
  
  /// @brief Angle value of this region i.e. the center of region 
  float angle() { return (high(2) + low(2)) / 2.0; }
  
  /// @brief scale value of this region i.e. the center of region 
  float scale() { return (high(3) + low(3)) / 2.0; }
  
  /// @brief rotation value, combination of rotation and scale, of this region i.e. the center of region 
  vec2 rotation() {
    float a = angle();
    float s = scale();
    return vec2(s * cos(a), s * sin(a));
  }
  
  float tdelta() {
    return 1.5 * max((high(0) - low(0)) / 2.0, (high(1) - low(1)) / 2.0);
  }
  float adelta() { return (high(2) - low(2)) / 2.0; }
  
  /// @brief scale max of region
  float smax() { return high(3); }
  /// @brief scale min of region
  float smin() { return low(3); }
  /// \f$ \Delta \f$ of scale calculated by \f${scale_{max} - scale_{min}} / 2 \f$
  float sdelta() { return (high(3) - low(3)) / 2.0; }
};

struct IMPair {
  short msource;
  short ipoint;
  IMPair() {}
  IMPair(int ms, int ip) {
    msource = ms;
    ipoint = ip;
  }
};

typedef narray<IMPair> Pairs;
typedef counted<Pairs> CPairs;
typedef narray<Msource> MsourceStack;
typedef narray<Ipoint> IpointStack;

class CRastS2D;

/// a state associates with a region, a bound and list of all matches under that bound
struct State {
  int depth;
  int generation;
  Region region;
  CPairs parent_matches;
  CPairs matches;
  float lbound;
  float ubound;

  void print(FILE *stream = stdout) {
    fprintf(stream, "<%d [%g:%g] (%g:%g %g:%g) %g:%g %g:%g>", depth, lbound,
            ubound, region.low(0), region.high(0), region.low(1),
            region.high(1), region.low(2), region.high(2), region.low(3),
            region.high(3));
  }

  void set(int depth, Region &oregion, CPairs omatches) {
    this->depth = depth;
    region = oregion;
    parent_matches = omatches;
  }

  void init(narray<Msource> &msources, narray<Ipoint> &ipoints, CRastS2D &env);
  
  /**
   * @brief evaluate a state under CRastS2D environment
   * @param[in] CRastS2D env  
   */
  void eval(CRastS2D &env);
};



/**
 * @brief main class for RastS2D but struct is used instead
 */
struct CRastS2D : RastS2D {
  rastUtils::vector<float> splitscale;

  bool final(Region &r, float delta) {
    for (int i = 0; i < r.low.length(); i++) {
      float v = (r.high(i) - r.low(i)) * double(splitscale(i));
      if (v > delta) return false;
    }
    return true;
  }

  /// @brief split the region into two lelf/right regions like binary split
  void split(Region &left, Region &right, Region &r) {
    int dim = r.low.length();
    int mi = -1;
    float mv = 0.0;
    for (int i = 0; i < dim; i++) {
      float v = (r.high(i) - r.low(i)) * splitscale(i);
      if (v < mv) continue;
      mv = v;
      mi = i;
    }
    float meanv = (r.high(mi) + r.low(mi)) / 2.0;
    left.low = r.low;
    left.high = r.high.with(mi, meanv);
    right.low = r.low.with(mi, meanv);
    right.high = r.high;
  }

  narray<Ipoint> ipoints;
  narray<float> iorient;
  narray<Msource> msources;
  narray<float> mbound;
  narray<float> mangle;
  narray<float> mabound;
  narray<bool> used;

  typedef counted<State> CState;
  heap<CState> queue;
  narray<CState> results;

  bool verbose;
  float tolerance;
  float min_q;
  int maxresults;
  rastUtils::vector<float> tlow;
  rastUtils::vector<float> thigh;
  int generation;
  bool use_lsq;
  bool unoriented;
  float eps;
  float aeps;
  float qtolerance;
  float model_total;
  bool eps_scales;
  float ieps;

  /**
   * @brief default constructor with default parameters
   * @details 
   * \n set the parameter space high/low threshold e.g. tlow, thigh, eps, aeps...
   */
  CRastS2D() {
    verbose = false;
    tolerance = 1e-3;
    min_q = 3.0;
    maxresults = 1;
    splitscale.set(1.0, 1.0, 500.0, 500.0);
    tlow.set(-1000.0, -1000.0, 0.0, 0.9);
    thigh.set(1000.0, 1000.0, 2 * M_PI, 1.1);
    generation = 1;
    use_lsq = false;
    unoriented = true;
    eps = 3.0;
    aeps = 0.05;
    qtolerance = 1e-4;
    eps_scales = 0;
    ieps = eps;
  }

  /**
   * @details priority = state->ubound + 1e-4 * state->lbound
   */
  double priority(CState state) {
    double priority = 1e30;
    priority = state->ubound + 1e-4 * state->lbound;
    if (priority >= state->ubound + 1) throw "error";
    return priority;
  }

  int n_nodes;
  int n_transforms;
  int n_distances;

  
  /**
   * @brief invokes matching 
   * @details
   * \n used has size equals all image points; At the beginning, none image point is used
   * \n 1) initially, parameter space for searching is [tlow, thigh] and all match pairs possible
   * \n 2) queue (which is a heap).insert(state, ubound); here ubound is priority of that state
   * \n 3) queue starts with 1 element i.e. initial state
   * \n 4) algorithm 1 is applied:
   * \n --> pop top priority element & eval (i.e. calc upper bound, matches) it
   * \n --> remove state with top priority
   * \n --> check stop condition to terminate   
   * \n --> if not terminate, split regions and eval them and put them to queue
   */
  void start_match() {
    n_nodes = 0;
    n_transforms = 0;
    n_distances = 0;
    results.clear();
    queue.clear();
    used.resize(ipoints.length());
    for (int i = 0; i < used.length(); i++) used[i] = false;
    
    CState initial_state;
    initial_state->init(msources, ipoints, *this);
    initial_state->region.low.copyfrom(tlow);
    initial_state->region.high.copyfrom(thigh);
    initial_state->eval(*this);
    initial_state->generation = generation;
    
    queue.insert(initial_state, initial_state->ubound);
    
    for (int iter = 0;; iter++) {
      if (queue.length() < 1) break;  // stop condition
      
      CState top;
      top = queue.extractMax();
      // top->print(); printf("\n");
      if (top->generation != generation) {
        top->eval(*this);
        top->generation = generation;
        queue.insert(top, priority(top));
        continue;
      }
      if (verbose && iter % 10000 == 0) {
        float q = results.length() > 0 ? results[0]->ubound : 0.0;
        fprintf(stderr, "# %10d result %6g queue %7d", iter, q,
                1 + queue.length());
        fprintf(stderr, "   ");
        top->print(stderr);
        fprintf(stderr, "\n");
      }
#if 0
		bool done = (1.0-qtolerance)*top->ubound<=top->lbound;
#else
      // use scale normalized termination conditions
      float subound = top->ubound / top->region.smin();
      float slbound = top->lbound / top->region.smax();
      bool done = (1.0 - qtolerance) * subound <= slbound;
#endif
      if (done || final(top->region, tolerance)) {
        results.push(top);
        Pairs &matches = top->matches;
        for (int i = 0; i < matches.length(); i++) {
          used[matches[i].ipoint] = true;
        }
        generation++;
        if (results.length() >= maxresults) return;
        continue;
      }
      Region subregions[2];
      CState substates[2];
      split(subregions[0], subregions[1], top->region);
      for (int i = 0; i < 2; i++) {
        substates[i]->set(top->depth + 1, subregions[i], top->matches);
        substates[i]->eval(*this);
        substates[i]->generation = generation;
        if (substates[i]->ubound < min_q) continue;
        queue.insert(substates[i], priority(substates[i]));
      }
    }
  }

  // defining the model and the image

  void clear_msources() { msources.clear(); }
  void add_mseg(float x, float y, float x1, float y1) {
    Msource &ms = msources.push();
    ms.set(vec2(x, y), vec2(x1, y1));
  }

  void clear_ipoints() { ipoints.clear(); }
  void add_iseg(float x, float y, float x1, float y1) {
    Ipoint &ip = ipoints.push();
    ip.set(vec2(x, y), vec2(x1, y1));
  }

  // setting match parameters

  void set_maxresults(int n) { maxresults = n; }
  void set_verbose(bool value) { verbose = value; }
  void set_eps(float eps, float aeps) {
    this->eps = eps;
    this->aeps = aeps;
  }
  void set_scale_eps(bool value, float ieps) {
    this->eps_scales = value;
    this->ieps = ieps;
  }
  void set_tolerance(float value) {
    if (value < 1e-3)
      throw "tolerance too small; would fail to converge occasionally";
    tolerance = value;
  }
  void set_min_q(float min_q) { this->min_q = min_q; }
  void set_xrange(float x0, float x1) {
    tlow(0) = x0;
    thigh(0) = x1;
  }
  void set_yrange(float y0, float y1) {
    tlow(1) = y0;
    thigh(1) = y1;
  }
  void set_arange(float a0, float a1) {
    tlow(2) = a0;
    thigh(2) = a1;
  }
  void set_srange(float s0, float s1) {
    tlow(3) = s0;
    thigh(3) = s1;
  }
  void set_lsq(bool value) { use_lsq = value; }
  void set_qtolerance(float value) { qtolerance = value; }
  void set_unoriented(bool value) { unoriented = value; }
  
  /// @brief another call to invoke match
  void match() { start_match(); }

  // reading out results

  int nresults() { return results.length(); }
  float ubound(int rank) { return results[rank]->ubound; }
  float lbound(int rank) { return results[rank]->lbound; }
  float translation(int rank, int dim) {
    return results[rank]->region.translation()[dim];
  }
  float angle(int rank) { return results[rank]->region.angle(); }
  float scale(int rank) { return results[rank]->region.scale(); }
};

/**
 * @details
 * @n lbound = 0. 
 * @n ubound = msources.length() 
 * @n matches = all matches possible -> exhaustive search
 */
void State::init(narray<Msource> &msources, narray<Ipoint> &ipoints,
                 CRastS2D &env) {
  depth = 0;
  region.low.set(0.0, 0.0, 0.0);
  region.high.set(0.0, 0.0, 0.0);
  lbound = 0.0;
  ubound = msources.length();
  Pairs &omatches = parent_matches;
  omatches.clear();
  
  // matches starts with exhaustive number of available matches i.e. all model points paired with all image points
  for (int i = 0; i < msources.length(); i++) {
    for (int j = 0; j < ipoints.length(); j++) {
      omatches.push(IMPair(i, j));
    }
  }

  float total = 0.0;
  for (int i = 0; i < msources.length(); i++) {
    total += msources[i].length();
  }
  
  // length of each model segment affects its weight
  for (int i = 0; i < msources.length(); i++) {
    msources[i].weight = msources[i].length() / total;
  }
  env.model_total = total;
}


/**
 * @details within or lsq (selected by option) will be used to check the matching of image point
 * @n       the matching point search algorithm is coded here   
 * @param[out] lbound, ubound, nmatches (referenced to matches) belonging to state
 * \n where \f$  \f$ 
 */
void State::eval(CRastS2D &env) {
  env.n_nodes++;
  MsourceStack &msources = env.msources;
  IpointStack &ipoints = env.ipoints;
  narray<bool> &used = env.used;

  Pairs &nmatches = matches;
  nmatches.clear();
  lbound = 0.0;
  ubound = 0.0;

  vec2 translation = region.translation();
  vec2 rotation = region.rotation();
  float tdelta = region.tdelta();
  float angle = region.angle();
  float adelta = region.adelta();
  float sdelta = region.sdelta();
  float smax = region.smax();

  Pairs &omatches = parent_matches;
  int n = omatches.length();
  float eps = env.eps;
  if (env.eps_scales) eps = min(smax * eps, env.ieps);
  float aeps = env.aeps;
  for (int i = 0; i < n;) {
    env.n_transforms++;
    int msource_index = omatches[i].msource;
    
    // projected segment using region parameters
    Msource &msource = msources[msource_index];   // model segment
    vec2 tmpoint0 = cmul(rotation, msource.p) + translation;  
    vec2 tmpoint1 = cmul(rotation, msource.q) + translation;
    Segment tmseg(tmpoint0, tmpoint1);  
    float tmlength = tmseg.length();
    float nmsource = max(norm(msource.p), norm(msource.q));
    float delta = tdelta + nmsource * smax * adelta + nmsource * sdelta;

    float tangle = angle + msource.a;
    float aloose = aeps + adelta;

    float llbound = 0.0;
    float lubound = 0.0;
    for (; i < n; i++) {
      if (omatches[i].msource != msource_index) break;
      env.n_distances++;
      int ipoint_index = omatches[i].ipoint;  // get the image point
      if (used[ipoint_index]) continue;
      Ipoint &ipoint = ipoints[ipoint_index];
      float adiff;
      
      // calculate angle diff
      if (env.unoriented)
        adiff = unoriented_angle_diff(ipoint.a, tangle);
      else
        adiff = angle_diff(ipoint.a, tangle);
      
      // if angle diff is too big, this ipoint is skipped
      if (adiff > aloose) continue;
      
      // check ipoint is within projected model segment tmseg
      if (!env.use_lsq) {
        float uq = tmseg.within(eps, delta, ipoint.p, ipoint.q);
        if (uq <= 0.0) continue;
        float lq = tmseg.within(eps, 0.0, ipoint.p, ipoint.q);
        lubound += uq;
        llbound += lq;
        nmatches.push(IMPair(msource_index, ipoint_index));
      } else {
        float uq = tmseg.lsq(eps, delta, ipoint.p, ipoint.q);
        if (uq <= 0.0) continue;
        float lq = tmseg.lsq(eps, 0.0, ipoint.p, ipoint.q);
        lubound += uq;
        llbound += lq;
        nmatches.push(IMPair(msource_index, ipoint_index));
      }
    }
    llbound = min(llbound, tmlength);
    lubound = min(lubound, tmlength);
#if 0
	    // fraction matched per segment
	    llbound = 100.0*llbound/tmlength*msource.weight;
	    lubound = 100.0*lubound/tmlength*msource.weight;
#endif
    lbound += llbound;
    ubound += lubound;
  }
#if 0
	// normalize matches by scale
	lbound = (lbound * 100.0) / env.model_total / region.smax();
	ubound = (ubound * 100.0) / env.model_total / region.smin();
#endif
}
}

RastS2D *makeRastS2D() { return new lumo_crasts2d::CRastS2D(); }

#ifdef TEST

using namespace lumo_crasts2d;

bool within(float x, float y, float eps) { return fabs(x - y) <= eps; }

struct seg2 {
  vec2 u, v;
  seg2() {}
  seg2(vec2 u, vec2 v) {
    this->u = u;
    this->v = v;
  }
  seg2 transform(vec2 rot, vec2 tr) {
    return seg2(cmul(rot, u) + tr, cmul(rot, v) + tr);
  }
  seg2 operator+(vec2 tr) { return seg2(u + tr, v + tr); }
  float angle() { return angleOf(v - u); }
  float length() { return (v - u).magnitude(); }
};

int main(int argc, char **argv) {
  srand48(1);
  for (int trial = 0; trial < 100; trial++) {
    bool trivial = false;
    printf("trial %d\n", trial);
    float eps = 3.0;
    float aeps = 100.0;
    vec2 tr;
    if (trivial) {
      tr = vec2(100.001, 100.001);
    } else {
      tr = vec2(urand(0.0, 512.0), urand(0.0, 512.0));
    }
    float arange = M_PI;
    float srange = 0.1;
    float alpha = urand(-arange, arange);
    float scale = urand(1.0 - srange, 1.0 + srange);
    vec2 rot = vec2(scale * cos(alpha), scale * sin(alpha));
    int nmodel = 20;
    narray<seg2> model;
    for (int i = 0; i < nmodel; i++) {
      vec2 u(urand(-100.0, 100.0), urand(-100.0, 100.0));
      vec2 v(urand(-100.0, 100.0), urand(-100.0, 100.0));
      model.push(seg2(u, v));
    }
    narray<seg2> image;
    for (int i = 0; i < model.length(); i++) {
      image.push(model[i].transform(rot, tr));
    }
    autodel<CRastS2D> rast = new CRastS2D();
    if (trivial) {
      rast->set_xrange(99, 100);
      rast->set_yrange(99, 100);
    } else {
      rast->set_xrange(0, 512);
      rast->set_yrange(0, 512);
    }
    rast->set_arange(-arange, arange);
    rast->set_srange(1.0 - srange, 1.0 + srange);
    rast->set_eps(eps, aeps, eps);
    rast->verbose = 2;
    for (int i = 0; i < image.length(); i++) {
      seg2 &s = image[i];
      rast->add_iseg(s.u[0], s.u[1], s.v[0], s.v[1]);
    }
    for (int i = 0; i < model.length(); i++) {
      seg2 &s = model[i];
      rast->add_mseg(s.u[0], s.u[1], s.v[0], s.v[1]);
    }
    rast->match();
    assert(rast->nresults() > 0);
    printf(
        "*** %3d: %8.4f %8.4f %8.4f %8.4f\n"
        "         %8.4f %8.4f %8.4f %8.4f\n"
        "    %8.4f %8.4f\n",
        trial, rast->translation(0, 0), rast->translation(0, 1), rast->angle(0),
        rast->scale(0), tr[0], tr[1], alpha, scale, rast->ubound(0),
        rast->lbound(0));
    assert(within(rast->translation(0, 0), tr[0], 2.0 * eps));
    assert(within(rast->translation(0, 1), tr[1], 2.0 * eps));
    // assert(within(rast->angle(0),alpha,0.05));
    // assert(within(rast->scale(0),scale,0.05));
  }
}
#endif
