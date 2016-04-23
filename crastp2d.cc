/* Copyright (c) 1990-1995 by Thomas M. Breuel */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "misc.h"
#include "narray.h"
#include "vec2.h"
using namespace colib;

#include "util.h"
#include "rast.h"

namespace lumo_crastp2d {

struct Ipoint {
  vec2 p;
  float a;
};

struct Msource {
  vec2 p;
  float a;
  float eps;
  float aeps;
};

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
  narray<float> low;
  narray<float> high;
  void operator=(Region &other) {
    copy(low, other.low);
    copy(high, other.high);
  }
  vec2 translation() {
    return vec2((high(0) + low(0)) / 2.0, (high(1) + low(1)) / 2.0);
  }
  float angle() { return (high(2) + low(2)) / 2.0; }
  float scale() { return (high(3) + low(3)) / 2.0; }
  
  /// @return vec2(s * cos(1), s * sin(a))
  vec2 rotation() {
    float a = angle();
    float s = scale();
    return vec2(s * cos(a), s * sin(a));
  }
  /**
   * @details
   * \n original equation: \f$ t^{\delta} = \sqrt{(t^d_x)^2 + (t^d_y)^2} \f$
   * \n approx. equation:
   * \n with \f$ t^d = max(t^d_x, t^d_y) \f$ 
   * \n \f$ \iff t^{\delta} = \sqrt{2}.\sqrt{(t^d)^2}  \f$ 
   * \n \f$ \iff t^{\delta} = 1.4|t^d|  \f$ 
   */
  float tdelta() {
    return 1.5 * max((high(0) - low(0)) / 2.0, (high(1) - low(1)) / 2.0);
  }
  float adelta() { return (high(2) - low(2)) / 2.0; }
  float smax() { return high(3); }
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

class CRastP2D;

struct State {
  int depth;
  int generation;
  Region region;
  CPairs parent_matches;
  CPairs matches;
  float lbound;
  float ubound;

  void print(FILE *stream) {
    fprintf(stream, "<%d u %g l %g low %g %g %g %g high %g %g %g %g>", depth,
            ubound, lbound, region.low(0), region.low(1), region.low(2),
            region.low(3), region.high(0), region.high(1), region.high(2),
            region.high(3));
  }

  void set(int depth, Region &oregion, CPairs omatches) {
    this->depth = depth;
    region = oregion;
    parent_matches = omatches;
  }

  void init(narray<Msource> &msources, narray<Ipoint> &ipoints) {
    depth = 0;
    ::set(region.low, 0.0, 0.0, 0.0);
    ::set(region.high, 0.0, 0.0, 0.0);
    lbound = 0.0;
    ubound = msources.length();
    Pairs &omatches = parent_matches;
    omatches.clear();
    for (int i = 0; i < msources.length(); i++) {
      for (int j = 0; j < ipoints.length(); j++) {
        omatches.push(IMPair(i, j));
      }
    }
  }
  void eval(CRastP2D &env);
};

struct CRastP2D : RastP2D {
  
  /// provides a mean to relatively compare between different units i.e. different dimensions of the search space
  narray<float> splitscale;
  
  /**
   * @brief    check terminate condition that all dimensions are smaller than delta threshold 
   * @details  v (magnitude of dimension i) = (r.high(i) - r.low(i)) * double(splitscale(i))  
   */
  bool final(Region &r, float delta) {
    for (int i = 0; i < r.low.length(); i++) {
      float v = (r.high(i) - r.low(i)) * double(splitscale(i));
      if (v > delta) return false;
    }
    return true;
  }

  /// @brief split the space at the largest dimension
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
    copy(left.low, r.low);
    copy(left.high, r.high);
    left.high(mi) = meanv;
    copy(right.low, r.low);
    right.low(mi) = meanv;
    copy(right.high, r.high);
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
  narray<float> tlow;
  narray<float> thigh;
  int generation;
  bool use_lsq;
  bool unoriented;

  CRastP2D() {
    verbose = false;
    tolerance = 1e-3;
    min_q = 3.0;
    maxresults = 1;
    ::set(splitscale, 1.0, 1.0, 700.0, 10.0);
//    ::set(tlow, -1000.0, -1000.0, 0.0, 0.9);
    ::set(tlow, -15.0, -15.0, 0.0, 0.95);
//    ::set(thigh, 1000.0, 1000.0, 2 * M_PI, 1.1);
    ::set(thigh, 15.0, 15.0, M_PI/4., 1.05);
    generation = 1;
    use_lsq = false;
    unoriented = true;
  }

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
   * @details the stop condition is (1) top.ubound = top.lbound or (2) "final" is reached.
   *          number matches defined in maxresults 
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
    initial_state->init(msources, ipoints);
    copy(initial_state->region.low, tlow);
    copy(initial_state->region.high, thigh);
    initial_state->eval(*this);
    initial_state->generation = generation;
    queue.insert(initial_state, initial_state->ubound);
    for (int iter = 0;; iter++) {
      if (queue.length() < 1) break;
      CState top;
      top = queue.extractMax();
      if (top->generation != generation) {
        top->eval(*this);
        top->generation = generation;
        queue.insert(top, priority(top));
        continue;
      }
      if (verbose && iter % 1000 == 0) {  // 10000
        float q = results.length() > 0 ? results[0]->ubound : 0.0;
        fprintf(stderr, "# %10d result %6g queue %7d", iter, q,
                1 + queue.length());
        fprintf(stderr, "   ");
        top->print(stderr);
        fprintf(stderr, "\n");
      }
      if (top->ubound == top->lbound || final(top->region, tolerance)) {
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
  void add_msource(float x, float y, float a, float eps, float aeps) {
    Msource &ms = msources.push();
    ms.p = vec2(x, y);
    ms.a = a;
    ms.eps = eps;
    ms.aeps = aeps;
  }

  void clear_ipoints() { ipoints.clear(); }
  void add_ipoint(float x, float y, float a) {
    Ipoint &ip = ipoints.push();
    ip.p = vec2(x, y);
    ip.a = a;
  }

  // setting match parameters

  void set_maxresults(int n) { maxresults = n; }
  void set_verbose(bool value) { verbose = value; }
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
  void set_unoriented(bool value) { unoriented = value; }
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
 * \n loose = \f$ \epsilon + \delta \f$ 
 * \n delta = delta(tdelta, adelta, sdelta) where tdelta, adelta, sdelta are calculated from region
 */
void State::eval(CRastP2D &env) {
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
  for (int i = 0; i < n;) {
    env.n_transforms++;
    int msource_index = omatches[i].msource;
    Msource &msource = msources[msource_index];
    vec2 tmpoint = cmul(rotation, msource.p) + translation;
    float eps2 = sqr(msource.eps);
    float nmsource = norm(msource.p);
    float strict = msource.eps;
    float delta = tdelta + nmsource * smax * adelta + nmsource * sdelta;
    float loose = strict + delta;

    float tangle = angle + msource.a;
    float astrict = msource.aeps;
    float aloose = astrict + adelta;

    float llbound = 0.0;
    float lubound = 0.0;

    for (; i < n && omatches[i].msource == msource_index; i++) {
      env.n_distances++;
      int ipoint_index = omatches[i].ipoint;
      if (used[ipoint_index]) continue;
      Ipoint &ipoint = ipoints[ipoint_index];
      float adiff;
      if (env.unoriented)
        adiff = unoriented_angle_diff(ipoint.a, tangle);
      else
        adiff = angle_diff(ipoint.a, tangle);
      if (adiff > aloose) continue;
      float err = rastUtils::distance(tmpoint, ipoint.p);
      if (err > loose) continue;
      if (env.use_lsq) {
        float ud = max(0.0F, err - delta);
        float uq = max(0.0F, 1.0F - ud * ud / eps2);
        lubound = max(lubound, uq);
        float ld = err;
        float lq = max(0.0F, 1.0F - ld * ld / eps2);
        llbound = max(llbound, lq);
      } else {
        if (err < strict && adiff < astrict) llbound = 1.0;
        lubound = 1.0;
      }
      nmatches.push(IMPair(msource_index, ipoint_index));
    }
    lbound += llbound;
    ubound += lubound;
  }
}
}

RastP2D *makeRastP2D() { return new lumo_crastp2d::CRastP2D(); }
