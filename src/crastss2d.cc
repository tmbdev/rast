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

namespace lumo_crastss2d {

struct Segment {
  vec2 p, q;
  float a;
  vec2 dir;
  float l0, l1;
  vec2 normal;
  float d;
  Segment() {}
  Segment(vec2 p, vec2 q) { set(p, q); }
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
  vec2 sub(float l) { return p + dir * l; }
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
  vector<float> low;
  vector<float> high;
  vec2 translation() {
    return vec2((high(0) + low(0)) / 2.0, (high(1) + low(1)) / 2.0);
  }
  float angle() { return (high(2) + low(2)) / 2.0; }
  float scale() { return (high(3) + low(3)) / 2.0; }
  vec2 rotation() {
    float a = angle();
    float s = scale();
    return vec2(s * cos(a), s * sin(a));
  }
  float tdelta() {
    return 1.5 * max((high(0) - low(0)) / 2.0, (high(1) - low(1)) / 2.0);
  }
  float adelta() { return (high(2) - low(2)) / 2.0; }
  float smax() { return high(3); }
  float sdelta() { return (high(3) - low(3)) / 2.0; }
};

struct IMPair {
#if 0
	short msource;
	short ipoint;
	float l;
#else
  unsigned char msource;
  unsigned char l;
  short ipoint;
#endif
  IMPair() {}
  IMPair(int ms, int ip, float l) {
    msource = ms;
    ipoint = ip;
    this->l = (unsigned char)l;
  }
};

typedef narray<IMPair> Pairs;
typedef counted<Pairs> CPairs;
typedef narray<Msource> MsourceStack;
typedef narray<Ipoint> IpointStack;

class CRastSS2D;

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

  void init(narray<Msource> &msources, narray<Ipoint> &ipoints, float sdist) {
    if (msources.length() > 255) throw "too many line segments in model";
    if (ipoints.length() > 32000) throw "too many line segments in image";
    depth = 0;
    region.low.set(0.0, 0.0, 0.0);
    region.high.set(0.0, 0.0, 0.0);
    lbound = 0.0;
    ubound = msources.length();
    Pairs &omatches = parent_matches;
    omatches.clear();
    for (int i = 0; i < msources.length(); i++) {
      for (float l = 0.0, ml = msources[i].length(); l < ml / 2; l += sdist) {
        // if(ml>=256.0) throw "model line segment too long";
        for (int j = 0; j < ipoints.length(); j++) {
          omatches.push(IMPair(i, j, l));
          omatches.push(IMPair(i, j, ml - l));
        }
      }
    }
  }
  void eval(CRastSS2D &env);
};

struct CRastSS2D : RastSS2D {
  vector<float> splitscale;

  bool final(Region &r, float delta) {
    for (int i = 0; i < r.low.length(); i++) {
      float v = (r.high(i) - r.low(i)) * double(splitscale(i));
      if (v > delta) return false;
    }
    return true;
  }

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
  vector<float> tlow;
  vector<float> thigh;
  int generation;
  bool use_lsq;
  bool unoriented;
  float eps;
  float aeps;
  float sdist;
  float qtolerance;

  CRastSS2D() {
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
    sdist = eps;
    qtolerance = 1e-4;
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

  void start_match() {
    n_nodes = 0;
    n_transforms = 0;
    n_distances = 0;
    results.clear();
    queue.clear();
    used.resize(ipoints.length());
    for (int i = 0; i < used.length(); i++) used[i] = false;
    CState initial_state;
    initial_state->init(msources, ipoints, sdist);
    initial_state->region.low.copyfrom(tlow);
    initial_state->region.high.copyfrom(thigh);
    initial_state->eval(*this);
    initial_state->generation = generation;
    queue.insert(initial_state, initial_state->ubound);
    for (int iter = 0;; iter++) {
      if (queue.length() < 1) break;
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
      if ((1.0 - qtolerance) * top->ubound <= top->lbound ||
          final(top->region, tolerance)) {
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
  void set_eps(float eps, float aeps, float sdist) {
    this->eps = eps;
    this->aeps = aeps;
    this->sdist = sdist;
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

void State::eval(CRastSS2D &env) {
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
  float aeps = env.aeps;
  for (int i = 0; i < n;) {
    env.n_transforms++;
    int msource_index = omatches[i].msource;
    Msource &msource = msources[msource_index];
    vec2 tmpoint0 = cmul(rotation, msource.p) + translation;
    vec2 tmpoint1 = cmul(rotation, msource.q) + translation;
    Segment tmseg(tmpoint0, tmpoint1);
    float nmsource = max(norm(msource.p), norm(msource.q));
    float delta = tdelta + nmsource * smax * adelta + nmsource * sdelta;

    float tangle = angle + msource.a;
    float aloose = aeps + adelta;

    for (; i < n;) {
      if (omatches[i].msource != msource_index) break;
      float start_l = omatches[i].l;
      float llbound = 0.0;
      float lubound = 0.0;
      for (; i < n; i++) {
        float l = omatches[i].l;
        if (l != start_l) break;
        vec2 tmpoint = tmseg.sub(l);
        // printf(".. %g %g\n",l,(tmseg.p-tmpoint).magnitude());
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
        if (!env.use_lsq) {
          float uq = ipoint.within(eps, delta, tmpoint);
          if (uq <= 0.0) continue;
          float lq = ipoint.within(eps, 0.0, tmpoint);
          lubound = max(lubound, uq);
          llbound = max(llbound, lq);
          nmatches.push(IMPair(msource_index, ipoint_index, l));
        } else {
          float uq = ipoint.lsq(eps, delta, tmpoint);
          if (uq <= 0.0) continue;
          float lq = ipoint.lsq(eps, 0.0, tmpoint);
          lubound = max(lubound, uq);
          llbound = max(llbound, lq);
          nmatches.push(IMPair(msource_index, ipoint_index, l));
        }
      }
      lbound += llbound;
      ubound += lubound;
    }
  }
}
}

RastSS2D *makeRastSS2D() { return new lumo_crastss2d::CRastSS2D(); }

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
    stack<seg2> model;
    for (int i = 0; i < nmodel; i++) {
      vec2 u(urand(-100.0, 100.0), urand(-100.0, 100.0));
      vec2 v(urand(-100.0, 100.0), urand(-100.0, 100.0));
      model.push(seg2(u, v));
    }
    stack<seg2> image;
    for (int i = 0; i < model.length(); i++) {
      image.push(model[i].transform(rot, tr));
    }
    autodel<CRastSS2D> rast = new CRastSS2D();
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
