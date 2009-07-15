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


namespace lumo_clinesp2d {


    inline float min(float a,float b,float c) {
	return ::min(a,::min(b,c));
    }

    inline float hside(float a) {
	return a<0?0:a;
    }

    struct Point {
	vec2 p;
	float a;
	float w;
    };

#ifdef TEST
    float test_angle,test_offset;
#endif

    struct LineRegion {
	float th0,ux0,uy0;
	float th1,ux1,uy1;
	float thm,uxm,uym;
	float r0,r1,rm;
	float rerr,therr;

	void print(FILE *stream=stdout) {
	    fprintf(stream,"<LineRegion %g %g   %g %g>",th0,th1,r0,r1);
	}

	void set(float th0,float th1,float r0,float r1) {
	    if(th1<=th0) throw "parameters (th)";
	    if(r1<=r0) throw "parameters (r)";
	    this->th0 = th0;
	    this->th1 = th1;
	    this->r0 = r0;
	    this->r1 = r1;
	    ux0 = cos(th0);
	    uy0 = sin(th0);
	    ux1 = cos(th1);
	    uy1 = sin(th1);
	    thm = ((th0+th1)/2);
	    uxm = cos(thm);
	    uym = sin(thm);
	    double factor = cos((th1-th0)/2);
	    this->rm = max(0.0,r0 * factor);
	    rerr = (r1-r0)/2.0;
	    therr = (th1-th0)/2.0;
	}

	float dist(float x,float y) {
	    float dot0 = ux0*x+uy0*y;
	    float dot1 = ux1*x+uy1*y;
	    float d0 = dot0-r1;
	    float d1 = dot1-r1;
	    float dotm = uxm*x+uym*y;
	    float d2 = -(dot0-r0);
	    float d3 = -(dotm-rm);
	    float d4 = -(dot1-r0);
	    float upper = 0.0;
	    if(d0>=0 && d1>=0) upper = ::min(d0,d1);
	    float lower = 0.0;
	    if(d2>=0 && d3>=0 && d4>=0) lower = min(d2,d3,d4);
	    return max(upper,lower);
	}

	// compute a lower bound on the difference between a line
	// in this line region and the given angle
	float adist(float a,bool unoriented) {
	    float diff = a-thm;
	    if(unoriented) {
		while(diff<-M_PI/4) diff += M_PI/2;
		while(diff>M_PI/4) diff -= M_PI/2;
		diff = fabs(diff);
		diff -= therr;
		if(diff<0) return 0.0;
		return diff;
	    } else {
		while(diff<-M_PI/2) diff += M_PI;
		while(diff>M_PI/2) diff -= M_PI;
		diff = fabs(diff);
		diff -= therr;
		if(diff<0) return 0.0;
		return diff;
	    }
	}

	void split(LineRegion *result,float thscale) {
	    float dr = r1-r0;
	    float da = thscale * (th1-th0);
	    if(da>dr) {
		float m = (th0+th1)/2;
		result[0].set(th0,m,r0,r1);
		result[1].set(m,th1,r0,r1);
	    } else {
		float m = (r0+r1)/2;
		result[0].set(th0,th1,r0,m);
		result[1].set(th0,th1,m,r1);
	    }
	}
    };

    typedef narray<int> Matches;
    typedef counted<Matches> CMatches;

    struct State {
	int generation;
	float weight;
	LineRegion region;
	CMatches matches;
    };

    typedef counted<State> CState;

    struct CLinesP2D : LinesP2D {
	float eps,aeps;
	float tol,atol;
	float minweight;
	float maxoffset;
	int maxresults;
	int generation;
	int verbose;
	bool unoriented;
	bool lsq;

	CLinesP2D() {
	    eps = 2.0; aeps = 0.05;
	    tol = 0.1; atol = 0.001;
	    minweight = 0.0;
	    maxoffset = 3000.0;
	    maxresults = 1;
	    generation = 0;
	    verbose = 0;
	    lsq = 1;
	    unoriented = 1;
	}

	narray<Point> points;
	narray<bool> used;

	heap<CState> queue;
	narray<CState> results;

	void filter(CState &state) {
	    LineRegion &region = state->region;
	    Matches &input = state->matches;
	    CMatches cresult;
	    Matches &result = cresult;
	    float weight = 0.0;
	    float eps2 = eps*eps;
	    float aeps2 = aeps*aeps;
	    for(int i=0;i<input.length();i++) {
		int index = input[i];
		if(used[index]) continue;
		Point p = points[index];
		float q = 1.0;
		float da = region.adist(p.a,unoriented);
		if(lsq) {
		    q *= max(0.0,1.0 - da*da/aeps2);
		} else {
		    if(da>aeps) q = 0.0;
		}
		if(q==0.0) continue;
		float d = region.dist(p.p[0],p.p[1]);
		if(lsq) {
		    q *= max(0.0,1.0 - d*d/eps2);
		} else {
		    if(d>eps) q = 0.0;
		}
		if(q==0.0) continue;
		q *= p.w;
		weight += q;
		result.push(index);
	    }
	    state->weight = weight;
	    state->matches = cresult;
	}

	void compute() {
	    if(verbose) fprintf(stderr,"[#segments %d]\n",points.length());
	    generation = 0;
	    used.resize(points.length());
	    fill(used,false);

#if 1
	    for(int i=0;i<8;i++) {
		CState initial;
		for(int j=0;j<points.length();j++) initial->matches->push(j);
		if(unoriented)
		    initial->region.set(i*M_PI/4,(i+1)*M_PI/4,0.0,maxoffset);
		else
		    initial->region.set(i*M_PI/4,(i+1)*M_PI/4,-maxoffset,maxoffset);
		filter(initial);
		queue.insert(initial,initial->weight);
	    }
#else
	    {
		CState initial;
		for(int j=0;j<points.length();j++) initial->matches->push(j);
		if(unoriented)
		    initial->region.set(0,2*M_PI,0.0,maxoffset);
		else
		    initial->region.set(0,2*M_PI,-maxoffset,maxoffset);
		filter(initial);
		queue.insert(initial,initial->weight);
	    }
#endif
	    
	    for(int iter=0;;iter++) {
		if(queue.length()<1) break;
		CState state;
		state = queue.extractMax();
		if(state->generation!=generation) {
		    filter(state);
		    state->generation = generation;
		    queue.insert(state,state->weight);
		    continue;
		}
		LineRegion &pregion = state->region;
		if(verbose>1) {
		    pregion.print(stderr); fprintf(stderr,"(%g)\n",state->weight);
		} else if(verbose>0) {
		    if(iter%1000==0 && queue.length()>0) {
			fprintf(stderr,"[%d %d %g %d]\n",
				iter,queue.length(),
				queue.topPriority(),results.length());
		    }
		}

		if(state->weight<minweight) continue;

		if(pregion.rerr<tol && pregion.therr<atol) {
		    results.push(state);
		    Matches &matches = state->matches;
		    for(int i=0;i<matches.length();i++) {
			used[matches[i]] = true;
			Point &p = points[matches[i]];
			if(verbose) printf("P %g %g %g\n",p.p[0],p.p[1],p.a);
		    }
		    generation++;
		    if(results.length()>=maxresults) break;
		    continue;
		}

		LineRegion regions[2];
		state->region.split(regions,tol/atol);
		for(int i=0;i<2;i++) {
		    CState child;
		    child->matches = state->matches;
		    child->weight = state->weight;
		    child->region = regions[i];
		    filter(child);
		    queue.insert(child,child->weight);
		}
	    }
	}

	void clear_ipoints() {
	    points.clear();
	}
	void add_ipoint(float x,float y,float a,float w=1.0) {
	    // we are only using "short" for image point indexes,
	    // so make sure we don't get too many
	    //if(points.length()>32000) throw "too many image points";
	    Point &p = points.push();
	    p.p = vec2(x,y);
	    p.a = a;
	    p.w = w;
	}

	void set_maxresults(int n) {
	    maxresults = n;
	}
	void set_minweight(float value) {
	    minweight = value;
	}
	void set_maxoffset(float value) {
	    maxoffset = value;
	}
	void set_lsq(bool value) {
	    lsq = value;
	}
	void set_unoriented(bool value) {
	    unoriented = value;
	}
	void set_error(float eps,float aeps) {
	    this->eps = eps;
	    this->aeps = aeps;
	}
	void set_tolerance(float tol,float atol) {
	    this->tol = tol;
	    this->atol = atol;
	}
	void set_breakpenalty(float eps,float cost) {
	    //throw "not implemented";
	}
	void set_verbose(int value) {
	    verbose = value;
	}

	int nresults() {
	    return results.length();
	}
	int nmatches(int i) {
	    return results[i]->matches->length();
	}
	float weight(int i) {
	    return results[i]->weight;
	}
	float angle(int i) {
	    return results[i]->region.thm;
	}
	float offset(int i) {
	    return results[i]->region.rm;
	}
    };

}

LinesP2D *makeLinesP2D() {
    return new lumo_clinesp2d::CLinesP2D();
}

#ifdef TEST
int main(int argc,char **argv) {
    autodel<CLinesP2D> lines = new CLinesP2D();
    lines->set_verbose(1);
    lines->set_error(1.0,8.0);
    float cx = 256, cy = 256;
    test_angle = atof(argv[1]);
    test_offset = -sin(test_angle)*cx + cos(test_angle)*cy;
    printf("test_angle=%g test_offset=%g\n",test_angle,test_offset);
    for(double l = -200.0;l<=200.5;l+=2.0) {
	float x = cx+l*cos(test_angle);
	float y = cy+l*sin(test_angle);
	lines->add_ipoint(x,y,test_angle,1.0);
    }
    lines->compute();
    printf("weight=%g   angle=%g offset=%g\n",lines->weight(0),lines->angle(0),lines->offset(0));
}
#endif
