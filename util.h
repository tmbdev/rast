#ifndef util_h__
#define util_h__

#include "narray.h"
#include <string.h>

using namespace colib;

template <class T, class S>
inline void set(narray<T> &a, S v0) {
  a.resize(1);
  a(0) = v0;
}

template <class T, class S>
inline void set(narray<T> &a, S v0, S v1) {
  a.resize(1);
  a(0) = v0;
  a(1) = v1;
}

template <class T, class S>
inline void set(narray<T> &a, S v0, S v1, S v2) {
  a.resize(3);
  a(0) = v0;
  a(1) = v1;
  a(2) = v2;
}

template <class T, class S>
inline void set(narray<T> &a, S v0, S v1, S v2, S v3) {
  a.resize(4);
  a(0) = v0;
  a(1) = v1;
  a(2) = v2;
  a(3) = v3;
}

inline vec2 normal(const vec2 &v) { return vec2(-v[1], v[0]); }

inline vec2 cmul(const vec2 &a, const vec2 &b) {
  return vec2(a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]);
}

inline vec2 cdiv(const vec2 &a, const vec2 &b) {
  double n = sqr(b[0]) + sqr(b[1]);
  return vec2((a[0] * b[0] + a[1] * b[1]) / n, (a[1] * b[0] - a[0] * b[1]) / n);
}

inline float angleOf(const vec2 &v) { return atan2(v[1], v[0]); }

inline float normangleOf(float a) {
  while (a < 0) a += 2 * M_PI;
  while (a >= 2 * M_PI) a -= 2 * M_PI;
  return a;
}

/// namespace is created to avoid conflict with std::distance
namespace rastUtils {

inline float distance(const vec2 &a, const vec2 &b) {
  return (a - b).magnitude();
}

}

inline float norm(vec2 &v) { return v.magnitude(); }

inline float normalize_orientation(float a) {
  while (a < 0) a += M_PI;
  while (a >= M_PI) a -= M_PI;
  return a;
}

inline float normalize_angle_centered(float a) {
  while (a < -M_PI) a += 2 * M_PI;
  while (a >= M_PI) a -= 2 * M_PI;
  return a;
}

inline double urand(double low, double high) {
  return drand48() * (high - low) + low;
}

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

inline int clamp(int i, int n) {
  if (i < 0) return 0;
  if (i >= n) return n - 1;
  return i;
}

inline vec2 randomUniformVectorFromCircle(float epsilon) {
  for (;;) {
    vec2 v(urand(-1.0, 1.0), urand(-1.0, 1.0));
    if (v.magnitude() < 1.0) return v * epsilon;
  }
}

#include "colib.h"

#if 0
template <class T>
struct counted {
    struct TC : T {
	int refcount_;
    };
    TC *p;
    counted() {
	p = 0;
    }
    counted(const counted<T> &other) {
	other.incref();
	p = other.p;
    }
    counted(counted<T> &other) {
	other.incref();
	p = other.p;
    }
    ~counted() {
	decref();
	p = 0;
    }
    void operator=(counted<T> &other) {
	other.incref();
	decref();
	p = other.p;
    }
    void operator=(const counted<T> &other) {
	other.incref();
	decref();
	p = other.p;
    }
    void operator*=(counted<T> &other) {
	other.incref();
	decref();
	p = other.p;
	other.drop();
    }
    void operator*=(const counted<T> &other) {
	other.incref();
	decref();
	p = other.p;
	other.drop();
    }
    bool allocated() {
	return !p;
    }
    void allocate() {
	p = new TC();
	p->refcount_ = 1;
    }
    operator bool() {
	return !!p;
    }
    void drop() {
	decref();
	p = 0;
    }
    T &operator *() {
	if(!p) allocate();
	return *(T*)p;
    }
    T *operator->() {
	if(!p) allocate();
	return (T*)p;
    }
    operator T&() {
	if(!p) allocate();
	return *(T*)p;
    }
    operator T*() {
	if(!p) allocate();
	return (T*)p;
    }
    void incref() const {
	check();
	if(p) {
	    if(p->refcount_>10000000) abort();
	    if(p->refcount_<0) abort();
	    p->refcount_++;
	}
    }
    void decref() const {
	check();
	if(p) {
	    if(--p->refcount_==0) delete p;
	    ((counted<T>*)this)->p = 0;
	}
    }
    void check() const {
	if(!p) return;
	if(p->refcount_>10000000) abort();
	if(p->refcount_<0) abort();
    }
};

template <class T,bool use_bitswap=true>
struct heap {
    struct Item {
	double priority;
	T object;
	bool operator>(Item &other) { return priority>other.priority; }
    };
    static void bitswap(Item &a,Item &b) {
	char buf[sizeof (Item)];
	memcpy(buf,&a,sizeof (Item));
	memcpy(&a,&b,sizeof (Item));
	memcpy(&b,buf,sizeof (Item));
    }
    narray<Item> data;
    int left(int i) { return 2*i; }
    int right(int i) { return 2*i+1; }
    int parent(int i) { return int(i/2); }
    Item &A(int i) { return data.at(i-1); }
    void clear() {
	data.clear();
    }
    void heapify(int i) {
	int heapsize = data.length();
	int l = left(i);
	int r = right(i);
	int largest = -1;
	if(l<=heapsize && A(l).priority>A(i).priority)
	    largest = l;
	else
	    largest = i;
	if(r<=heapsize && A(r).priority > A(largest).priority)
	    largest = r;
	if(largest!=i) {
	    if(use_bitswap) bitswap(A(i),A(largest));
	    else swap(A(i),A(largest));
	    heapify(largest);
	}
    }
    float topPriority() {
	if(data.length()<1) throw "heap empty";
	return data.at(0).priority;
    }
    T &top() {
	if(data.length()<1) throw "heap empty";
	return data.at(0).object;
    }
    T &extractMax() {
	if(data.length()<1) throw "heap empty";
	if(use_bitswap) bitswap(A(1),A(data.length()));
	else swap(A(1),A(data.length()));
	T &result = data.pop().object;
	heapify(1);
	return result;
    }
    void insert(T &object,float priority) {
	Item &item = data.push();
	item.object = object;
	item.priority = priority;
	int i = data.length();
	while(i>1 && priority>A(parent(i)).priority) {
	    A(i) = A(parent(i));
	    i = parent(i);
	}
	A(i).object = object;
	A(i).priority = priority;
    }
    int length() {
	return data.length();
    }
};
#endif

/**
 * @brief namespace rastUtils is created to discriminate "vector" used in RAST and std::vector
 */
namespace rastUtils {
template <class T>
struct vector {
  narray<T> v;
  T &operator()(int i) { return v(i); }
  void set(T a) {
    v.resize(1);
    v(0) = a;
  }
  void set(T a, T b) {
    v.resize(2);
    v(0) = a;
    v(1) = b;
  }
  void set(T a, T b, T c) {
    v.resize(3);
    v(0) = a;
    v(1) = b;
    v(2) = c;
  }
  void set(T a, T b, T c, T d) {
    v.resize(4);
    v(0) = a;
    v(1) = b;
    v(2) = c;
    v(3) = d;
  }
  int length() { return v.length(); }
  vector<T> with(int index, T value) {
    vector<T> result;
    result = *this;
    result(index) = value;
    return result;
  }
  void copyfrom(vector<T> &source) { copy(v, source.v); }
  vector() {}
  vector(const vector<T> &source) { copy(v, source.v); }
  vector(vector<T> &source) { copy(v, source.v); }
  void operator=(const vector<T> &source) { copy(v, source.v); }
  void operator=(vector<T> &source) { copy(v, source.v); }
};
}

inline int mkseed() {
  int seed = igetenv("seed", -1);
  if (seed != -1) return seed;
  return 0;
}

#include <sys/time.h>

inline double now() {
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1e-6 * tv.tv_usec;
}

#endif
