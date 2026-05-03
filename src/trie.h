/* Copyright (c) 1990-1995 by Thomas M. Breuel */

#ifndef trie__
#define trie__

#include "util.h"
#include "narray.h"
#include "narray-util.h"
using namespace colib;

#if 0
inline int clamp(int x, int n) {
  if (x < 0) return 0;
  if (x >= n) return n-1;
  return x;
}
#endif

template <class T>
struct Trie1 {
  struct Item {
    T key;
    float x;
  };
  float eps;
  narray<narray<Item *> > buckets;
  void init(float eps, int w) {
    this->eps = eps;
    buckets.resize(int(w / eps) + 1);
  }
  void add(float x, T key) {
    int i = int(x / eps);
    Item &item = buckets(i).push();
    item.x = x;
    item.key = key;
  }
  void query(narray<T> &keys, float x0, float x1) {
    for (int i = clamp(x0 / eps, buckets.dim(0)),
             n = clamp(x1 / eps, buckets.dim(0));
         i <= n; i++) {
      for (int k = 0, r = buckets(i).length(); k < r; k++) {
        Item &item = buckets(i)[k];
        if (item.x >= x0 && item.x < x1) {
          keys.push(item.key);
        }
      }
    }
  }
};

template <class T>
struct Trie2 {
  struct Item {
    T key;
    float x, y;
  };
  float eps;
  narray<narray<Item> > buckets;
  int xoffset, yoffset;
  void init(float eps, int xmax, int ymax, int xmin = 0, int ymin = 0) {
    int w = xmax - xmin;
    int h = ymax - ymin;
    this->xoffset = xmin;
    this->yoffset = ymin;
    this->eps = eps;
    buckets.resize(int(w / eps) + 1, int(h / eps) + 1);
  }
  void add(float x, float y, T key) {
    x -= xoffset;
    y -= yoffset;
    int i = int(x / eps);
    int j = int(y / eps);
    Item &item = buckets(i, j).push();
    item.x = x;
    item.y = y;
    item.key = key;
  }
  void query(narray<T> &keys, float x0, float y0, float x1, float y1) {
    x0 -= xoffset;
    y0 -= yoffset;
    x1 -= xoffset;
    y1 -= yoffset;
    for (int i = clamp(int(x0 / eps), buckets.dim(0)),
             n = clamp(int(x1 / eps), buckets.dim(0));
         i <= n; i++) {
      for (int j = clamp(int(y0 / eps), buckets.dim(1)),
               m = clamp(int(y1 / eps), buckets.dim(1));
           j <= m; j++) {
        for (int k = 0, r = buckets(i, j).length(); k < r; k++) {
          Item &item = buckets(i, j)[k];
          if (item.x >= x0 && item.x < x1 && item.y >= y0 && item.y < y1) {
            keys.push(item.key);
          }
        }
      }
    }
  }
};

#endif
