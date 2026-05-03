// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#ifndef RAST_SRC_TRIE_H_
#define RAST_SRC_TRIE_H_

#include <vector>

#include "util.h"

template <class T>
struct Trie1 {
  struct Item {
    T key;
    float x;
  };
  float eps{0};
  std::vector<std::vector<Item>> buckets;

  void init(float eps_, int w) {
    this->eps = eps_;
    buckets.assign(int(w / eps_) + 1, {});
  }
  void add(float x, T key) {
    int i = int(x / eps);
    buckets[i].push_back({key, x});
  }
  void query(std::vector<T> &keys, float x0, float x1) {
    const int n = clamp(int(x1 / eps), int(buckets.size()));
    for (int i = clamp(int(x0 / eps), int(buckets.size())); i <= n; i++) {
      for (const auto &item : buckets[i]) {
        if (item.x >= x0 && item.x < x1) keys.push_back(item.key);
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
  float eps{0};
  int xoffset{0};
  int yoffset{0};
  int dimX{0};
  int dimY{0};
  std::vector<std::vector<Item>> buckets;  // flat: buckets[i * dimY + j]

  void init(float eps_, int xmax, int ymax, int xmin = 0, int ymin = 0) {
    int w = xmax - xmin;
    int h = ymax - ymin;
    this->xoffset = xmin;
    this->yoffset = ymin;
    this->eps = eps_;
    this->dimX = int(w / eps_) + 1;
    this->dimY = int(h / eps_) + 1;
    buckets.assign(dimX * dimY, {});
  }
  void add(float x, float y, T key) {
    x -= xoffset;
    y -= yoffset;
    int i = int(x / eps);
    int j = int(y / eps);
    buckets[i * dimY + j].push_back({key, x, y});
  }
  void query(std::vector<T> &keys, float x0, float y0, float x1, float y1) {
    x0 -= xoffset;
    y0 -= yoffset;
    x1 -= xoffset;
    y1 -= yoffset;
    const int iEnd = clamp(int(x1 / eps), dimX);
    const int jEnd = clamp(int(y1 / eps), dimY);
    for (int i = clamp(int(x0 / eps), dimX); i <= iEnd; i++) {
      for (int j = clamp(int(y0 / eps), dimY); j <= jEnd; j++) {
        for (const auto &item : buckets[i * dimY + j]) {
          if (item.x >= x0 && item.x < x1 && item.y >= y0 && item.y < y1) {
            keys.push_back(item.key);
          }
        }
      }
    }
  }
};

#endif  // RAST_SRC_TRIE_H_
