// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#ifndef RAST_SRC_HEAP_H_
#define RAST_SRC_HEAP_H_

#include <cassert>
#include <cstddef>
#include <queue>
#include <utility>
#include <vector>

// Max-heap priority queue, thin wrapper over std::priority_queue.
// Each entry carries an explicit priority alongside the value.

template <class T>
struct Heap {
  struct Item {
    double priority;
    T value;
    bool operator<(const Item &other) const { return priority < other.priority; }
  };

  std::priority_queue<Item> q;

  std::size_t length() const { return q.size(); }
  bool empty() const { return q.empty(); }
  void clear() { q = std::priority_queue<Item>{}; }

  void insert(T value, double priority) { q.push({priority, std::move(value)}); }

  T extractMax() {
    assert(!q.empty());
    T result = q.top().value;
    q.pop();
    return result;
  }

  double topPriority() const {
    assert(!q.empty());
    return q.top().priority;
  }
  const T &top() const {
    assert(!q.empty());
    return q.top().value;
  }
};

#endif  // RAST_SRC_HEAP_H_
