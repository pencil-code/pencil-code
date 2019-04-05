#pragma once

template <class T> 
static inline const T max(const T& a, const T& b) {
  return a > b ? a : b;
}

template <class T> 
static inline const T min(const T& a, const T& b) {
  return a < b ? a : b;
}
/*
template <class T> 
static inline const T min(const T& a, T b) {
  return a < b ? a : b;
}

template <class T> 
static inline const T min(T a, const T& b) {
  return a < b ? a : b;
}

template <class T> 
static inline const T min(T a, T b) {
  return a < b ? a : b;
}
*/
template <class T> 
static inline const T sum(const T& a, const T& b) {
  return a + b;
}
