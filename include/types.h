#ifndef TYPES_H
#define TYPES_H
#include <array>
#include <oopsconfig.h>

// A pair of ordered pairs.
template<class T> using pair2 = std::array<std::array<T, 2>, 2>;

// A general enumerator to describe error codes.
enum Result{
  SUCCESS,
  FAILURE,
  BAD_ALLOC,
  OUT_OF_BOUNDS,
  INVALID_STAGE,
  UNRECOGNIZED_PARAMS,
  FIELD_EXISTS,
  UNRECOGNIZED_FIELD,
  UNINITIALIZED,
  UNEQUAL_SPACING,
};

enum Boundary{
  NONE,
  LEFT,
  RIGHT,
  UP,
  DOWN
};

// If we're using g++, we have access to its restrict keyword, which can be useful
// for optimization. To keep the code from breaking, we still define it for other
// compilers even though it won't do anything.
#ifndef __GNUG__
#define RESTRICT
#else
#define RESTRICT __restrict__
#endif

#endif
