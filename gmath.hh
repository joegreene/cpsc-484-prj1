
//
// gmath.hh
//
// graphics math library.
//
// CPSC 484, CSU Fullerton, Spring 2016, Prof. Kevin Wortman
// Project 1
//
// Name:
//   Kyle Terrien
//   Adam Beck
//   Joe Greene
//
// In case it ever matters, this file is hereby placed under the MIT
// License:
//
// Copyright (c) 2016, Kevin Wortman
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#pragma once

#include <algorithm> // for swap()
#include <cassert>   // for assert()
#include <cmath>     // for acos(), sqrt()
#include <memory>    // for shared_ptr
#include <iomanip>   // for setw
#include <iostream>  // for cout

namespace gmath {

  // Data type for an angle; default to double for compatibility with
  // <cmath> trigonometry functions.
  typedef double angle_type;

  // Multiplicative error tolerated by the approximate_equal function
  // below.
  const double APPROXIMATE_EQUAL_EPSILON = .001; // 0.1 percent

  // Return true if left and right are approximately equal, up to a
  // multiplicative factor equal to APPROXIMATE_EQUAL_EPSILON.
  //
  // This function exists so that floating-point values can be compare
  // for equality, without being confused by the small rounding errors
  // that are inevitable with floating-point arithmetic.
  template <typename SCALAR>
  bool approximate_equal(SCALAR left, SCALAR right) {
    if (left == right) {
      // If the very strict == operator says the values are equal,
      // they're definitely equal.
      return true;
    }
    if (right == 0) {
      // Make sure right is nonzero, to avoid a potential
      // divide-by-zero in the next statement. Note that, if both
      // arguments were zero, we would've already returned true, so if
      // we get to this line, at least one argument must be nonzero.
      std::swap(left, right);
    }
    assert(right != 0);
    // Compute the relative error between the arguments.
    SCALAR delta = (left - right) / right;
    // Make sure delta is non-negative. (We don't use the abs()
    // function because it may not be defined for whatever data type
    // SCALAR happens to be.)
    if (delta < 0) {
      delta = -delta;
    }
    assert(delta >= 0);
    return (delta <= APPROXIMATE_EQUAL_EPSILON);
  }

  // Vector<SCALAR, DIMENSION> represents a mathematical vector of
  // DIMENSION length, where each base element is of type
  // SCALAR. These vectors' dimension (length) is fixed; they are not
  // resizeable the way the std::vector data structure is. Making the
  // DIMENSION a const template parameter means that the compiler can
  // type-check for incorrect vector or matrix dimensions at
  // compile-time.
  template <typename SCALAR, const int DIMENSION>
  class Vector {
  public:
    
    // Alias for this very type.
    typedef Vector<SCALAR, DIMENSION> same_type;

    // Alias for a shared_ptr to same_type.
    typedef std::shared_ptr<same_type> ptr_type;

  private:
    SCALAR _elements[DIMENSION];
  
  public:
    // Initialize all elements to initializer, 0 by default.
    Vector(SCALAR initializer = 0) {
      // TODO: implement this function
    }

    // Copy constructor.
    Vector(const same_type& other) {
      // TODO: implement this function
    }

    // Equality comparison operator; uses approximate_equal to compare
    // individual elements.
    bool operator== (const same_type& right) const {
      for (int i = 0; i < DIMENSION; ++i) {
	if (!approximate_equal(_elements[i], right._elements[i])) {
	  return false;
	}
      }
      return true;
    }

    // Equality comparison operator with a ptr_type.
    bool operator== (ptr_type right) const {
      return (*this) == (*right);
    }

    
    // Inequality comparison operator; uses approximate_equal to
    // compare individual elements.
    bool operator!= (const same_type& right) const {
      return !(*this == right);
    }

    // Inequality comparison operator with a ptr_type.
    bool operator!= (ptr_type right) const {
      return (*this) != (*right);
    }

    // Dereference operator, to obtain a const reference to an
    // element.
    const SCALAR& operator[] (int index) const {
      assert(is_index(index));
      // TODO: implement this function
    }
  
    // Dereference operator, to obtain a non-const reference to an
    // element.
    SCALAR& operator[] (int index) {
      assert(is_index(index));
      // TODO: implement this function
    }

    // Assignment from a scalar. Overwrites each element with s.
    same_type& operator= (SCALAR s) {
      // TODO: implement this function
      return (*this);
    }

    // Assignment from another vector (copy).
    same_type& operator= (const same_type& right) {
      // TODO: implement this function
      return (*this);
    }

    // Assignment from a pointer to another vector (copy).
    same_type& operator= (ptr_type right) {
      return (*this) = (*right);
    }

    // Addition with a vector.
    ptr_type operator+ (const same_type& right) const {
      // TODO: implement this function
    }

    // Addition with a pointer to a vector.
    ptr_type operator+ (const ptr_type right) const {
      return (*this) + (*right);
    }

    // Subtraction with a vector.
    ptr_type operator- (const same_type& right) const {
      // TODO: implement this function
    }

    // Subtraction with a pointer to a vector.
    ptr_type operator- (const ptr_type right) const {
      return (*this) - (*right);
    }

    // Negation operator. Negates (toggles sign of) each element of
    // the vector.
    ptr_type operator- () const {
      // TODO: implement this function

      // hint: reuse other overloaded operators, instead of writing
      // this from scratch
    }

    // Multiply by a scalar.
    ptr_type operator* (SCALAR s) const {
      // TODO: implement this function
    }

    // Multiply by another vector (dot product).
    SCALAR operator* (const same_type& right) const {
      // TODO: implement this function
    }

    // Multiply by a pointer to another vector (dot product).
    SCALAR operator* (const ptr_type right) const {
      return (*this) * (*right);
    }

    // Division by a scalar.
    ptr_type operator/ (SCALAR s) const {
      // TODO: implement this function
    }

    // Return true iff i is a valid element index.
    bool is_index(int i) const {
      return ((i >= 0) && (i < DIMENSION));
    }

    // Return true iff the homogeneous element indicates this is a
    // point. (I.e. if the last element is nonzero.)
    bool is_homogeneous_point() const {
      assert(DIMENSION > 1);
      // TODO: implement this function
    }
    
    // Return true iff the homogeneous element indicates this is a
    // translation vector. (I.e. if the last element is zero.)
    bool is_homogeneous_translation() const {
      assert(DIMENSION > 1);
      // TODO: implement this function
    }

    // Return true iff all elements are zero.
    bool is_zero() const {
      // TODO: implement this function
    }

    // Return the angle, in radians, between this vector and
    // another. This function assumes that the two vectors each
    // represent a point relative to the origin.
    angle_type angle(const same_type& right) const {
      // TODO: implement this function

      // hint: reuse other overloaded operators, instead of writing
      // this from scratch
    }

    // Compute the cross product between this and another vector. The
    // cross product is only defined for certain dimensions; this
    // function asserts that DIMENSION is 3.
    ptr_type cross(const same_type& right) const {
      assert(DIMENSION == 3);
      // TODO: implement this function
    }

    // Return the dimension of this vector.
    int dimension() const {
      return DIMENSION;
    }

    // Return the scalar distance between this and another
    // vector. This function assumes that each vector represents a
    // point. The returned value is always a non-negative scalar.
    SCALAR distance(const same_type& p) const {
      // TODO: implement this function

      // hint: reuse other overloaded operators, instead of writing
      // this from scratch
    }

    // Return the magnitude of this vector, i.e. for vector v, the
    // quantity denoted ||v|| .
    SCALAR magnitude() const {
      // TODO: implement this function

      // hint: reuse other overloaded operators, instead of writing
      // this from scratch
    }

    // Return a pointer to a normalized version of this vector,
    // i.e. another vector with identical distance, and magnitude
    // exactly 1.
    ptr_type normalized() const {
      assert(!is_zero());

      // TODO: implement this function

      // hint: reuse other overloaded operators, instead of writing
      // this from scratch
    }

    // Utility function to print a representation of this vector to
    // cout. This is intende for debugging purposes.
    void print() {
      std::cout << "<";
      for (int i = 0; i < DIMENSION; ++i) {
	if (i > 0) {
	  std::cout << ", ";
	}
	std::cout << _elements[i];
      }
      std::cout << ">";       
    }

    // Return the homogeneous element of this vector (i.e. the last
    // element).
    SCALAR w() const {
      // TODO: implement this function
    }
  };
  
  // Matrix<SCALAR, HEIGHT, WIDTH> represents a mathematical vector of
  // HEIGHT x WIDTH dimension, where each base element is of type
  // SCALAR.
  template <typename SCALAR, const int HEIGHT, const int WIDTH>
  class Matrix {
  public:

    // Alias for this very type.
    typedef Matrix<SCALAR, HEIGHT, WIDTH> same_type;

    // Alias for a shared_ptr to same_type.
    typedef std::shared_ptr<same_type> ptr_type;

    // Data type of one column of this matrix.
    typedef Vector<SCALAR, HEIGHT> column_type;

    // Data type of one row of this matrix.
    typedef Vector<SCALAR, WIDTH> row_type;

  private:
    row_type _rows[HEIGHT];

  public:

    // Initialize all elements to initializer, 0 by default.
    Matrix(SCALAR initializer = 0) {
      // TODO: implement this function

      // hint: reuse other overloaded operators, instead of writing
      // this from scratch
    }

    // Copy constructor.
    Matrix(const same_type& other) {
      // TODO: implement this function

      // hint: reuse other overloaded operators, instead of writing
      // this from scratch
    }

    // Equality comparison operator; uses approximate_equal to compare
    // individual elements.
    bool operator== (const same_type& right) const {
      // TODO: implement this function
    }

    // Equality comparison operator with a ptr_type.
    bool operator== (ptr_type right) const {
      return (*this) == (*right);
    }

    // Inequality comparison operator; uses approximate_equal to
    // compare individual elements.
    bool operator!= (const same_type& right) const {
      return !( (*this) == right );
    }

    // Inequality comparison operator with a ptr_type.
    bool operator!= (ptr_type right) const {
      return (*this) != (*right);
    }
    
    // Dereference operator, to obtain a const reference to a row.
    const row_type& operator[] (int row) const {
      assert(is_row(row));
      return _rows[row];
    }

    // Dereference operator, to obtain a non-const reference to row.
    row_type& operator[] (int row) {
      assert(is_row(row));
      return _rows[row];
    }

    // Assignment from a scalar. Overwrites each element with s.
    same_type& operator= (SCALAR s) {
      // TODO: implement this function
    }

    // Assignment from another vector (copy).
    same_type& operator= (const same_type& right) {
      // TODO: implement this function
    }

    // Addition with a matrix.
    ptr_type operator+ (const same_type& right) const {
      // TODO: implement this function
    }

    // Addition with a pointer to a matrix.
    ptr_type operator+ (ptr_type right) const {
      return (*this) + (*right);
    }

    // Subtraction with a vector.
    ptr_type operator- (const same_type& right) const {
      // TODO: implement this function
    }

    // Subtraction with a pointer to a vector.
    ptr_type operator- (ptr_type right) const {
      // TODO: implement this function

      // hint: reuse other overloaded operators, instead of writing
      // this from scratch
    }

    // Negation operator.
    ptr_type operator- () const {
      // TODO: implement this function

      // hint: reuse other overloaded operators, instead of writing
      // this from scratch
    }

    // Multiplication by a scalar.
    ptr_type operator* (SCALAR s) const {
      // TODO: implement this function
    }

    // Multiplication by a vector. Note the dimensions of the input
    // and output vectors.
    std::shared_ptr<column_type> operator* (const row_type& v) const {
      // TODO: implement this function
    }
    
    // Multiplication by a pointer to a vector. Note the dimensions of
    // the input and output vectors.
    std::shared_ptr<column_type> operator* (std::shared_ptr<row_type> v) const {
      return (*this) * (*v);
    }

    // Multiplication by a matrix. This and the operand matrix must
    // both be square and have identical dimensions.
    //
    // Note that in mathematics, matrices are not necessarily required
    // to be square in order to be multiplied. However in graphics we
    // will only ever multiply square matrices, and making this
    // assumption makes this function significantly easier to declare
    // and implement.
    ptr_type operator* (const same_type& right) const {
      assert(is_square());
      // TODO: implement this function
    }

    // Multiplication by a pointer to a matrix.
    ptr_type operator* (ptr_type right) const {
      return (*this) * (*right);
    }

    // Return true iff j is a valid column index.
    bool is_column(int j) const {
      return ((j >= 0) && (j < WIDTH));
    }

    // Return true iff i is a valid row index.
    bool is_row(int i) const {
      return ((i >= 0) && (i < HEIGHT));
    }

    // Return true iff this is a square matrix.
    bool is_square() const {
      return (WIDTH == HEIGHT);
    }

    // Return the transpose of this matrix.
    std::shared_ptr<Matrix<SCALAR, WIDTH, HEIGHT> > transpose() const {
      // TODO: implement this function
    }

    // Return the height of this matrix.
    int height() const {
      return HEIGHT;
    }

    // Utility function to print a representation of this matrix to
    // cout. This is intende for debugging purposes.
    void print(int element_width = 6) const {
      for (int i = 0; i < HEIGHT; ++i) {
	if (i == 0) {
	  std::cout << "┌";
	} else if (i < (HEIGHT - 1)) {
	  std::cout << "│";
	} else {
	  std::cout << "└";
	}

	for (int j = 0; j < WIDTH; ++j) {
	  std::cout << ' ' << std::right << std::setw(element_width) << _rows[i][j];
	}

	std::cout << ' ';

	if (i == 0) {
	  std::cout << "┐";
	} else if (i < (HEIGHT - 1)) {
	  std::cout << "│";
	} else {
	  std::cout << "┘";
	}

	std::cout << std::endl;
      }
    }

    // Return the width of this matrix.
    int width() const {
      return WIDTH;
    }

  };

  // Compute the determinant of a matrix. We only define this function
  // for 2x2 and 3x3 square matrices.
  
  template <typename SCALAR>
  SCALAR determinant(const Matrix<SCALAR, 2, 2>& m) {
    // TODO: implement this function

    // hint: reuse other overloaded operators, instead of writing
    // this from scratch
  }

  template <typename SCALAR>
  SCALAR determinant(const Matrix<SCALAR, 3, 3>& m) {
    // TODO: implement this function

    // hint: reuse other overloaded operators, instead of writing
    // this from scratch
  }

  // Compute the inverse of a matrix. We only define this function for
  // 2x2 and 3x3 square matrices.

  template <typename SCALAR>
  std::shared_ptr<Matrix<SCALAR, 2, 2> > inverse(const Matrix<SCALAR, 2, 2>& m) {
    // TODO: implement this function

    // hint: reuse other overloaded operators, instead of writing
    // this from scratch
  }

  template <typename SCALAR>
  std::shared_ptr<Matrix<SCALAR, 3, 3> > inverse(const Matrix<SCALAR, 3, 3>& m) {
    // TODO: implement this function

    // hint: reuse other overloaded operators, instead of writing
    // this from scratch
  }
  
}
