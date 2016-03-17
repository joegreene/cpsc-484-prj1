//
// gmath_test.cc
//
// graphics math library
// (unit tests)
//
// CPSC 484, CSU Fullerton, Spring 2016, Prof. Kevin Wortman
// Project 1
//
// You shouldn't need to modify this file.
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

#include <cassert>
#include <cmath>

#include "gmath.hh"

using namespace gmath;
using namespace std;

int main() {

  // approximate_equal
  assert(approximate_equal(0.0, 0.0));
  assert(!approximate_equal(1.0, 2.0));
  assert(!approximate_equal(2.0, 1.0));
  assert(approximate_equal(1.0, 1.0001));
  assert(approximate_equal(1.0, 0.9999));

  // Vector
  {
    Vector<double, 3> a, b(0.0), c(7.1), d(7.1), e(c), diff;
    Vector<double, 3>::ptr_type p(new Vector<double, 3>(e));

    // constructors
    assert(a == b);
    assert(c == d);
    assert(e == c);

    // ==, !=
    assert(c == d);
    assert(e == p);
    assert(b != c);
    assert(b != p);

    // []
    assert(0.0 == a[0]);
    assert(0.0 == a[1]);
    assert(0.0 == a[2]);
    assert(7.1 == c[0]);
    assert(7.1 == c[1]);
    assert(7.1 == c[2]);
    diff[0] = 7;
    diff[1] = 8;
    diff[2] = 9;
    assert(7 == diff[0]);
    assert(8 == diff[1]);
    assert(9 == diff[2]);

    // =
    b = 9;
    assert(9 == b[0]);
    assert(9 == b[1]);
    assert(9 == b[2]);
    b = diff;
    assert(b == diff);
    b = p;
    assert(b == p);

    // +
    {
      auto sum1 = c + diff;
      assert( (7.1 + 7) == (*sum1)[0] );
      assert( (7.1 + 8) == (*sum1)[1] );
      assert( (7.1 + 9) == (*sum1)[2] );
      auto sum2 = diff + e;
      assert(*sum1 == sum2);
    }

    // - (subtraction)
    {
      auto dif1 = c - diff;
      assert( (7.1 - 7) == (*dif1)[0] );
      assert( (7.1 - 8) == (*dif1)[1] );
      assert( (7.1 - 9) == (*dif1)[2] );
      auto dif2 = diff - p;
      assert( (7 - 7.1) == (*dif2)[0] );
      assert( (8 - 7.1) == (*dif2)[1] );
      assert( (9 - 7.1) == (*dif2)[2] );
    }

    // - (negation)
    {
      auto neg = -diff;
      assert(-7 == (*neg)[0]);
      assert(-8 == (*neg)[1]);
      assert(-9 == (*neg)[2]);
    }

    // *
    {
      auto product1 = diff * 3;
      assert(7*3 == (*product1)[0]);
      assert(8*3 == (*product1)[1]);
      assert(9*3 == (*product1)[2]);
      auto dot1 = diff * c;
      assert( (7*7.1 + 8*7.1 + 9*7.1) == dot1 );
      auto dot2 = diff * p;
      assert(dot1 == dot2);
    }

    // / (division by scalar)
    {
      auto quotient = diff / 2.0;
      assert( (7.0 / 2.0) == (*quotient)[0] );
      assert( (8.0 / 2.0) == (*quotient)[1] );
      assert( (9.0 / 2.0) == (*quotient)[2] );
    }

    // is_index
    assert(!a.is_index(-1));
    assert(a.is_index(0));
    assert(a.is_index(1));
    assert(a.is_index(2));
    assert(!a.is_index(3));

    // is_homogeneous_point
    assert(!a.is_homogeneous_point());
    assert(b.is_homogeneous_point());

    // is_homogeneous_translation
    assert(a.is_homogeneous_translation());
    assert(!b.is_homogeneous_translation());

    // is_zero
    assert(a.is_zero());
    assert(!b.is_zero());
    assert(!c.is_zero());
    assert(!d.is_zero());
    assert(!e.is_zero());
    assert(!diff.is_zero());
    assert(!p->is_zero());

    // angle
    {
      // Example from https://www.ltcconline.net/greenl/courses/107/vectors/dotcros.htm
      Vector<double, 3> v, w;
      v[0] = 2;
      v[1] = 3;
      v[2] = 1;
      w[0] = 4;
      w[1] = 1;
      w[2] = 2;
      assert(approximate_equal(acos(13.0 / (sqrt(14.0) * sqrt(21.0))),
                               v.angle(w)));
    }

    // cross
    {
      // Example 1 from http://mathinsight.org/cross_product_examples
      Vector<int, 3> a, b;
      a[0] = 3;
      a[1] = -3;
      a[2] = 1;
      b[0] = 4;
      b[1] = 9;
      b[2] = 2;
      auto axb = a.cross(b);
      assert(-15 == (*axb)[0]);
      assert(-2 == (*axb)[1]);
      assert(39 == (*axb)[2]);
    }
    {
      // check higher-dimension behavior
      Vector<int, 3> a(1), b(1); // initialize all elements to 1
      a[0] = 3;
      a[1] = -3;
      a[2] = 1;
      b[0] = 4;
      b[1] = 9;
      b[2] = 2;
      auto axb = a.cross(b);
      assert(-15 == (*axb)[0]);
      assert(-2 == (*axb)[1]);
      assert(39 == (*axb)[2]);
    }

    // dimension
    assert(3 == a.dimension());

    // distance
    assert(0 == a.distance(a));
    assert(0 == diff.distance(diff));
    {
      double dist = sqrt(pow(7.1-7, 2) + pow(7.1-8, 2) + pow(7.1-9, 2));
      assert(dist == b.distance(diff));
      assert(dist == diff.distance(b));
    }

    // magnitude
    assert(0 == a.magnitude());
    assert(approximate_equal(sqrt(7*7 + 8*8 + 9*9),
                             diff.magnitude()));

    // normalized
    {
      auto norm = diff.normalized();
      auto mag = diff.magnitude();
      assert( (7 / mag) == (*norm)[0] );
      assert( (8 / mag) == (*norm)[1] );
      assert( (9 / mag) == (*norm)[2] );
    }

    // w
    assert(0 == a.w());
    assert(9 == diff.w());
  }

  // Matrix
  /*
  {
    Matrix<double, 2, 2> zero(0), // all zeroes
      ones(1), // all ones
      id, // identity matrix
      a; // distinct made-up values

    // initialize id
    id[0][0] = 1;
    id[1][1] = 1;

    // initialize a
    //
    // | 3 -2 |
    // | 1  4 |
    //
    a[0][0] = 3;
    a[0][1] = -2;
    a[1][0] = 1;
    a[1][1] = 4;

    Matrix<double, 2, 2>::ptr_type p_a(new Matrix<double, 2, 2>(a));

    // default constructor
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        assert(0 == zero[i][j]);
        assert(1 == ones[i][j]);
      }
    }

    // copy constructor
    {
      auto copy(a);
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          assert(copy[i][j] == a[i][j]);
        }
      }
    }

    // ==
    assert(zero == zero);
    assert(a == p_a);
    assert(!(zero == ones));
    assert(!(zero == p_a));

    // !=
    assert(zero != ones);
    assert(zero != p_a);
    assert(!(zero != zero));
    assert(!(a != p_a));

    // []
    // (already tested by initializations above)

    // =
    {
      // assign from matrix
      auto result = a;
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          assert(result[i][j] == a[i][j]);
        }
      }

      // assign from scalar
      result = 7;
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          assert(result[i][j] == 7);
        }
      }
    }

    // +
    {
      auto result = id + a;
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          assert( (*result)[i][j] == (id[i][j] + a[i][j]) );
        }
      }

      result = id + p_a;
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          assert( (*result)[i][j] == (id[i][j] + a[i][j]) );
        }
      }
    }

    // - (subtraction)
    {
      auto result = id - a;
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          assert( (*result)[i][j] == (id[i][j] - a[i][j]) );
        }
      }

      result = id - p_a;
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          assert( (*result)[i][j] == (id[i][j] - a[i][j]) );
        }
      }
    }

    // - (negation)
    {
      auto result = -a;
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          assert( (*result)[i][j] == (-a[i][j]) );
        }
      }
    }

    // * (scalar multiply)
    {
      auto result = a * 7;
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          assert( (*result)[i][j] == (a[i][j] * 7) );
        }
      }
    }

    // * (vector multiply)
    {
      Vector<double, 2> v;
      v[0] = 7;
      v[1] = 8;
      auto result = a * v;
      assert((*result)[0] == (3*7-2*8));
      assert((*result)[1] == (1*7+4*8));
      std::shared_ptr<Vector<double, 2> > ptr(new Vector<double, 2>(v));
      auto result2 = a * ptr;
      assert((*result) == (*result2));
    }

    // * (matrix multiply)
    {
      auto result = id * id;
      assert(*result == id);

      result = id * a;
      assert(*result == a);

      result = a * zero;
      assert(*result == zero);

      result = a * a;
      assert(7 == (*result)[0][0]);
      assert(-14 == (*result)[0][1]);
      assert(7 == (*result)[1][0]);
      assert(14 == (*result)[1][1]);

      auto ptr_result = a * p_a;
      assert( (*result) == (*ptr_result) );
    }

    // is_column
    assert(!id.is_column(-1));
    assert(id.is_column(0));
    assert(id.is_column(1));
    assert(!id.is_column(2));

    // is_row
    assert(!id.is_row(-1));
    assert(id.is_row(0));
    assert(id.is_row(1));
    assert(!id.is_row(2));

    // is_square
    // height
    // width
    assert(id.is_square());
    {
      assert(id.is_square());

      Matrix<double, 2, 3> wide;
      assert(!wide.is_square());
      assert(2 == wide.height());
      assert(3 == wide.width());

      Matrix<double, 3, 2> narrow;
      assert(!narrow.is_square());
      assert(3 == narrow.height());
      assert(2 == narrow.width());
    }

    // transpose
    {
      auto at = a.transpose();
      assert(3 == (*at)[0][0]);
      assert(1 == (*at)[0][1]);
      assert(-2 == (*at)[1][0]);
      assert(4 == (*at)[1][1]);
    }

    // 2D determinant and inverse
    // example from http://www.mathwords.com/i/inverse_of_a_matrix.htm
    {
      Matrix<double, 2, 2> A;
      A[0][0] = 4;
      A[0][1] = 3;
      A[1][0] = 3;
      A[1][1] = 2;
      assert(!approximate_equal(0.0, determinant(A)));
      assert(approximate_equal(8.0-9.0, determinant(A)));
      auto A_inv = inverse(A);
      assert(-2 == (*A_inv)[0][0]);
      assert(3 == (*A_inv)[0][1]);
      assert(3 == (*A_inv)[1][0]);
      assert(-4 == (*A_inv)[1][1]);
    }

    // 3D determinant and inverse
    {
      Matrix<double, 3, 3> m;
      m[0][0] = 9;
      m[0][1] = 3;
      m[0][2] = 5;
      m[1][0] = -6;
      m[1][1] = -9;
      m[1][2] = 7;
      m[2][0] = -1;
      m[2][1] = -8;
      m[2][2] = 1;
      assert(approximate_equal(615.0, determinant(m)));
      auto inv = inverse(m);
      assert(approximate_equal(0.0764228, (*inv)[0][0]));
      assert(approximate_equal(-0.0699187, (*inv)[0][1]));
      assert(approximate_equal(0.107317, (*inv)[0][2]));
      assert(approximate_equal(-0.00162602, (*inv)[1][0]));
      assert(approximate_equal(0.0227642, (*inv)[1][1]));
      assert(approximate_equal(-0.15122, (*inv)[1][2]));
      assert(approximate_equal(0.0634146, (*inv)[2][0]));
      assert(approximate_equal(0.112195, (*inv)[2][1]));
      assert(approximate_equal(-0.102439, (*inv)[2][2]));
    }

    // Show off that we can perform linear algebra operations in a way
    // that is convenient and resembles plain C++ arithmetic.
    auto mr = *(a*ones) + a*a;
  }
  */

  return 0;
}

// vim: et ts=2 sw=2 :
