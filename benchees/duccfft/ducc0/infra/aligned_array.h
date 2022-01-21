/*
 *  This file is part of the MR utility library.
 *
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This code is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/** \file ducc0/infra/aligned_array.h
 *
 * \copyright Copyright (C) 2019-2021 Max-Planck-Society
 * \author Martin Reinecke
 */

#ifndef DUCC0_ALIGNED_ARRAY_H
#define DUCC0_ALIGNED_ARRAY_H

#include <cstdlib>
#include <new>

namespace ducc0 {

namespace detail_aligned_array {

using namespace std;

/// Bare bones array class.
/** Mostly useful for uninitialized temporary buffers.
 *  \note Since this class operates on raw memory, it should only be used with
 *        POD types, and even then only with caution! */
template<typename T, size_t alignment=alignof(T)> class array_base
  {
  private:
    T *p;
    size_t sz;

    static T *ralloc(size_t num)
      {
      if constexpr(alignment<=alignof(max_align_t))
        {
        void *res = malloc(num*sizeof(T));
        if (!res) throw bad_alloc();
        return reinterpret_cast<T *>(res);
        }
      else
        {
        if (num==0) return nullptr;
// FIXME: let's not use aligned_alloc on Apple for the moment,
// it's only supported from 10.15 on...
#if ((__cplusplus >= 201703L) && (!defined(__APPLE__)))
        // aligned_alloc requires the allocated size to be a multiple of the
        // requested alignment, so increase size if necessary
        void *res = aligned_alloc(alignment,((num*sizeof(T)+alignment-1)/alignment)*alignment);
        if (!res) throw bad_alloc();
#else // portable emulation
        void *ptr = malloc(num*sizeof(T)+alignment);
        if (!ptr) throw bad_alloc();
        void *res = reinterpret_cast<void *>((reinterpret_cast<size_t>(ptr) & ~(size_t(alignment-1))) + alignment);
        (reinterpret_cast<void**>(res))[-1] = ptr;
#endif
        return reinterpret_cast<T *>(res);
        }
      }
    static void dealloc(T *ptr)
      {
      if constexpr(alignment<=alignof(max_align_t))
        free(ptr);
      else
#if ((__cplusplus >= 201703L) && (!defined(__APPLE__)))
        free(ptr);
#else
        if (ptr) free((reinterpret_cast<void**>(ptr))[-1]);
#endif
      }

  public:
    /// Creates a zero-sized array with no associated memory.
    array_base() : p(nullptr), sz(0) {}
    /// Creates an array with \a n entries.
    /** \note Memory is not initialized! */
    array_base(size_t n) : p(ralloc(n)), sz(n) {}
    array_base(array_base &&other)
      : p(other.p), sz(other.sz)
      { other.p=nullptr; other.sz=0; }
    ~array_base() { dealloc(p); }

    /// If \a n is different from the currnt size, resizes the array to hold
    /// \a n elements.
    /** \note No data content is copied, the new array is uninitialized! */
    void resize(size_t n)
      {
      if (n==sz) return;
      dealloc(p);
      p = ralloc(n);
      sz = n;
      }

    /// Returns a writeable reference to the element at index \a idx.
    T &operator[](size_t idx) { return p[idx]; }
    /// Returns a read-only reference to the element at index \a idx.
    const T &operator[](size_t idx) const { return p[idx]; }

    /// Returns a writeable pointer to the array data.
    T *data() { return p; }
    /// Returns a read-only pointer to the array data.
    const T *data() const { return p; }

    /// Returns the size of the array.
    size_t size() const { return sz; }
  };

template<typename T> using quick_array = array_base<T>;
template<typename T> using aligned_array = array_base<T,64>;

}

using detail_aligned_array::aligned_array;
using detail_aligned_array::quick_array;

}

#endif

