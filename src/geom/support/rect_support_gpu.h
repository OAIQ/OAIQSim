// -----------------------------------------------------------------------
// ************************ OAIQSIM TOOLBOX ******************************
// file: rect_support.h
// purpose: Contains class definitions and implementations for a 3D "box" 
//          support set class.
// Author:  Nick Henscheid
// Date:    3-2019
// Contact: nhenscheid@math.arizona.edu
// References: 
// This software is in the public domain, furnished "as is", without 
// technical support, and with no warranty, express or implied, as to its 
// usefulness for any purpose.
// ------------------------------------------------------------------------
#ifndef __RAY_TRACING_UTIL_H__
#define __RAY_TRACING_UTIL_H__

#include <assert.h>

class Vector3 {
  public:
    __device__ Vector3() { };
    __device__ Vector3(float x, float y, float z) { d[0] = x; d[1] = y; d[2] = z; }
    __device__ Vector3(const Vector3 &v)
      { d[0] = v.d[0]; d[1] = v.d[1]; d[2] = v.d[2]; }

    __device__ float x() const { return d[0]; }
    __device__ float y() const { return d[1]; }
    __device__ float z() const { return d[2]; }

    __device__ float operator[](int i) const { return d[i]; }
    
    __device__ float length() const
      { return sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]); }
    __device__ void normalize() {
      float temp = length();
      if (temp == 0.0)
        return;	// 0 length vector
      // multiply by 1/magnitude
      temp = 1 / temp;
      d[0] *= temp;
      d[1] *= temp;
      d[2] *= temp;
    }

    /////////////////////////////////////////////////////////
    // Overloaded operators
    /////////////////////////////////////////////////////////
  
    __device__ Vector3 operator+(const Vector3 &op2) const {   // vector addition
      return Vector3(d[0] + op2.d[0], d[1] + op2.d[1], d[2] + op2.d[2]);
    }
    __device__ Vector3 operator-(const Vector3 &op2) const {   // vector subtraction
      return Vector3(d[0] - op2.d[0], d[1] - op2.d[1], d[2] - op2.d[2]);
    }
    __device__ Vector3 operator-() const {                    // unary minus
      return Vector3(-d[0], -d[1], -d[2]);
    }
    __device__ Vector3 operator*(float s) const {            // scalar multiplication
      return Vector3(d[0] * s, d[1] * s, d[2] * s);
    }
    __device__ void operator*=(float s) {
      d[0] *= s;
      d[1] *= s;
      d[2] *= s;
    }
    __device__ Vector3 operator/(float s) const {            // scalar division
      return Vector3(d[0] / s, d[1] / s, d[2] / s);
    }
    __device__ float operator*(const Vector3 &op2) const {   // dot product
      return d[0] * op2.d[0] + d[1] * op2.d[1] + d[2] * op2.d[2];
    }
    __device__ Vector3 operator^(const Vector3 &op2) const {   // cross product
      return Vector3(d[1] * op2.d[2] - d[2] * op2.d[1], d[2] * op2.d[0] - d[0] * op2.d[2],
                    d[0] * op2.d[1] - d[1] * op2.d[0]);
    }
    __device__ bool operator==(const Vector3 &op2) const {
      return (d[0] == op2.d[0] && d[1] == op2.d[1] && d[2] == op2.d[2]);
    }
    __device__ bool operator!=(const Vector3 &op2) const {
      return (d[0] != op2.d[0] || d[1] != op2.d[1] || d[2] != op2.d[2]);
    }
    __device__ bool operator<(const Vector3 &op2) const {
      return (d[0] < op2.d[0] && d[1] < op2.d[1] && d[2] < op2.d[2]);
    }
    __device__ bool operator<=(const Vector3 &op2) const {
      return (d[0] <= op2.d[0] && d[1] <= op2.d[1] && d[2] <= op2.d[2]);
    }
  
  private:
    float d[3];
};


/*
 * Ray class, for use with the optimized ray-box intersection test
 * described in:
 *
 *      Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley
 *      "An Efficient and Robust Ray-Box Intersection Algorithm"
 *      Journal of graphics tools, 10(1):49-54, 2005
 * 
 */

class Ray {
  public:
    __device__ Ray() { }
    __device__ Ray(Vector3 o, Vector3 d) {
      origin = o;
      direction = d;
      inv_direction = Vector3(1/d.x(), 1/d.y(), 1/d.z());
      sign[0] = (inv_direction.x() < 0);
      sign[1] = (inv_direction.y() < 0);
      sign[2] = (inv_direction.z() < 0);
    }
    __device__ Ray(const Ray &r) {
      origin = r.origin;
      direction = r.direction;
      inv_direction = r.inv_direction;
      sign[0] = r.sign[0]; sign[1] = r.sign[1]; sign[2] = r.sign[2];
    }

    Vector3 origin;
    Vector3 direction;
    Vector3 inv_direction;
    int sign[3];
};

/*
 * Axis-aligned bounding box class, for use with the optimized ray-box
 * intersection test described in:
 *
 *      Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley
 *      "An Efficient and Robust Ray-Box Intersection Algorithm"
 *      Journal of graphics tools, 10(1):49-54, 2005
 *
 */
class Box {
  public:
    __device__ Box() { }
    __device__ Box(const Vector3 &min, const Vector3 &max) {
      assert(min < max);
      parameters[0] = min;
      parameters[1] = max;
    }
    // (t0, t1) is the interval for valid hits
    //bool intersect(const Ray &, float t0, float t1) const;
    __device__ float intersect(const Ray &r, float t0, float t1) const {
      float tmin, tmax, tymin, tymax, tzmin, tzmax;

      tmin = (parameters[r.sign[0]].x() - r.origin.x()) * r.inv_direction.x();
      tmax = (parameters[1-r.sign[0]].x() - r.origin.x()) * r.inv_direction.x();
      tymin = (parameters[r.sign[1]].y() - r.origin.y()) * r.inv_direction.y();
      tymax = (parameters[1-r.sign[1]].y() - r.origin.y()) * r.inv_direction.y();
      if ( (tmin > tymax) || (tymin > tmax) ) 
        //return false;
        return 0.0;
      if (tymin > tmin)
        tmin = tymin;
      if (tymax < tmax)
        tmax = tymax;
      tzmin = (parameters[r.sign[2]].z() - r.origin.z()) * r.inv_direction.z();
      tzmax = (parameters[1-r.sign[2]].z() - r.origin.z()) * r.inv_direction.z();
      if ( (tmin > tzmax) || (tzmin > tmax) ) 
        //return false;
        return 0.0;
      if (tzmin > tmin)
        tmin = tzmin;
      if (tzmax < tmax)
        tmax = tzmax;
      //return ( (tmin < t1) && (tmax > t0) );
        return ((tmin < t1)&&(tmax > t0))*(tmax-tmin);
    }

    // corners
    Vector3 parameters[2];
};


#endif