#ifndef TAG_INTERSECTION_H_
#define TAG_INTERSECTION_H_

#include <algorithm>

#include <TooN/helpers.h>

namespace tag {

/// @defgroup intersection 3D intersections
/// This group contains a set of functions to compute intersections
/// between various 3D structures

/// computes the point of intersection between a plane (given as normal and constant)
/// and a line (given as two points on the line). It does not do any
/// normalization of plane normal or distance between points! therefore the
/// test for plane and line being parallel might depend on the magnitude
/// of these vectors.
/// @ingroup intersection
inline bool intersect_plane_line( const TooN::Vector<3> & normal, const double d, const TooN::Vector<3> & p1, const TooN::Vector<3> & p2, TooN::Vector<3> & i){
    const double EPSILON = 0.000001;

    TooN::Vector<3> dir = p2 - p1;
    double c = (normal * dir);
    if( fabs(c) < EPSILON )
        return false;
    double t = (d - normal * p1) / c;
    i = p1 + t * dir;
    return true;
}

/// computes the intersection between a line and a triangle. 
/// @ingroup intersection
inline bool intersect_triangle(const TooN::Vector<3> & orig, const TooN::Vector<3> & dir,
                        const TooN::Vector<3> & vert0, const TooN::Vector<3> & vert1, const TooN::Vector<3> & vert2,
                        double & t, double & u, double & v)
{
   const double EPSILON = 0.000001;

   double det,inv_det;

   // find vectors for two edges sharing vert0
   TooN::Vector<3> edge1 = vert1 - vert0;
   TooN::Vector<3> edge2 = vert2 - vert0;

   // begin calculating determinant - also used to calculate U parameter
   TooN::Vector<3> pvec = dir ^ edge2;

   // if determinant is near zero, ray lies in plane of triangle
   det = edge1 * pvec;

   if (det > -EPSILON && det < EPSILON)
     return false;
   inv_det = 1.0 / det;

   // calculate distance from vert0 to ray origin
   TooN::Vector<3> tvec = orig - vert0;

   // calculate U parameter and test bounds
   u = (tvec * pvec) * inv_det;
   if (u < 0.0 || u > 1.0)
     return false;

   // prepare to test V parameter
   TooN::Vector<3> qvec = tvec ^ edge1;

   // calculate V parameter and test bounds
   v = (dir * qvec) * inv_det;
   if (v < 0.0 || u + v > 1.0)
     return false;

   // calculate t, ray intersects triangle
   t = (edge2 * qvec) * inv_det;
   return true;
}

/// computes the intersection between a line and a backfaced culled triangle.
/// @ingroup intersection
inline bool intersect_culled_triangle(const TooN::Vector<3> & orig, const TooN::Vector<3> & dir,
                               const TooN::Vector<3> & vert0, const TooN::Vector<3> & vert1, const TooN::Vector<3> & vert2,
                               double & t, double & u, double & v)
{
   const double EPSILON = 0.000001;

   double det,inv_det;

   // find vectors for two edges sharing vert0
   TooN::Vector<3> edge1 = vert1 - vert0;
   TooN::Vector<3> edge2 = vert2 - vert0;

   // begin calculating determinant - also used to calculate U parameter
   TooN::Vector<3> pvec = dir ^ edge2;

   // if determinant is near zero, ray lies in plane of triangle
   det = edge1 * pvec;

   if (det < EPSILON)
      return false;

   // calculate distance from vert0 to ray origin
   TooN::Vector<3> tvec = orig - vert0;

   // calculate U parameter and test bounds
   u = tvec * pvec;
   if (u < 0.0 || u > det)
      return false;

   // prepare to test V parameter
   TooN::Vector<3> qvec = tvec ^ edge1;

    // calculate V parameter and test bounds
   v = dir * qvec;
   if (v < 0.0 || u + v > det)
      return false;

   // calculate t, scale parameters, ray intersects triangle
   t = edge2 * qvec;
   inv_det = 1.0 / det;
   t *= inv_det;
   u *= inv_det;
   v *= inv_det;
   return true;
}

/// computes the intersection between two triangles.
/// @ingroup intersection
inline bool intersect_triangles( const TooN::Vector<3> & v1, const TooN::Vector<3> & v2, const TooN::Vector<3> & v3,
                                 const TooN::Vector<3> & w1, const TooN::Vector<3> & w2, const TooN::Vector<3> & w3,
                                 TooN::Vector<3> & p1, TooN::Vector<3> & p2 ){
    const double EPSILON = 0.000001;

    const TooN::Vector<3> * tv[3]; // = {&v1, &v2, &v3};
    tv[0] = &v1; tv[1] = &v2; tv[2] = &v3;
    const TooN::Vector<3> * tw[3]; //= {&w1, &w2, &w3};
    tw[0] = &w1; tw[1] = &w2; tw[2] = &w3;

    // normal vector of plane of triangle v
    TooN::Vector<3> nv = (v2 - v1) ^ ( v3 - v1 );
    double t1w = nv * w1;
    double t2w = nv * w2;
    double t3w = nv * w3;
    // all of triangle w on one side of plane of v ?
    if( (t1w < -EPSILON && t2w < -EPSILON && t3w < -EPSILON) ||
        (t1w >  EPSILON && t2w >  EPSILON && t3w >  EPSILON) ) {
        return false;
    }

    // normal vector of plane of triangle w
    TooN::Vector<3> nw = (w2 - w1) ^ ( w3 - w1 );
    double t1v = nw * v1;
    double t2v = nw * v2;
    double t3v = nw * v3;
    // all of triangle v on one side of plane of w ?
    if( (t1v < -EPSILON && t2v < -EPSILON && t3v < -EPSILON) ||
        (t1v >  EPSILON && t2v >  EPSILON && t3v >  EPSILON) ) {
        return false;
    }

    TooN::normalize(nv);
    TooN::normalize(nw);
    // direction of line of intersection of planes
    TooN::Vector<3> d = nv ^ nw;
    // are the supporting planes almost parallel ? -> not dealing with this case
    if( d*d < EPSILON ){
        return false;
    }

    // find out which one is alone
    int iv, iw;
    if( t1v * t2v > 0 )
        iv = 2;
    else if( t1v * t3v > 0 )
        iv = 1;
    else
        iv = 0;
    if( t1w * t2w > 0 )
        iw = 2;
    else if( t1w * t3w > 0 )
        iw = 1;
    else
        iw = 0;

    // compute interval points by intersecting with respective planes
    // we know that they intersect, therefore no test for failure
    double dw = w1 * nw;
    TooN::Vector<3> iv1, iv2;
    intersect_plane_line( nw, dw, *tv[iv], *tv[(iv+1)%3], iv1 );
    intersect_plane_line( nw, dw, *tv[iv], *tv[(iv+2)%3], iv2 );
    double dv = v1 * nv;
    TooN::Vector<3> iw1, iw2;
    intersect_plane_line( nv, dv, *tw[iw], *tw[(iw+1)%3], iw1 );
    intersect_plane_line( nv, dv, *tw[iw], *tw[(iw+2)%3], iw2 );

    // project onto line
    double tiv1 = d * iv1,
           tiv2 = d * iv2;
    double tiw1 = d * iw1,
           tiw2 = d * iw2;

    // calculate interval extensions
    double tvmin = std::min(tiv1, tiv2),
           tvmax = std::max(tiv1, tiv2);
    double twmin = std::min(tiw1, tiw2),
           twmax = std::max(tiw1, tiw2);
    double intmin = std::max(tvmin, twmin),
           intmax = std::min(tvmax, twmax);

    // no intersection!
    if( intmin > intmax )
        return false;

    // determine intersection points
    if( intmin == tiv1 )
        p1 = iv1;
    else if( intmin == tiv2)
        p1 = iv2;
    else if( intmin == tiw1 )
        p1 = iw1;
    else if( intmin == tiw2 )
        p1 = iw2;
    if( intmax == tiv1 )
        p2 = iv1;
    else if( intmax == tiv2)
        p2 = iv2;
    else if( intmax == tiw1 )
        p2 = iw1;
    else if( intmax == tiw2 )
        p2 = iw2;

    return true;
}

}

#endif
