#include "rdr/accel.h"

#include "rdr/canary.h"
#include "rdr/interaction.h"
#include "rdr/math_aliases.h"
#include "rdr/platform.h"
#include "rdr/shape.h"

RDR_NAMESPACE_BEGIN

/* ===================================================================== *
 *
 * AABB Implementations
 *
 * ===================================================================== */

bool AABB::isOverlap(const AABB &other) const {
  return ((other.low_bnd[0] >= this->low_bnd[0] &&
              other.low_bnd[0] <= this->upper_bnd[0]) ||
             (this->low_bnd[0] >= other.low_bnd[0] &&
                 this->low_bnd[0] <= other.upper_bnd[0])) &&
         ((other.low_bnd[1] >= this->low_bnd[1] &&
              other.low_bnd[1] <= this->upper_bnd[1]) ||
             (this->low_bnd[1] >= other.low_bnd[1] &&
                 this->low_bnd[1] <= other.upper_bnd[1])) &&
         ((other.low_bnd[2] >= this->low_bnd[2] &&
              other.low_bnd[2] <= this->upper_bnd[2]) ||
             (this->low_bnd[2] >= other.low_bnd[2] &&
                 this->low_bnd[2] <= other.upper_bnd[2]));
}

bool AABB::intersect(const Ray &ray, Float *t_in, Float *t_out) const {
  // TODO(HW3): implement ray intersection with AABB.
  // ray distance for two intersection points are returned by pointers.
  //
  // This method should modify t_in and t_out as the "time"
  // when the ray enters and exits the AABB respectively.
  //
  // And return true if there is an intersection, false otherwise.
  //
  // Useful Functions:
  // @see Ray::safe_inverse_direction
  //    for getting the inverse direction of the ray.
  // @see Min/Max/ReduceMin/ReduceMax
  //    for vector min/max operations.
  
  // Use the slab method for ray-AABB intersection
  // For each axis, compute the intersection times with the two planes
  Vec3f t0 = (low_bnd - ray.origin) * ray.safe_inverse_direction;
  Vec3f t1 = (upper_bnd - ray.origin) * ray.safe_inverse_direction;
  
  // Make sure t0 has the near intersection and t1 has the far intersection
  Vec3f t_near = Min(t0, t1);
  Vec3f t_far = Max(t0, t1);
  
  // Find the largest t_near and smallest t_far across all axes
  Float t_enter = ReduceMax(t_near);
  Float t_exit = ReduceMin(t_far);
  
  // Check if there's a valid intersection
  // The ray intersects the AABB if t_enter <= t_exit
  if (t_enter > t_exit) {
    return false;
  }
  
  // Set the output parameters
  *t_in = t_enter;
  *t_out = t_exit;
  
  return true;
}

/* ===================================================================== *
 *
 * Accelerator Implementations
 *
 * ===================================================================== */

bool TriangleIntersect(Ray &ray, const uint32_t &triangle_index,
    const ref<TriangleMeshResource> &mesh, SurfaceInteraction &interaction) {
  using InternalScalarType = Double;
  using InternalVecType    = Vec<InternalScalarType, 3>;

  AssertAllValid(ray.direction, ray.origin);
  AssertAllNormalized(ray.direction);

  const auto &vertices = mesh->vertices;
  const Vec3u v_idx(&mesh->v_indices[3 * triangle_index]);
  assert(v_idx.x < mesh->vertices.size());
  assert(v_idx.y < mesh->vertices.size());
  assert(v_idx.z < mesh->vertices.size());

  InternalVecType dir = Cast<InternalScalarType>(ray.direction);
  InternalVecType v0  = Cast<InternalScalarType>(vertices[v_idx[0]]);
  InternalVecType v1  = Cast<InternalScalarType>(vertices[v_idx[1]]);
  InternalVecType v2  = Cast<InternalScalarType>(vertices[v_idx[2]]);

  // TODO(HW3): implement ray-triangle intersection test.
  // You should compute the u, v, t as InternalScalarType
  //
  //   InternalScalarType u = ...;
  //   InternalScalarType v = ...;
  //   InternalScalarType t = ...;
  //
  // And exit early with `return false` if there is no intersection.
  //
  // The intersection points is denoted as:
  // (1 - u - v) * v0 + u * v1 + v * v2 == ray.origin + t * ray.direction
  // where the left side is the barycentric interpolation of the triangle
  // vertices, and the right side is the parametric equation of the ray.
  //
  // You should also make sure that:
  // u >= 0, v >= 0, u + v <= 1, and, ray.t_min <= t <= ray.t_max
  //
  // Useful Functions:
  // You can use @see Cross and @see Dot for determinant calculations.

  // Möller–Trumbore intersection algorithm
  // Solve: origin + t*dir = (1-u-v)*v0 + u*v1 + v*v2
  // Rearranged: origin + t*dir = v0 + u*(v1-v0) + v*(v2-v0)
  // So: t*dir = (v0 - origin) + u*(v1-v0) + v*(v2-v0)
  // Or: -t*dir + u*e1 + v*e2 = origin - v0
  // Matrix form: [-dir, e1, e2] * [t, u, v]^T = origin - v0
  
  InternalVecType origin = Cast<InternalScalarType>(ray.origin);
  InternalVecType e1 = v1 - v0;
  InternalVecType e2 = v2 - v0;
  InternalVecType s = origin - v0;
  
  // Calculate pvec = dir × e2
  InternalVecType pvec = Cross(dir, e2);
  
  // Calculate determinant = e1 · (dir × e2) = e1 · pvec
  InternalScalarType det = Dot(e1, pvec);
  
  // If determinant is near zero, ray is parallel to triangle
  if (abs(det) < InternalScalarType(1e-10)) {
    return false;
  }
  
  InternalScalarType inv_det = InternalScalarType(1) / det;
  
  // Calculate u = (s · pvec) / det
  InternalScalarType u = Dot(s, pvec) * inv_det;
  
  // Check if u is in valid range [0, 1]
  if (u < InternalScalarType(0) || u > InternalScalarType(1)) {
    return false;
  }
  
  // Calculate qvec = s × e1
  InternalVecType qvec = Cross(s, e1);
  
  // Calculate v = (dir · qvec) / det
  InternalScalarType v = Dot(dir, qvec) * inv_det;
  
  // Check if v is in valid range and u+v <= 1
  if (v < InternalScalarType(0) || u + v > InternalScalarType(1)) {
    return false;
  }
  
  // Calculate t = (e2 · qvec) / det
  InternalScalarType t = Dot(e2, qvec) * inv_det;
  
  // Check if intersection is within ray's valid range
  if (t < static_cast<InternalScalarType>(ray.t_min) || 
      t > static_cast<InternalScalarType>(ray.t_max)) {
    return false;
  }

  // We will reach here if there is an intersection

  CalculateTriangleDifferentials(interaction,
      {static_cast<Float>(1 - u - v), static_cast<Float>(u),
          static_cast<Float>(v)},
      mesh, triangle_index);
  AssertNear(interaction.p, ray(t));
  assert(ray.withinTimeRange(t));
  ray.setTimeMax(t);
  return true;
}

void Accel::setTriangleMesh(const ref<TriangleMeshResource> &mesh) {
  // Build the bounding box
  AABB bound(Vec3f(Float_INF, Float_INF, Float_INF),
      Vec3f(Float_MINUS_INF, Float_MINUS_INF, Float_MINUS_INF));
  for (auto &vertex : mesh->vertices) {
    bound.low_bnd   = Min(bound.low_bnd, vertex);
    bound.upper_bnd = Max(bound.upper_bnd, vertex);
  }

  this->mesh  = mesh;   // set the pointer
  this->bound = bound;  // set the bounding box
}

void Accel::build() {}

AABB Accel::getBound() const {
  return bound;
}

bool Accel::intersect(Ray &ray, SurfaceInteraction &interaction) const {
  bool success = false;
  for (int i = 0; i < mesh->v_indices.size() / 3; i++)
    success |= TriangleIntersect(ray, i, mesh, interaction);
  return success;
}

RDR_NAMESPACE_END
