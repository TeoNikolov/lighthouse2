#include "core_settings.h"

namespace lh2core {
	bool Mesh::IntersectTriangle(const Ray& ray, const float3 &v0, const float3 &v1,
		const float3 &v2, float& tHit, float& u, float& v) {
		float3 e1, e2, h, s, q;
		float a, f;
		e1 = v1 - v0;
		e2 = v2 - v0;
		h = cross(ray.d, e2);
		a = dot(e1, h);
		if (a > -0.00001 && a < 0.00001) {
			return false;
		}
		f = 1 / a;
		s = ray.o - v0;
		u = f * (dot(s, h));
		if (u < 0.0 || u > 1.0) {
			return false;
		}
		q = cross(s, e1);
		v = f * dot(ray.d, q);
		if (v < 0.0 || u + v > 1.0) {
			return false;
		}
		// at this stage we can compute t to find out where
		// the intersection point is on the line
		tHit = f * dot(e2, q);
		if (tHit > 0.00001) {// ray intersection
			return true;
		}
		else {// this means that there is a line intersection but not a ray intersection
			return false;
		}
	}

	bool Mesh::Intersect(const Ray& ray, HitInfo* hitInfo) const {
		float t = INFINITY;
		float2 uv = make_float2(-1, -1);
		CoreTri *intersectedTriangle;
		bool intersection = false;

		for (int i = 0; i < vcount / 3; i++) {
			CoreTri *triangle = &triangles[i];
			float ttemp = 0;
			float u = -1;
			float v = -1;
			if (IntersectTriangle(ray, triangle->vertex0, triangle->vertex1, triangle->vertex2, ttemp, u, v)) {
				intersection = true;
				if (ttemp < t) {
					intersectedTriangle = triangle;
					t = ttemp;
					uv.x = u;
					uv.y = v;
				}
			}
		}

		if (intersection) {
			hitInfo->tHit = t;
			hitInfo->intersects = true;
			hitInfo->hitNormal = float3({ intersectedTriangle->Nx, intersectedTriangle->Ny, intersectedTriangle->Nz });
			hitInfo->intersection = ray.o + ray.d * t;
			hitInfo->triangle = intersectedTriangle;
			hitInfo->uv = uv;
		}

		return intersection;
	}
}