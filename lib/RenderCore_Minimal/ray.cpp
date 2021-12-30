#include "core_settings.h"

namespace lh2core {
	bool Ray::IntersectTriangle(const float3& v0, const float3& v1,
		const float3& v2, float& tHit, float& u, float& v) const {
		float3 e1, e2, h, s, q;
		float a, f;
		e1 = v1 - v0;
		e2 = v2 - v0;
		h = cross(d, e2);
		a = dot(e1, h);
		if (a > -0.00001 && a < 0.00001) {
			return false;
		}
		f = 1 / a;
		s = o - v0;
		u = f * (dot(s, h));
		if (u < 0.0 || u > 1.0) {
			return false;
		}
		q = cross(s, e1);
		v = f * dot(d, q);
		if (v < 0.0 || u + v > 1.0) {
			return false;
		}
		// at this stage we can compute t to find out where
		// the intersection point is on the line
		tHit = f * dot(e2, q);
		if (tHit > 0.00001 && tHit < tMax) {// ray intersection
			return true;
		} else {// this means that there is a line intersection but not a ray intersection
			return false;
		}
	}
}