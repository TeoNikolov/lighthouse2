#include "core_settings.h"

namespace lh2core {
	Frustum::Frustum(Ray& r0, Ray& r1, Ray& r2, Ray& r3) {
		cornerRays[0] = r0;
		cornerRays[1] = r1;
		cornerRays[2] = r2;
		cornerRays[3] = r3;
		normals[0] = -normalize(cross(r0.d, r1.d));
		normals[1] = -normalize(cross(r1.d, r2.d));
		normals[2] = -normalize(cross(r2.d, r3.d));
		normals[3] = -normalize(cross(r3.d, r0.d));
		b[0] = dot(normals[0], cornerRays[0].o);
		b[1] = dot(normals[1], cornerRays[1].o);
		b[2] = dot(normals[2], cornerRays[2].o);
		b[3] = dot(normals[3], cornerRays[3].o);
	}

	bool Frustum::Intersect(AABB& aabb) {
		int out, in;
		for (int i = 0; i < 4; i++) {
			out = 0;
			in = 0;
			for (int j = 0; j < 8 && (in == 0 || out == 0); j++) {
				if ((dot(normals[i], aabb.GetCornerCoords(j)) - b[i]) > 0) {
					out++;
				} else {
					in++;
				}
			}
			if (!in) {
				return false;
			} // else if (out) the box intersects
		}
		return true;
	}

	bool Frustum::Intersect(CoreTri& tri) {
		for (int i = 0; i < 4; i++) {
			bool v0 = (dot(normals[i], tri.vertex0) - b[i]) > 0;
			bool v1 = (dot(normals[i], tri.vertex1) - b[i]) > 0;
			bool v2 = (dot(normals[i], tri.vertex2) - b[i]) > 0;

			if (v0 && v1 && v2) {
				return false;
			}
		}
		return true;
	}

	int RayPacket::GetLastHit(AABB& aabb, int index) {
		return 1;
	}

	int RayPacket::PartRays(Frustum& frustum, AABB& aabb, int* indices, int lastIndex) {
		if (!frustum.Intersect(aabb)) return 0;
		int index = 0;
		for (int i = 0; i < lastIndex; i++) {
			float t0, t1;
			if (aabb.Intersect(rays[indices[i]], t0, t1)) std::swap(indices[index++], indices[i]);
		}
		return index;
	}
}