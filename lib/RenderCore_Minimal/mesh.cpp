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

	void Mesh::ConstructBVH() {
		if (VERBOSEBVHCONSTRUCTION) {
			printf("Mesh %d: Constructing BVH...\n", meshIdx);
		}
		bvh->mesh = this;
		bvh->Build();
	}

	bool Mesh::Intersect(RayPacket& rayPacket, Frustum& frustum, HitInfo* hitInfo) {
		return bvh->Traverse(rayPacket, frustum, hitInfo, hitInfo->traverseDepth);
	}

	bool Mesh::Intersect(const Ray& ray, HitInfo* hitInfo) const {
#if USEBVH
		return bvh->Intersect(ray, hitInfo, hitInfo->traverseDepth);
#else
		float t = INFINITY;
		float2 uv = make_float2(-1, -1);
		CoreTri* intersectedTriangle;
		bool intersection = false;

		for (int i = 0; i < triangleCount; i++) {
			CoreTri* triangle = &triangles[i];
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
#endif
	}

	void MeshInstance::UpdateTransform(const mat4& transform) {
		this->transform = transform;
		this->inverseTransform = this->transform.Inverted();
		this->directionTransform = this->transform;
		this->directionTransform.cell[3] =
		this->directionTransform.cell[7] =
		this->directionTransform.cell[11] = 0.0;
		this->inverseDirectionTransform = directionTransform.Inverted();
	}
}