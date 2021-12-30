#include "core_settings.h"

namespace lh2core {
	void AABB::Expand(float3 point) {
		pmin.x = min(pmin.x, point.x);
		pmin.y = min(pmin.y, point.y);
		pmin.z = min(pmin.z, point.z);
		pmax.x = max(pmax.x, point.x);
		pmax.y = max(pmax.y, point.y);
		pmax.z = max(pmax.z, point.z);
	}

	bool AABB::IsContained(AABB& aabb) {
		bool mask = true;
		mask &= (!(aabb.pmin.x < pmin.x) & 0x1) << 0;
		mask &= (!(aabb.pmin.y < pmin.y) & 0x1) << 1;
		mask &= (!(aabb.pmin.z < pmin.z) & 0x1) << 2;
		mask &= (!(aabb.pmax.x > pmax.x) & 0x1) << 3;
		mask &= (!(aabb.pmax.y > pmax.y) & 0x1) << 4;
		mask &= (!(aabb.pmax.z > pmax.z) & 0x1) << 5;
		return !mask;
	}

	AABB AABB::Merge(AABB& aabb) {
		AABB r = AABB();
		r.pmin.x = min(pmin.x, aabb.pmin.x);
		r.pmin.y = min(pmin.y, aabb.pmin.y);
		r.pmin.z = min(pmin.z, aabb.pmin.z);
		r.pmax.x = max(pmax.x, aabb.pmax.x);
		r.pmax.y = max(pmax.y, aabb.pmax.y);
		r.pmax.z = max(pmax.z, aabb.pmax.z);
		return r;
	}

	// 0-3 = near corners CCW, 4-7 = far corners CCW
	float3 AABB::GetCornerCoords(uint cornerId) {
		switch (cornerId) {
		default:
			throw runtime_error("Invalid cornerID: Must be between 0-7 inclusive.");
		case 0:
			return make_float3(pmin.x, pmin.y, pmin.z);
		case 1:
			return make_float3(pmax.x, pmin.y, pmin.z);
		case 2:
			return make_float3(pmax.x, pmax.y, pmin.z);
		case 3:
			return make_float3(pmin.x, pmax.y, pmin.z);
		case 4:
			return make_float3(pmin.x, pmin.y, pmax.z);
		case 5:
			return make_float3(pmax.x, pmin.y, pmax.z);
		case 6:
			return make_float3(pmax.x, pmax.y, pmax.z);
		case 7:
			return make_float3(pmin.x, pmax.y, pmax.z);
		}
	}

	AABB AABB::TransformedAABB(mat4& transform) {
		float3 corners[] = {
			GetCornerCoords(0),
			GetCornerCoords(1),
			GetCornerCoords(2),
			GetCornerCoords(3),
			GetCornerCoords(4),
			GetCornerCoords(5),
			GetCornerCoords(6),
			GetCornerCoords(7) };
		AABB aabb = AABB();
		aabb.pmin = make_float3(INFINITY, INFINITY, INFINITY);
		aabb.pmax = make_float3(-INFINITY, -INFINITY, -INFINITY);
		for (int i = 0; i < 8; i++) {
			aabb.Expand(make_float3(transform * make_float4(corners[i], 1.0)));
		}
		return aabb;
	}

	float AABB::ComputeArea() {
		float3 diff = make_float3(pmax.x - pmin.x, pmax.y - pmin.y, pmax.z - pmin.z);
		return (diff.x * diff.y + diff.x * diff.z + diff.y * diff.z) * 2.0f;
	}

	bool AABB::Intersect(const Ray& ray, float& hitt0, float& hitt1) {
		float t0 = 0;
		float t1 = ray.tMax;
		float invDir, tNear, tFar;
		invDir = 1.0f / ray.d.x;
		tNear = (pmin.x - ray.o.x) * invDir;
		tFar = (pmax.x - ray.o.x) * invDir;
		if (tNear > tFar) std::swap(tNear, tFar);
		t0 = tNear > t0 ? tNear : t0;
		t1 = tFar < t1 ? tFar : t1;
		if (t0 > t1) return false;
		invDir = 1.0f / ray.d.y;
		tNear = (pmin.y - ray.o.y) * invDir;
		tFar = (pmax.y - ray.o.y) * invDir;
		if (tNear > tFar) std::swap(tNear, tFar);
		t0 = tNear > t0 ? tNear : t0;
		t1 = tFar < t1 ? tFar : t1;
		if (t0 > t1) return false;
		invDir = 1.0f / ray.d.z;
		tNear = (pmin.z - ray.o.z) * invDir;
		tFar = (pmax.z - ray.o.z) * invDir;
		if (tNear > tFar) std::swap(tNear, tFar);
		t0 = tNear > t0 ? tNear : t0;
		t1 = tFar < t1 ? tFar : t1;
		if (t0 > t1) return false;
		hitt0 = t0;
		hitt1 = t1;
		return true;
	}

	float3 AABB::Diagonal() {
		return pmax - pmin;
	}

	uint AABB::GetLongestAxis() {
		float3 d = Diagonal();
		if (d.x > d.y && d.x > d.z) {
			return 0; // x
		}
		else if (d.y > d.z) {
			return 1; // y
		}
		else {
			return 2; // z
		}
	}
}