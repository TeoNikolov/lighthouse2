#pragma once

namespace lh2core
{
	const float3 zeroFloat = float3({ 0.0, 0.0, 0.0 });
	struct HitInfo {
		HitInfo() { Reset(); }
		void Reset() {
			intersects = false;
			intersection = zeroFloat;
			hitNormal = zeroFloat;
			tHit = INFINITY;
			uv = make_float2(-1, -1);
		}
		bool intersects;
		float3 intersection;
		CoreTri* triangle;
		float3 hitNormal;
		float2 uv;
		float tHit; // distance from ray to intersection
	};
}