#pragma once

namespace lh2core
{
	struct Ray
	{
		Ray() { Reset(); }
		Ray( float3 origin, float3 direction, float tMax, int depth ) :
			o( origin ), d( direction ), tMax(tMax) {}
		float3 o; // 12 bytes
		mutable float tMax; // 4 bytes
		int dummy;
		float3 d; // 12 bytes

		void Reset() {
			o = float3({ 0.0, 0.0, 0.0 });
			d = float3({ 0.0, 0.0, 1.0 });
			tMax = INFINITY;
		}
		bool IntersectTriangle(const float3& v0, const float3& v1,
			const float3& v2, float& tHit, float& u, float& v) const;
		// get point along the ray 
		float3 operator()( float t ) const { return o + d * t; }
	};
}