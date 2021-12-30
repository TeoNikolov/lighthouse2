#pragma once

#define RECURSIONLIMIT 4

namespace lh2core
{
	struct Ray
	{
		Ray() { Reset(); }
		Ray( float3 origin, float3 direction, float tMax, int depth ) :
			o( origin ), d( direction ), tMax(tMax) {}
		float3 o; // 12 bytes
		float3 d; // 12 bytes
		mutable float tMax; // 4 bytes
		mutable bool inDielectric; // 4 bytes IOR of ray origin
		
		void Reset() {
			o = float3({ 0.0, 0.0, 0.0 });
			d = float3({ 0.0, 0.0, 1.0 });
			inDielectric = false;
			tMax = INFINITY;
		}
		// get point along the ray 
		float3 operator()( float t ) const { return o + d * t; }
	};
}