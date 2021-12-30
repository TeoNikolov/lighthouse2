#pragma once

namespace lh2core
{

class Mesh
{
public:
	static bool IntersectTriangle(const Ray &ray,const float3 &v0, const float3 &v1, const float3 &v2,
		float &tHit, float &u, float &v);
	bool Intersect(const Ray &ray, HitInfo *hitInfo) const;

	float4* vertices = 0;							// vertex data received via SetGeometry
	int vcount = 0;									// vertex count
	CoreTri* triangles = 0;							// 'fat' triangle data
};

} // namespace Tmpl8