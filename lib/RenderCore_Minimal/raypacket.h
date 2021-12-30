#pragma once

namespace lh2core {
	class Frustum {
	public:
		Frustum(Ray& r0, Ray& r1, Ray& r2, Ray& r3);
		bool Intersect(AABB& aabb);
		bool Intersect(CoreTri& tri);
		Ray cornerRays[4];
		float3 normals[4];
		float b[4];
	};

	class RayPacket
	{
	public:
		int GetLastHit(AABB& aabb, int index);
		int PartRays(Frustum& frustum, AABB& aabb, int* indices, int lastIndex);
		Ray rays[PACKETSIZE];
	};
}