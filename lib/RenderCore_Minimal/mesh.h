#pragma once

namespace lh2core
{

class Mesh
{
public:
	static bool IntersectTriangle(const Ray &ray,const float3 &v0, const float3 &v1, const float3 &v2,
		float &tHit, float &u, float &v);
	bool Intersect(RayPacket &ray, Frustum& frustum, HitInfo *hitInfo);
	bool Intersect(const Ray &ray, HitInfo *hitInfo) const;
	void ConstructBVH();
	float4* vertices = 0;							// vertex data received via SetGeometry
	//int vcount = 0;									// vertex count
	uint triangleCount = 0;
	uint primitiveIndex = -1;
	uint meshIdx = -1;
	uint lightIdx = -1; // when mesh is an area light, this idx links the mesh with a
	CoreTri* triangles = 0;							// 'fat' triangle data
	float3* centroids = 0; // centroids of the triangles
	BVH* bvh = new BVH();
};

class MeshInstance {
public:
	MeshInstance(int index, Mesh* mesh, mat4 transform) :
		index(index), mesh(mesh), transform(transform), inverseTransform(transform.Inverted()) {}
	void UpdateTransform(const mat4& transform);
	int index = -1;
	Mesh* mesh;
	mat4 directionTransform;
	mat4 inverseDirectionTransform;
	mat4 transform;
	mat4 inverseTransform;
};

} // namespace Tmpl8