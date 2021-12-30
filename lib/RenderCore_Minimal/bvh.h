#pragma once

namespace lh2core {
	enum BVHBuildType {
		Midpoint,
		Binned1,
		Binned3
	};
	
	enum BVHUpdateType {
		Refit,
		Rebuild
	};

	class Mesh;
	class MeshInstance;

	struct BVHNode {
	public:
		void CalculateBounds(CoreTri* primitives, uint* indices, int firstPrim, int primCount);
		bool RefitBounds(CoreTri* primitives, uint* indices, BVHNode** pool); // true if parent needs refitting
		void Subdivide(CoreTri* primitives, float3* centroids, uint* indices,
			BVHNode** pool, uint& poolPtr, int firstPrim, int primCount, uint depth);
		uint Partition(CoreTri* primitives, float3* centroids, uint* indices, int firstPrim, int primCount, uint depth);
		bool Intersect(CoreTri* primitives, uint* indices, BVHNode** pool, const Ray& r, HitInfo* hitInfo, uint& traverseCounter);
		bool IntersectPrimitives(CoreTri* primitives, uint* indices, const Ray& ray, HitInfo* hitInfo);

		// data members
		AABB bounds;
		int leftFirst, count;
	};

	class TopBVHNode {
	public:
		bool Intersect(Ray& r, HitInfo* hitInfo, uint& traverseCounter);
		// given an array of AABBs, find a node for which a parent AABB has the smallest SA
		// 'index' is the index of THIS node in the nodes array (-1 if it is not contained in the array)
		// data members
		AABB bounds;
		bool isLeaf;
		TopBVHNode* left, *right;
		MeshInstance* instance;
	};

	class BVH {
	public:
		void Update();
		void Build();
		void Refit();
		bool Intersect(const Ray& ray, HitInfo* hitInfo, uint &traverseCounter);
		bool Traverse(RayPacket& ray, Frustum& frustum, HitInfo* hitInfo, uint &traverseCounter);
		void* operator new(size_t i)
		{
			return _mm_malloc(i, 128);
		}
		void operator delete(void* p)
		{
			_mm_free(p);
		}
		Mesh* mesh;
		BVHNode* root;
		BVHUpdateType updateType = BVHUpdateType::Refit;
	private:
		uint* indices;
		alignas(128) BVHNode** pool;
		uint poolPtr;
		bool initialized = false;
	};

	class TopBVH {
	public:
		void Build(vector<MeshInstance*>& instances);
		bool Intersect(Ray& ray, HitInfo* hitInfo, uint& traverseCounter);
		int FindBestMatch(vector<TopBVHNode*> nodeList, TopBVHNode* node, int index, AABB& matchAABB);
	private:
		TopBVHNode* root;
	};

	class StackNode {
	public:
		BVHNode* node;
		int index;
	};
}