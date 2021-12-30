#pragma once

namespace lh2core {
	class AABB {
	public:
		uint GetLongestAxis();
		bool Intersect(const Ray& ray, float& hitt0, float& hitt1);
		bool IsContained(AABB& aabb);
		float3 Diagonal();
		float ComputeArea();
		void Expand(float3 point);
		float3 GetCornerCoords(uint cornerId);
		AABB TransformedAABB(mat4& transform);
		AABB Merge(AABB& aabb);
		// data variables
		float3 pmin;
		float3 pmax;
	};
}