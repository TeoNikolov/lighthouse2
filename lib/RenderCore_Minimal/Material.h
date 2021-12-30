#pragma once

namespace lh2core {
	enum MaterialType {
		DIFFUSE,
		MIRROR,
		DIELECTRIC
	};

	class Material {
		public:
		// constructor / destructor
		Material() = default;
		float3 GetTextureColor(float2 uv);
		// data members
		MaterialType type = DIFFUSE;
		float3 diffuse = float3({ 1.0, 1.0, 1.0 });		// diffuse material color
		float3 absorption = float3({ 0.0, 0.0, 0.0 });		// diffuse material color
		float specularity = 0.0;
		float ior = 1.0;
		bool emissive = false;
		bool textured = false;
		Texture* texture = 0;			// texture
	};
}