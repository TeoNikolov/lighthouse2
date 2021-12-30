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
		float3 diffuse = WHITE;		// diffuse material color
		float3 absorption = BLACK;		// diffuse material color
		float3 specular = WHITE;
		float specularity = 0.0;
		float ior = 1.0;
		float gloss = 1.0f; // alpha factor
		bool emissive = false;
		bool textured = false;
		Texture* texture = 0;			// texture
	};
}