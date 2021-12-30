#pragma once

namespace lh2core {
	class Sky {
	public:
		float3 GetColor(Ray& ray);
		// data variables
		float3* pixels;
		uint width;
		uint height;
		mat4 worldToLight;
	};
}