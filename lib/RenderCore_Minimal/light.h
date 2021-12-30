#pragma once

namespace lh2core {

	enum class LightType {
		PointLight,
		SpotLight,
		DirectionalLight,
		AreaLight
	};

	class Light {
	public:
		Light(LightType type) : type(type) {};
		float3 position = float3({ 0.0f, 0.0f, 0.0f });
		float3 radiance = float3({ 100.0f, 100.0f, 100.0f });
		LightType type;
	};

	class PointLight : public Light {
	public:
		PointLight() : Light(LightType::PointLight) {};
	};

	class SpotLight : public Light {
	public:
		SpotLight() : Light(LightType::SpotLight) {};
		float3 direction = float3({0.0f, 0.0f, 1.0f});
		float cosInner = 1.0f;
		float cosOuter = 1.0f;
	};

	class DirectionalLight : public Light {
	public:
		DirectionalLight() : Light(LightType::DirectionalLight) {};
	};

	class AreaLight : public Light {
	public:
		AreaLight() : Light(LightType::AreaLight) {};
		float area = 0.0f;
		float3 v1v0;
		float3 v2v0;
		float3 v0, v1, v2;
		vector<CoreLightTri> tris;
	};
}