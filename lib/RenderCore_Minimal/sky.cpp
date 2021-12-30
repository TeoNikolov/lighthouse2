#include "core_settings.h"

namespace lh2core {
	float3 Sky::GetColor(Ray& ray) {
		double u = 1.0 + atan2(ray.d.x, -ray.d.z) * INVPI;
		double v = acos(ray.d.y) * INVPI;
		uint index = (uint)(u * width) + (uint)(v * height) * width;
		return pixels[index];
	}
}