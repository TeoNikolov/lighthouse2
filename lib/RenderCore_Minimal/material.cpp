#include "core_settings.h"

namespace lh2core {
	float3 Material::GetTextureColor(float2 uv) {
		int x = (int)texture->width * uv.x;
		int y = (int)texture->height * uv.y;
		return texture->pixels[y * texture->width + x];
	}
}