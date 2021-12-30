#pragma once

namespace lh2core {
	class Texture
	{
	public:
		// constructor / destructor
		Texture() = default;
		Texture(int w, int h) : width(w), height(h) { pixels = new float3[w * h]; }
		//Texture(int w, int h) : width(w), height(h) { pixels = (uint*)MALLOC64(w * h * sizeof(uint)); }
		~Texture() { }
		//~Texture() { FREE64(pixels); }
		// data members
		int width = 0, height = 0;
		float3* pixels;
		//uint* pixels = 0;
	};
}