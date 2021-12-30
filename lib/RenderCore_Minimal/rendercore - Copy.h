/* rendercore.h - Copyright 2019 Utrecht University

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

	   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#pragma once

namespace lh2core
{

//  +-----------------------------------------------------------------------------+
//  |  RenderCore                                                                 |
//  |  Encapsulates device code.                                            LH2'19|
//  +-----------------------------------------------------------------------------+
class RenderCore : public CoreAPI_Base
{
public:
	// methods
	void Init();
	void SetTarget( GLTexture* target, const uint spp );
	void SetGeometry( const int meshIdx, const float4* vertexData, const int vertexCount, const int triangleCount, const CoreTri* triangles, const uint* alphaFlags = 0 );
	void Render( const ViewPyramid& view, const Convergence converge );
	CoreStats GetCoreStats() const override;
	void Shutdown();

	// unimplemented for the minimal core
	inline void SetProbePos(const int2 pos) override {};
	inline void Setting(const char* name, float value) override {};
	inline void SetSkyData(const float3* pixels, const uint width, const uint height, const mat4& worldToLight) override {};
	
	// implemented for the minimal core
	void SetTextures(const CoreTexDesc* tex, const int textureCount);
	void SetMaterials(CoreMaterial* mat, const int materialCount);
	void SetLights(const CoreLightTri* areaLights, const int areaLightCount,
		const CorePointLight* pointLights, const int pointLightCount,
		const CoreSpotLight* spotLights, const int spotLightCount,
		const CoreDirectionalLight* directionalLights, const int directionalLightCount);
	void SetInstance(const int instanceIdx, const int modelIdx, const mat4& transform);
	void UpdateToplevel() override;

	// custom methods
	bool TraceShadow(const Ray);
	float3 DiffuseReflection(float3 normal);
	float3 WhittedTrace(Ray& ray, const int depth, int& traverseDepth);
	float3 WhittedDirectIllumination(const float3 intersection, const float3 hitNormal);
	float3 PathTrace(Ray& ray, const int depth, int& traverseDepth);
	float3 PathDirectIllumination(const float3 intersection, const float3 hitNormal);
	void Trace(RayPacket& rayPacket, Frustum& frustum, const int depth, int* traverseDepth, float3* colors);
	bool NearestIntersection(Ray& ray, HitInfo* hitInfo);
	bool NearestIntersection(RayPacket& rayPacket, Frustum& frustum, HitInfo* hitInfo);

	CoreStats coreStats;							// rendering statistics
private:
	// data members
	Bitmap* screen = 0;								// temporary storage of RenderCore output; will be copied to render target
	int targetTextureID = 0;						// ID of the target OpenGL texture
	vector<Mesh*> meshes;							// mesh data storage
	vector<MeshInstance*> instances;
	vector<Light*> lights;
	vector<Material*> mats;
	vector<Texture*> texs;
	TopBVH* topBVH;
};

} // namespace lh2core

// EOF