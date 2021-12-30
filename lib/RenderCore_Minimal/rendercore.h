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
class RenderCore
{
public:
	// methods
	void Init();
	void SetTarget( GLTexture* target );
	void SetGeometry( const int meshIdx, const float4* vertexData, const int vertexCount, const int triangleCount, const CoreTri* triangles, const uint* alphaFlags = 0 );
	void SetInstance(const int instanceIdx, const int modelIdx, const mat4& transform = mat4::Identity());
	void SetMaterials(CoreMaterial* mat, const CoreMaterialEx* matEx, const int materialCount);
	void SetTextures(const CoreTexDesc* tex, const int textureCount);
	void Render( const ViewPyramid& view, const Convergence converge );
	void Shutdown();
	void SetLights(const CoreLightTri* areaLights, const int areaLightCount,
		const CorePointLight* pointLights, const int pointLightCount,
		const CoreSpotLight* spotLights, const int spotLightCount,
		const CoreDirectionalLight* directionalLights, const int directionalLightCount);
	void UpdateToplevel();
	bool TraceShadow(const Ray);
	float3 Trace(Ray& ray, const int depth, int& traverseDepth);
	void Trace(RayPacket& rayPacket, Frustum& frustum, const int depth, int* traverseDepth, float3* colors);
	float3 DirectIllumination(const float3 intersection, const float3 hitNormal);
	bool NearestIntersection(Ray& ray, HitInfo* hitInfo);
	bool NearestIntersection(RayPacket& rayPacket, Frustum& frustum, HitInfo* hitInfo);
	// internal methods
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
public:
	CoreStats coreStats;							// rendering statistics
};

} // namespace lh2core

// EOF