/* rendercore.cpp - Copyright 2019 Utrecht University

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

#include "core_settings.h"

using namespace lh2core;
int lightsCount, matsCount, texsCount;
Ray ray = Ray();

RayPacket rayPacket = RayPacket();
HitInfo tempHitInfoArray[PACKETSIZE];
HitInfo hitInfoArray[PACKETSIZE];
int ttdArray[PACKETSIZE];
int tdArray[PACKETSIZE];
float3 radiance[PACKETSIZE];

bool doLinesPerFrame = false;
int linesPerFrame = 4;
int lineCounter = 0;
int linesElapsed = 0;
bool useAccumulator = false;
float3 accumulator[SCRWIDTH * SCRHEIGHT];
bool geometryInitialized = false;

// analysis
double lastNoiseVal = 0.0f;
double start; // multi-thread timing
double elapsedTime = 0;
double bvhConstructionTime = 0;
uint frames = 1;
//  +-----------------------------------------------------------------------------+
//  |  RenderCore::Init                                                           |
//  |  Initialization.                                                      LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::Init()
{
	topBVH = new TopBVH();
	for (int i = 0; i < PACKETSIZE; i++) {
		tempHitInfoArray[i] = HitInfo();
		hitInfoArray[i] = HitInfo();
		radiance[i] = make_float3(0, 0, 0);
	}
	memset(ttdArray, 0, sizeof(int) * PACKETSIZE);
	memset(tdArray, 0, sizeof(int) * PACKETSIZE);
}

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::SetTarget                                                      |
//  |  Set the OpenGL texture that serves as the render target.             LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::SetTarget( GLTexture* target, const uint spp )
{
	// synchronize OpenGL viewport
	targetTextureID = target->ID;
	if (screen != 0 && target->width == screen->width && target->height == screen->height) return; // nothing changed
	delete screen;
	screen = new Bitmap( target->width, target->height );
}

void RenderCore::SetLights(const CoreLightTri* areaLightTris, const int areaLightCount,
	const CorePointLight* pointLights, const int pointLightCount,
	const CoreSpotLight* spotLights, const int spotLightCount,
	const CoreDirectionalLight* directionalLights, const int directionalLightCount) {
	printf("Light Count: area %d, point %d, spot %d, directional %d\n", areaLightCount, pointLightCount, spotLightCount, directionalLightCount);
	
	for (int i = 0; i < pointLightCount; i++) {
		PointLight* l = new PointLight();
		l->radiance = pointLights[i].radiance;
		l->position = pointLights[i].position;
		lights.push_back(l);
	}
	for (int i = 0; i < spotLightCount; i++) {
		SpotLight* l = new SpotLight();
		l->radiance = spotLights[i].radiance;
		l->position = spotLights[i].position;
		l->direction = spotLights[i].direction;
		l->cosInner = spotLights[i].cosInner;
		l->cosOuter = spotLights[i].cosOuter;
		if (l->cosOuter > l->cosInner) {
			printf("Warning: spotlight had an outer radius smaller than the inner radius. Clamping to inner radius.");
			l->cosOuter = l->cosInner;
		}
		lights.push_back(l);
	}
	for (int i = 0; i < directionalLightCount; i++) {
		DirectionalLight* l = new DirectionalLight();
		l->radiance = directionalLights[i].radiance;
		l->position = directionalLights[i].direction;
		lights.push_back(l);
	}
	vector<tuple<int, int>> helper;
	for (int i = 0; i < areaLightCount; i++) {
		useAccumulator = true;
		int lightIdx = -1;
		for (int j = 0; j < helper.size(); j++) {
			if (get<0>(helper[j]) == areaLightTris[i].instIdx) {
				lightIdx = get<1>(helper[j]);
				break;
			}
		}
		AreaLight* l;
		if (lightIdx == -1) { // light not stored yet			
			printf("Storing area light...\n");
			l = new AreaLight();
			tuple<int, int> t = make_tuple(areaLightTris[i].instIdx, lights.size());
			helper.push_back(t);
			lights.push_back(l);
			l->normal = areaLightTris[i].N;
			l->v1v0 = areaLightTris[i].vertex1 - areaLightTris[i].vertex0;
			l->v2v0 = areaLightTris[i].vertex2 - areaLightTris[i].vertex0;
			l->v0 = areaLightTris[i].vertex0;
			l->v1 = areaLightTris[i].vertex1;
			l->v2 = areaLightTris[i].vertex2;
		} else { // light already stored
			printf("Area light already stored (idx == %d)...\n", lightIdx);
			l = (AreaLight*)lights[lightIdx];
		}
		l->area += areaLightTris[i].area;
		l->tris.push_back(areaLightTris[i]);
		l->radiance = areaLightTris[i].radiance;
		printf("TriIdx == %d\n", areaLightTris[i].triIdx);
		printf("Area light instIdx: %d, Total light area %f\n", areaLightTris[i].instIdx, l->area);
	}
}

void RenderCore::SetSkyData(const float3* pixels, const uint width, const uint height, const mat4& worldToLight) {
	sky = new Sky();
	sky->width = width;
	sky->height = height;
	sky->worldToLight = worldToLight;
	sky->pixels = new float3[width * height + width * 2];
	memcpy(sky->pixels, pixels, width * height * sizeof(float3));
	sky->pixels[width * height + width * 2] = sky->pixels[0];
}

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::SetGeometry                                                    |
//  |  Set the geometry data for a model.                                   LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::SetGeometry( const int meshIdx, const float4* vertexData, const int vertexCount, const int triangleCount, const CoreTri* triangleData, const uint* alphaFlags )
{
	//if (meshIdx == 6) return; // remove geometry for skybox testing
	if (geometryInitialized) {
		Mesh* mesh = meshes[meshIdx];
		if (triangleCount != mesh->triangleCount) runtime_error("Updated mesh has different TRICount. Fatal error.");
		memcpy(mesh->vertices, vertexData, vertexCount * sizeof(float4));
		memcpy(mesh->triangles, triangleData, (vertexCount / 3) * sizeof(CoreTri));
		for (int i = 0; i < triangleCount; i++) {
			mesh->centroids[i] =
				(triangleData[i].vertex0 +
					triangleData[i].vertex1 +
					triangleData[i].vertex2) / 3.0f;
		}
		double timestart = omp_get_wtime();
		mesh->bvh->Update();
		int a = 0;
#if VERBOSELOGGING
		printf("Mesh %d: Geometry refitted in %.3fs\n", meshIdx, (omp_get_wtime() - timestart));
#endif
	} else {
#if VERBOSELOGGING
		printf("==============================================\n");
		printf("Mesh %d: Processing mesh with %d triangles...\n", meshIdx, triangleCount);
#endif
		Mesh* newMesh = new Mesh();
		// copy the supplied vertices; we cannot assume that the render system does not modify
		// the original data after we leave this function.
		newMesh->vertices = new float4[vertexCount];
		//newMesh.vcount = vertexCount;
		newMesh->triangleCount = triangleCount;
		memcpy(newMesh->vertices, vertexData, vertexCount * sizeof(float4));
		// copy the supplied 'fat triangles'
		newMesh->triangles = new CoreTri[triangleCount];
		newMesh->centroids = new float3[triangleCount];
		for (int i = 0; i < triangleCount; i++) {
			newMesh->centroids[i] =
				(triangleData[i].vertex0 +
					triangleData[i].vertex1 +
					triangleData[i].vertex2) / 3.0f;
		}
		newMesh->meshIdx = meshIdx;
		memcpy(newMesh->triangles, triangleData, (vertexCount / 3) * sizeof(CoreTri));
		meshes.push_back(newMesh);

		if (USEBVH) {
			double bvhTime = omp_get_wtime();
			newMesh->ConstructBVH();
			bvhTime = omp_get_wtime() - bvhTime;
			bvhConstructionTime += bvhTime;
			switch (meshIdx) {
			case 0:
				newMesh->bvh->updateType = BVHUpdateType::Refit;
				break;
			case 1:
				newMesh->bvh->updateType = BVHUpdateType::Rebuild;
				break;
			}
//#if VERBOSELOGGING
			printf("Mesh %d: BVH constructed in %fms | Total thus far %fms\n",
				meshIdx, bvhTime * 1000.0f, bvhConstructionTime * 1000.0f);
//#endif
		}
	}
}

// assume instances arrive ordered and only the last instance may have -1 modelIdx
void RenderCore::SetInstance(const int instanceIdx, const int modelIdx, const mat4& transform) {
	if (modelIdx > -1) {
		if (instances.size() < instanceIdx + 1) { // does instance exist ? no : yes
			printf("Instance %d: Storing with modelIdx %d\n", instanceIdx, modelIdx);
			instances.push_back(new MeshInstance(instanceIdx, meshes[modelIdx], transform));
			instances[instanceIdx]->UpdateTransform(transform);
		} else {
			instances[instanceIdx]->UpdateTransform(transform);
		}
	}
}

void RenderCore::SetMaterials( CoreMaterial* mat, const int materialCount ) {
	// copy the supplied array of materials
	for (int i = 0; i < materialCount; i++)
	{
		printf("Material %d: Loading...\n", i);
		Material* m;
		if (i < mats.size()) m = mats[i];
		else mats.push_back(m = new Material());
		m->texture = 0;
		int texID = mat[i].color.textureID;
		if (texID == -1)
		{
			float r = mat[i].color.value.x / 255.0f, g = mat[i].color.value.y / 255.0f, b = mat[i].color.value.z / 255.0f;
			m->diffuse.x = mat[i].color.value.x;
			m->diffuse.y = mat[i].color.value.y;
			m->diffuse.z = mat[i].color.value.z;
			m->absorption.x = mat[i].absorption.value.x;
			m->absorption.y = mat[i].absorption.value.y;
			m->absorption.z = mat[i].absorption.value.z;
			m->specular = m->diffuse;
			m->specularity = mat[i].specular.value;
			m->ior = mat[i].eta.value;
			m->gloss = 5.0f;
			//m->specularity = (float)((mat[i].parameters.x & 0xff0000) >> 16);
			//m->ior = (float)(mat[i].parameters.w & 0xff);
			if (m->diffuse.x > 1.0 || m->diffuse.y > 1.0 || m->diffuse.z > 1.0) {
				m->emissive = true;
			}
			// assign material type
			if (m->specularity > 0.5) {
				m->type = MIRROR;
			} else if (mat[i].transmission.value > EPSILON) {
			//} else if (mat[i].flags & 0x1) {
				m->type = DIELECTRIC;
			} else {
				m->type = DIFFUSE;
			}

			switch (SCENETYPE) {
			case SceneType::WHITTED: {
				switch (i) {
				case 1: // cone material
					m->specularity = 1.0;
					m->type = MIRROR;
					m->diffuse.x = 1.0;
					m->diffuse.y = 1.0;
					m->diffuse.z = 1.0;
					break;
				}
				break;
			} case SceneType::SPOTLIGHT: {
				break;
			} case SceneType::DIRECTIONALLIGHT: {
				break;
			} case SceneType::AREALIGHT: {
				break;
			} case SceneType::PATHTRACER: {
				switch (i) {
				case 3: // cube material
					m->absorption.x = 0.0f;
					m->absorption.y = 0.0f;
					m->absorption.z = 0.02f;
					m->diffuse.x = 1.0;
					m->diffuse.y = 1.0;
					m->diffuse.z = 1.0;
					m->ior = 1.2;
					m->type = DIELECTRIC;
					break;
				case 5: // cone material
					m->specularity = 1.0;
					m->type = MIRROR;
					m->diffuse.x = 1.0;
					m->diffuse.y = 1.0;
					m->diffuse.z = 1.0;
					break;
				case 6: // sphere material
					break;
					m->absorption.x = 0.0f;
					m->absorption.y = 0.0f;
					m->absorption.z = 0.0f;
					m->diffuse.x = 1.0;
					m->diffuse.y = 0.75;
					m->diffuse.z = 1.0;
					m->ior = 1.08;
					//m->type = DIELECTRIC;
					break;
				case 2:
					m->specular = float3({0.95f, 0.64f, 0.54f});
					break;
				}
				break;
			} case SceneType::DIELECTRIC: {
				switch (i) {
				case 0: // cube material
					m->absorption.x = 0.25f;
					m->absorption.y = 0.0f;
					m->absorption.z = 0.25f;
					m->diffuse.x = 1.0;
					m->diffuse.y = 1.0;
					m->diffuse.z = 1.0;
					m->ior = 1.05;
					m->type = DIELECTRIC;
					break;
				case 1: // cone material
					m->absorption.x = 0.0f;
					m->absorption.y = 0.0f;
					m->absorption.z = 0.0f;
					m->diffuse.x = 1.0;
					m->diffuse.y = 1.0;
					m->diffuse.z = 1.0;
					m->ior = 1.02;
					m->type = DIELECTRIC;
					break;
				}
				break;
			}
			}

			/*printf("=== %d ===\n", i);
			printf("Absorption RGB %f, %f, %f | ", m->absorption.x, m->absorption.y, m->absorption.z);
			printf("Spec %f | ", m->specularity);
			printf("IOR %f | ", m->ior);
			printf("Diffuse RGB %f, %f, %f\n", m->diffuse.x, m->diffuse.y, m->diffuse.z);*/
		} else {
			printf("Material %d: Is textured.\n", i);
			m->textured = true;
			m->texture = texs[texID];
			//m->texture->width = mat[i].texwidth0; // we know this only now, so set it properly
			//m->texture->height = mat[i].texheight0;
		}
	}

}

void RenderCore::SetTextures(const CoreTexDesc* tex, const int textures) {
//texs = new CoreTexDesc[textureCount];
	//memcpy(texs, tex, textureCount * sizeof(CoreTexDesc));
	//texsCount = textureCount;

	// copy the supplied array of texture descriptors
	for (int i = 0; i < textures; i++)
	{
		printf("Setting texture... %d\n", i);
		Texture* t;
		if (i < texs.size()) t = texs[i];
		else texs.push_back(t = new Texture());
		t->pixels = new float3[tex[i].pixelCount];
		t->width = tex[i].width;
		t->height = tex[i].height;
		for (uint p = 0; p < tex[i].pixelCount; p++) {
			t->pixels[p] = make_float3(
				(float)(tex[i].idata[p].x) / 255.0f,
				(float)(tex[i].idata[p].y) / 255.0f,
				(float)(tex[i].idata[p].z) / 255.0f);
		}
		//t->pixels = (uint*)MALLOC64(tex[i].pixelCount * sizeof(uint));
		//if (tex[i].idata) memcpy(t->pixels, tex[i].idata, tex[i].pixelCount * sizeof(uint));
		//else memcpy(t->pixels, 0, tex[i].pixelCount * sizeof(uint) /* assume integer textures */);
		// Note: texture width and height are not known yet, will be set when we get the materials.
	}
}

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::Render                                                         |
//  |  Produce one image.                                                   LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::Render( const ViewPyramid& view, const Convergence converge )
{
	geometryInitialized = true;
	if (converge) {
		printf("converge\n");
		memset(accumulator, 0, sizeof(accumulator));
		linesElapsed = 0;
		frames = 1;
	}

	float3 scrdX = (view.p2 - view.p1) / SCRWIDTH; // screen plane horizontal delta
	float3 scrdY = (view.p3 - view.p1) / SCRHEIGHT; // screen plane vertical delta
	float3 scrXD = normalize(view.p2 - view.p1); // the horizontal direction of the plane
	float3 scrYD = normalize(view.p3 - view.p1); // the vertical direction of the plane (negative Y)

	screen->Clear();
	start = omp_get_wtime();

#if USEPACKETS
	//float3* radiance = new float3[PACKETSIZE];
	//for (int i = 0; i < PACKETSIZE; i++) {
	//	radiance[i] = float3({ 0, 0, 0 });
	//}
	//int* td = new int[PACKETSIZE];
	//memset(td, 0, sizeof(int) * PACKETSIZE);
	//int sqrtps = (int)sqrtf(PACKETSIZE);
	//rayPacket.rays[0].o = view.pos;
	//rayPacket.rays[1].o = view.pos;
	//rayPacket.rays[2].o = view.pos;
	//rayPacket.rays[3].o = view.pos;
	//rayPacket.rays[0].d = normalize(view.p1 - view.pos);
	//rayPacket.rays[1].d = normalize(view.p2 - view.pos);
	//rayPacket.rays[2].d = normalize(view.p2 + view.p3 - view.p1 - view.pos);
	//rayPacket.rays[3].d = normalize(view.p3 - view.pos);
	//Frustum frustum = Frustum(
	//		rayPacket.rays[0],
	//		rayPacket.rays[1],
	//		rayPacket.rays[2],
	//		rayPacket.rays[3]);
	//		Trace(rayPacket, frustum, 0, td, radiance);
	//		int xp = 0;
	//		int yp = 0;
	//		for (int k = 0; k < PACKETSIZE; k++) {
	//			switch (k) {
	//			case 1:
	//				xp = SCRWIDTH - 1;
	//				break;
	//			case 2:
	//				xp = SCRWIDTH - 1;
	//				yp = SCRHEIGHT - 1;
	//				break;
	//			case 3:
	//				yp = SCRHEIGHT - 1;
	//				break;
	//			}
	//			uint r = (uint)(radiance[k].x * 255.0);
	//			uint g = (uint)(radiance[k].y * 255.0);
	//			uint b = (uint)(radiance[k].z * 255.0);
	//			if (r > 255) r = 255;
	//			if (g > 255) g = 255;
	//			if (b > 255) b = 255;
	//			screen->Plot(xp, yp, r + (g << 8) + (b << 16));
	//		}

	for (int i = 0; i < PACKETSIZE; i++) {
		radiance[i] = float3({ 0, 0, 0 });
		tempHitInfoArray[i] = HitInfo();
		hitInfoArray[i] = HitInfo();
	}
	memset(tdArray, 0, sizeof(int) * PACKETSIZE);
	int sqrtps = (int)sqrtf(PACKETSIZE);
	#pragma omp parallel for schedule(dynamic) private(rayPacket, ttdArray, tdArray, hitInfoArray, tempHitInfoArray, radiance) num_threads(8)	
	for (int j = 0; j < SCRHEIGHT / sqrtps; j++) {
		for (int i = 0; i < SCRWIDTH / sqrtps; i++) {
			int xsq = i * sqrtps;
			int ysq = j * sqrtps;
			for (int k = 0; k < PACKETSIZE; k++) {
				float3 scrX = scrdX * (xsq + k % sqrtps);
				float3 scrY = scrdY * (ysq + k / sqrtps);
				rayPacket.rays[k].o = view.pos;
				rayPacket.rays[k].d = normalize(view.p1 + scrX + scrY - view.pos);
			}
			Frustum frustum = Frustum(
				rayPacket.rays[0],
				rayPacket.rays[sqrtps - 1],
				rayPacket.rays[PACKETSIZE - 1],
				rayPacket.rays[PACKETSIZE - sqrtps]);
			Trace(rayPacket, frustum, 0, tdArray, radiance);

			for (int k = 0; k < PACKETSIZE; k++) {
				int xPixel = xsq + k % sqrtps;
				int yPixel = ysq + k / sqrtps;
				if (useAccumulator) {
					accumulator[yPixel * SCRWIDTH + xPixel] += radiance[k];
					radiance[k] = accumulator[yPixel * SCRWIDTH + xPixel] / frames;
				}
#if TRAVERSALDEPTHMODULATION
				radiance[k].y += tdArray[k] * TDMODULATIONRATE;
#endif
				uint r = (uint)(radiance[k].x * 255.0);
				uint g = (uint)(radiance[k].y * 255.0);
				uint b = (uint)(radiance[k].z * 255.0);
				if (r > 255) r = 255;
				if (g > 255) g = 255;
				if (b > 255) b = 255;
				screen->Plot(xPixel, yPixel, r + (g << 8) + (b << 16));
			}
		}
	}
#else
#pragma omp parallel for schedule(dynamic) private(ray) num_threads(8)
	for (int j = 0; j < SCRHEIGHT; j++) {
		for (int i = 0; i < SCRWIDTH; i++) {
			int td = 0;
			float3 radiance;
			float3 scrXY = scrdX * i + scrdY * j;
			if (ENABLE_DOF) {
				// sample aperture
				float2 apertureXY = RejectionSampleDisk(view.aperture * 50.0f);
				ray.o = view.pos + scrXD * apertureXY.x + scrYD * apertureXY.y;
				ray.d = normalize(view.p1 + scrXY - ray.o);
			} else {
				ray.o = view.pos;
				ray.d = normalize(view.p1 + scrXY - ray.o);
			}

			if (SCENETYPE == SceneType::PATHTRACER) {
				accumulator[j * SCRWIDTH + i] += PathTrace(ray, 0, td, true, WHITE);
				radiance = accumulator[j * SCRWIDTH + i] / (float)frames;
			} else {
				radiance = WhittedTrace(ray, 0, td);
				if (useAccumulator) {
					accumulator[j * SCRWIDTH + i] += radiance;
					radiance = accumulator[j * SCRWIDTH + i] / (float)frames;
				}
			}
			
#if TRAVERSALDEPTHMODULATION
			radiance.y += td * TDMODULATIONRATE;
#endif
			uint r = (uint)(radiance.x * 255.0);
			uint g = (uint)(radiance.y * 255.0);
			uint b = (uint)(radiance.z * 255.0);
			if (r > 255) r = 255;
			if (g > 255) g = 255;
			if (b > 255) b = 255;
			screen->Plot(i, j, r + (g << 8) + (b << 16));
		}
	}
#endif

	double elapsed = omp_get_wtime() - start;
	elapsedTime += elapsed;
	float totalval = 0.0f;
	float noiseval = 0.0f;
	for (int i = 0; i < SCRWIDTH * SCRHEIGHT; i++) {
		totalval += (accumulator[i].x + accumulator[i].y + accumulator[i].z);
	}
	for (int i = 0; i < SCRWIDTH * SCRHEIGHT; i++) {
		float3 top, right, bottom, left;
		float3 current = accumulator[i];
		float3 noisediff = BLACK;
		int divisor = 0;
		if (i >= SCRWIDTH) {
			top = accumulator[i - SCRWIDTH];
			noisediff += fabs(top - current);
			divisor++;
		}
		if (i % SCRWIDTH < SCRWIDTH - 1) {
			right = accumulator[i + 1];
			noisediff += fabs(right - current);
			divisor++;
		}
		if (i < (SCRHEIGHT - 1) * SCRWIDTH) {
			bottom = accumulator[i + SCRWIDTH];
			noisediff += fabs(bottom - current);
			divisor++;
		}
		if (i % SCRWIDTH > 0) {
			left = accumulator[i - 1];
			noisediff += fabs(left - current);
			divisor++;
		}
		noiseval += noisediff.x + noisediff.y + noisediff.z / ((float) divisor);
	}
	noiseval /= frames;
	if (COLLECTDATA) {
		if (frames > DATACOUNT) {
			printf("F");
		}
		else {
			printf("%d\n", (int)(noiseval));
		}
	} else {
		printf("#%d: T<elapsed, total>: %.3fs, %.3fs RGB <total, diff, %% diff>: %d, %d, %.2f\n",
			frames, elapsed, elapsedTime, (int)(totalval / frames),
			(int)(noiseval), (lastNoiseVal - noiseval) / lastNoiseVal * 100.0f);
	}
	lastNoiseVal = noiseval;
	frames++;
	// copy pixel buffer to OpenGL render target texture
	glBindTexture( GL_TEXTURE_2D, targetTextureID );
	glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, screen->width, screen->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, screen->pixels );
}

float2 RenderCore::RejectionSampleDisk( float radius ) {
	while (true) {
		float rx = ((float)rand() / (float)RAND_MAX) * 2.0f - 1.0f;
		float ry = ((float)rand() / (float)RAND_MAX) * 2.0f - 1.0f;
		float dist2 = rx * rx + ry * ry;
		// terminate if sample not in sphere
		if (dist2 > 1.0) continue;
		return float2({ rx * radius, ry * radius });
	}
}

float3 RenderCore::CosineWeightedDiffuseReflection(float3 normal) {
	// random point on hemisphere
	float r0 = (float)rand() / (float)RAND_MAX;
	float r1 = (float)rand() / (float)RAND_MAX;
	float r = sqrtf( r0 );
	float theta = 2.0f * PI * r1;
	float x = r * cosf( theta );
	float y = r * sinf( theta );
	float3 xyz = float3({x, sqrtf( 1.0f - r0 ), y});
	// convert to tangent space
	float3 Nt;
	if (fabs(normal.x) > fabs(normal.y)) {
		Nt = normalize(float3({normal.z, 0, -normal.x}));
	} else {
		Nt = normalize(float3({0, -normal.z, normal.y}));
	}
	float3 Nb = -cross(normal, Nt);
	float3 result = float3({ xyz.x * Nb.x + xyz.y * normal.x + xyz.z * Nt.x,
					xyz.x * Nb.y + xyz.y * normal.y + xyz.z * Nt.y,
					xyz.x * Nb.z + xyz.y * normal.z + xyz.z * Nt.z });
	return result;
}

float3 RenderCore::DiffuseReflection(float3 normal) {
	while (true) {
		float rx = ((float)rand() / (float)RAND_MAX) * 2.0f - 1.0f;
		float ry = ((float)rand() / (float)RAND_MAX) * 2.0f - 1.0f;
		float rz = ((float)rand() / (float)RAND_MAX) * 2.0f - 1.0f;
		float3 d = float3({ rx, ry, rz });
		float dist2 = dot(d, d);
		// terminate if sample not in sphere
		if (dist2 > 1.0) continue;
		d = normalize(d);
		bool inHalfSpace = dot(d, normal) > 0.0;
		return d * (1.0f - 2.0f * !inHalfSpace);
	}
}

float3 RenderCore::PathTrace(Ray& ray, int depth, int& traverseDepth, bool lastSpecular, float3 pweight) {
	// return black if ray recursion limit is reached

#if (1)
	HitInfo hitInfo = HitInfo();
	float3 T = WHITE, E = BLACK;
	while (1) {
		if (depth > RECURSIONLIMIT) {
			break;
		}
		// terminate if ray left the scene
		float3 oldBRDF;
		float3 oldHitNormal = hitInfo.hitNormal;
		float oldTHit = hitInfo.tHit;
		hitInfo.Reset();
		int ttd = 0;
		if (!NearestIntersection(ray, &hitInfo)) {
			if (ENABLE_SKYBOX) {
				E += T * sky->GetColor(ray);
			}
			break;
		}
		traverseDepth = hitInfo.traverseDepth;
		Material* m = mats[hitInfo.triangle->material];
		float3 hitDiffuse = m->diffuse;
		if (m->textured) {
			hitDiffuse = m->GetTextureColor(hitInfo.uv);
		}
		// terminate if we hit a light source
		if (m->emissive) {
			if (lastSpecular) {
				if (dot(hitInfo.hitNormal, -ray.d) > 0.0f) {
					E += T * hitDiffuse;
				}
			} else {
				if (ENABLE_MULTIPLEIMPORTANCESAMPLING) {
					AreaLight* light = (AreaLight*)lights[hitInfo.triangle->ltriIdx / 2];
					float cos_o = dot(-ray.d, light->normal);
					float cos_i = dot(ray.d, oldHitNormal);
					if ((cos_o > 0) && (cos_i > 0)) {
						float misPDF, brdfPDF, lightPDF;
						float dist2 = dot(ray.d * hitInfo.tHit, ray.d * hitInfo.tHit); // distance to light
						float solidAngle = (cos_o * light->area) / dist2;
						lightPDF = 1.0f / solidAngle;
						if (ENABLE_IMPORTANCESAMPLING) {
							brdfPDF = cos_i / PI;
						} else {
							brdfPDF = 1 / (2.0f * PI);
						}
						misPDF = lightPDF + brdfPDF;
						E += T * lights.size() * light->radiance * (cos_i / misPDF) * oldBRDF;
					}
				}
			}
			break;
		}

#ifndef loalefa
		// surface interaction
		if (m->type == MIRROR) {
			// continue in fixed direction
			lastSpecular = true;
			ray.o = hitInfo.intersection + hitInfo.hitNormal * EPSILON;
			ray.d = reflect(ray.d, hitInfo.hitNormal);
			T *= hitDiffuse;
			continue;
		} else if (m->type == DIELECTRIC) {
			lastSpecular = true;
			// determine if inside medium and appropriate IORs
			bool inDielectric = dot(hitInfo.hitNormal, ray.d) > 0.0;
			const int boolInv = 1 + -2 * inDielectric;
			float n1 = 1.0f * !inDielectric + m->ior * inDielectric;
			float n2 = 1.0f * inDielectric + m->ior * !inDielectric;
			float n12 = n1 / n2;
			float t1 = abs(dot(ray.d, hitInfo.hitNormal));
			float k = 1 - (n12 * n12) * (1 - t1 * t1);

			// total internal reflection
			if (k < 0.0) {
				ray.o = hitInfo.intersection + hitInfo.hitNormal * EPSILON * boolInv;
				ray.d = reflect(ray.d, hitInfo.hitNormal);
				T *= hitDiffuse;
			} else {
				// determine probability of particle being reflective (Fresnel)
				float3 Tr = n12 * ray.d + hitInfo.hitNormal * (n12 * t1 - sqrt(k));
				float t2 = abs(dot(Tr, hitInfo.hitNormal));
				float n1t1 = n1 * t1;
				float n1t2 = n1 * t2;
				float n2t1 = n2 * t1;
				float n2t2 = n2 * t2;
				float Fs = sqr((n1t1 - n2t2) / (n1t1 + n2t2));
				float Fp = sqr((n1t2 - n2t1) / (n1t2 + n2t1));
				float Fr = (Fs + Fp) / 2.0f;
				float rprob = max(EPSILON, (float)rand() / (float)RAND_MAX);

				if (rprob < Fr) {
					// particle is reflective
					ray.o = hitInfo.intersection + hitInfo.hitNormal * EPSILON * boolInv;
					ray.d = reflect(ray.d, hitInfo.hitNormal);
					T *= hitDiffuse;
				} else {
					float beerDistance = length(hitInfo.intersection - ray.o) * inDielectric;
					// particle is transmissive
					ray.o = hitInfo.intersection - hitInfo.hitNormal * EPSILON * boolInv;
					ray.d = Tr;
					// Beer's Law

					float beerR = exp(-m->absorption.x * beerDistance);
					float beerG = exp(-m->absorption.y * beerDistance);
					float beerB = exp(-m->absorption.z * beerDistance);
					T *= float3({ beerR, beerG, beerB });
					T *= hitDiffuse;
				}
			}
			continue;
		} else {
			lastSpecular = false;
		}
#endif
		float hemiPDF;
		float3 R, Ei;
		float3 BRDF = hitDiffuse / PI;
		// sample a random light source
		if (ENABLE_MICROFACETS) {
			E += T * PathDirectIllumination(
				hitInfo.intersection, hitInfo.hitNormal, ray.d, m->specular, hitDiffuse, m->gloss, true);
		} else {
			E += T * PathDirectIllumination(
				hitInfo.intersection, hitInfo.hitNormal, ray.d, m->specular, hitDiffuse, m->gloss, true) * BRDF;
		}
		// continue random walk
		if (ENABLE_IMPORTANCESAMPLING) {
			R = CosineWeightedDiffuseReflection(hitInfo.hitNormal);
			hemiPDF = max(dot(R, hitInfo.hitNormal) / PI, EPSILON);
		}
		else {
			R = DiffuseReflection(hitInfo.hitNormal);
			hemiPDF = 1.0f / (2.0f * PI);
		}
		if (ENABLE_MICROFACETS) {
			float3 V = -ray.d;
			float3 H = normalize(V + R);
			float NH = max(EPSILON, dot(hitInfo.hitNormal, H));
			float VH = max(EPSILON, dot(V, H));
			float brdfD = ((m->gloss + 2.0f) / (2.0f * PI)) * pow(NH, m->gloss);
			float brdfG = min(1.0f, min(
				(2.0f * NH * dot(hitInfo.hitNormal, V)) / VH,
				(2.0f * NH * dot(hitInfo.hitNormal, R)) / VH
			));
			float3 brdfF = m->specular + (1 - m->specular) * pow((1 - (dot(H, R))), 5.0f);
			BRDF = brdfF * brdfG * brdfD /
				(4.0f * max(EPSILON, dot(hitInfo.hitNormal, R)) * max(EPSILON, dot(hitInfo.hitNormal, V)));
		}
		ray.o = hitInfo.intersection + hitInfo.hitNormal * EPSILON;
		ray.d = R;
		// update throughput
		if (ENABLE_RUSSIANROULETTE) {
			//float q = clamp(max({ hitDiffuse.x, hitDiffuse.y, hitDiffuse.z }), 0.1f, 0.9f);
			//float q = clamp(1.0f - max({ pweight.x, pweight.y, pweight.z }), 0.05f, 0.95f);
			float q = clamp(1.0f - max({ T.x, T.y, T.z }), 0.05f, 0.95f);
			if (((float)rand() / (float)RAND_MAX) < q) {
				break;
			} else {
				float invq = 1.0f - q;
				T *= (dot(R, hitInfo.hitNormal) / (invq * hemiPDF)) * BRDF;
			}
		} else {
			T *= (dot(R, hitInfo.hitNormal) / hemiPDF) * BRDF;
		}
		oldBRDF = BRDF;
		depth++;
	}
	return E;
#else
	if (depth > RECURSIONLIMIT) {
		return BLACK;
	}
	// terminate if ray left the scene
	HitInfo hitInfo = HitInfo();
	if (!NearestIntersection(ray, &hitInfo)) {
		if (ENABLE_SKYBOX) {
			return sky->GetColor(ray);
		} else {
			return BLACK;
		}
	}
	int ttd = 0;
	traverseDepth = hitInfo.traverseDepth;
	Material* m = mats[hitInfo.triangle->material];
	float3 hitDiffuse = m->diffuse;
	// terminate if we hit a light source
	if (m->emissive) {
		if (lastSpecular) {
			if (dot(hitInfo.hitNormal, -ray.d) > 0.0f) {
				return hitDiffuse;
			}
		}
		return BLACK;
		if (ENABLE_MULTIPLEIMPORTANCESAMPLING) {
			AreaLight* light = (AreaLight*)lights[hitInfo.triangle->ltriIdx / 2];
			float cos_o = dot(-ray.d, light->normal);
			float cos_i = dot(ray.d, hitInfo.hitNormal);
			if ((cos_o > 0) && (cos_i > 0)) {
				float dist2 = dot(ray.d, ray.d); // distance to light
				float dist = sqrtf(dist2);
				float solidAngle = (cos_o * light->area) / dist2;
				float misPDF;
				float lightPDF = 1.0f / solidAngle;
				float brdfPDF = cos_i / PI;
				misPDF = lightPDF + brdfPDF;
				return lights.size() * light->radiance * (cos_i / misPDF) / brdfPDF;
			}
		}
		return BLACK;
	}
	if (m->textured) {
		hitDiffuse = m->GetTextureColor(hitInfo.uv);
	}

#ifndef loalefa
	// surface interaction
	if (m->type == MIRROR) {
		// continue in fixed direction
		ray.o = hitInfo.intersection + hitInfo.hitNormal * EPSILON;
		ray.d = reflect(ray.d, hitInfo.hitNormal);
		return hitDiffuse * PathTrace(ray, depth, traverseDepth, true, pweight);
	}
	if (m->type == DIELECTRIC) {
		// determine if inside medium and appropriate IORs
		bool inDielectric = dot(hitInfo.hitNormal, ray.d) > 0.0;
		const int boolInv = 1 + -2 * inDielectric;
		float n1 = 1.0f * !inDielectric + m->ior * inDielectric;
		float n2 = 1.0f * inDielectric + m->ior * !inDielectric;
		float n12 = n1 / n2;
		float t1 = abs(dot(ray.d, hitInfo.hitNormal));
		float k = 1 - (n12 * n12) * (1 - t1 * t1);

		// total internal reflection
		if (k < 0.0) {
			ray.o = hitInfo.intersection + hitInfo.hitNormal * EPSILON * boolInv;
			ray.d = reflect(ray.d, hitInfo.hitNormal);
			return hitDiffuse * PathTrace(ray, depth, ttd, true, pweight);
		} else {
			// determine probability of particle being reflective (Fresnel)
			float3 T = n12 * ray.d + hitInfo.hitNormal * (n12 * t1 - sqrt(k));
			float t2 = abs(dot(T, hitInfo.hitNormal));
			float n1t1 = n1 * t1;
			float n1t2 = n1 * t2;
			float n2t1 = n2 * t1;
			float n2t2 = n2 * t2;
			float Fs = sqr((n1t1 - n2t2) / (n1t1 + n2t2));
			float Fp = sqr((n1t2 - n2t1) / (n1t2 + n2t1));
			float Fr = (Fs + Fp) / 2.0f;
			float rprob = max(EPSILON, (float)rand() / (float)RAND_MAX);

			if ( rprob < Fr ) {
				// particle is reflective
				ray.o = hitInfo.intersection + hitInfo.hitNormal * EPSILON * boolInv;
				ray.d = reflect(ray.d, hitInfo.hitNormal);
				return hitDiffuse * PathTrace(ray, depth, ttd, true, pweight);
			} else {
				float beerDistance = length(hitInfo.intersection - ray.o) * inDielectric;
				// particle is transmissive
				ray.o = hitInfo.intersection - hitInfo.hitNormal * EPSILON * boolInv;
				ray.d = T;
				// Beer's Law
				float3 c = PathTrace(ray, depth, ttd, true, pweight);
				float beerR = exp(-m->absorption.x * beerDistance);
				float beerG = exp(-m->absorption.y * beerDistance);
				float beerB = exp(-m->absorption.z * beerDistance);
				c *= float3({ beerR, beerG, beerB });
				return hitDiffuse * c;
			}
		}
	}
#endif

	float misPDF;
	float3 R, Ei, BRDF = hitDiffuse / PI;
	// randomly sample a random light source
	float3 Ld = PathDirectIllumination(hitInfo.intersection, hitInfo.hitNormal, false) * BRDF;
	// continue random walk
	if (ENABLE_IMPORTANCESAMPLING) {
		R = CosineWeightedDiffuseReflection(hitInfo.hitNormal);
		misPDF = max(dot(R, hitInfo.hitNormal) / PI, EPSILON);
	} else {
		R = DiffuseReflection(hitInfo.hitNormal);
		misPDF = 1.0f / (2.0f * PI);
	}
	pweight *= BRDF / misPDF;
	ray.o = hitInfo.intersection + hitInfo.hitNormal * EPSILON;
	ray.d = R;
	// update throughput
	if (ENABLE_RUSSIANROULETTE) {
		float q = clamp(1.0f - max({ pweight.x, pweight.y, pweight.z}), 0.05f, 0.95f);
		if (((float)rand() / (float)RAND_MAX) < q) {
			Ei = BLACK;
		} else {
			float invq = 1.0f - q;
			pweight /= invq;
			Ei = PathTrace(ray, depth + 1, traverseDepth, false, pweight) * dot(R, hitInfo.hitNormal) / (invq * misPDF); // with IS
		}
	} else {
		Ei = PathTrace(ray, depth + 1, traverseDepth, false, pweight) * dot(R, hitInfo.hitNormal) / misPDF; // with IS
	}
	return BRDF * Ei + Ld;
#endif
}

float3 RenderCore::PathDirectIllumination(const float3 intersection, const float3 hitNormal, float3 viewDir,
	float3 specular, float3 diffuse, float alpha, bool testMIS) {
	float3 illumination = BLACK;
	HitInfo tHitInfo = HitInfo();
	Ray shadowRay = Ray();

	int lightIdx = rand() % lights.size();
	AreaLight* l = (AreaLight*)lights[lightIdx];
	// construct vector to random point on light
	float3 xy;
	int stratumID = frames % (STRATACOUNTX * STRATACOUNTY);
	float r1 = ((float)rand() / (float)RAND_MAX);
	float r2 = ((float)rand() / (float)RAND_MAX);
	if (ENABLE_STRATIFICATION) {
		float stratum_width = 1.0f / STRATACOUNTX;
		float stratum_height = 1.0f / STRATACOUNTY;
		float stratum_x = (stratumID % STRATACOUNTX) * stratum_width;
		float stratum_y = (stratumID / STRATACOUNTY) * stratum_height;
		float xpos = stratum_x + stratum_width * r1;
		float ypos = stratum_y + stratum_height * r2;
		xy = l->v0 + l->v1v0 * xpos + l->v2v0 * ypos;
	} else {
		xy = l->v0 + l->v1v0 * r1 + l->v2v0 * r2;
	}
	float3 L = xy - intersection;
	float dist2 = dot(L, L); // distance to light
	float dist = sqrtf(dist2);
	L = L / dist; // direction to light
	float cos_o = dot(-L, l->normal);
	float cos_i = dot(L, hitNormal);
	if ((cos_o > 0) && (cos_i > 0)) {
		// light is not behind surface point, trace shadow ray
		shadowRay.o = intersection + hitNormal * EPSILON;
		shadowRay.d = L;
		shadowRay.tMax = dist - 2 * EPSILON;
		// terminate if light is occluded
		if (NearestIntersection(shadowRay, &tHitInfo)) return BLACK;
		float misPDF, brdfPDF, lightPDF;
		float solidAngle = (cos_o * l->area) / dist2;
		lightPDF = 1.0f / solidAngle;
		if (ENABLE_MULTIPLEIMPORTANCESAMPLING && testMIS) {
			if (ENABLE_IMPORTANCESAMPLING) {
				brdfPDF = cos_i / PI;
			} else {
				brdfPDF = 1 / (2.0f * PI);
			}
			misPDF = lightPDF + brdfPDF;
		} else {
			misPDF = lightPDF;
		}
		if (ENABLE_MICROFACETS) {
			float3 V = -viewDir;
			float3 H = normalize(V + L);
			float NH = max(EPSILON, dot(hitNormal, H));
			float VH = max(EPSILON, dot(V, H));
			float brdfD = ((alpha + 2.0f) / (2.0f * PI)) * pow(NH, alpha);
			float brdfG = min(1.0f, min(
				(2.0f * NH * dot(hitNormal, V)) / VH,
				(2.0f * NH * dot(hitNormal, L)) / VH
			));
			float3 brdfF = specular + (1.0f - specular) * pow((1.0f - (dot(H, L))), 5.0f);
			float3 BRDF = brdfF * brdfG * brdfD /
				(4.0f * max(EPSILON, dot(hitNormal, L)) * max(EPSILON, dot(hitNormal, V)));
			illumination = lights.size() * l->radiance * (cos_i / misPDF) * BRDF;
		} else {
			illumination = lights.size() * l->radiance * (cos_i / misPDF);	
		}

	}
	return illumination;
}

float3 RenderCore::WhittedTrace(Ray& ray, const int depth, int& traverseDepth) {
	// return black if ray recursion limit is reached
	if (depth > RECURSIONLIMIT) {
		return BLACK;
	}
	HitInfo hitInfo = HitInfo();
	NearestIntersection(ray, &hitInfo);
	traverseDepth = hitInfo.traverseDepth;
	int ttd = 0;
	// return black if no intersection
	if (!hitInfo.intersects) {
		return BLACK;
	}
	Material* m = mats[hitInfo.triangle->material];
	float3 hitColor = m->diffuse;
	if (m->textured) {
		hitColor = m->GetTextureColor(hitInfo.uv);
	}
	if (m->type == DIFFUSE) {
		if (m->emissive) { // woah, we hit an emissive surface
			if (dot(ray.d, hitInfo.hitNormal) > 0.0) { // WOAH, we can't see the surface!
				return float3({ 0.0, 0.0, 0.0 });
			}
			return hitColor;
		}
		return hitColor * WhittedDirectIllumination(hitInfo.intersection, hitInfo.hitNormal);
	} else if (m->type == MIRROR) {
		Ray r = Ray();
		r.o = hitInfo.intersection + hitInfo.hitNormal * EPSILON;
		r.d = reflect(ray.d, hitInfo.hitNormal);
		return hitColor * WhittedTrace(r, depth + 1, ttd);
	} else if (m->type == DIELECTRIC) {
		bool inDielectric = dot(hitInfo.hitNormal, ray.d) > 0.0;
		const int boolInv = 1 + -2 * inDielectric;
		float n1 = 1.0f * !inDielectric + m->ior * inDielectric;
		float n2 = 1.0f * inDielectric + m->ior * !inDielectric;
		float beerDistance = length(hitInfo.intersection - ray.o) * inDielectric;
		float n12 = n1 / n2;
		float t1 = abs(dot(ray.d, hitInfo.hitNormal));
		float k = 1 - (n12 * n12) * (1 - t1 * t1);
		float3 R = reflect(ray.d, hitInfo.hitNormal);
		if (k < 0.0f) {
			// Total Internal Reflection
			Ray r = Ray();
			r.o = hitInfo.intersection + hitInfo.hitNormal * EPSILON * boolInv;
			r.d = R;
			return hitColor * WhittedTrace(r, depth + 1, ttd);
		} else {
			float3 T = n12 * ray.d + hitInfo.hitNormal * (n12 * t1 - sqrt(k));
			Ray r = Ray();
			Ray t = Ray();
			r.o = hitInfo.intersection + hitInfo.hitNormal * EPSILON * boolInv;
			t.o = hitInfo.intersection - hitInfo.hitNormal * EPSILON * boolInv;
			r.d = R;
			t.d = T;

			// Fresnel
			float t2 = abs(dot(T, hitInfo.hitNormal));
			float n1t1 = n1 * t1;
			float n1t2 = n1 * t2;
			float n2t1 = n2 * t1;
			float n2t2 = n2 * t2;
			float Fs = sqr((n1t1 - n2t2) / (n1t1 + n2t2));
			float Fp = sqr((n1t2 - n2t1) / (n1t2 + n2t1));
			float Fr = (Fs + Fp) / 2.0f;

			// Beer's Law
			float3 c = (WhittedTrace(t, depth + 1, ttd) * (1.0 - Fr) + WhittedTrace(r, depth + 1, ttd) * Fr);
			float beerR = exp(-m->absorption.x * beerDistance);
			float beerG = exp(-m->absorption.y * beerDistance);
			float beerB = exp(-m->absorption.z * beerDistance);
			c *= float3({ beerR, beerG, beerB });
			return hitColor * c;
		}
	} else {
		throw runtime_error("Invalid material type: How did you even get here?");
	}
}

void RenderCore::Trace(RayPacket& rayPacket, Frustum& frustum, const int depth, int* traverseDepth, float3* colors) {
	// return black if ray recursion limit is reached

	for (int i = 0; i < PACKETSIZE; i++) colors[i] = float3({0, 0, 0});
	if (depth > RECURSIONLIMIT) {
		return;
	}
	for (int i = 0; i < PACKETSIZE; i++) hitInfoArray[i].Reset();
	bool packetIntersects = NearestIntersection(rayPacket, frustum, hitInfoArray);
	for (int i = 0; i < PACKETSIZE; i++) traverseDepth[i] = hitInfoArray[i].traverseDepth;
	memset(ttdArray, 0, sizeof(int) * PACKETSIZE);

	// return black if no intersection
	if (!packetIntersects) {
		return;
	}
	Material** m = new Material*[PACKETSIZE];
	for (int i = 0; i < PACKETSIZE; i++) {
		m[i] = mats[hitInfoArray[i].triangle->material];
		colors[i] = m[i]->diffuse;
		if (m[i]->textured) {
			colors[i] = m[i]->GetTextureColor(hitInfoArray[i].uv);
		}
	}
	for (int i = 0; i < PACKETSIZE; i++) {
		if (!hitInfoArray[i].intersects) {
			continue;
		}
		if (m[i]->type == DIFFUSE) {
			if (m[i]->emissive) { // woah, we hit an emissive surface
				if (dot(rayPacket.rays[i].d, hitInfoArray[i].hitNormal) > 0.0) { // WOAH, we can't see the surface!
					continue;
				}
			}
			colors[i] = colors[i] * WhittedDirectIllumination(hitInfoArray[i].intersection, hitInfoArray[i].hitNormal);
		} else if (m[i]->type == MIRROR) {
			Ray r = Ray();
			r.o = hitInfoArray[i].intersection + hitInfoArray[i].hitNormal * EPSILON;
			r.d = reflect(rayPacket.rays[i].d, hitInfoArray[i].hitNormal);
			colors[i] = colors[i] * WhittedTrace(r, depth + 1, ttdArray[i]);
		} else if (m[i]->type == DIELECTRIC) {
			bool inDielectric = dot(hitInfoArray[i].hitNormal, rayPacket.rays[i].d) > 0.0;
			const int boolInv = 1 + -2 * inDielectric;
			float n1 = 1.0f * !inDielectric + m[i]->ior * inDielectric;
			float n2 = 1.0f * inDielectric + m[i]->ior * !inDielectric;
			float beerDistance = length(hitInfoArray[i].intersection - rayPacket.rays[i].o) * inDielectric;
			float n12 = n1 / n2;
			float t1 = abs(dot(rayPacket.rays[i].d, hitInfoArray[i].hitNormal));
			float k = 1 - (n12 * n12) * (1 - t1 * t1);
			float3 R = reflect(rayPacket.rays[i].d, hitInfoArray[i].hitNormal);
			if (k < 0.0f) {
				// Total Internal Reflection
				Ray r = Ray();
				r.o = hitInfoArray[i].intersection + hitInfoArray[i].hitNormal * EPSILON * boolInv;
				r.d = R;
				colors[i] = colors[i] * WhittedTrace(r, depth + 1, ttdArray[i]);
			} else {
				float3 T = n12 * rayPacket.rays[i].d + hitInfoArray[i].hitNormal * (n12 * t1 - sqrt(k));
				Ray r = Ray();
				Ray t = Ray();
				r.o = hitInfoArray[i].intersection + hitInfoArray[i].hitNormal * EPSILON * boolInv;
				t.o = hitInfoArray[i].intersection - hitInfoArray[i].hitNormal * EPSILON * boolInv;
				r.d = R;
				t.d = T;

				// Fresnel
				float t2 = abs(dot(T, hitInfoArray[i].hitNormal));
				float n1t1 = n1 * t1;
				float n1t2 = n1 * t2;
				float n2t1 = n2 * t1;
				float n2t2 = n2 * t2;
				float Fs = sqr((n1t1 - n2t2) / (n1t1 + n2t2));
				float Fp = sqr((n1t2 - n2t1) / (n1t2 + n2t1));
				float Fr = (Fs + Fp) / 2.0f;

				// Beer's Law
				float3 c = (WhittedTrace(t, depth + 1, ttdArray[i]) * (1.0 - Fr) + WhittedTrace(r, depth + 1, ttdArray[i]) * Fr);
				float beerR = exp(-m[i]->absorption.x * beerDistance);
				float beerG = exp(-m[i]->absorption.y * beerDistance);
				float beerB = exp(-m[i]->absorption.z * beerDistance);
				c *= float3({ beerR, beerG, beerB });
				colors[i] = colors[i] * c;
			}
		} else {
			throw runtime_error("Invalid material type: How did you even get here?");
		}
	}
}

bool RenderCore::NearestIntersection(RayPacket& rayPacket, Frustum& frustum, HitInfo* hitInfo) {
	bool intersects = false;
	for (Mesh* mesh : meshes) {
		for (int i = 0; i < PACKETSIZE; i++) tempHitInfoArray[i].Reset();
		if (mesh->Intersect(rayPacket, frustum, tempHitInfoArray)) {
			for (int i = 0; i < PACKETSIZE; i++) {
				if (tempHitInfoArray[i].tHit < hitInfo[i].tHit) {
					intersects = true;
					hitInfo[i] = tempHitInfoArray[i];
				}
			}
		}
		for (int i = 0; i < PACKETSIZE; i++) hitInfo[i].traverseDepth += tempHitInfoArray[i].traverseDepth;
	}
	return intersects;
}

bool RenderCore::NearestIntersection(Ray& ray, HitInfo* hitInfo) {
#if USETOPBVH
	return (topBVH->Intersect(ray, hitInfo, hitInfo->traverseDepth));
#else
	HitInfo tempHitInfo = HitInfo();
	for (Mesh* mesh : meshes) {
		tempHitInfo.Reset();
		if (mesh->Intersect(ray, &tempHitInfo)) {
			if (tempHitInfo.tHit < hitInfo->tHit) {
				*hitInfo = tempHitInfo;
				hitInfo->lightIdx = mesh->lightIdx;
			}
		}
		hitInfo->traverseDepth += tempHitInfo.traverseDepth;
	}
	return hitInfo->intersects;
#endif
}

float3 RenderCore::WhittedDirectIllumination(const float3 intersection, const float3 hitNormal) {
	float3 illumination = float3({ 0.0, 0.0, 0.0 });
	HitInfo tHitInfo = HitInfo();
	Ray shadowRay = Ray();
	for (int i = 0; i < lights.size(); i++) {
		tHitInfo.Reset();
		switch (lights[i]->type) {
		case LightType::PointLight: {
			PointLight* l = (PointLight*)lights[i];
			float3 L = l->position - intersection;
			float D2 = dot(L, L); // distance to light
			L = L / sqrtf(D2); // direction to light
			float a = dot(L, hitNormal);
			if (a < 0.0) { // is light behind the face?
				continue;
			}
			shadowRay.o = intersection + L * EPSILON;
			shadowRay.d = L;
			NearestIntersection(shadowRay, &tHitInfo); // is the light occluded?
			if (tHitInfo.tHit * tHitInfo.tHit < D2) {
				continue;
			}
			illumination += l->radiance * a / D2;
			continue;
		}
		case LightType::SpotLight:
		{
			SpotLight* l = (SpotLight*)lights[i];
			float3 L = l->position - intersection;
			float D2 = dot(L, L); // distance to light
			L = L / sqrtf(D2); // direction to light
			shadowRay.o = intersection + L * EPSILON;
			shadowRay.d = L;
			NearestIntersection(shadowRay, &tHitInfo); // is the light occluded?
			if (tHitInfo.tHit * tHitInfo.tHit < D2) {
				continue;
			}
			float a = dot(L, hitNormal);
			if (a < 0.0f) { // is light behind the face?
				continue;
			}
			float b = dot(-L, l->direction); // cos between light dir and dir from light to intersection
			if (b < l->cosOuter) { // is intersection outside spotlight cone?
				continue;
			}
			if (b > l->cosInner) { // is intersection in inner radius?
				illumination += l->radiance * a / D2;
			}
			else { // intersection is in outer radius
				float baseDiff = l->cosOuter - l->cosInner;
				float cosDiff = l->cosOuter - b;
				illumination += (l->radiance * a * cosDiff) / (D2 * baseDiff);
			}
			continue;
		}
		case LightType::DirectionalLight:
		{
			DirectionalLight* l = (DirectionalLight*)lights[i];
			float3 L = normalize(-l->position); // direction to light
			float a = dot(L, hitNormal);
			if (a < 0.0) { // is light behind the face?
				continue;
			}
			shadowRay.o = intersection + L * EPSILON;
			shadowRay.d = L;
			if (NearestIntersection(shadowRay, &tHitInfo)) { // is the light occluded?
				continue;
			}
			illumination += l->radiance * a;
			continue;
		}
		case LightType::AreaLight: {
			AreaLight* l = (AreaLight*)lights[i];
			//CoreLightTri* tri = &l->tris[0];
			float r1 = ((float)rand() / (float)RAND_MAX);
			float r2 = ((float)rand() / (float)RAND_MAX);
			float3 xy = l->v0 + l->v1v0 * r1 + l->v2v0 * r2;
			float3 L = xy - intersection;
			float D2 = dot(L, L); // distance to light
			L = L / sqrtf(D2); // direction to light
			float c = dot(l->tris[0].N, L);
			if (c >= 0.0) { // is intersection behind light?
				continue;
			}
			float a = dot(L, hitNormal);
			if (a < 0.0) { // is light behind the face?
				continue;
			}
			shadowRay.o = intersection + L * EPSILON;
			shadowRay.d = L;
			NearestIntersection(shadowRay, &tHitInfo);
			tHitInfo.tHit += 2 * EPSILON;
			if (tHitInfo.tHit * tHitInfo.tHit < D2) { // is the light occluded?
				continue;
			}
			illumination += l->radiance * a * -c / D2;
			continue;
		}
		}
	}
	return illumination;
}

void RenderCore::UpdateToplevel() {
#if USETOPBVH
	double st = omp_get_wtime();
	topBVH->Build(instances);
#if VERBOSELOGGING
	printf("TopBVH: Rebuilt in %.5fms\n", (omp_get_wtime() - st) * 1000.0);
#endif
#endif
}

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::GetCoreStats                                                   |
//  |  Get a copy of the counters.                                          LH2'19|
//  +-----------------------------------------------------------------------------+
CoreStats RenderCore::GetCoreStats() const
{
	return coreStats;
}

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::Shutdown                                                       |
//  |  Free all resources.                                                  LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::Shutdown()
{
	delete screen;
}

// EOF