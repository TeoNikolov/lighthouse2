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
const Mesh* meshes;
vector<Light*> lights;
vector<Material*> mats;
vector<Texture*> texs;
int lightsCount, matsCount, texsCount;
Ray ray = Ray();
bool doLinesPerFrame = false;
int linesPerFrame = 4;
int lineCounter = 0;
int linesElapsed = 0;
double start; // multi-thread timing
bool useAccumulator = false;
float3 accumulator[SCRWIDTH * SCRHEIGHT];

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::Init                                                           |
//  |  Initialization.                                                      LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::Init()
{
	if (TESTTYPE == TestType::AREALIGHT) {
		useAccumulator = true;
	}
	// initialize core
}

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::SetTarget                                                      |
//  |  Set the OpenGL texture that serves as the render target.             LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::SetTarget( GLTexture* target )
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
			l->v1v0 = areaLightTris[i].vertex1 - areaLightTris[i].vertex0;
			l->v2v0 = areaLightTris[i].vertex2 - areaLightTris[i].vertex0;
			l->v0 = areaLightTris[i].vertex0;
			l->v1 = areaLightTris[i].vertex1;
			l->v2 = areaLightTris[i].vertex2;
		} else { // light already stored
			l = (AreaLight*)lights[lightIdx];
		}
		l->area += areaLightTris[i].area;
		l->tris.push_back(areaLightTris[i]);
		l->radiance = areaLightTris[i].radiance;
		printf("Area light instIdx: %d, Total light area %f\n", areaLightTris[i].instIdx, l->area);
	}
	printf("helper size: %d\n", helper.size());
}

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::SetGeometry                                                    |
//  |  Set the geometry data for a model.                                   LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::SetGeometry( const int meshIdx, const float4* vertexData, const int vertexCount, const int triangleCount, const CoreTri* triangleData, const uint* alphaFlags )
{
	Mesh newMesh;
	// copy the supplied vertices; we cannot assume that the render system does not modify
	// the original data after we leave this function.
	newMesh.vertices = new float4[vertexCount];
	newMesh.vcount = vertexCount;
	memcpy( newMesh.vertices, vertexData, vertexCount * sizeof( float4 ) );
	// copy the supplied 'fat triangles'
	newMesh.triangles = new CoreTri[vertexCount / 3];
	memcpy( newMesh.triangles, triangleData, (vertexCount / 3) * sizeof( CoreTri ) );
	meshes.push_back( newMesh );
	/*for (int i = 0; i < newMesh.vcount / 3; i++) {
		printf("%d: %f, %f, %f | %f, %f, %f\n", i,
			newMesh.triangles[i].u0,
			newMesh.triangles[i].u1,
			newMesh.triangles[i].u2,
			newMesh.triangles[i].v0,
			newMesh.triangles[i].v1,
			newMesh.triangles[i].v2
		);
	}
	printf("MeshIdx %d, Tricount %d\n", meshIdx, triangleCount);*/
}

void RenderCore::SetInstance(const int instanceIdx, const int modelIdx, const mat4& transform) {
	printf("InstanceIdx %d, modelIdx %d\n", instanceIdx, modelIdx);
}

void RenderCore::SetMaterials(CoreMaterial* mat, const CoreMaterialEx* matEx, const int materialCount) {
	// copy the supplied array of materials
	for (int i = 0; i < materialCount; i++)
	{
		printf("Loading Material... %d\n", i);
		Material* m;
		if (i < mats.size()) m = mats[i];
		else mats.push_back(m = new Material());
		m->texture = 0;
		int texID = matEx[i].texture[TEXTURE0];
		if (texID == -1)
		{
			float r = mat[i].diffuse_r, g = mat[i].diffuse_g, b = mat[i].diffuse_b;
			m->diffuse.x = (float)mat[i].diffuse_r;
			m->diffuse.y = (float)mat[i].diffuse_g;
			m->diffuse.z = (float)mat[i].diffuse_b;
			m->absorption.x = 1.0f - (float)mat[i].transmittance_r;
			m->absorption.y = 1.0f - (float)mat[i].transmittance_g;
			m->absorption.z = 1.0f - (float)mat[i].transmittance_b;
			m->specularity = mat[i].specular();
			m->ior = mat[i].eta();
			//m->specularity = (float)((mat[i].parameters.x & 0xff0000) >> 16);
			//m->ior = (float)(mat[i].parameters.w & 0xff);
			if (m->diffuse.x > 1.0 || m->diffuse.y > 1.0 || m->diffuse.z > 1.0) {
				m->emissive = true;
			}
			// assign material type
			if (m->specularity > 0.5) {
				m->type = MIRROR;
			} else if (mat[i].flags & 0x1) {
				m->type = DIELECTRIC;
			} else {
				m->type = DIFFUSE;
			}

			switch (TESTTYPE) {
			case TestType::WHITTED: {
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
			}
			case TestType::SPOTLIGHT: {
				break;
			}
			case TestType::DIRECTIONALLIGHT: {
				break;
			}
			case TestType::AREALIGHT: {
				break;
			}
			case TestType::DIELECTRIC: {
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
			}
			}

			/*printf("=== %d ===\n", i);
			printf("Absorption RGB %f, %f, %f | ", m->absorption.x, m->absorption.y, m->absorption.z);
			printf("Spec %f | ", m->specularity);
			printf("IOR %f | ", m->ior);
			printf("Diffuse RGB %f, %f, %f\n", m->diffuse.x, m->diffuse.y, m->diffuse.z);*/
		}
		else
		{
			printf("Material has texture...\n", i);
			m->textured = true;
			m->texture = texs[texID];
			m->texture->width = mat[i].texwidth0; // we know this only now, so set it properly
			m->texture->height = mat[i].texheight0;
		}
	}

}

void RenderCore::SetTextures(const CoreTexDesc* tex, const int textureCount) {
	//texs = new CoreTexDesc[textureCount];
	//memcpy(texs, tex, textureCount * sizeof(CoreTexDesc));
	//texsCount = textureCount;

	// copy the supplied array of texture descriptors
	for (int i = 0; i < textureCount; i++)
	{
		printf("Setting texture... %d\n", i);
		Texture* t;
		if (i < texs.size()) t = texs[i];
		else texs.push_back(t = new Texture());
		t->pixels = new float3[tex[i].pixelCount];
		for (int p = 0; p < tex[i].pixelCount; p++) {
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
	if (useAccumulator) {
		if (converge) {
			memset(accumulator, 0, sizeof(accumulator));
			linesElapsed = 0;
		}
	}
	float3 scrdX = (view.p2 - view.p1) / SCRWIDTH; // screen plane horizontal delta
	float3 scrdY = (view.p3 - view.p1) / SCRHEIGHT; // screen plane vertical delta
	if (doLinesPerFrame) {
		// not accumulator-friendly with 'omp parallel'
		for (int l = 0; l < linesPerFrame; l++) {
			for (int i = 0; i < SCRWIDTH; i++) {
				float3 scrXY = scrdX * i + scrdY * lineCounter;
				ray.o = view.pos;
				ray.d = normalize(view.p1 + scrXY - view.pos);
				float3 radiance = Trace(ray, 0);
				if (useAccumulator) {
					accumulator[lineCounter * SCRWIDTH + i] += radiance;
					radiance = accumulator[lineCounter * SCRWIDTH + i] / (1 + linesElapsed / SCRHEIGHT);
				}
				uint r = (uint)(radiance.x * 255.0);
				uint g = (uint)(radiance.y * 255.0);
				uint b = (uint)(radiance.z * 255.0);
				if (r > 255) r = 255;
				if (g > 255) g = 255;
				if (b > 255) b = 255;
				screen->Plot(i, lineCounter, r + (g << 8) + (b << 16));
			}
			lineCounter++;
			linesElapsed++;
			if (lineCounter >= SCRHEIGHT) {
				printf("elapsed: %f\n", omp_get_wtime() - start);
				lineCounter = 0;
				start = omp_get_wtime();
			}
		}
	} else {
		screen->Clear();
		start = omp_get_wtime();
		#pragma omp parallel for collapse(2) schedule(dynamic) private(ray) num_threads(8)
		for (int j = 0; j < SCRHEIGHT; j++) {
			for (int i = 0; i < SCRWIDTH; i++) {
				float3 scrXY = scrdX * i + scrdY * j;
				ray.o = view.pos;
				ray.d = normalize(view.p1 + scrXY - view.pos);
				float3 radiance = Trace(ray, 0);
				if (useAccumulator) {
					accumulator[j * SCRWIDTH + i] += radiance;
					radiance = accumulator[j * SCRWIDTH + i] / (1.0f + (float)linesElapsed / SCRHEIGHT);
				}
				uint r = (uint)(radiance.x * 255.0);
				uint g = (uint)(radiance.y * 255.0);
				uint b = (uint)(radiance.z * 255.0);
				if (r > 255) r = 255;
				if (g > 255) g = 255;
				if (b > 255) b = 255;
				screen->Plot(i, j, r + (g << 8) + (b << 16));
			}
		}
		linesElapsed += SCRHEIGHT;
		printf("elapsed: %f\n", omp_get_wtime() - start);
	}

	// copy pixel buffer to OpenGL render target texture
	glBindTexture( GL_TEXTURE_2D, targetTextureID );
	glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, screen->width, screen->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, screen->pixels );
}

float3 RenderCore::Trace(const Ray& ray, const int depth) {
	// return black if ray recursion limit is reached
	if (depth == RECURSIONLIMIT) {
		return float3({ 0.0, 0.0, 0.0 });
	}
	HitInfo hitInfo = HitInfo();
	NearestIntersection(ray, &hitInfo);
	// return black if no intersection
	if (!hitInfo.intersects) {
		return float3({ 0.0, 0.0, 0.0 });
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
		return hitColor * DirectIllumination(hitInfo.intersection, hitInfo.hitNormal);
	}
	if (m->type == MIRROR) {
		Ray r = Ray();
		r.o = hitInfo.intersection + hitInfo.hitNormal * EPSILON;
		r.d = reflect(ray.d, hitInfo.hitNormal);
		return hitColor * Trace(r, depth + 1);
	}
	if (m->type == DIELECTRIC) {
		float n1 = 1.0f;
		float n2 = m->ior;
		float beerDistance = 0.0f;
		if (ray.inDielectric) {
			// exiting material
			n1 = m->ior;
			n2 = 1.0f;
			beerDistance = length(hitInfo.intersection - ray.o);
		}
		float n12 = n1 / n2;
		//float t1 = abs(dot(make_float3(-ray.d.x, -ray.d.y, -ray.d.z), hitInfo.hitNormal));
		float t1 = abs(dot(ray.d, hitInfo.hitNormal));
		float k = 1 - (n12 * n12) * (1 - t1 * t1);
		float3 R = reflect(ray.d, hitInfo.hitNormal);
		if (k < 0.0f) {
			// Total Internal Reflection
			Ray r = Ray();
			r.o = hitInfo.intersection + hitInfo.hitNormal * EPSILON * (1 + -2 * ray.inDielectric);
			r.d = R;
			r.inDielectric = ray.inDielectric;
			return hitColor * Trace(r, depth + 1);
		} else {
			float3 T = n12 * ray.d + hitInfo.hitNormal * (n12 * t1 - sqrt(k));
			Ray r = Ray();
			Ray t = Ray();
			r.o = hitInfo.intersection + hitInfo.hitNormal * EPSILON * (1 + -2 * ray.inDielectric);
			t.o = hitInfo.intersection - hitInfo.hitNormal * EPSILON * (1 + -2 * ray.inDielectric);
			r.d = R;
			t.d = T;
			r.inDielectric = ray.inDielectric;
			t.inDielectric = !ray.inDielectric;

			// Fresnel
			//float t2 = abs(dot(T, -hitInfo.hitNormal));
			float t2 = abs(dot(T, hitInfo.hitNormal));
			float n1t1 = n1 * t1;
			float n1t2 = n1 * t2;
			float n2t1 = n2 * t1;
			float n2t2 = n2 * t2;
			float Fs = sqr((n1t1 - n2t2) / (n1t1 + n2t2));
			float Fp = sqr((n1t2 - n2t1) / (n1t2 + n2t1));
			float Fr = (Fs + Fp) / 2.0f;

			// Beer's Law
			float3 c = (Trace(t, depth + 1) * (1.0 - Fr) + Trace(r, depth + 1) * Fr);
			float beerR = exp(-m->absorption.x * beerDistance);
			float beerG = exp(-m->absorption.y * beerDistance);
			float beerB = exp(-m->absorption.z * beerDistance);
			c *= float3({beerR, beerG, beerB});
			return hitColor * c;
		}
	}
}

bool RenderCore::NearestIntersection(const Ray ray, HitInfo* hitInfo) {
	HitInfo tempHitInfo = HitInfo();
	for (Mesh& mesh : meshes) {
		tempHitInfo.Reset();
		if (mesh.Intersect(ray, &tempHitInfo)) {
			if (tempHitInfo.tHit < hitInfo->tHit) {
				*hitInfo = tempHitInfo;
			}
		}
	}
	return hitInfo->intersects;
}

float3 RenderCore::DirectIllumination(const float3 intersection, const float3 hitNormal) {
	float3 illumination = float3({0.0, 0.0, 0.0});
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
			} else { // intersection is in outer radius
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
			//float3 x = (l->tris[0].vertex1 - l->tris[0].vertex0) * r1;
			//float3 y = (l->tris[0].vertex2 - l->tris[0].vertex0) * r2;
			//float3 xy = l->tris[0].vertex0 + x + y;
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

//  +-----------------------------------------------------------------------------+
//  |  RenderCore::Shutdown                                                       |
//  |  Free all resources.                                                  LH2'19|
//  +-----------------------------------------------------------------------------+
void RenderCore::Shutdown()
{
	delete screen;
}

// EOF