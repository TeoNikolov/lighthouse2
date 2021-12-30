/* common_settings.h - Copyright 2019 Utrecht University

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

	   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   The settings and classes in this file are global:
   - available in host and device code
   - the same for each core.
   Settings that can be configured per core can be found in core_settings.h.
*/

#pragma once

// global settings
#define CACHEIMAGES					// imported images will be saved to bin files (faster)

enum SceneType {
	NONE = 0,
	BVH = 1,
	BVH100K = 2,
	RAYPACKET = 3,
	WHITTED = 4,
	SPOTLIGHT = 5,
	DIRECTIONALLIGHT = 6,
	AREALIGHT = 7,
	DIELECTRIC = 8,
	PATHTRACER = 9
};

// global settings
#define CACHEIMAGES					// imported images will be saved to bin files (faster)
// #define ZIPIMGBINS				// cached images will be zipped (slower but smaller)

// default screen size
#define SCRWIDTH			720
#define SCRHEIGHT			406

// settings
#define COUNTMASK 0x3FFFFFFF
#define AXISMASK 0xC0000000
#define SCENETYPE SceneType::PATHTRACER
#define DEFAULTCAMERAPOSITION float3({0, 0, 11.994759})
#define DEFAULTCAMERADIRECTION float3({-0.0029208283, -0.021947961, -0.99975485})
#define DIRECTIONALCAMERAPOSITION float3({0, 5, 11.994759})
#define DIELECTRICCAMERAPOSITION float3({8.0302057, 0.43783143, -6.3139606})
#define DIELECTRICCAMERADIRECTION float3({-0.57642096, -0.01068913, 0.817083})
#define FOXCAMERAPOSITION float3({-203.01425, 118.17175, 39.328934})
#define FOXCAMERADIRECTION float3({0.80994076, -0.32763427, -0.48646861})
#define PATHTRACERCAMERAPOSITION float3({3.6539237, 0.055381246, 12.382888})
#define PATHTRACERCAMERADIRECTION float3({-0.44229344, -0.015886446, -0.89672971})
#define PATHTRACERCAMERAPOSITION2 float3({-1.0222399, 0.35080683, 3.4248352})
#define PATHTRACERCAMERADIRECTION2 float3({0.55486554, 0.18313919, -0.81153208})
#define PATHTRACERCAMERAPOSITION3 float3({2.8449929, 7.8905797, 6.0796051})
#define PATHTRACERCAMERADIRECTION3 float3({-0.29400155, -0.74716669, -0.5960747})
// pathtrace.gltf scene
// front-view
#define PATHTRACERCAMERAPOSITION4 float3({-1.3425646, 3.1349559, 3.680275})
#define PATHTRACERCAMERADIRECTION4 float3({0.72637802, -0.21603858, -0.65245867})
// back-view
#define PATHTRACERCAMERAPOSITION5 float3({9.2037716, 1.4102788, -4.4751797})
#define PATHTRACERCAMERADIRECTION5 float3({-0.89656603, -0.11002588, 0.42902642})
#define DEFAULTFOV 90.0f
#define MULTITHREAD TRUE
#define USEBVH TRUE
#define USETOPBVH FALSE
#define USEPACKETS FALSE
#define PACKETSIZE 16


#define RECURSIONLIMIT 1024 // keep large for uncapped russian roulette

// path tracer variables
#define ENABLE_SKYBOX TRUE
#define ENABLE_MICROFACETS TRUE
#define ENABLE_IMPORTANCESAMPLING TRUE
#define ENABLE_MULTIPLEIMPORTANCESAMPLING TRUE
#define ENABLE_STRATIFICATION TRUE
#define ENABLE_RUSSIANROULETTE TRUE
#define ENABLE_DOF TRUE
#define STRATACOUNTX 4
#define STRATACOUNTY 4

// debug
#define COLLECTDATA FALSE // changes printing format per frame
#define DATACOUNT 128 // same as number of frames
#define TRAVERSALDEPTHMODULATION FALSE
#define TDMODULATIONRATE 1.0f / 24.0f
#define VERBOSELOGGING FALSE
#define VERBOSEBVHCONSTRUCTION FALSE
#define VERBOSEBVHTRAVERSAL FALSE // Warning: rapes the FPS

// skydome defines
// #define IBL						// calculate pdf and cdf for ibl renderer
// #define TESTSKY					// red/green/blue area lights for debugging
#define IBLWIDTH			512
#define IBLHEIGHT			256
#define IBLWBITS			9
#define IBLHBITS			8

// low level settings
#define PI					3.14159265358979323846264f
#define INVPI				0.31830988618379067153777f
#define INV2PI				0.15915494309189533576888f
#define TWOPI				6.28318530717958647692528f
#define SQRT_PI_INV			0.56418958355f
#define LARGE_FLOAT			1e34f
#define EPSILON				0.0001f
#define MINROUGHNESS		0.0001f	// minimal GGX roughness
#define BLACK				make_float3( 0 )
#define WHITE				make_float3( 1 )
#define MIPLEVELCOUNT		5

// file format versions
#define BINTEXFILEVERSION	0x10001001

// tools

// nan chasing
#ifndef __OPENCLCC__
#define FIXNAN_FLOAT3(a)	{if(!isfinite(a.x+a.y+a.z))a=make_float3(0);}
#define FIXNAN_FLOAT4(a)	{if(!isfinite(a.x+a.y+a.z+a.w))a=make_float4(0);}
#define REPORTNAN_FLOAT3(a)	{if(!isfinite(a.x+a.y+a.z))printf("getting NaNs here!");}
#define REPORTNAN_FLOAT4(a)	{if(!isfinite(a.x+a.y+a.z))printf("getting NaNs here!");}
#else
#define FIXNAN_FLOAT3(a)	{if(isnan(a.x+a.y+a.z))a=make_float3(0);}
#define FIXNAN_FLOAT4(a)	{if(isnan(a.x+a.y+a.z+a.w))a=make_float4(0);}
#define REPORTNAN_FLOAT3(a)	{if(isnan(a.x+a.y+a.z))printf("getting NaNs here!");}
#define REPORTNAN_FLOAT4(a)	{if(isnan(a.x+a.y+a.z))printf("getting NaNs here!");}
#endif

// EOF