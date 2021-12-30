/* main.cpp - Copyright 2019 Utrecht University

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

#include "platform.h"
#include "rendersystem.h"
#include <bitset>

static RenderAPI* renderer = 0;
static GLTexture* renderTarget = 0;
static Shader* shader = 0;
static uint scrwidth = 0, scrheight = 0, scrspp = 1;
static bool camMoved = false, spaceDown = false, hasFocus = true, running = true, animPaused = false;
static std::bitset<1024> keystates;
static std::bitset<8> mbstates;
static string materialFile;

// material editing
HostMaterial currentMaterial;
int currentMaterialID = -1;
static CoreStats coreStats;

#include "main_tools.h"

//  +-----------------------------------------------------------------------------+
//  |  PrepareScene                                                               |
//  |  Initialize a scene.                                                  LH2'19|
//  +-----------------------------------------------------------------------------+
void PrepareScene()
{
#if 1
	if (SCENETYPE != SceneType::NONE) {
		Camera* camera = renderer->GetCamera();
		camera->position = DEFAULTCAMERAPOSITION;
		camera->direction = DEFAULTCAMERADIRECTION;
	}
	switch (SCENETYPE) {
	case SceneType::NONE: {
		materialFile = string("data/pica/pica_materials.xml");
		renderer->AddScene("default_scene.gltf", "data/pica/", mat4::Translate(0, -10.2f, 0));
		int rootNode = renderer->FindNode("RootNode (gltf orientation matrix)");
		renderer->SetNodeTransform(rootNode, mat4::RotateX(-PI / 2));
		renderer->AddPointLight(float3({ 5.0, 7.0, 5.0 }), float3({ 200, 200, 200 }));
		break;
	}
	case SceneType::PATHTRACER: {
		materialFile = string("data/pica/pica_materials.xml");
		int lightMat = renderer->AddMaterial(make_float3(7, 7, 7));
		int quad = renderer->AddQuad(float3({ 0, -1, 0 }), float3({ 0, 10, 0 }), 5.0f, 5.0f, lightMat);
		//renderer->AddScene("default_scene.gltf", "data/pica/");
		renderer->AddScene("pathtracing.gltf", "data/pica/");
		//int quad0 = renderer->AddQuad(float3({ 0, -1, 0 }), float3({-4.5, 10, 0}), 2.0f, 2.0f, lightMat);
		//int quad1 = renderer->AddQuad(float3({ 0, -1, 0 }), float3({0, 10, -4.5}), 2.0f, 2.0f, lightMat);
		//int quad2 = renderer->AddQuad(float3({ 0, -1, 0 }), float3({4.5, 10, 0}), 2.0f, 2.0f, lightMat);
		//int quad3 = renderer->AddQuad(float3({ 0, -1, 0 }), float3({0, 10, 4.5}), 2.0f, 2.0f, lightMat);
		renderer->AddInstance(quad);
		//renderer->AddInstance(quad0);
		//renderer->AddInstance(quad1);
		//renderer->AddInstance(quad2);
		//renderer->AddInstance(quad3);
		Camera* camera = renderer->GetCamera();
		camera->position = PATHTRACERCAMERAPOSITION4;
		camera->direction = PATHTRACERCAMERADIRECTION4;

		auto scene = renderer->GetScene();
		auto sky = new HostSkyDome();
		sky->Load( "data/sky_15.hdr" );
		// Compensate for different evaluation in PBRT
		sky->worldToLight = mat4::RotateX( -PI / 2 );
		scene->SetSkyDome( sky );
		break;
	}
	case SceneType::RAYPACKET: {
		materialFile = string("data/pica/pica_materials.xml");
		renderer->AddPointLight(float3({ 0, 5, 5 }), float3({ 200, 200, 200 }));
		renderer->AddScene("default_scene.gltf", "data/pica/");
		Camera* camera = renderer->GetCamera();
		camera->position = DEFAULTCAMERAPOSITION;
		camera->direction = DEFAULTCAMERADIRECTION;
		break;
	}
	case SceneType::BVH: {
		materialFile = string("data/pica/pica_materials.xml");
		renderer->AddPointLight(float3({ -200, 200, 75 }), float3({ 40000, 40000, 40000 }));
		renderer->AddScene("scene.gltf", "data/pica/fox/", mat4::Translate(0, 0, 0));
		renderer->AddScene("scene.gltf", "data/pica/fox/", mat4::Translate(-50, 0, -150));
		renderer->AddScene("scene.gltf", "data/pica/fox/", mat4::Translate(50, 0, -150));
		renderer->AddScene("untitled.gltf", "data/pica/", mat4::Translate(0, 75, 125));
		int gun = renderer->AddMesh("hp.obj", "data/verybigmeshes/gun2/", 4);
		renderer->AddInstance(gun, mat4::Translate(0, 70, 15) * mat4::RotateY(PI));
		Camera* camera = renderer->GetCamera();
		camera->position = FOXCAMERAPOSITION;
		camera->direction = FOXCAMERADIRECTION;
		break;
	}
	case SceneType::WHITTED: {
		materialFile = string("data/pica/pica_materials.xml");
		renderer->AddScene("default_scene.gltf", "data/pica/", mat4::Translate(0, -10.2f, 0));
		int rootNode = renderer->FindNode("RootNode (gltf orientation matrix)");
		renderer->SetNodeTransform(rootNode, mat4::RotateX(-PI / 2));
		renderer->AddPointLight(float3({ 5.0, 7.0, 5.0 }), float3({ 200, 200, 200 }));
		int cube = renderer->AddMesh("cube_textured.obj", "data/pica/", 2);
		mat4 transform = mat4();
		transform.Identity();
		transform.cell[0] = 4.0;
		int instance1 = renderer->AddInstance(cube, transform);
		break;
	}
	case SceneType::DIELECTRIC: {
		materialFile = string("data/pica/pica_materials.xml");
		renderer->AddScene("default_scene.gltf", "data/pica/", mat4::Translate(0, -10.2f, 0));
		int rootNode = renderer->FindNode("RootNode (gltf orientation matrix)");
		renderer->SetNodeTransform(rootNode, mat4::RotateX(-PI / 2));
		renderer->AddPointLight(float3({ 0.0, 8.0, 15.0 }), float3({ 200, 200, 200 }));
		Camera* camera = renderer->GetCamera();
		camera->position = DIELECTRICCAMERAPOSITION;
		camera->direction = DIELECTRICCAMERADIRECTION;
		camera->FOV = DEFAULTFOV;
		break;
	}
	case SceneType::SPOTLIGHT: {
		materialFile = string("data/pica/pica_materials.xml");
		renderer->AddScene("default_scene.gltf", "data/pica/", mat4::Translate(0, -10.2f, 0));
		int rootNode = renderer->FindNode("RootNode (gltf orientation matrix)");
		renderer->SetNodeTransform(rootNode, mat4::RotateX(-PI / 2));
		renderer->AddSpotLight(float3({ 5, 8, 5 }), normalize(float3({ 0, -1, -1 })),
			0.8, 0.5, float3({ 100, 100, 100 }));
		renderer->AddSpotLight(float3({ -5, 8, 5 }), normalize(float3({ 0, -1, -1 })),
			0.9, 0.8, float3({ 150, 150, 150 }));
		break;
	}
	case SceneType::DIRECTIONALLIGHT: {
		materialFile = string("data/pica/pica_materials.xml");
		renderer->AddScene("directional_scene.gltf", "data/pica/", mat4::Translate(0, -10.2f, 0));
		int rootNode = renderer->FindNode("RootNode (gltf orientation matrix)");
		renderer->SetNodeTransform(rootNode, mat4::RotateX(-PI / 2));
		renderer->AddDirectionalLight(normalize(float3({ -1.0, -2.0, -0.5 })), float3({ 1.0, 1.0, 1.0 }));
		Camera* camera = renderer->GetCamera();
		camera->position = DIRECTIONALCAMERAPOSITION;
		break;
	}
	case SceneType::AREALIGHT: {
		materialFile = string("data/pica/pica_materials.xml");
		renderer->AddScene("default_scene.gltf", "data/pica/", mat4::Translate(0, -10.2f, 0));
		int rootNode = renderer->FindNode("RootNode (gltf orientation matrix)");
		renderer->SetNodeTransform(rootNode, mat4::RotateX(-PI / 2));
		int lightMat = renderer->AddMaterial(make_float3(200, 200, 200));
		int lightQuad = renderer->AddQuad(normalize(make_float3(0, 1.25, -1)), make_float3(0, -3, 9), 2.0f, 4.0f, lightMat);
		int lightInst = renderer->AddInstance(lightQuad);
		break;
	}
	}
#endif
#if 1
	// overhead light, use regular PT
	//int lightMat = renderer->AddMaterial( make_float3( 50, 50, 45 ) );
	//int lightQuad = renderer->AddQuad(normalize(make_float3(0, 1.5, -1)), make_float3(0, -3, 9), 3.0f, 3.0f, lightMat);
#else
	// difficult light; use BDPT
	//int lightMat = renderer->AddMaterial( make_float3( 500, 500, 400 ) );
	//int lightQuad = renderer->AddQuad( make_float3( 0.15188693, -0.32204545, 0.93446094 ), make_float3( -12.938412, -5.0068984, -25.725601 ), 1.9f, 1.9f, lightMat );
#endif
	//int lightInst = renderer->AddInstance( lightQuad );

	// renderer->AddPointLight( make_float3( 20, 26, 20 ), make_float3( 1000, 1000, 1000 ) );
	// optional animated models
	// renderer->AddScene( "CesiumMan.glb", "data/", mat4::Translate( 0, -2, -9 ) );
	// renderer->AddScene( "project_polly.glb", "data/", mat4::Translate( 4.5f, -5.45f, -5.2f ) * mat4::Scale( 2 ) );
	// load changed materials
	renderer->DeserializeMaterials( materialFile.c_str() );
}

//  +-----------------------------------------------------------------------------+
//  |  HandleInput                                                                |
//  |  Process user input.                                                  LH2'19|
//  +-----------------------------------------------------------------------------+
bool HandleInput( float frameTime )
{
	// handle keyboard input
	float tspd = (keystates[GLFW_KEY_LEFT_SHIFT] ? 50.0f : 5.0f) * frameTime, rspd = 2.5f * frameTime;
	bool changed = false;
	Camera *camera = renderer->GetCamera();
	if (keystates[GLFW_KEY_A]) { changed = true; camera->TranslateRelative( make_float3( -tspd, 0, 0 ) ); }
	if (keystates[GLFW_KEY_D]) { changed = true; camera->TranslateRelative( make_float3( tspd, 0, 0 ) ); }
	if (keystates[GLFW_KEY_W]) { changed = true; camera->TranslateRelative( make_float3( 0, 0, tspd ) ); }
	if (keystates[GLFW_KEY_S]) { changed = true; camera->TranslateRelative( make_float3( 0, 0, -tspd ) ); }
	if (keystates[GLFW_KEY_R]) { changed = true; camera->TranslateRelative( make_float3( 0, tspd, 0 ) ); }
	if (keystates[GLFW_KEY_F]) { changed = true; camera->TranslateRelative( make_float3( 0, -tspd, 0 ) ); }
	if (keystates[GLFW_KEY_B]) changed = true; // force restart
	if (keystates[GLFW_KEY_UP]) { changed = true; camera->TranslateTarget( make_float3( 0, -rspd, 0 ) ); }
	if (keystates[GLFW_KEY_DOWN]) { changed = true; camera->TranslateTarget( make_float3( 0, rspd, 0 ) ); }
	if (keystates[GLFW_KEY_LEFT]) { changed = true; camera->TranslateTarget( make_float3( -rspd, 0, 0 ) ); }
	if (keystates[GLFW_KEY_RIGHT]) { changed = true; camera->TranslateTarget( make_float3( rspd, 0, 0 ) ); }
	if (!keystates[GLFW_KEY_SPACE]) spaceDown = false; else { if (!spaceDown) animPaused = !animPaused, changed = true; spaceDown = true; }
	// process left button click
	if (mbstates[GLFW_MOUSE_BUTTON_1] && keystates[GLFW_KEY_LEFT_SHIFT])
	{
		int selectedMaterialID = renderer->GetTriangleMaterialID( coreStats.probedInstid, coreStats.probedTriid );
		if (selectedMaterialID != -1)
		{
			currentMaterial = *renderer->GetMaterial( selectedMaterialID );
			currentMaterialID = selectedMaterialID;
			currentMaterial.Changed(); // update checksum so we can track changes
		}
		camera->focalDistance = coreStats.probedDist;
		changed = true;
	}
	// let the main loop know if the camera should update
	return changed;
}

//  +-----------------------------------------------------------------------------+
//  |  HandleMaterialChange                                                       |
//  |  Update a scene material based on AntTweakBar.                        LH2'19|
//  +-----------------------------------------------------------------------------+
bool HandleMaterialChange()
{
	if (currentMaterial.Changed() && currentMaterialID != -1)
	{
		// put it back
		*renderer->GetMaterial( currentMaterialID ) = currentMaterial;
		renderer->GetMaterial( currentMaterialID )->MarkAsDirty();
		return true;
	}
	return false;
}

//  +-----------------------------------------------------------------------------+
//  |  main                                                                       |
//  |  Application entry point.                                             LH2'19|
//  +-----------------------------------------------------------------------------+
int main()
{
	// initialize OpenGL and ImGui
	InitGLFW();
	InitImGui();

	// initialize renderer: pick one
	// renderer = RenderAPI::CreateRenderAPI( "RenderCore_Optix7filter" );		// OPTIX7 core, with filtering (static scenes only for now)
	// renderer = RenderAPI::CreateRenderAPI( "RenderCore_Optix7" );			// OPTIX7 core, best for RTX devices
	// renderer = RenderAPI::CreateRenderAPI( "RenderCore_OptixPrime_B" );			// OPTIX PRIME, best for pre-RTX CUDA devices
	// renderer = RenderAPI::CreateRenderAPI( "RenderCore_PrimeRef" );			// REFERENCE, for image validation
	 //renderer = RenderAPI::CreateRenderAPI( "RenderCore_SoftRasterizer" );	// RASTERIZER, your only option if not on NVidia
	 renderer = RenderAPI::CreateRenderAPI( "RenderCore_Minimal" );			// MINIMAL example, to get you started on your own core
	// renderer = RenderAPI::CreateRenderAPI( "RenderCore_Vulkan_RT" );			// Meir's Vulkan / RTX core
	// renderer = RenderAPI::CreateRenderAPI( "RenderCore_OptixPrime_BDPT" );	// Peter's OptixPrime / BDPT core
	// renderer = RenderAPI::CreateRenderAPI( "RenderCore_OptixPrime_PBRT" );	// Marijn's PBRT core
	// renderer = RenderAPI::CreateRenderAPI( "RenderCore_Optix7Guiding" );		// OPTIX7 core with path guiding for next event estimation (under construction)
	
	renderer->DeserializeCamera( "camera.xml" );
	// initialize scene
	PrepareScene();
	// set initial window size
	ReshapeWindowCallback( 0, SCRWIDTH, SCRHEIGHT );
	// enter main loop
	Timer timer;
	timer.reset();
	float deltaTime = 0;
	while (!glfwWindowShouldClose( window ))
	{
		// detect camera changes
		camMoved = renderer->GetCamera()->Changed();
		deltaTime = timer.elapsed();
		if (HandleInput( deltaTime )) camMoved = true;
		// handle material changes
		if (HandleMaterialChange()) camMoved = true;
		// poll events, may affect probepos so needs to happen between HandleInput and Render
		glfwPollEvents();
		// update animations
		if (!animPaused) for (int i = 0; i < renderer->AnimationCount(); i++)
		{
			renderer->UpdateAnimation( i, deltaTime );
			camMoved = true; // will remain false if scene has no animations
		}
		renderer->SynchronizeSceneData();
		// render
		timer.reset();
		renderer->Render( camMoved ? Restart : Converge );
		// postprocess
		shader->Bind();
		shader->SetInputTexture( 0, "color", renderTarget );
		shader->SetInputMatrix( "view", mat4::Identity() );
		shader->SetFloat( "contrast", renderer->GetCamera()->contrast );
		shader->SetFloat( "brightness", renderer->GetCamera()->brightness );
		shader->SetFloat( "gamma", renderer->GetCamera()->gamma );
		shader->SetInt( "method", renderer->GetCamera()->tonemapper );
		DrawQuad();
		shader->Unbind();
		// gui
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		ImGui::Begin( "Render statistics", 0 );
		coreStats = renderer->GetCoreStats();
		SystemStats systemStats = renderer->GetSystemStats();
		ImGui::Text( "Frame time:   %6.2fms", coreStats.renderTime * 1000 );
		ImGui::Text( "Scene update: %6.2fms", systemStats.sceneUpdateTime * 1000 );
		ImGui::Text( "Primary rays: %6.2fms", coreStats.traceTime0 * 1000 );
		ImGui::Text( "Secondary:    %6.2fms", coreStats.traceTime1 * 1000 );
		ImGui::Text( "Deep rays:    %6.2fms", coreStats.traceTimeX * 1000 );
		ImGui::Text( "Shadow rays:  %6.2fms", coreStats.shadowTraceTime * 1000 );
		ImGui::Text( "Shading time: %6.2fms", coreStats.shadeTime * 1000 );
		ImGui::Text( "Filter time:  %6.2fms", coreStats.filterTime * 1000 );
		ImGui::Text( "# primary:    %6ik (%6.1fM/s)", coreStats.primaryRayCount / 1000, coreStats.primaryRayCount / (max( 1.0f, coreStats.traceTime0 * 1000000 )) );
		ImGui::Text( "# secondary:  %6ik (%6.1fM/s)", coreStats.bounce1RayCount / 1000, coreStats.bounce1RayCount / (max( 1.0f, coreStats.traceTime1 * 1000000 )) );
		ImGui::Text( "# deep rays:  %6ik (%6.1fM/s)", coreStats.deepRayCount / 1000, coreStats.deepRayCount / (max( 1.0f, coreStats.traceTimeX * 1000000 )) );
		ImGui::Text( "# shadw rays: %6ik (%6.1fM/s)", coreStats.totalShadowRays / 1000, coreStats.totalShadowRays / (max( 1.0f, coreStats.shadowTraceTime * 1000000 )) );
		ImGui::End();
		ImGui::Begin( "Camera parameters", 0 );
		float3 camPos = renderer->GetCamera()->position;
		float3 camDir = renderer->GetCamera()->direction;
		ImGui::Text( "position: %5.2f, %5.2f, %5.2f", camPos.x, camPos.y, camPos.z );
		ImGui::Text( "viewdir:  %5.2f, %5.2f, %5.2f", camDir.x, camDir.y, camDir.z );
		ImGui::SliderFloat( "FOV", &renderer->GetCamera()->FOV, 10, 90 );
		ImGui::SliderFloat( "aperture", &renderer->GetCamera()->aperture, 0.001, 0.025f );
		ImGui::SliderFloat( "focaldistance", &renderer->GetCamera()->focalDistance, 0.1, 20.0f );
		ImGui::SliderFloat( "distortion", &renderer->GetCamera()->distortion, 0, 0.5f );
		ImGui::Combo( "tonemap", &renderer->GetCamera()->tonemapper, "clamp\0reinhard\0reinhard ext\0reinhard lum\0reinhard jodie\0uncharted2\0\0" );
		ImGui::SliderFloat( "brightness", &renderer->GetCamera()->brightness, 0, 0.5f );
		ImGui::SliderFloat( "contrast", &renderer->GetCamera()->contrast, 0, 0.5f );
		ImGui::SliderFloat( "gamma", &renderer->GetCamera()->gamma, 1, 2.5f );
		ImGui::End();
		ImGui::Begin( "Material parameters", 0 );
		ImGui::Text( "name:    %s", currentMaterial.name.c_str() );
		ImGui::ColorEdit3( "color", (float*)&currentMaterial.color() );
		ImGui::ColorEdit3( "absorption", (float*)&currentMaterial.absorption() );
		ImGui::SliderFloat( "metallic", &currentMaterial.metallic(), 0, 1 );
		ImGui::SliderFloat( "subsurface", &currentMaterial.subsurface(), 0, 1 );
		ImGui::SliderFloat( "specular", &currentMaterial.specular(), 0, 1 );
		ImGui::SliderFloat( "roughness", &currentMaterial.roughness(), 0, 1 );
		ImGui::SliderFloat( "specularTint", &currentMaterial.specularTint(), 0, 1 );
		ImGui::SliderFloat( "anisotropic", &currentMaterial.anisotropic(), 0, 1 );
		ImGui::SliderFloat( "sheen", &currentMaterial.sheen(), 0, 1 );
		ImGui::SliderFloat( "sheenTint", &currentMaterial.sheenTint(), 0, 1 );
		ImGui::SliderFloat( "clearcoat", &currentMaterial.clearcoat(), 0, 1 );
		ImGui::SliderFloat( "clearcoatGloss", &currentMaterial.clearcoatGloss(), 0, 1 );
		ImGui::SliderFloat( "transmission", &currentMaterial.transmission(), 0, 1 );
		ImGui::SliderFloat( "eta (1/ior)", &currentMaterial.eta(), 0.25f, 1.0f );
		ImGui::End();
		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData( ImGui::GetDrawData() );
		// finalize
		glfwSwapBuffers( window );
		// terminate
		if (!running) break;
	}
	// save material changes
	renderer->SerializeMaterials( materialFile.c_str() );
	// clean up
	renderer->SerializeCamera( "camera.xml" );
	renderer->Shutdown();
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
	glfwDestroyWindow( window );
	glfwTerminate();
	return 0;
}

// EOF