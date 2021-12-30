#include "core_settings.h"

namespace lh2core {
	const float _third = 1.0f / 3.0f;

	void BVHNode::CalculateBounds(CoreTri* primitives, uint* indices, int firstPrim, int primCount) {
		bounds = AABB();
		bounds.pmin = make_float3(INFINITY, INFINITY, INFINITY);
		bounds.pmax = make_float3(-INFINITY, -INFINITY, -INFINITY);
		for (int i = 0; i < primCount; i++) {
			float3 v0 = primitives[indices[firstPrim + i]].vertex0;
			float3 v1 = primitives[indices[firstPrim + i]].vertex1;
			float3 v2 = primitives[indices[firstPrim + i]].vertex2;
			// min x
			bounds.pmin.x = min(bounds.pmin.x, v0.x);
			bounds.pmin.x = min(bounds.pmin.x, v1.x);
			bounds.pmin.x = min(bounds.pmin.x, v2.x);
			// min y
			bounds.pmin.y = min(bounds.pmin.y, v0.y);
			bounds.pmin.y = min(bounds.pmin.y, v1.y);
			bounds.pmin.y = min(bounds.pmin.y, v2.y);
			// min z
			bounds.pmin.z = min(bounds.pmin.z, v0.z);
			bounds.pmin.z = min(bounds.pmin.z, v1.z);
			bounds.pmin.z = min(bounds.pmin.z, v2.z);
			// max x
			bounds.pmax.x = max(bounds.pmax.x, v0.x);
			bounds.pmax.x = max(bounds.pmax.x, v1.x);
			bounds.pmax.x = max(bounds.pmax.x, v2.x);
			// max y
			bounds.pmax.y = max(bounds.pmax.y, v0.y);
			bounds.pmax.y = max(bounds.pmax.y, v1.y);
			bounds.pmax.y = max(bounds.pmax.y, v2.y);
			// max z
			bounds.pmax.z = max(bounds.pmax.z, v0.z);
			bounds.pmax.z = max(bounds.pmax.z, v1.z);
			bounds.pmax.z = max(bounds.pmax.z, v2.z);
		}
#if VERBOSEBVHCONSTRUCTION
		printf("Node: Computed bounds %.2f, %.2f, %.2f x %.2f, %.2f, %.2f\n",
			bounds.pMin.x, bounds.pMin.y, bounds.pMin.z, bounds.pMax.x, bounds.pMax.y, bounds.pMax.z);
#endif
	}

	bool BVHNode::RefitBounds(CoreTri* primitives, uint* indices, BVHNode** pool) {
		if (count > 0) {
			CalculateBounds(primitives, indices, leftFirst, count);
		} else {
			//bool leftcontained = pool[leftFirst]->bounds.IsContained(bounds);
			//bool rightcontained = pool[leftFirst + 1]->bounds.IsContained(bounds);
			//if (leftcontained && rightcontained) return false; // does parent need refitting?
			bounds = pool[leftFirst]->bounds.Merge(pool[leftFirst + 1]->bounds);
		}
		return true;
	}

	bool BVHNode::IntersectPrimitives(CoreTri* primitives, uint* indices, const Ray& ray, HitInfo* hitInfo) {
		CoreTri* intersectedTriangle = nullptr;
		float u, v;
		float t = INFINITY;
		float2 uv = make_float2(-1, -1);
		for (int i = 0; i < count; i++) {
			CoreTri* tri = &primitives[indices[leftFirst + i]];
			float ttemp = 0;
			if (ray.IntersectTriangle(tri->vertex0, tri->vertex1, tri->vertex2, ttemp, u, v)) {
				if (ttemp < t) {
					intersectedTriangle = tri;
					t = ttemp;
					uv.x = u;
					uv.y = v;
				}
			}
		}

		if (intersectedTriangle) {
			hitInfo->tHit = t;
			hitInfo->intersects = true;
			hitInfo->hitNormal = float3({ intersectedTriangle->Nx, intersectedTriangle->Ny, intersectedTriangle->Nz });
			hitInfo->intersection = ray.o + ray.d * t;
			hitInfo->triangle = intersectedTriangle;
			float w0 = 1 - (uv.x + uv.y);
			float w1 = uv.x;
			float w2 = 1 - w0 - w1;
			uv.x = intersectedTriangle->u0 * w0 + intersectedTriangle->u1 * w1 + intersectedTriangle->u2 * w2;
			uv.y = intersectedTriangle->v0 * w0 + intersectedTriangle->v1 * w1 + intersectedTriangle->v2 * w2;
			hitInfo->uv = uv;
			return true;
		}
		return false;
	}

	bool BVHNode::Intersect(CoreTri* primitives, uint* indices, BVHNode** pool, const Ray& r, HitInfo* hitInfo, uint& traverseCounter) {
		float t0, t1;
		if (!bounds.Intersect(r, t0, t1)) return false;
		traverseCounter++;
		if (count > 0) {
#if VERBOSEBVHTRAVERSAL
			printf("Node: Is a leaf.\n");
#endif
			return IntersectPrimitives(primitives, indices, r, hitInfo);
		}
		else {
			// determine traversal order (distance to children)
			BVHNode* nodeNear = pool[leftFirst];
			BVHNode* nodeFar = pool[leftFirst + 1];
			float nhitt0, nhitt1, fhitt0, fhitt1;
			nodeNear->bounds.Intersect(r, nhitt0, nhitt1);
			nodeFar->bounds.Intersect(r, fhitt0, fhitt1);
			if (fhitt0 < nhitt0) {
				std::swap(nodeNear, nodeFar);
				std::swap(nhitt0, fhitt0);
				std::swap(nhitt1, fhitt1);
			}

			// traverse
			HitInfo h1, h2;
			nodeNear->Intersect(primitives, indices, pool, r, &h1, traverseCounter);
			// is intersected primitive closer than the 'far' bounding box?
			if (h1.intersects && h1.tHit < fhitt0) {
				// then we should not intersect 'far' bounding box
				*hitInfo = h1;
				return true;
			}
			nodeFar->Intersect(primitives, indices, pool, r, &h2, traverseCounter);
			if (h1.intersects && h1.tHit < h2.tHit) {
				// intersection found in 'near' node
				*hitInfo = h1;
				return true;
			}
			else if (h2.intersects && h2.tHit < h1.tHit) {
				// intersection found in 'far' node
				*hitInfo = h2;
				return true;
			}
			else {
				// no intersection
				return false;
			}
		}
	}

	static uint nodeCount = 0;

	void BVHNode::Subdivide(CoreTri* primitives, float3* centroids, uint* indices,
		BVHNode** pool, uint& poolPtr, int firstPrim, int primCount, uint depth) {
#if VERBOSEBVHCONSTRUCTION
		printf("Node %d: Subdividing node with 'first' %d and 'count' %d...\n", depth, leftFirst, count);
#endif

		if (primCount < 3) {
#if VERBOSEBVHCONSTRUCTION
			printf("Node %d: Primitive count below threshold. Done.\n", depth);
#endif
			count = primCount;
			leftFirst = firstPrim;
			return;
		}
		leftFirst = poolPtr;
		BVHNode* leftNode = pool[poolPtr++];
		BVHNode* rightNode = pool[poolPtr++];
		uint splitIndex = Partition(primitives, centroids, indices, firstPrim, primCount, depth);
		//leftNode->count = splitIndex - firstPrim;
		//rightNode->count = count - leftNode->count;
		int leftCount = splitIndex - firstPrim;
		int rightCount = primCount - leftCount;
		leftNode->count = 0;
		rightNode->count = 0;
		leftNode->CalculateBounds(primitives, indices, firstPrim, leftCount);
		rightNode->CalculateBounds(primitives, indices, splitIndex, rightCount);
		leftNode->Subdivide(primitives, centroids, indices, pool, poolPtr, firstPrim, leftCount, depth + 1);
		rightNode->Subdivide(primitives, centroids, indices, pool, poolPtr, splitIndex, rightCount, depth + 1);
#if VERBOSEBVHCONSTRUCTION
		nodeCount++;
#endif
	}

	uint BVHNode::Partition(CoreTri* primitives, float3* centroids, uint* indices, int firstPrim, int primCount, uint depth) {
#if VERBOSEBVHCONSTRUCTION
		printf("Node %d: Partitioning...\n", depth);
#endif
		AABB centroidBounds = AABB();
		centroidBounds.pmin = make_float3(INFINITY, INFINITY, INFINITY);
		centroidBounds.pmax = make_float3(-INFINITY, -INFINITY, -INFINITY);
		for (int i = 0; i < primCount; i++) {
			centroidBounds.Expand(centroids[indices[firstPrim + i]]);
		}

		uint axis = centroidBounds.GetLongestAxis();
#if VERBOSEBVHCONSTRUCTION
		printf("Node %d: Split axis %d.\n", depth, axis);
#endif
		uint j = firstPrim - 1;
		float p;
		if (axis == 0) { // x-split
			p = (centroidBounds.pmin.x + centroidBounds.pmax.x) * 0.5f;
			for (int i = 0; i < primCount; i++) {
				CoreTri* tri = &primitives[indices[firstPrim + i]];
				if (centroids[indices[firstPrim + i]].x < p) {
					j++;
					swap(indices[j], indices[firstPrim + i]);
				}
			}
		}
		else if (axis == 1) { // y-split
			p = (centroidBounds.pmin.y + centroidBounds.pmax.y) * 0.5f;
			for (int i = 0; i < primCount; i++) {
				CoreTri* tri = &primitives[indices[firstPrim + i]];
				if (centroids[indices[firstPrim + i]].y < p) {
					j++;
					swap(indices[j], indices[firstPrim + i]);
				}
			}
		}
		else { // z-split
			p = (centroidBounds.pmin.z + centroidBounds.pmax.z) * 0.5f;
			for (int i = 0; i < primCount; i++) {
				CoreTri* tri = &primitives[indices[firstPrim + i]];
				if (centroids[indices[firstPrim + i]].z < p) {
					j++;
					swap(indices[j], indices[firstPrim + i]);
				}
			}
		}
#if VERBOSEBVHCONSTRUCTION
		printf("Node %d: Got split index %d.\n", depth, j + 1);
		printf("Node %d: Indices [", depth);
		for (int i = 0; i < primCount; i++) {
			printf("%d,", indices[leftFirst + i]);
		}
		printf("]\n");
#endif
		return j + 1;
	}

	bool BVH::Intersect(const Ray& ray, HitInfo* hitInfo, uint& traverseCounter) {
		return root->Intersect(mesh->triangles, indices, pool, ray, hitInfo, traverseCounter);
	}

	void BVH::Update() {
		switch (updateType) {
		case BVHUpdateType::Refit:
			Refit();
			break;
		case BVHUpdateType::Rebuild:
			Build();
			break;
		}
	}

	void BVH::Build() {
		if (!initialized) {
			indices = new uint[mesh->triangleCount];
			pool = new BVHNode * [mesh->triangleCount * 2];
		}
		for (uint i = 0; i < mesh->triangleCount; i++) indices[i] = i;
		for (uint i = 0; i < mesh->triangleCount * 2; i++) pool[i] = new BVHNode();
		root = pool[0];
		poolPtr = 2;

		//subdivide root node
		//root->count = mesh->triangleCount;
		root->count = 0;
		root->CalculateBounds(mesh->triangles, indices, 0, mesh->triangleCount);
		root->Subdivide(mesh->triangles, mesh->centroids, indices, pool, poolPtr, 0, mesh->triangleCount, 0);
		initialized = true;
#if VERBOSEBVHCONSTRUCTION
		printf("BVH: total %d nodes.\n", nodeCount);
#endif
	}

	void TopBVH::Build(vector<MeshInstance*>& instances) {
		vector<TopBVHNode*> nodes;
		for (int i = 0; i < instances.size(); i++) {
			TopBVHNode* node = new TopBVHNode();
			node->isLeaf = true;
			node->instance = instances[i];
			node->bounds = node->instance->mesh->bvh->root->bounds.TransformedAABB(node->instance->transform);
			nodes.push_back(node);
		}
		// Note: 1. Find the two elements in the list for which the AABB has the smallest SA
		uint indexA, indexB, indexC;
		AABB taabb = AABB();
		indexA = 0;
		TopBVHNode* A = nodes[indexA];
		indexB = FindBestMatch(nodes, A, indexA, taabb);
		TopBVHNode* B = nodes[indexB];
		while (nodes.size() > 1) {
			indexC = FindBestMatch(nodes, B, indexB, taabb);
			TopBVHNode* C = nodes[indexC];

			if (indexA == indexC) {
				nodes.erase(find(nodes.begin(), nodes.end(), A));
				nodes.erase(find(nodes.begin(), nodes.end(), B));
				TopBVHNode* newnode = new TopBVHNode();
				newnode->left = A;
				newnode->right = B;
				newnode->bounds = taabb;
				newnode->isLeaf = false;
				nodes.push_back(newnode);
				indexA = nodes.size() - 1;
				A = newnode;
				indexB = FindBestMatch(nodes, A, indexA, taabb);
				B = nodes[indexB];
			} else {
				indexA = indexB;
				indexB = indexC;
				A = B;
				B = C;
			}
		}
		this->root = nodes[0];
	}

	// given an array of AABBs, find a node for which a parent AABB has the smallest SA
// returns the index of the match
	int TopBVH::FindBestMatch(vector<TopBVHNode*> nodeList, TopBVHNode* node, int index, AABB& matchAABB) {
		int matchIndex = 0;
		float matchArea = INFINITY;
		AABB aabb = AABB();
		for (int i = 0; i < nodeList.size(); i++) {
			if (i == index) continue;
			aabb = nodeList[i]->bounds.Merge(node->bounds);
			float area = aabb.ComputeArea();
			if (area < matchArea) { // this is a better match
				matchIndex = i;
				matchArea = area;
				matchAABB = aabb;
			}

		}
		return matchIndex;
	};

	bool TopBVHNode::Intersect(Ray& ray, HitInfo* hitInfo, uint& traverseCounter) {
		float t0, t1;
		if (!bounds.Intersect(ray, t0, t1)) return false;
		traverseCounter++;
		if (isLeaf) {
			ray.o = make_float3(instance->inverseTransform * make_float4(ray.o, 1.0));
			ray.d = make_float3(instance->inverseDirectionTransform * make_float4(ray.d, 0.0));
			bool result = instance->mesh->bvh->Intersect(ray, hitInfo, traverseCounter);
			hitInfo->hitNormal = make_float3(instance->directionTransform * make_float4(hitInfo->hitNormal, 0.0));
			hitInfo->intersection = make_float3(instance->transform * make_float4(hitInfo->intersection, 1.0));
			ray.o = make_float3(instance->transform * make_float4(ray.o, 1.0));
			ray.d = make_float3(instance->directionTransform * make_float4(ray.d, 0.0));
			return result;
		} else {
			float nhitt0, nhitt1, fhitt0, fhitt1;
			TopBVHNode* nodeNear = left;
			TopBVHNode* nodeFar = right;
			nodeNear->bounds.Intersect(ray, nhitt0, nhitt1);
			nodeFar->bounds.Intersect(ray, fhitt0, fhitt1);
			if (fhitt0 < nhitt0) {
				std::swap(nodeNear, nodeFar);
				std::swap(nhitt0, fhitt0);
				std::swap(nhitt1, fhitt1);
			}
			HitInfo h1, h2;
			nodeNear->Intersect(ray, &h1, traverseCounter);
			// is intersected primitive closer than the 'far' bounding box?
			if (h1.intersects && h1.tHit < fhitt0) {
				// then we should not intersect 'far' bounding box
				*hitInfo = h1;
				return true;
			}
			nodeFar->Intersect(ray, &h2, traverseCounter);
			if (h1.intersects && h1.tHit < h2.tHit) {
				// intersection found in 'near' node
				*hitInfo = h1;
				return true;
			} else if (h2.intersects && h2.tHit < h1.tHit) {
				// intersection found in 'far' node
				*hitInfo = h2;
				return true;
			} else {
				// no intersection
				return false;
			}
		}
	}

	bool TopBVH::Intersect(Ray& ray, HitInfo* hitInfo, uint& traverseCounter) {
		return root->Intersect(ray, hitInfo, traverseCounter);
	}

	void BVH::Refit() {
		for (int i = poolPtr - 1; i >= 2; i--) {
			pool[i]->RefitBounds(mesh->triangles, indices, pool);
		}
		pool[0]->RefitBounds(mesh->triangles, indices, pool);
	}

	bool BVH::Traverse(RayPacket& rayPacket, Frustum& frustum, HitInfo* hitInfo, uint& traverseCounter) {
		bool intersection = false;
		BVHNode* curNode = root;
		vector<StackNode> traversalStack;
		int rayIndices[PACKETSIZE];
		for (int i = 0; i < PACKETSIZE; i++) rayIndices[i] = i;
		int lastActiveIndex = PACKETSIZE;

		float2 negfloat = make_float2(-1, -1);
		CoreTri* intersectedTriangle[PACKETSIZE];
		float u[PACKETSIZE], v[PACKETSIZE], t[PACKETSIZE], ttemp[PACKETSIZE];
		float2 uv[PACKETSIZE];
		for (int i = 0; i < PACKETSIZE; i++) {
			t[i] = INFINITY;
			intersectedTriangle[i] = nullptr;
			uv[i] = negfloat;
		}

		while (true) {
			lastActiveIndex = rayPacket.PartRays(frustum, curNode->bounds, rayIndices, lastActiveIndex);
			if (lastActiveIndex > 0) {
				if (curNode->count == 0) { // is innner node
					StackNode stackNode = StackNode();

					float nhitt0, nhitt1, fhitt0, fhitt1;
					BVHNode* nodeNear = pool[curNode->leftFirst];
					BVHNode* nodeFar = pool[curNode->leftFirst + 1];
					nodeNear->bounds.Intersect(rayPacket.rays[0], nhitt0, nhitt1);
					nodeFar->bounds.Intersect(rayPacket.rays[0], fhitt0, fhitt1);
					if (fhitt0 < nhitt0) {
						std::swap(nodeNear, nodeFar);
						std::swap(nhitt0, fhitt0);
						std::swap(nhitt1, fhitt1);
					}

					stackNode.node = nodeFar;
					stackNode.index = lastActiveIndex;
					traversalStack.push_back(stackNode);
					curNode = nodeNear;
					continue;
				} else {
					// indexe
					for (int i = 0; i < curNode->count; i++) {
						CoreTri* tri = &mesh->triangles[indices[curNode->leftFirst + i]];
						if (frustum.Intersect(*tri)) {
							memset(ttemp, 0, sizeof(float) * PACKETSIZE);
							for (int j = 0; j < lastActiveIndex; j++) {
								int currIdx = rayIndices[j];

								if (rayPacket.rays[currIdx].IntersectTriangle(tri->vertex0, tri->vertex1, tri->vertex2, t[currIdx], u[currIdx], v[currIdx])) {
									if (ttemp[currIdx] < t[currIdx]) {
										intersectedTriangle[currIdx] = tri;
										t[currIdx] = ttemp[currIdx];
										uv[currIdx].x = u[currIdx];
										uv[currIdx].y = v[currIdx];
									}
								}
							}
						}
					}

					for (int j = 0; j < lastActiveIndex; j++) {
						int currIdx = rayIndices[j];
						if (intersectedTriangle[currIdx]) {
							HitInfo* currHitInfo = hitInfo + currIdx;
							currHitInfo->tHit = t[currIdx];
							currHitInfo->intersects = true;
							currHitInfo->hitNormal = float3({
								intersectedTriangle[currIdx]->Nx,
								intersectedTriangle[currIdx]->Ny,
								intersectedTriangle[currIdx]->Nz });
							currHitInfo->intersection = rayPacket.rays[currIdx].o + rayPacket.rays[currIdx].d * t[currIdx];
							currHitInfo->triangle = intersectedTriangle[currIdx];
							float w0 = 1 - (uv[currIdx].x + uv[currIdx].y);
							float w1 = uv[currIdx].x;
							float w2 = 1 - w0 - w1;
							uv[currIdx].x = intersectedTriangle[currIdx]->u0 * w0 + intersectedTriangle[currIdx]->u1 * w1 + intersectedTriangle[currIdx]->u2 * w2;
							uv[currIdx].y = intersectedTriangle[currIdx]->v0 * w0 + intersectedTriangle[currIdx]->v1 * w1 + intersectedTriangle[currIdx]->v2 * w2;
							currHitInfo->uv = uv[currIdx];
							intersection = true;
						}
					}
				}
			}
			if (traversalStack.size() == 0) break;
			StackNode stackNode = traversalStack.back();
			traversalStack.pop_back();
			curNode = stackNode.node;
			lastActiveIndex = stackNode.index;
		}

		return intersection;
	}
}