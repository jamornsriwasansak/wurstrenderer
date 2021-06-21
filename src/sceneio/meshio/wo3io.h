#pragma once

#include "common/wurst.h"

#include "geometry/trianglemesh.h"

struct Wo3InOut
{
	struct Wo3Vertex
	{
		float mPos[3];
		float mNormal[3];
		float mUv[2];
	};

	struct Wo3Triangle
	{
		uint32_t mIndices[3];
		uint32_t mSimdPadding; // >:( whyyyyyy
	};

	static shared_ptr<TriangleMeshGeometry> LoadFromWo3(const filesystem::path & objFilepath)
	{
		assert(objFilepath.is_absolute());
		uint64_t numVertices;
		uint64_t numTriangles;
		std::vector<Wo3Vertex> wo3Vertices;
		std::vector<Wo3Triangle> wo3Triangles;

		std::ifstream is(objFilepath, std::ios::binary | std::ios::in);
		assert(is.is_open());

		// read #vertices 
		{
			size_t extracted = is.read(reinterpret_cast<char *>(&numVertices), sizeof(uint64_t)).gcount();
			assert(extracted == sizeof(uint64_t));
		}

		// read wo3vertex
		{
			wo3Vertices.resize(numVertices);
			size_t extracted = is.read(reinterpret_cast<char *>(&wo3Vertices[0]), sizeof(Wo3Vertex) * numVertices).gcount();
			assert(extracted == sizeof(Wo3Vertex) * numVertices);
		}

		// read #triangles
		{
			size_t extracted = is.read(reinterpret_cast<char *>(&numTriangles), sizeof(uint64_t)).gcount();
			assert(extracted == sizeof(uint64_t));
		}

		// read wo3triangle
		{
			wo3Triangles.resize(numTriangles);
			size_t extracted = is.read(reinterpret_cast<char *>(&wo3Triangles[0]), sizeof(Wo3Triangle) * numTriangles).gcount();
			assert(extracted == sizeof(Wo3Triangle) * numTriangles);
		}

		// finish reading, close the file manually, this helps me prevent bugs.
		is.close();

		// organize it in triangle geometry
		shared_ptr<TriangleMeshGeometry> result = make_shared<TriangleMeshGeometry>();
		result->mPositions.reserve(numVertices);
		result->mNormals.reserve(numVertices);
		result->mTextureCoords.reserve(numVertices);
		for (const Wo3Vertex & vertex : wo3Vertices)
		{
			result->mPositions.emplace_back(Vec3(vertex.mPos[0], vertex.mPos[1], vertex.mPos[2]));
			result->mNormals.emplace_back(Vec3(vertex.mNormal[0], vertex.mNormal[1], vertex.mNormal[2]));
			result->mTextureCoords.emplace_back(Vec2(vertex.mUv[0], vertex.mUv[1]));
		}
		result->mTriangles.reserve(numTriangles);
		for (const Wo3Triangle & triangle : wo3Triangles)
		{
			result->mTriangles.emplace_back(Uvec3(triangle.mIndices[0], triangle.mIndices[1], triangle.mIndices[2]));
		}

		return result;
	}
};