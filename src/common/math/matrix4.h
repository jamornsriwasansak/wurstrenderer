#pragma once

#include <array>
#include "vecmath.h"

struct Matrix4 
{
	friend std::ostream & operator<<(std::ostream & out, const Matrix4 & v)
	{
		out << "[";
		for (size_t i = 0; i < 3; i++) { out << v.cols[i] << ", "; }
		out << v.cols[3] << "]";
		return out;
	}

	static Matrix4 Transpose(const Matrix4 & m)
	{
		Matrix4 result;
		for (Uint r = 0; r < 4; r++) for (Uint c = 0; c < 4; c++) { result.cols[r][c] = m.cols[c][r]; }
		return result;
	}

	static Matrix4 Scale(const Vec3 & v)
	{
		Matrix4 result;
		for (Uint i = 0; i < 3; i++) result.d[i][i] = v[i];
		return result;
	}

	static Matrix4 Translate(const Vec3 & v)
	{
		Matrix4 result;
		for (Uint i = 0; i < 3; i++) result.cols[3][i] = v[i];
		return result;
	}

	static Matrix4 Zero()
	{
		Matrix4 result;
		for (int i = 0; i < 16; i++) { result.c[i] = 0.0; }
		return result;
	}

	static Matrix4 FromQuaternion(const double w, const double x, const double y, const double z)
	{
		const double xx = x * x;
		const double xy = x * y;
		const double xz = x * z;
		const double xw = x * w;

		const double yy = y * y;
		const double yz = y * z;
		const double yw = y * w;

		const double zz = z * z;
		const double zw = z * w;

		assert(Math::IsApprox(xx + yy + zz + w*w, 1.0));

		return Matrix4(1.0 - 2.0 * (yy + zz),	2.0 * (xy - zw),		2.0 * (xz + yw),		0.0,
					   2.0 * (xy + zw),			1.0 - 2.0 * (xx + zz),	2.0 * (yz - xw),		0.0,
					   2.0 * (xz - yw),			2.0 * (yz + xw),		1.0 - 2.0 * (xx + yy),	0.0,
					   0.0,						0.0,					0.0,					1.0);
	}

	// an answer from https://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix
	static Matrix4 Inverse(const Matrix4 & m)
	{
		// TODO:: SIMD!! OTL
		double A2323 = m.d[2][2] * m.d[3][3] - m.d[3][2] * m.d[2][3];
		double A1323 = m.d[1][2] * m.d[3][3] - m.d[3][2] * m.d[1][3];
		double A1223 = m.d[1][2] * m.d[2][3] - m.d[2][2] * m.d[1][3];
		double A0323 = m.d[0][2] * m.d[3][3] - m.d[3][2] * m.d[0][3];
		double A0223 = m.d[0][2] * m.d[2][3] - m.d[2][2] * m.d[0][3];
		double A0123 = m.d[0][2] * m.d[1][3] - m.d[1][2] * m.d[0][3];
		double A2313 = m.d[2][1] * m.d[3][3] - m.d[3][1] * m.d[2][3];
		double A1313 = m.d[1][1] * m.d[3][3] - m.d[3][1] * m.d[1][3];
		double A1213 = m.d[1][1] * m.d[2][3] - m.d[2][1] * m.d[1][3];
		double A2312 = m.d[2][1] * m.d[3][2] - m.d[3][1] * m.d[2][2];
		double A1312 = m.d[1][1] * m.d[3][2] - m.d[3][1] * m.d[1][2];
		double A1212 = m.d[1][1] * m.d[2][2] - m.d[2][1] * m.d[1][2];
		double A0313 = m.d[0][1] * m.d[3][3] - m.d[3][1] * m.d[0][3];
		double A0213 = m.d[0][1] * m.d[2][3] - m.d[2][1] * m.d[0][3];
		double A0312 = m.d[0][1] * m.d[3][2] - m.d[3][1] * m.d[0][2];
		double A0212 = m.d[0][1] * m.d[2][2] - m.d[2][1] * m.d[0][2];
		double A0113 = m.d[0][1] * m.d[1][3] - m.d[1][1] * m.d[0][3];
		double A0112 = m.d[0][1] * m.d[1][2] - m.d[1][1] * m.d[0][2];

		double det = m.d[0][0] * (m.d[1][1] * A2323 - m.d[2][1] * A1323 + m.d[3][1] * A1223)
			- m.d[1][0] * (m.d[0][1] * A2323 - m.d[2][1] * A0323 + m.d[3][1] * A0223)
			+ m.d[2][0] * (m.d[0][1] * A1323 - m.d[1][1] * A0323 + m.d[3][1] * A0123)
			- m.d[3][0] * (m.d[0][1] * A1223 - m.d[1][1] * A0223 + m.d[2][1] * A0123);
		det = 1 / det;

		return Matrix4(
		   det * (m.d[1][1] * A2323 - m.d[2][1] * A1323 + m.d[3][1] * A1223), // m00
		   det * -(m.d[1][0] * A2323 - m.d[2][0] * A1323 + m.d[3][0] * A1223), // m01
		   det * (m.d[1][0] * A2313 - m.d[2][0] * A1313 + m.d[3][0] * A1213), // m02
		   det * -(m.d[1][0] * A2312 - m.d[2][0] * A1312 + m.d[3][0] * A1212), // m03
		   det * -(m.d[0][1] * A2323 - m.d[2][1] * A0323 + m.d[3][1] * A0223), // m10
		   det * (m.d[0][0] * A2323 - m.d[2][0] * A0323 + m.d[3][0] * A0223), // m11
		   det * -(m.d[0][0] * A2313 - m.d[2][0] * A0313 + m.d[3][0] * A0213), // m12
		   det * (m.d[0][0] * A2312 - m.d[2][0] * A0312 + m.d[3][0] * A0212), // m13
		   det * (m.d[0][1] * A1323 - m.d[1][1] * A0323 + m.d[3][1] * A0123), // m20
		   det * -(m.d[0][0] * A1323 - m.d[1][0] * A0323 + m.d[3][0] * A0123), // m21
		   det * (m.d[0][0] * A1313 - m.d[1][0] * A0313 + m.d[3][0] * A0113), // m22
		   det * -(m.d[0][0] * A1312 - m.d[1][0] * A0312 + m.d[3][0] * A0112), // m23
		   det * -(m.d[0][1] * A1223 - m.d[1][1] * A0223 + m.d[2][1] * A0123), // m30
		   det * (m.d[0][0] * A1223 - m.d[1][0] * A0223 + m.d[2][0] * A0123), // m31
		   det * -(m.d[0][0] * A1213 - m.d[1][0] * A0213 + m.d[2][0] * A0113), // m32
		   det * (m.d[0][0] * A1212 - m.d[1][0] * A0212 + m.d[2][0] * A0112) // m33
		);
	}

	// copied and transpose from tungsten renderer https://github.com/tunabrain/tungsten/
	// necessary for reading tungsten file format
	static Matrix4 RotateYXZ(const Vec3 & v)
	{
		double c[] = { std::cos(v[0]), std::cos(v[1]), std::cos(v[2]) };
		double s[] = { std::sin(v[0]), std::sin(v[1]), std::sin(v[2]) };

		return Matrix4(
			c[1] * c[2] - s[1] * s[0] * s[2], // m00
			-c[1] * s[2] - s[1] * s[0] * c[2], // m01
			-s[1] * c[0], // m02
			0.0f, // m03
			c[0] * s[2], // m10
			c[0] * c[2], // m11
			-s[0], // m12
			0.0f, // m13
			s[1] * c[2] + c[1] * s[0] * s[2], // m20
			-s[1] * s[2] + c[1] * s[0] * c[2], // m21
			c[1] * c[0], // m22
			0.0f, // m23
			0.0f, // m30
			0.0f, // m31
			0.0f, // m32
			1.0f // m33
		);
	}

	// create identity matrix
	Matrix4()
	{
		for (Uint i = 0; i < 16; i++)
			c[i] = (i % 5 == 0) ? 1.0 : 0.0;
	}

	Matrix4(const Vec4 & col0, const Vec4 & col1, const Vec4 & col2, const Vec4 & col3): cols{col0, col1, col2, col3} {}

	// by column constructor
	Matrix4(const double s00, const double s01, const double s02, const double s03,
			const double s10, const double s11, const double s12, const double s13,
			const double s20, const double s21, const double s22, const double s23,
			const double s30, const double s31, const double s32, const double s33): 
		c{ s00, s10, s20, s30, s01, s11, s21, s31, s02, s12, s22, s32, s03, s13, s23, s33 }
	{}

	Matrix4(const Matrix4 & m)
	{
		for (Uint i = 0; i < 4; i++)
			cols[i] = m.cols[i];
	}

	Matrix4 operator*(const Matrix4 & m) const
	{
		Matrix4 result;
		// TODO:: simd
		for (Uint col = 0; col < 4; col++)
			for (Uint row = 0; row < 4; row++)
			{
				result.d[col][row] = 0;
				for (Uint k = 0; k < 4; k++) result.d[col][row] += d[k][row] * m.d[col][k];
			}
		return result;
	}

	Vec3 transformPoint(const Vec3 & v) const
	{
		Vec4 hResult;
		for (Uint i = 0; i < 3; i++) { hResult += v[i] * cols[i]; }
		hResult += cols[3];
		return Vec3(hResult[0], hResult[1], hResult[2]);
	}

	Vec3 transformDirection(const Vec3 & v) const
	{
		Vec4 hResult;
		for (Uint i = 0; i < 3; i++) { hResult += v[i] * cols[i]; }
		return Vec3(hResult[0], hResult[1], hResult[2]);
	}

	Vec4 transform(const Vec4 & v) const
	{
		Vec4 hResult;
		for (Uint i = 0; i < 4; i++) { hResult += v[i] * cols[i]; }
		return hResult;
	}

	std::array<float, 16> getFloats() const 
	{
		std::array<float, 16> result;
		for (Uint i = 0; i < 16; i++) { result[i] = static_cast<float>(c[i]); }
		return result;
	}

	union
	{
		Vec4 cols[4];
		double c[16];
		double d[4][4];
	};
};