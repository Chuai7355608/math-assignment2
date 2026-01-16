#pragma once
#define TEST_VIEW_PERSPECTIVE
#define USE_LH
#include <cmath>
#include <cstdint>

class mat4
{
public:
	vec4 m[4];   //create a vec4 with 4 index , a square

	mat4()    //initialize identity matrix
	{
		m[0] = vec4(1.0f, 0.0f, 0.0f, 0.0f);
		m[1] = vec4(0.0f, 1.0f, 0.0f, 0.0f);
		m[2] = vec4(0.0f, 0.0f, 1.0f, 0.0f);
		m[3] = vec4(0.0f, 0.0f, 0.0f, 1.0f);
	}


	//initialize matrix with specific vec4 rows
	mat4(vec4 const& r0, vec4 const& r1, vec4 const& r2, vec4 const& r3)
	{
		m[0] = vec4(r0.x, r1.x, r2.x, r3.x);
		m[1] = vec4(r0.y, r1.y, r2.y, r3.y);
		m[2] = vec4(r0.z, r1.z, r2.z, r3.z);
		m[3] = vec4(r0.w, r1.w, r2.w, r3.w);
	}

	mat4(mat4 const& m2) //make a copy of m2
	{
		m[0] = m2.m[0];
		m[1] = m2.m[1];
		m[2] = m2.m[2];
		m[3] = m2.m[3];
	}


	mat4& operator=(mat4 const& rhs)   //same with two vec above,just change values into rows
	{
		m[0] = rhs.m[0];
		m[1] = rhs.m[1];
		m[2] = rhs.m[2];
		m[3] = rhs.m[3];
		return *this;
	}

	//rhs as right hand mat, use equation of matrix product to solve this
	mat4 operator*(mat4 const& rhs) const
	{
		mat4 result;
		for (int row = 0; row < 4; row++)
		{
			for (int col = 0; col < 4; col++)
			{
				float num =
					m[row][0] * rhs.m[0][col] +
					m[row][1] * rhs.m[1][col] +
					m[row][2] * rhs.m[2][col] +
					m[row][3] * rhs.m[3][col];
				result.m[row][col] = num;
			}
		}
		return result;  //return a new mat4 
	}

	//use every row of mat to product with a 1*4 vec
	//the result will be a 1*4 vec in the end
	vec4 operator*(vec4 const& rhs) const
	{
		return vec4
		(
			dot(m[0], rhs),
			dot(m[1], rhs),
			dot(m[2], rhs),
			dot(m[3], rhs)
		);
	}


	bool operator==(mat4 const& rhs) const
	{
		for (int i = 0; i < 4; i++)
		{
			if (m[i] != rhs.m[i])   //check that if every row of vectors are equal,use function write above(vec4)
				return false;     //if now then return false    
		}
		return true;
	}

	bool operator!=(mat4 const& rhs) const
	{
		return !(*this == rhs);        //if mat not equal then return false
	}

	vec4& operator[](uint32_t const i) // range [0,3]
	{
		return m[i % 4];   //i Row index (0-3, modulo 4 for safety)
	}

	vec4 const& operator[](uint32_t const i) const // range [0,3]
	{
		return m[i % 4];   //i Row index (0-3, modulo 4 for safety)
	}
};

inline float determinant(mat4 const& m) // find the determinant for mat4
{
	//first find every single value in every position for later calculation
	const float& m00 = m[0][0], & m01 = m[0][1], & m02 = m[0][2], & m03 = m[0][3];
	const float& m10 = m[1][0], & m11 = m[1][1], & m12 = m[1][2], & m13 = m[1][3];
	const float& m20 = m[2][0], & m21 = m[2][1], & m22 = m[2][2], & m23 = m[2][3];
	const float& m30 = m[3][0], & m31 = m[3][1], & m32 = m[3][2], & m33 = m[3][3];

	//use lambda function here to calculate cofactor matrix,we also need to use it for inverse
	auto det3 = []
	(float a00, float a01, float a02,
		float a10, float a11, float a12,
		float a20, float a21, float a22) -> float
		{
			return a00 * a11 * a22 + a01 * a12 * a20 + a02 * a10 * a21  //cross product
				- a02 * a11 * a20 - a01 * a10 * a22 - a00 * a12 * a21;
		};

	//calculation here is based on equation for 4*4 matrix determinant
	//first row times cofactor matrix and the sign before is based on addition of row and col num
	float det = m00 * det3(m11, m12, m13, m21, m22, m23, m31, m32, m33)
		- m01 * det3(m10, m12, m13, m20, m22, m23, m30, m32, m33)
		+ m02 * det3(m10, m11, m13, m20, m21, m23, m30, m31, m33)
		- m03 * det3(m10, m11, m12, m20, m21, m22, m30, m31, m32);
	return det;
}

inline mat4 transpose(mat4 const& m)  //exchange rows and cols
{
	mat4 result;
	for (int row = 0; row < 4; row++)  //for loop every rows and cols
	{
		for (int col = 0; col < 4; col++)
		{
			result.m[row][col] = m.m[col][row];
		}
	}
	return result;
}

inline mat4 inverse(mat4 const& m)
{
	//use the equation for inverse matrix
	//inverse matraix = 1/determinant * adjugate matrix
	//adjugate matrix = transpose of cofactor matrix

	float det = determinant(m);

	auto det3 = []
	(float a00, float a01, float a02,
		float a10, float a11, float a12,
		float a20, float a21, float a22) -> float
		{
			return a00 * a11 * a22 + a01 * a12 * a20 + a02 * a10 * a21
				- a02 * a11 * a20 - a01 * a10 * a22 - a00 * a12 * a21;
		};

	const float& m00 = m[0][0], & m01 = m[0][1], & m02 = m[0][2], & m03 = m[0][3];  //find every single value in every position
	const float& m10 = m[1][0], & m11 = m[1][1], & m12 = m[1][2], & m13 = m[1][3];
	const float& m20 = m[2][0], & m21 = m[2][1], & m22 = m[2][2], & m23 = m[2][3];
	const float& m30 = m[3][0], & m31 = m[3][1], & m32 = m[3][2], & m33 = m[3][3];

	mat4 cofactor;   //create a function for cofactor matrix

	cofactor[0][0] = det3(m11, m12, m13, m21, m22, m23, m31, m32, m33);  // calculate every position in 4*4 cofactor matrix
	cofactor[0][1] = -det3(m10, m12, m13, m20, m22, m23, m30, m32, m33);
	cofactor[0][2] = det3(m10, m11, m13, m20, m21, m23, m30, m31, m33);
	cofactor[0][3] = -det3(m10, m11, m12, m20, m21, m22, m30, m31, m32);

	cofactor[1][0] = -det3(m01, m02, m03, m21, m22, m23, m31, m32, m33);
	cofactor[1][1] = det3(m00, m02, m03, m20, m22, m23, m30, m32, m33);
	cofactor[1][2] = -det3(m00, m01, m03, m20, m21, m23, m30, m31, m33);
	cofactor[1][3] = det3(m00, m01, m02, m20, m21, m22, m30, m31, m32);

	cofactor[2][0] = det3(m01, m02, m03, m11, m12, m13, m31, m32, m33);
	cofactor[2][1] = -det3(m00, m02, m03, m10, m12, m13, m30, m32, m33);
	cofactor[2][2] = det3(m00, m01, m03, m10, m11, m13, m30, m31, m33);
	cofactor[2][3] = -det3(m00, m01, m02, m10, m11, m12, m30, m31, m32);

	cofactor[3][0] = -det3(m01, m02, m03, m11, m12, m13, m21, m22, m23);
	cofactor[3][1] = det3(m00, m02, m03, m10, m12, m13, m20, m22, m23);
	cofactor[3][2] = -det3(m00, m01, m03, m10, m11, m13, m20, m21, m23);
	cofactor[3][3] = det3(m00, m01, m02, m10, m11, m12, m20, m21, m22);

	float inv_det = 1.0f / det;
	mat4 adj = transpose(cofactor);    //the adjugate matrix

	mat4 result;

	for (int row = 0; row < 4; row++)
	{
		result[row] = adj[row] * inv_det; // use the equation above for every row and assign it into final result
	}

	return result;
}


inline mat4 rotationx(float const rad) // radians
{
	float s = sin(rad);      //use the equation for rotationx
	float c = cos(rad);
	mat4 rx;
	rx.m[1][1] = c;
	rx.m[1][2] = -s;
	rx.m[2][1] = s;
	rx.m[2][2] = c;

	return rx;
}

inline mat4 rotationy(float const rad)
{
	float s = sin(rad);   //use the equation for rotationy
	float c = cos(rad);
	mat4 ry;

	ry.m[0][0] = c;
	ry.m[0][2] = s;
	ry.m[2][0] = -s;
	ry.m[2][2] = c;

	return ry;
}

inline mat4 rotationz(float const rad)
{
	float s = sin(rad);  // //use the equation for rotationz
	float c = cos(rad);
	mat4 rz;

	rz.m[0][0] = c;
	rz.m[0][1] = -s;
	rz.m[1][0] = s;
	rz.m[1][1] = c;

	return rz;
}


//use random vector as axis(two method)

inline mat4 rotationaxis(vec3 const& v, float const rad)
{
	//vec3 axis = normalize(v);
	//float nx = axis.x, ny = axis.y, nz = axis.z;

	//if (nx == 1 && ny == 0 && nz == 0)     //check if it is already on a axis
	//{
	//	return rotationx(rad);
	//}
	//if (nx == 0 && ny == 1 && nz == 0)
	//{
	//	return rotationy(rad);
	//}
	//if (nx == 0 && ny == 0 && nz == 1)
	//{
	//	return rotationz(rad);
	//}

	////if not then start rotate it to z axis and then calculate and then rotate back
	//float r = sqrt(nx * nx + nz * nz); //find the length that reflect to xz plane
	//float F1 = 0.f;
	//if (r > 0)  //if exist
	//{
	//	F1 = -atan2f(nx, nz);  //use negative sign to turn x value into 0
	//}
	//float F2 = 0.f;
	//if (length(axis) > 0)
	//{
	//	F2 = -atan2f(ny, r);   //use negative sign to turn x value into 0
	//}

	//mat4 rot_y = rotationy(F1);  //rotate y first, let the vector go into zy plane 
	//mat4 rot_x = rotationx(F2);  //then rotate x mkae the vector go to the same position with z axis
	//mat4 M = rot_x * rot_y;

	//mat4 inv = transpose(M);  //because of it s rotate equation, inverse equals to transpose 

	//mat4 rot_z = rotationz(rad); //rotate rad angle

	//mat4 result = inv * rot_z * M;
	//return result;
	vec3 x;
	vec3 y;
	vec3 z;
	x = vec3(1.0f, 0.0f, 0.0f);
	y = vec3(0.0f, 1.0f, 0.0f);
	z = vec3(0.0f, 0.0f, 1.0f);

	x = x* cosf(rad) + cross(v, x) * sinf(rad) + v * dot(v, x) * (1 - cosf(rad));
	y = y* cosf(rad) + cross(v, y) * sinf(rad) + v * dot(v, y) * (1 - cosf(rad));
	z = z* cosf(rad) + cross(v, z) * sinf(rad) + v * dot(v, z) * (1 - cosf(rad));

	return mat4(vec4(x.x, x.y, x.z, 0), vec4(y.x, y.y, y.z, 0), vec4(z.x, z.y, z.z, 0), vec4(0, 0, 0, 1));
}

inline mat4 perspective(float const fovy, float const aspect, float const near, float const far) // fovy expressed in radians
{
	mat4 proj;

	float half_fovy = fovy / 2;
	float inv_tan_half_fovy = 1.f / tan(half_fovy);
	float inv_aspect_tan = inv_tan_half_fovy / aspect;

	//just put every value into the right position then it is equation for perspective
	proj[0][0] = inv_aspect_tan;
	proj[1][1] = inv_tan_half_fovy;
	proj[2][2] = (far + near) / (near - far);
	proj[2][3] = (2 * far * near) / near - far;
	proj[3][2] = -1;
	proj[3][3] = 0;

	return proj;
}
inline mat4 lookat(vec3 const& eye, vec3 const& at, vec3 const& up)
{
	mat4 view;

	// Calculate camera axes
	vec3 zaxis = normalize(eye - at);
	vec3 xaxis = normalize(cross(up, zaxis));
	vec3 yaxis = normalize(cross(zaxis, xaxis));

	// Set rotation components (orientation)
	view[0][0] = xaxis.x;
	view[0][1] = xaxis.y;
	view[0][2] = xaxis.z;
	view[1][0] = yaxis.x;
	view[1][1] = yaxis.y;
	view[1][2] = yaxis.z;
	view[2][0] = zaxis.x;
	view[2][1] = zaxis.y;
	view[2][2] = zaxis.z;

	// Set translation components (position)
	view[0][3] = -xaxis.x * eye.x - xaxis.y * eye.y - xaxis.z * eye.z;
	view[1][3] = -yaxis.x * eye.x - yaxis.y * eye.y - yaxis.z * eye.z;
	view[2][3] = -zaxis.x * eye.x - zaxis.y * eye.y - zaxis.z * eye.z;

	return transpose(view);
}