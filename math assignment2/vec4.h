#pragma once
#define TEST_VIEW_PERSPECTIVE
#define USE_LH
#include <cmath>
#include <cstdint>


class vec4
{
public:

	float x;          // values for 4Dvec
	float y;
	float z;
	float w;

	vec4() :x(0.0f), y(0.0f), z(0.0f), w(0.0f) {};   //initialize zero length vector
	vec4(const float init_x, const float init_y, const float init_z, const float init_w) :x(init_x), y(init_y), z(init_z), w(init_w) {};
	//assign values

	vec4(vec4 const& v) :x(v.x), y(v.y), z(v.z), w(v.w) {};  //make a copy of vec4


	vec4& operator=(vec4 const& rhs)   //same with vec3 just one more value w
	{
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
		w = rhs.w;
		return *this;
	}
	vec4 operator-() const
	{
		return vec4(-x, -y, -z, -w);
	}

	vec4 operator+(vec4 const& rhs) const
	{
		return vec4(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w);
	}
	vec4& operator+=(vec4 const& rhs)
	{
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
		w += rhs.w;
		return *this;
	}
	vec4 operator-(vec4 const& rhs) const
	{
		return vec4(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w);
	}
	vec4& operator-=(vec4 const& rhs)
	{
		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
		w -= rhs.w;
		return *this;
	}
	vec4& operator*=(float const scalar)
	{
		x *= scalar;
		y *= scalar;
		z *= scalar;
		w *= scalar;
		return *this;
	}
	vec4 operator*(float const scalar) const
	{
		return vec4(x * scalar, y * scalar, z * scalar, w * scalar);
	}
	bool operator==(vec4 const& rhs) const
	{
		return (x == rhs.x) && (y == rhs.y) && (z == rhs.z) && (w == rhs.w);
	}
	bool operator!=(vec4 const& rhs) const
	{
		return !(*this == rhs);
	}
	float& operator[](uint32_t const i) // range [0,3]  //same with vec3 just one more value w
	{
		if (i == 1) return y;
		if (i == 0) return x;
		if (i == 2) return z;
		return w;
	}
	float const& operator[](uint32_t const i) const // range [0,3]
	{
		if (i == 1) return y;
		if (i == 0) return x;
		if (i == 2) return z;
		return w;
	}
};

inline float dot(vec4 const& a, vec4 const& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

inline float length(vec4 const& v)
{
	return std::sqrt(dot(v, v));
}

inline vec4 normalize(vec4 const& v)
{
	float len = length(v);
	return vec4(v.x / len, v.y / len, v.z / len, v.w / len);
}