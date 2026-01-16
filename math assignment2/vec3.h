#pragma once
#define TEST_VIEW_PERSPECTIVE
#define USE_LH
#include <cmath>
#include <cstdint>

class vec3              //class for 3D vector
{
public:
	float x;
	float y;
	float z;

	vec3() :x(0.0f), y(0.0f), z(0.0f) {};     //initialize zero length vector

	vec3(const float init_x, const float init_y, const float init_z) :x(init_x), y(init_y), z(init_z) {};
	//initialize with specific coordinates

	vec3(vec3 const& v) :x(v.x), y(v.y), z(v.z) {};  //make a copy of another vector

	vec3& operator=(vec3 const& rhs)      //assign coordinates from another vector
	{
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
		return *this;        //return a pointer to self
	}

	vec3 operator-() const         //negate all vector values
	{
		return vec3(-x, -y, -z);      //create a new vec3 and then return the value
	}

	vec3 operator+(vec3 const& rhs) const       //addition vectors
	{
		return vec3(x + rhs.x, y + rhs.y, z + rhs.z);       //return a new vec3 with this + target vector
	}

	vec3& operator+=(vec3 const& rhs)      // add and then assign
	{
		x += rhs.x;            //every value add target value
		y += rhs.y;
		z += rhs.z;
		return *this;          //return pointer to self
	}

	vec3 operator-(vec3 const& rhs) const     //subtraction vectors
	{
		return vec3(x - rhs.x, y - rhs.y, z - rhs.z);     //return a new vec3 with subtrection values
	}

	vec3& operator-=(vec3 const& rhs)       //subtract and then assgin
	{
		x -= rhs.x;             //every value subtract target value
		y -= rhs.y;
		z -= rhs.z;
		return *this;          //return a pointer to self
	}

	vec3& operator*=(float const scalar)     // multiply by scalar and then assign
	{
		x *= scalar;             //every value multiply by scalar
		y *= scalar;
		z *= scalar;
		return *this;           //return a pointer to self
	}

	vec3 operator*(float const scalar) const           // multipllication vectors
	{
		return vec3(x * scalar, y * scalar, z * scalar);     //return a new vec with multiplied numbers
	}

	bool operator==(vec3 const& rhs) const   // check if two vectors equals
	{
		return (x == rhs.x) && (y == rhs.y) && (z == rhs.z);    //check every values are the same with target values
	}

	bool operator!=(vec3 const& rhs) const    // check if two vectors are not equal
	{
		return !(*this == rhs);   //use == write above to check if it s not equal
	}

	float& operator[](uint32_t const i)// range [0,2]
	{
		if (i == 1) return y;         //eference to the component at specified index
		if (i == 0) return x;
		return z;
	}

	float const& operator[](uint32_t const i) const// range [0,2]
	{
		if (i == 1) return y;          //const reference to the component at specified index
		if (i == 0) return x;
		return z;
	}
};

inline float dot(vec3 const& a, vec3 const& b)     //dot production
{
	return a.x * b.x + a.y * b.y + a.z * b.z;       //return the number of addition of every same position product
}

inline float length(vec3 const& v)      //length
{
	return std::sqrt(dot(v, v));       //use the equation to get the length of tow dot or vector
}

inline vec3 cross(vec3 const& a, vec3 const& b)     //cross production
{
	return vec3(a.y * b.z - b.y * a.z, b.x * a.z - a.x * b.z, a.x * b.y - b.x * a.y);    //use the equation of cross product
}

inline vec3 normalize(vec3 const& v)   //normalize vector
{
	float len = length(v);
	return vec3(v.x / len, v.y / len, v.z / len);   //return new vec3 with original value divide by length
}

