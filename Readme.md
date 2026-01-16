# Vector and Matrix Library

## vec3 Class

### 1. Member Variables
- `float x`: x-coordinate.
- `float y`: y-coordinate.
- `float z`: z-coordinate.

### 2. Constructors
- `vec3()`: Default constructor, initializes to the zero vector (0.0, 0.0, 0.0).
- `vec3(const float init_x, const float init_y, const float init_z)`: Parameterized constructor, initializes the vector with specified x, y, z values.
- `vec3(vec3 const& v)`: Copy constructor, creates a copy of the given vector.

### 3. Operator Overloads
- `vec3& operator=(vec3 const& rhs)`: Assignment operator, assigns the values from another vector to the current vector.
- `vec3 operator-() const`: Negation operator, returns a new vector with all components negated.
- `vec3 operator+(vec3 const& rhs) const`: Addition operator, returns the result of adding two vectors.
- `vec3& operator+=(vec3 const& rhs)`: Addition assignment operator, adds another vector to the current vector.
- `vec3 operator-(vec3 const& rhs) const`: Subtraction operator, returns the result of subtracting two vectors.
- `vec3& operator-=(vec3 const& rhs)`: Subtraction assignment operator, subtracts another vector from the current vector.
- `vec3& operator*=(float const scalar)`: Scalar multiplication assignment operator, multiplies each component of the current vector by a scalar.
- `vec3 operator*(float const scalar) const`: Scalar multiplication operator, returns the result of multiplying a vector by a scalar.
- `bool operator==(vec3 const& rhs) const`: Equality operator, checks if two vectors are equal.
- `bool operator!=(vec3 const& rhs) const`: Inequality operator, checks if two vectors are not equal.
- `float& operator[](uint32_t const i)`: Index operator, returns a reference to the component at the specified index (read/write).
- `float const& operator[](uint32_t const i) const`: Index operator, returns a const reference to the component at the specified index (read-only).

### 4. External Functions
- `float dot(vec3 const& a, vec3 const& b)`: Dot product function, calculates the dot product of two vectors.
- `float length(vec3 const& v)`: Length function, calculates the length (magnitude) of a vector.
- `vec3 cross(vec3 const& a, vec3 const& b)`: Cross product function, calculates the cross product of two vectors.
- `vec3 normalize(vec3 const& v)`: Normalize function, returns a vector with the same direction as the given vector but with a length of 1.

## vec4 Class

### 1. Member Variables
- `float x`: x-coordinate.
- `float y`: y-coordinate.
- `float z`: z-coordinate.
- `float w`: w-coordinate.

### 2. Constructors
- `vec4()`: Default constructor, initializes to the zero vector (0.0, 0.0, 0.0, 0.0).
- `vec4(const float init_x, const float init_y, const float init_z, const float init_w)`: Parameterized constructor, initializes the vector with specified x, y, z, w values.
- `vec4(vec4 const& v)`: Copy constructor, creates a copy of the given vector.

### 3. Operator Overloads
- `vec4& operator=(vec4 const& rhs)`: Assignment operator, assigns the values from another vector to the current vector.
- `vec4 operator-() const`: Negation operator, returns a new vector with all components negated.
- `vec4 operator+(vec4 const& rhs) const`: Addition operator, returns the result of adding two vectors.
- `vec4& operator+=(vec4 const& rhs)`: Addition assignment operator, adds another vector to the current vector.
- `vec4 operator-(vec4 const& rhs) const`: Subtraction operator, returns the result of subtracting two vectors.
- `vec4& operator-=(vec4 const& rhs)`: Subtraction assignment operator, subtracts another vector from the current vector.
- `vec4& operator*=(float const scalar)`: Scalar multiplication assignment operator, multiplies each component of the current vector by a scalar.
- `vec4 operator*(float const scalar) const`: Scalar multiplication operator, returns the result of multiplying a vector by a scalar.
- `bool operator==(vec4 const& rhs) const`: Equality operator, checks if two vectors are equal.
- `bool operator!=(vec4 const& rhs) const`: Inequality operator, checks if two vectors are not equal.
- `float& operator[](uint32_t const i)`: Index operator, returns a reference to the component at the specified index (read/write).
- `float const& operator[](uint32_t const i) const`: Index operator, returns a const reference to the component at the specified index (read-only).

### 4. External Functions
- `float dot(vec4 const& a, vec4 const& b)`: Dot product function, calculates the dot product of two vectors.
- `float length(vec4 const& v)`: Length function, calculates the length (magnitude) of a vector.
- `vec4 normalize(vec4 const& v)`: Normalize function, returns a vector with the same direction as the given vector but with a length of 1.

## mat4 Class

### 1. Member Variables
- `vec4 m[4]`: An array containing four `vec4` objects, representing the four rows of the matrix.

### 2. Constructors
- `mat4()`: Default constructor, initializes to the identity matrix.
- `mat4(vec4 const& r0, vec4 const& r1, vec4 const& r2, vec4 const& r3)`: Parameterized constructor, initializes the matrix with the specified four `vec4` objects as rows.
- `mat4(mat4 const& m2)`: Copy constructor, creates a copy of the given matrix.

### 3. Operator Overloads
- `mat4& operator=(mat4 const& rhs)`: Assignment operator, assigns the values from another matrix to the current matrix.
- `mat4 operator*(mat4 const& rhs) const`: Matrix multiplication operator, returns the result of multiplying two matrices.
- `vec4 operator*(vec4 const& rhs) const`: Matrix-vector multiplication operator, returns the result of multiplying a matrix by a vector.
- `bool operator==(mat4 const& rhs) const`: Equality operator, checks if two matrices are equal.
- `bool operator!=(mat4 const& rhs) const`: Inequality operator, checks if two matrices are not equal.
- `vec4& operator[](uint32_t const i)`: Index operator, returns a reference to the row at the specified index (read/write).
- `vec4 const& operator[](uint32_t const i) const`: Index operator, returns a const reference to the row at the specified index (read-only).

### 4. External Functions
- `float determinant(mat4 const& m)`: Determinant function, calculates the determinant of the matrix.
- `mat4 inverse(mat4 const& m)`: Inverse function, calculates the inverse of the matrix.
- `mat4 transpose(mat4 const& m)`: Transpose function, returns the transpose of the matrix.
- `mat4 rotationx(float const rad)`: Rotation matrix around the X-axis.
- `mat4 rotationy(float const rad)`: Rotation matrix around the Y-axis.
- `mat4 rotationz(float const rad)`: Rotation matrix around the Z-axis.
- `mat4 rotationaxis(vec3 const& v, float const rad)`: Rotation matrix around an arbitrary axis, `v` is the rotation axis vector.
- `mat4 perspective(float const fovy, float const aspect, float const near, float const far)`: Perspective projection matrix, `fovy` is the vertical field of view in radians, `aspect` is the aspect ratio, `near` is the near plane distance, `far` is the far plane distance.
- `mat4 lookat(vec3 const& eye, vec3 const& at, vec3 const& up)`: View matrix, `eye` is the camera position, `at` is the target position, `up` is the up direction.