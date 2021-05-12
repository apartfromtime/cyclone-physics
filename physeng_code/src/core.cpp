/*
 * Implementation file for core functions in the library.
 * 
 * Part of the Cyclone physics system.
 * 
 * Copyright (c) Icosagon 2003. All Rights Reserved.
 *
 * This software is distributed under licence. Use of this software
 * implies agreement with all terms and conditions of the accompanying
 * software licence.
 */


#include <cyclone/core.h>

using namespace cyclone;

/**
 * Holds the value for energy under which a body will be put to
 * sleep. This is a global value for the whole solution. By
 * default it is 0.3, which is fine for simulation when gravity is
 * about 20 units per second squared, masses are about one, and
 * other forces are around that of gravity. It may need tweaking
 * if your simulation is drastically different to this.
 */
real_t sleepEpsilon = ( ( real_t )0.3f );

/*
 * Functions to change sleepEpsilon.
 */
void cyclone::SetSleepEpsilon(real_t value)
{
    sleepEpsilon = value;
}

real_t cyclone::GetSleepEpsilon(void)
{
    return sleepEpsilon;
}


//-----------------------------------------------------------------------------
// VECTOR 3
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// adds two 3d vectors
//-----------------------------------------------------------------------------
vec3_t cyclone::Vec3Add(vec3_t a, vec3_t b)
{
    vec3_t v;

    v.x = ( a.x + b.x );
    v.y = ( a.y + b.y );
    v.z = ( a.z + b.z );

    return v;
}

//-----------------------------------------------------------------------------
// subtracts two 3d vectors
//-----------------------------------------------------------------------------
vec3_t cyclone::Vec3Subtract(vec3_t a, vec3_t b)
{
    vec3_t v;

    v.x = ( a.x - b.x );
    v.y = ( a.y - b.y );
    v.z = ( a.z - b.z );

    return v;
}

//-----------------------------------------------------------------------------
// calculates the component-wise product of two 3d vectors ( Multiply )
//-----------------------------------------------------------------------------
vec3_t cyclone::Vec3ComponentProduct(vec3_t a, vec3_t b)
{
    vec3_t v;

    v.x = a.x * b.x;
    v.y = a.y * b.y;
    v.z = a.z * b.z;

    return v;
}

//-----------------------------------------------------------------------------
// scales two 3d vectors
//-----------------------------------------------------------------------------
vec3_t cyclone::Vec3Scale(vec3_t a, real_t s)
{
    vec3_t v = { a.x * s, a.y * s, a.z * s };

    return v;
}

//-----------------------------------------------------------------------------
// Calculates the vector product of two 3d vectors ( CrossProduct )
//-----------------------------------------------------------------------------
vec3_t cyclone::Vec3VectorProduct(vec3_t a, vec3_t b)
{
    vec3_t v;

    v.x = a.y * b.z - a.z * b.y;
    v.y = a.z * b.x - a.x * b.z;
    v.z = a.x * b.y - a.y * b.x;

    return v;
}

//-----------------------------------------------------------------------------
// calculates the scalar product of two 3d vectors ( DotProduct )
//-----------------------------------------------------------------------------
real_t cyclone::Vec3ScalarProduct(vec3_t a, vec3_t b)
{
    return ( ( a.x * b.x ) + ( a.y * b.y ) + ( a.z * b.z ) );
}

//-----------------------------------------------------------------------------
// gets the magnitude of a vector
//-----------------------------------------------------------------------------
real_t cyclone::Vec3Magnitude(vec3_t v)
{
    return R_sqrt( ( v.x * v.x ) + ( v.y * v.y ) + ( v.z * v.z ) );
}

//-----------------------------------------------------------------------------
// gets the squared magnitude of a vector
//-----------------------------------------------------------------------------
real_t cyclone::Vec3MagnitudeSqr(vec3_t v)
{
    return ( ( v.x * v.x ) + ( v.y * v.y ) + ( v.z * v.z ) );   
}

//-----------------------------------------------------------------------------
// turns a non-zero vector into a vector of unit length
//-----------------------------------------------------------------------------
vec3_t cyclone::Vec3Normalise(vec3_t v)
{
    vec3_t a;
    real_t l;

    a = v;
    l = cyclone::Vec3Magnitude( a );

    if ( l > 0.0f ) {

        a.x *= 1.0f / l;
        a.y *= 1.0f / l;
        a.z *= 1.0f / l;
    }

    return a;
}

//-----------------------------------------------------------------------------
// zero all the components of the vector
//-----------------------------------------------------------------------------
vec3_t cyclone::Vec3Clear(void)
{
    vec3_t v = { 0.0f };
    return v;
}

//-----------------------------------------------------------------------------
// flips all the components of the vector
//-----------------------------------------------------------------------------
vec3_t cyclone::Vec3Invert(vec3_t v)
{
    vec3_t a;

    a.x = -v.x;
    a.y = -v.y;
    a.z = -v.z;

    return a;
}


//-----------------------------------------------------------------------------
// QUATERNION
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// zero all the components of the quaternion
//-----------------------------------------------------------------------------
quat_t cyclone::QuatClear(void)
{
    quat_t q = { 1.0f, 0.0f, 0.0f, 0.0f };
    return q;
}

//-----------------------------------------------------------------------------
// normalises the quaternion to unit length, making it a valid orientation
// quaternion
//-----------------------------------------------------------------------------
quat_t cyclone::QuatNormalise(quat_t q)
{
    quat_t a = { 0.0f };

    a = q;

    real_t d = ( q.r * q.r ) + ( q.i * q.i ) + ( q.j * q.j ) +
        ( q.k * q.k );

    /* Check for zero length quaternion, and use the no-rotation quaternion
    in that case */
    if ( d != 0.0f ) {
        
        d = 1.0f / R_sqrt( d );
        a.r *= d;
        a.i *= d;
        a.j *= d;
        a.k *= d;

    } else {
        a.r = 1.0f;
    }

    return a;
}

//-----------------------------------------------------------------------------
// multiplies two quaternions
//-----------------------------------------------------------------------------
quat_t cyclone::QuatMultiply(quat_t a, quat_t b)
{
    quat_t q = { 0.0f };

    q.r = a.r * b.r - a.i * b.i - a.j * b.j - a.k * b.k;
    q.i = a.r * b.i + a.i * b.r + a.j * b.k - a.k * b.j;
    q.j = a.r * b.j + a.j * b.r + a.k * b.i - a.i * b.k;
    q.k = a.r * b.k + a.k * b.r + a.i * b.j - a.j * b.i;

    return q;
}

//-----------------------------------------------------------------------------
// adds the given vector to quaternion, scaled by the given amount.
// this is used to update the orientation quaternion by a rotation and time
//-----------------------------------------------------------------------------
quat_t cyclone::QuatAddScaledVector(quat_t q, vec3_t v, real_t s)
{
    quat_t a = { 0.0f };
    quat_t b = { 0.0f };

    a.r = 0.0f;
    a.i = v.x * s;
    a.j = v.y * s;
    a.k = v.z * s;

    a = QuatMultiply( q, a );

    b = q;
    b.r += a.r * 0.5f;
    b.i += a.i * 0.5f;
    b.j += a.j * 0.5f;
    b.k += a.k * 0.5f;

    return b;
}

//-----------------------------------------------------------------------------
// rotates quaternion about given vector
//-----------------------------------------------------------------------------
quat_t cyclone::QuatRotateByVector(quat_t q, vec3_t v)
{
    quat_t a = { 0.0f };

    a.r = 0.0f;
    a.i = v.x;
    a.j = v.y;
    a.k = v.z;

    return QuatMultiply( q, a );
}


//-----------------------------------------------------------------------------
// MATRIX 4
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// sets the matrix to be a identity a matrix
//-----------------------------------------------------------------------------
mat4_t cyclone::Mat4Identity(void)
{
    mat4_t m = { 0.0f };

    m.n[ 1] = m.n[ 2] = m.n[ 3] = m.n[ 4] =
    m.n[ 6] = m.n[ 7] = m.n[ 8] = m.n[ 9] = m.n[11] = 0;
    m.n[ 0] = m.n[ 5] = m.n[10] = m.n[15] = 1;

    return m;
}

//-----------------------------------------------------------------------------
// sets the matrix to be a diagonal matrix with the given coefficients
//-----------------------------------------------------------------------------
mat4_t cyclone::Mat4SetDiagonal(real_t x, real_t y, real_t z)
{
    mat4_t m = { 0.0f };
    
    m.x0 = x;
    m.y1 = y;
    m.z2 = z;
    m.w3 = 1.0f;

    return m;
}

//-----------------------------------------------------------------------------
// returns a matrix which is this matrix multiplied by the given other matrix
//-----------------------------------------------------------------------------
mat4_t cyclone::Mat4Multiply(mat4_t a, mat4_t b)
{
    mat4_t m = { 0.0f };

    m.n[ 0] = ( b.n[ 0] * a.n[ 0] ) + ( b.n[ 4] * a.n[ 1] ) +
        ( b.n[ 8] * a.n[ 2] );
    m.n[ 4] = ( b.n[ 0] * a.n[ 4] ) + ( b.n[ 4] * a.n[ 5] ) +
        ( b.n[ 8] * a.n[ 6] );
    m.n[ 8] = ( b.n[ 0] * a.n[ 8] ) + ( b.n[ 4] * a.n[ 9] ) +
        ( b.n[ 8] * a.n[10] );

    m.n[ 1] = ( b.n[ 1] * a.n[ 0] ) + ( b.n[ 5] * a.n[ 1] ) +
        ( b.n[ 9] * a.n[ 2] );
    m.n[ 5] = ( b.n[ 1] * a.n[ 4] ) + ( b.n[ 5] * a.n[ 5] ) +
        ( b.n[ 9] * a.n[ 6] );
    m.n[ 9] = ( b.n[ 1] * a.n[ 8] ) + ( b.n[ 5] * a.n[ 9] ) +
        ( b.n[ 9] * a.n[10] );

    m.n[ 2] = ( b.n[ 2] * a.n[ 0] ) + ( b.n[ 6] * a.n[ 1] ) +
        ( b.n[10] * a.n[ 2] );
    m.n[ 6] = ( b.n[ 2] * a.n[ 4] ) + ( b.n[ 6] * a.n[ 5] ) +
        ( b.n[10] * a.n[ 6] );
    m.n[10] = ( b.n[ 2] * a.n[ 8] ) + ( b.n[ 6] * a.n[ 9] ) +
        ( b.n[10] * a.n[10] );

    m.n[ 3] = ( b.n[ 3] * a.n[ 0] ) + ( b.n[ 7] * a.n[ 1] ) +
        ( b.n[11] * a.n[ 2] ) + a.n[ 3];
    m.n[ 7] = ( b.n[ 3] * a.n[ 4] ) + ( b.n[ 7] * a.n[ 5] ) +
        ( b.n[11] * a.n[ 6] ) + a.n[ 7];
    m.n[11] = ( b.n[ 3] * a.n[ 8] ) + ( b.n[ 7] * a.n[ 9] ) +
        ( b.n[11] * a.n[10] ) + a.n[11];

    return m;
}

//-----------------------------------------------------------------------------
// transform the given vector by a matrix
//-----------------------------------------------------------------------------
vec3_t cyclone::Mat4Transform(vec3_t v, mat4_t m)
{
    vec3_t a = { 0.0f };

    a.x = v.x * m.n[ 0] + v.y * m.n[ 1] + v.z * m.n[ 2] + m.n[ 3];
    a.y = v.x * m.n[ 4] + v.y * m.n[ 5] + v.z * m.n[ 6] + m.n[ 7];
    a.z = v.x * m.n[ 8] + v.y * m.n[ 9] + v.z * m.n[10] + m.n[11];

    return a;
}

//-----------------------------------------------------------------------------
// transform the given vector by the transformational inverse of this matrix 
//-----------------------------------------------------------------------------
vec3_t cyclone::Mat4TransformInverse(vec3_t v, mat4_t m)
{
    vec3_t a = { 0.0f };
    vec3_t b = { 0.0f };

    b = v;

    b.x -= m.n[ 3];
    b.y -= m.n[ 7];
    b.z -= m.n[11];

    a.x = b.x * m.n[ 0] + b.y * m.n[ 4] + b.z * m.n[ 8];
    a.y = b.x * m.n[ 1] + b.y * m.n[ 5] + b.z * m.n[ 9];
    a.z = b.x * m.n[ 2] + b.y * m.n[ 6] + b.z * m.n[10];

    return a;
}

//-----------------------------------------------------------------------------
// returns the determinant of the matrix
//-----------------------------------------------------------------------------
real_t cyclone::Mat4Determinant(mat4_t m)
{
    real_t d = 0.0f;

    d = m.n[ 8] * m.n[ 5] * m.n[ 2] + m.n[ 4] * m.n[ 9] * m.n[ 2] +
        m.n[ 8] * m.n[ 1] * m.n[ 6] - m.n[ 0] * m.n[ 9] * m.n[ 6] -
        m.n[ 4] * m.n[ 1] * m.n[10] + m.n[ 0] * m.n[ 5] * m.n[10];

    return d;
}

//-----------------------------------------------------------------------------
// returns a new matrix containing the inverse of this matrix
//-----------------------------------------------------------------------------
mat4_t cyclone::Mat4Inverse(mat4_t m)
{
    mat4_t a = { 0.0f };
    real_t d = 0.0f;            /* make sure the determinant is non-zero */

    d = Mat4Determinant( m );

    if ( d == 0.0f ) {
        return m;
    }
    
    d = 1.0f / d;

    a.n[ 0] = (-m.n[ 9] * m.n[ 6] + m.n[ 5] * m.n[10] ) * d;
    a.n[ 4] = ( m.n[ 8] * m.n[ 6] - m.n[ 4] * m.n[10] ) * d;
    a.n[ 8] = (-m.n[ 8] * m.n[ 5] + m.n[ 4] * m.n[ 9] * m.n[15] ) * d;

    a.n[ 1] = ( m.n[ 9] * m.n[ 2] - m.n[ 1] * m.n[10] ) * d;
    a.n[ 5] = (-m.n[ 8] * m.n[ 2] + m.n[ 0] * m.n[10] ) * d;
    a.n[ 9] = ( m.n[ 8] * m.n[ 1] - m.n[ 0] * m.n[ 9] * m.n[15] ) * d;

    a.n[ 2] = (-m.n[ 5] * m.n[ 2] + m.n[ 1] * m.n[ 6] * m.n[15] ) * d;
    a.n[ 6] = ( m.n[ 4] * m.n[ 2] - m.n[ 0] * m.n[ 6] * m.n[15] ) * d;
    a.n[10] = (-m.n[ 4] * m.n[ 1] + m.n[ 0] * m.n[ 5] * m.n[15] ) * d;

    a.n[ 3] = ( m.n[ 9] * m.n[ 6] * m.n[ 3] -
        m.n[ 5] * m.n[10] * m.n[ 3] -
        m.n[ 9] * m.n[ 2] * m.n[ 7] +
        m.n[ 1] * m.n[10] * m.n[ 7] +
        m.n[ 5] * m.n[ 2] * m.n[11] -
        m.n[ 1] * m.n[ 6] * m.n[11] ) * d;
    a.n[ 6] = (-m.n[ 8] * m.n[ 6] * m.n[ 3] +
        m.n[ 4] * m.n[10] * m.n[ 3] +
        m.n[ 8] * m.n[ 2] * m.n[ 7] -
        m.n[ 0] * m.n[10] * m.n[ 7] -
        m.n[ 4] * m.n[ 2] * m.n[11] +
        m.n[ 0] * m.n[ 6] * m.n[11] ) * d;
    a.n[11] = ( m.n[ 8] * m.n[ 5] * m.n[ 3] -
        m.n[ 4] * m.n[ 9] * m.n[ 3] -
        m.n[ 8] * m.n[ 1] * m.n[ 7] +
        m.n[ 0] * m.n[ 9] * m.n[ 7] +
        m.n[ 4] * m.n[ 1] * m.n[11] -
        m.n[ 0] * m.n[ 5] * m.n[11] ) * d;

    return a;
}

//-----------------------------------------------------------------------------
// transform the given direction vector by this matrix
//-----------------------------------------------------------------------------
vec3_t cyclone::Mat4TransformDirection(vec3_t v, mat4_t m)
{
    vec3_t a = { 0.0f };

    a.x = v.x * m.n[ 0] + v.y * m.n[ 1] + v.z * m.n[ 2];
    a.y = v.x * m.n[ 4] + v.y * m.n[ 5] + v.z * m.n[ 6];
    a.z = v.x * m.n[ 8] + v.y * m.n[ 9] + v.z * m.n[10];

    return a;
}

//-----------------------------------------------------------------------------
// transform the given direction vector by the transformational inverse of this
// matrix
//-----------------------------------------------------------------------------
vec3_t cyclone::Mat4TransformInverseDirection(vec3_t v, mat4_t m)
{
    vec3_t a = { 0.0f };

    a.x = v.x * m.n[ 0] + v.y * m.n[ 4] + v.z * m.n[ 8];
    a.y = v.x * m.n[ 1] + v.y * m.n[ 5] + v.z * m.n[ 9];
    a.z = v.x * m.n[ 2] + v.y * m.n[ 6] + v.z * m.n[10];

    return a;
}

//-----------------------------------------------------------------------------
// gets a vector representing one axis (i.e. one column) in the matrix
//-----------------------------------------------------------------------------
vec3_t cyclone::Mat4AxisVector(mat4_t m, int i)
{
    vec3_t v = { 0.0f, 0.0f, 0.0f };

    if ( i < 0 || i > 3 ) {
        return v;
    }

    v.x = m.n[i+0];
    v.y = m.n[i+4];
    v.z = m.n[i+8];

    return v;
}

//-----------------------------------------------------------------------------
// sets this matrix to be the rotation matrix corresponding to the given
// quaternion
//-----------------------------------------------------------------------------
mat4_t cyclone::Mat4SetOrientationAndPos(quat_t q, vec3_t v)
{
    mat4_t m = { 0.0f };

    m.n[ 0] = 1 - ( 2 * q.j * q.j + 2 * q.k * q.k );
    m.n[ 1] = 2 * q.i * q.j + 2 * q.k * q.r;
    m.n[ 2] = 2 * q.i * q.k - 2 * q.j * q.r;
    m.n[ 3] = v.x;

    m.n[ 4] = 2 * q.i * q.j - 2 * q.k * q.r;
    m.n[ 5] = 1 - ( 2 * q.i * q.i  + 2 * q.k * q.k );
    m.n[ 6] = 2 * q.j * q.k + 2 * q.i * q.r;
    m.n[ 7] = v.y;

    m.n[ 8] = 2 * q.i * q.k + 2 * q.j * q.r;
    m.n[ 9] = 2 * q.j * q.k - 2 * q.i * q.r;
    m.n[10] = 1 - ( 2 * q.i * q.i  + 2 * q.j * q.j );
    m.n[11] = v.z;

    return m;
}

//-----------------------------------------------------------------------------
// fills the given array with this transform matrix, so it is usable as an
// open-gl transform matrix. OpenGL uses a column major format, so that the
// values are transposed as they are written
//-----------------------------------------------------------------------------
mat4_t cyclone::Mat4FillGLArray(mat4_t m)
{
    mat4_t a = { 0.0f };

    a.n[ 0] = m.n[ 0];
    a.n[ 1] = m.n[ 4];
    a.n[ 2] = m.n[ 8];
    a.n[ 3] = 0.0f;

    a.n[ 4] = m.n[ 1];
    a.n[ 5] = m.n[ 5];
    a.n[ 6] = m.n[ 9];
    a.n[ 7] = 0.0f;

    a.n[ 8] = m.n[ 2];
    a.n[ 9] = m.n[ 6];
    a.n[10] = m.n[10];
    a.n[11] = 0.0f;

    a.n[12] = m.n[ 3];
    a.n[13] = m.n[ 7];
    a.n[14] = m.n[11];
    a.n[15] = 1.0f;

    return a;
}


//-----------------------------------------------------------------------------
// MATRIX 3
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
mat3_t cyclone::Mat3Identity(void)
{
    mat3_t m = { 0.0f };

    m.n[0] = 1.0f;
    m.n[1] = 0.0f;
    m.n[2] = 0.0f;
    m.n[3] = 0.0f;
    m.n[4] = 1.0f;
    m.n[5] = 0.0f;
    m.n[6] = 0.0f;
    m.n[7] = 0.0f;
    m.n[8] = 1.0f;

    return m;
}

//-----------------------------------------------------------------------------
// sets the matrix to be a diagonal matrix with the given values along the
// leading diagonal
//-----------------------------------------------------------------------------
mat3_t cyclone::Mat3SetDiagonal(real_t x, real_t y, real_t z)
{
    mat3_t m = { 0.0f };

    m.x0 = x;
    m.y1 = y;
    m.z2 = z;

    return m;
}

//-----------------------------------------------------------------------------
// sets the value of the matrix from inertia tensor values
//-----------------------------------------------------------------------------
mat3_t cyclone::Mat3SetInertiaTensorCoeffs(real_t ix, real_t iy, real_t iz,
    real_t ixy, real_t ixz, real_t iyz)
{
    mat3_t m = { 0.0f };

    m.n[0] = ix;
    m.n[1] = m.n[3] = -ixy;
    m.n[2] = m.n[6] = -ixz;
    m.n[3] = iy;
    m.n[4] = m.n[7] = -iyz;
    m.n[5] = iz;

    return m;
}

//-----------------------------------------------------------------------------
// sets the value of the matrix as an inertia tensor of a rectangular block
// aligned with the body's coordinate  system with the given axis half-sizes
// and mass
//-----------------------------------------------------------------------------
mat3_t cyclone::Mat3SetBlockInertiaTensor(vec3_t halfSizes, real_t mass)
{
    vec3_t squares = Vec3ComponentProduct( halfSizes, halfSizes );

    mat3_t m = cyclone::Mat3SetInertiaTensorCoeffs( 0.3f * mass * ( squares.y +
        squares.z ), 0.3f * mass * ( squares.x + squares.z ), 0.3f * mass *
        ( squares.x + squares.y ) );

    return m;
}

//-----------------------------------------------------------------------------
// sets the matrix to be a skew symmetric matrix based on the given vector. the
// skew symmetric matrix is the equivalent of the vector product. So if a,b are
// vectors. a x b = A_s b where A_s is the skew symmetric form of a
//-----------------------------------------------------------------------------
mat3_t cyclone::Mat3SetSkewSymmetric(vec3_t v)
{
    mat3_t m = { 0.0f };

    m.n[0] = m.n[4] = m.n[8] = 0.0f;
    m.n[1] = -v.z;
    m.n[2] =  v.y;
    m.n[3] =  v.z;
    m.n[5] = -v.x;
    m.n[6] = -v.y;
    m.n[7] =  v.x;

    return m;
}

//-----------------------------------------------------------------------------
// sets the matrix values from the given three vector components. these are
// arranged as the three columns of the vector
//-----------------------------------------------------------------------------
mat3_t cyclone::Mat3SetComponents(vec3_t compOne, vec3_t compTwo,
    vec3_t compThree)
{
    mat3_t m = { 0.0f };

    m.n[0] = compOne.x;
    m.n[1] = compTwo.x;
    m.n[2] = compThree.x;
    m.n[3] = compOne.y;
    m.n[4] = compTwo.y;
    m.n[5] = compThree.y;
    m.n[6] = compOne.z;
    m.n[7] = compTwo.z;
    m.n[8] = compThree.z;

    return m;
}

//-----------------------------------------------------------------------------
// transform the given vector by this matrix
//-----------------------------------------------------------------------------
vec3_t cyclone::Mat3Transform(vec3_t v, mat3_t m)
{
    vec3_t a = { 0.0f };

    a.x = v.x * m.n[0] + v.y * m.n[1] + v.z * m.n[2];
    a.y = v.x * m.n[3] + v.y * m.n[4] + v.z * m.n[5];
    a.z = v.x * m.n[6] + v.y * m.n[7] + v.z * m.n[8];

    return a;
}

//-----------------------------------------------------------------------------
// transform the given vector by the transpose of this matrix
//-----------------------------------------------------------------------------
vec3_t cyclone::Mat3TransformTranspose(vec3_t v, mat3_t m)
{
    vec3_t a = { 0.0f };

    a.x = v.x * m.n[0] + v.y * m.n[3] + v.z * m.n[6];
    a.y = v.x * m.n[1] + v.y * m.n[4] + v.z * m.n[7];
    a.z = v.x * m.n[2] + v.y * m.n[5] + v.z * m.n[8];

    return a;
}

//-----------------------------------------------------------------------------
// gets a vector representing one row in the matrix
//-----------------------------------------------------------------------------
vec3_t cyclone::Mat3RowVector(mat3_t m, int i)
{
    vec3_t v = { m.n[i*3+0], m.n[i*3+1], m.n[i*3+2] };
    return v;
}

//-----------------------------------------------------------------------------
// gets a vector representing one axis (i.e. one column) in the matrix
//-----------------------------------------------------------------------------
vec3_t cyclone::Mat3AxisVector(mat3_t m, int i)
{
    vec3_t v = { m.n[i+0], m.n[i+3], m.n[i+6] };
    return v;
}

//-----------------------------------------------------------------------------
// sets the matrix to be the inverse of the given matrix
//-----------------------------------------------------------------------------
mat3_t cyclone::Mat3Inverse(mat3_t m)
{
    mat3_t a = { 0.0f };

    real_t t04 = m.n[0] * m.n[4];
    real_t t06 = m.n[0] * m.n[5];
    real_t t08 = m.n[1] * m.n[3];
    real_t t10 = m.n[2] * m.n[3];
    real_t t12 = m.n[1] * m.n[6];
    real_t t14 = m.n[2] * m.n[6];

    // calculate the determinant
    real_t t16 = ( t04 * m.n[8] - t06 * m.n[7] - t08 * m.n[8] +
        t10 * m.n[7] + t12 * m.n[5] - t14 * m.n[4] );

    // make sure the determinant is non-zero
    if ( t16 == 0.0f ) {
        return m;
    }
    
    real_t t17 = 1 / t16;

    a.n[0] =  ( m.n[4] * m.n[8] - m.n[5] * m.n[7] ) * t17;
    a.n[1] = -( m.n[1] * m.n[8] - m.n[2] * m.n[7] ) * t17;
    a.n[2] =  ( m.n[1] * m.n[5] - m.n[2] * m.n[4] ) * t17;
    a.n[3] = -( m.n[3] * m.n[8] - m.n[5] * m.n[6] ) * t17;
    a.n[4] =  ( m.n[0] * m.n[8] - t14 ) * t17;
    a.n[5] = -( t06 - t10 ) * t17;
    a.n[6] =  ( m.n[3] * m.n[7] - m.n[4] * m.n[6] ) * t17;
    a.n[7] = -( m.n[0] * m.n[7] - t12) * t17;
    a.n[8] =  ( t04 - t08 ) * t17;

    return a;
}

//-----------------------------------------------------------------------------
// sets the matrix to be the transpose of the given matrix
//-----------------------------------------------------------------------------
mat3_t cyclone::Mat3Transpose(mat3_t m)
{
    mat3_t a = { 0.0f };

    a.n[0] = m.n[0];
    a.n[1] = m.n[3];
    a.n[2] = m.n[6];
    a.n[3] = m.n[1];
    a.n[4] = m.n[4];
    a.n[5] = m.n[7];
    a.n[6] = m.n[2];
    a.n[7] = m.n[5];
    a.n[8] = m.n[8];

    return a;
}

//-----------------------------------------------------------------------------
// returns a matrix which is this matrix multiplied by the given other matrix
//-----------------------------------------------------------------------------
mat3_t cyclone::Mat3Multiply(mat3_t a, mat3_t b)
{
    mat3_t m = { 0.0f };

    m.n[0] = a.n[0] * b.n[0] + a.n[1] * b.n[3] + a.n[2] * b.n[6];
    m.n[1] = a.n[0] * b.n[1] + a.n[1] * b.n[4] + a.n[2] * b.n[7];
    m.n[2] = a.n[0] * b.n[2] + a.n[1] * b.n[5] + a.n[2] * b.n[8];

    m.n[3] = a.n[3] * b.n[0] + a.n[4] * b.n[3] + a.n[5] * b.n[6];
    m.n[4] = a.n[3] * b.n[1] + a.n[4] * b.n[4] + a.n[5] * b.n[7];
    m.n[5] = a.n[3] * b.n[2] + a.n[4] * b.n[5] + a.n[5] * b.n[8];

    m.n[6] = a.n[6] * b.n[0] + a.n[7] * b.n[3] + a.n[8] * b.n[6];
    m.n[7] = a.n[6] * b.n[1] + a.n[7] * b.n[4] + a.n[8] * b.n[7];
    m.n[8] = a.n[6] * b.n[2] + a.n[7] * b.n[5] + a.n[8] * b.n[8];

    return m;
}

//-----------------------------------------------------------------------------
// multiplies this matrix in place by the given scalar
//-----------------------------------------------------------------------------
mat3_t cyclone::Mat3Scale(mat3_t m, real_t scalar)
{
    mat3_t a = { 0.0f };

    a = m;

    a.n[0] *= scalar; a.n[1] *= scalar; a.n[2] *= scalar;
    a.n[3] *= scalar; a.n[4] *= scalar; a.n[5] *= scalar;
    a.n[6] *= scalar; a.n[7] *= scalar; a.n[8] *= scalar;

    return a;
}

//-----------------------------------------------------------------------------
// does a component-wise addition of this matrix and the given matrix
//-----------------------------------------------------------------------------
mat3_t cyclone::Mat3Add(mat3_t a, mat3_t b)
{
    mat3_t m = { 0.0f };

    m = a;

    m.n[0] += b.n[0];
    m.n[1] += b.n[1];
    m.n[2] += b.n[2];
    
    m.n[3] += b.n[3];
    m.n[4] += b.n[4];
    m.n[5] += b.n[5];
    
    m.n[6] += b.n[6];
    m.n[7] += b.n[7];
    m.n[8] += b.n[8];

    return m;
}

//-----------------------------------------------------------------------------
// sets this matrix to be the rotation matrix corresponding to the given
// quaternion
//-----------------------------------------------------------------------------
mat3_t cyclone::Mat3SetOrientation(quat_t q)
{
    mat3_t m = { 0.0f };

    m.n[0] = 1 - (2 * q.j * q.j + 2 * q.k * q.k);
    m.n[1] = 2 * q.i * q.j + 2 * q.k * q.r;
    m.n[2] = 2 * q.i * q.k - 2 * q.j * q.r;
    m.n[3] = 2 * q.i * q.j - 2 * q.k * q.r;
    m.n[4] = 1 - (2 * q.i * q.i  + 2 * q.k * q.k);
    m.n[5] = 2 * q.j * q.k + 2 * q.i * q.r;
    m.n[6] = 2 * q.i * q.k + 2 * q.j * q.r;
    m.n[7] = 2 * q.j * q.k - 2 * q.i * q.r;
    m.n[8] = 1 - (2 * q.i*q.i  + 2 * q.j * q.j);

    return m;
}

//-----------------------------------------------------------------------------
// interpolates a couple of matrices
//-----------------------------------------------------------------------------
mat3_t cyclone::Mat3LinearInterpolate(mat3_t a, mat3_t b, real_t prop)
{
    mat3_t m = { 0.0f };

    for (unsigned i = 0; i < 3; i++) {
        m.n[i*3+0] = a.n[i*3+0] * ( 1 - prop ) + b.n[i*3+0] * prop;
        m.n[i*3+1] = a.n[i*3+1] * ( 1 - prop ) + b.n[i*3+1] * prop;
        m.n[i*3+2] = a.n[i*3+2] * ( 1 - prop ) + b.n[i*3+2] * prop;
    }

    return m;
}

const Vector3 Vector3::GRAVITY = Vector3(0, -9.81f, 0);
const Vector3 Vector3::HIGH_GRAVITY = Vector3(0, -19.62f, 0);
const Vector3 Vector3::UP = Vector3(0, 1, 0);
const Vector3 Vector3::RIGHT = Vector3(1, 0, 0);
const Vector3 Vector3::OUT = Vector3(0, 0, 1);
const Vector3 Vector3::X = Vector3(1, 0, 0);
const Vector3 Vector3::Y = Vector3(0, 1, 0);
const Vector3 Vector3::Z = Vector3(0, 0, 1);

///>Matrix4Inverse
real Matrix4::getDeterminant() const
{
    return data[8]*data[5]*data[2]+
        data[4]*data[9]*data[2]+
        data[8]*data[1]*data[6]-
        data[0]*data[9]*data[6]-
        data[4]*data[1]*data[10]+
        data[0]*data[5]*data[10];
}

void Matrix4::setInverse(const Matrix4 &m)
{
    // Make sure the determinant is non-zero.
    real det = getDeterminant();
    if (det == 0) return;
    det = ((real)1.0)/det;

    data[0] = (-m.data[9]*m.data[6]+m.data[5]*m.data[10])*det;
    data[4] = (m.data[8]*m.data[6]-m.data[4]*m.data[10])*det;
    data[8] = (-m.data[8]*m.data[5]+m.data[4]*m.data[9]* m.data[15])*det;

    data[1] = (m.data[9]*m.data[2]-m.data[1]*m.data[10])*det;
    data[5] = (-m.data[8]*m.data[2]+m.data[0]*m.data[10])*det;
    data[9] = (m.data[8]*m.data[1]-m.data[0]*m.data[9]* m.data[15])*det;

    data[2] = (-m.data[5]*m.data[2]+m.data[1]*m.data[6]* m.data[15])*det;
    data[6] = (+m.data[4]*m.data[2]-m.data[0]*m.data[6]* m.data[15])*det;
    data[10] = (-m.data[4]*m.data[1]+m.data[0]*m.data[5]* m.data[15])*det;

    data[3] = (m.data[9]*m.data[6]*m.data[3]
               -m.data[5]*m.data[10]*m.data[3]
               -m.data[9]*m.data[2]*m.data[7]
               +m.data[1]*m.data[10]*m.data[7]
               +m.data[5]*m.data[2]*m.data[11]
               -m.data[1]*m.data[6]*m.data[11])*det;
    data[7] = (-m.data[8]*m.data[6]*m.data[3]
               +m.data[4]*m.data[10]*m.data[3]
               +m.data[8]*m.data[2]*m.data[7]
               -m.data[0]*m.data[10]*m.data[7]
               -m.data[4]*m.data[2]*m.data[11]
               +m.data[0]*m.data[6]*m.data[11])*det;
    data[11] =(m.data[8]*m.data[5]*m.data[3]
               -m.data[4]*m.data[9]*m.data[3]
               -m.data[8]*m.data[1]*m.data[7]
               +m.data[0]*m.data[9]*m.data[7]
               +m.data[4]*m.data[1]*m.data[11]
               -m.data[0]*m.data[5]*m.data[11])*det;
}

Matrix3 Matrix3::linearInterpolate(const Matrix3& a, const Matrix3& b, real prop)
{
	Matrix3 result;
	for (unsigned i = 0; i < 9; i++) {
		result.data[i] = a.data[i] * (1-prop) + b.data[i] * prop;
	}
	return result;
}

