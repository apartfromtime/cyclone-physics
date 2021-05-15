/*
 * Implementation file for the rigid body class.
 *
 * Part of the Cyclone physics system.
 *
 * Copyright (c) Icosagon 2003. All Rights Reserved.
 *
 * This software is distributed under licence. Use of this software
 * implies agreement with all terms and conditions of the accompanying
 * software licence.
 */


#include <cyclone/body.h>
#include <memory.h>
#include <assert.h>

using namespace cyclone;


/*
 * --------------------------------------------------------------------------
 * INTERNAL OR HELPER FUNCTIONS:
 * --------------------------------------------------------------------------
 */

/**
 * Internal function that checks the validity of an inverse inertia tensor.
 */
static inline void _checkInverseInertiaTensor(const mat3_t & iitWorld)
{
    // TODO: Perform a validity check in an assert.
}

///>BodyDerivedTransformIT
/**
 * Internal function to do an intertia tensor transform by a quaternion.
 * Note that the implementation of this function was created by an
 * automated code-generator and optimizer.
 */
static inline void _transformInertiaTensor(mat3_t & iitWorld, const quat_t & q,
    const mat3_t & iitBody, const mat4_t & rotmat)
{
    real_t t4 = rotmat.n[0]*iitBody.n[0]+
        rotmat.n[1]*iitBody.n[3]+
        rotmat.n[2]*iitBody.n[6];
    real_t t9 = rotmat.n[0]*iitBody.n[1]+
        rotmat.n[1]*iitBody.n[4]+
        rotmat.n[2]*iitBody.n[7];
    real_t t14 = rotmat.n[0]*iitBody.n[2]+
        rotmat.n[1]*iitBody.n[5]+
        rotmat.n[2]*iitBody.n[8];
    real_t t28 = rotmat.n[4]*iitBody.n[0]+
        rotmat.n[5]*iitBody.n[3]+
        rotmat.n[6]*iitBody.n[6];
    real_t t33 = rotmat.n[4]*iitBody.n[1]+
        rotmat.n[5]*iitBody.n[4]+
        rotmat.n[6]*iitBody.n[7];
    real_t t38 = rotmat.n[4]*iitBody.n[2]+
        rotmat.n[5]*iitBody.n[5]+
        rotmat.n[6]*iitBody.n[8];
    real_t t52 = rotmat.n[8]*iitBody.n[0]+
        rotmat.n[9]*iitBody.n[3]+
        rotmat.n[10]*iitBody.n[6];
    real_t t57 = rotmat.n[8]*iitBody.n[1]+
        rotmat.n[9]*iitBody.n[4]+
        rotmat.n[10]*iitBody.n[7];
    real_t t62 = rotmat.n[8]*iitBody.n[2]+
        rotmat.n[9]*iitBody.n[5]+
        rotmat.n[10]*iitBody.n[8];

    iitWorld.n[0] = t4*rotmat.n[0]+
        t9*rotmat.n[1]+
        t14*rotmat.n[2];
    iitWorld.n[1] = t4*rotmat.n[4]+
        t9*rotmat.n[5]+
        t14*rotmat.n[6];
    iitWorld.n[2] = t4*rotmat.n[8]+
        t9*rotmat.n[9]+
        t14*rotmat.n[10];
    iitWorld.n[3] = t28*rotmat.n[0]+
        t33*rotmat.n[1]+
        t38*rotmat.n[2];
    iitWorld.n[4] = t28*rotmat.n[4]+
        t33*rotmat.n[5]+
        t38*rotmat.n[6];
    iitWorld.n[5] = t28*rotmat.n[8]+
        t33*rotmat.n[9]+
        t38*rotmat.n[10];
    iitWorld.n[6] = t52*rotmat.n[0]+
        t57*rotmat.n[1]+
        t62*rotmat.n[2];
    iitWorld.n[7] = t52*rotmat.n[4]+
        t57*rotmat.n[5]+
        t62*rotmat.n[6];
    iitWorld.n[8] = t52*rotmat.n[8]+
        t57*rotmat.n[9]+
        t62*rotmat.n[10];
}
///<BodyDerivedTransformIT

///>BodyDerivedTransform
/**
 * Inline function that creates a transform matrix from a
 * position and orientation.
 */
static inline void _calculateTransformMatrix(mat4_t & transformMatrix,
    const vec3_t & position, const quat_t & orientation)
{
    transformMatrix.n[ 0] = 1-2*orientation.j*orientation.j-
        2*orientation.k*orientation.k;
    transformMatrix.n[ 1] = 2*orientation.i*orientation.j -
        2*orientation.r*orientation.k;
    transformMatrix.n[ 2] = 2*orientation.i*orientation.k +
        2*orientation.r*orientation.j;
    transformMatrix.n[ 3] = position.x;

    transformMatrix.n[ 4] = 2*orientation.i*orientation.j +
        2*orientation.r*orientation.k;
    transformMatrix.n[ 5] = 1-2*orientation.i*orientation.i-
        2*orientation.k*orientation.k;
    transformMatrix.n[ 6] = 2*orientation.j*orientation.k -
        2*orientation.r*orientation.i;
    transformMatrix.n[ 7] = position.y;

    transformMatrix.n[ 8] = 2*orientation.i*orientation.k -
        2*orientation.r*orientation.j;
    transformMatrix.n[ 9] = 2*orientation.j*orientation.k +
        2*orientation.r*orientation.i;
    transformMatrix.n[10] = 1-2*orientation.i*orientation.i-
        2*orientation.j*orientation.j;
    transformMatrix.n[11] = position.z;
}
///<BodyDerivedTransform

/*
 * --------------------------------------------------------------------------
 * FUNCTIONS DECLARED IN HEADER:
 * --------------------------------------------------------------------------
 */
///>BodyDerivedBase
void RigidBody::calculateDerivedData()
{
///<BodyDerivedBase
    orientation = QuatNormalise( orientation );

///>BodyDerivedTransform
    // Calculate the transform matrix for the body.
    _calculateTransformMatrix( transformMatrix, position, orientation );
///<BodyDerivedTransform

///>BodyDerivedTransformIT
    // Calculate the inertiaTensor in world space.
    _transformInertiaTensor(inverseInertiaTensorWorld,
        orientation,
        inverseInertiaTensor,
        transformMatrix);
///<BodyDerivedTransformIT

///>BodyDerivedBase
}
///<BodyDerivedBase

///>BodyIntegrateBase
void RigidBody::integrate(real_t duration)
{
///<BodyIntegrateBase
    if (!isAwake) return;

///>BodyIntegrate10
///>LastUpdateAcceleration
    // Calculate linear acceleration from force inputs.
    lastFrameAcceleration = Vec3Add( acceleration,
        Vec3Scale( forceAccum, inverseMass ) );
///<LastUpdateAcceleration

    // Calculate angular acceleration from torque inputs.
    vec3_t angularAcceleration = Mat3Transform( torqueAccum,
        inverseInertiaTensorWorld );

    // Adjust velocities
    // Update linear velocity from both acceleration and impulse.
    velocity = Vec3Add( velocity,
        Vec3Scale( lastFrameAcceleration, duration ) );

    // Update angular velocity from both acceleration and impulse.
    rotation = Vec3Add( rotation, Vec3Scale( angularAcceleration, duration ) );

    // Impose drag.
    velocity = Vec3Scale( velocity, R_pow( linearDamping, duration ) );
    rotation = Vec3Scale( rotation, R_pow( angularDamping, duration ) );

    // Adjust positions
    // Update linear position.
    position = Vec3Add( position, Vec3Scale( velocity, duration ) );

    // Update angular position.
    vec3_t rot = { rotation.x, rotation.y, rotation.z };
    orientation = QuatAddScaledVector( orientation, rot, duration);

    // Impose drag.
    velocity = Vec3Scale( velocity, R_pow( linearDamping, duration ) );
    rotation = Vec3Scale( rotation, R_pow( angularDamping, duration ) );

    // Normalise the orientation, and update the matrices with the new
    // position and orientation
    calculateDerivedData();

///>ClearAccumulators
    // Clear accumulators.
    clearAccumulators();
///<ClearAccumulators
///<BodyIntegrate10

    // Update the kinetic energy store, and possibly put the body to
    // sleep.
    if (canSleep) {
        real_t currentMotion = Vec3ScalarProduct( velocity, velocity ) +
            Vec3ScalarProduct( rotation, rotation );

        real_t bias = R_pow(0.5, duration);
        motion = bias*motion + (1-bias)*currentMotion;

        if ( motion < GetSleepEpsilon() ) {
            setAwake( false );
        } else if ( motion > 10 * GetSleepEpsilon() ) {
            motion = 10 * GetSleepEpsilon();
        }
    }
///>BodyIntegrateBase
}
///<BodyIntegrateBase

void RigidBody::setMass(const real_t mass)
{
    assert(mass != 0);
    RigidBody::inverseMass = ( ( real_t )1.0f) / mass;
}

real_t RigidBody::getMass(void) const
{
    if ( inverseMass == 0 ) {
        return REAL_MAX;
    } else {
        return ( ( real_t )1.0f) / inverseMass;
    }
}

void RigidBody::setInverseMass(const real_t inverseMass)
{
    RigidBody::inverseMass = inverseMass;
}

real_t RigidBody::getInverseMass(void) const
{
    return inverseMass;
}

bool RigidBody::hasFiniteMass(void) const
{
    return inverseMass >= 0.0f;
}

///>SetInertiaTensor
void RigidBody::setInertiaTensor(const mat3_t & inertiaTensor)
{
    inverseInertiaTensor = Mat3Inverse( inertiaTensor );
///<SetInertiaTensor
///>SetInertiaTensor
    _checkInverseInertiaTensor(inverseInertiaTensor);
}
///<SetInertiaTensor

void RigidBody::getInertiaTensor(mat3_t * inertiaTensor) const
{
    *inertiaTensor = Mat3Inverse( inverseInertiaTensor );
}

mat3_t RigidBody::getInertiaTensor(void) const
{
    mat3_t it;
    getInertiaTensor( &it );
    return it;
}

void RigidBody::getInertiaTensorWorld(mat3_t * inertiaTensor) const
{
    *inertiaTensor = Mat3Inverse( inverseInertiaTensorWorld );
}

mat3_t RigidBody::getInertiaTensorWorld(void) const
{
    mat3_t it;
    getInertiaTensorWorld( &it );
    return it;
}

void RigidBody::setInverseInertiaTensor(const mat3_t & inverseInertiaTensor)
{
    _checkInverseInertiaTensor(inverseInertiaTensor);
    RigidBody::inverseInertiaTensor = inverseInertiaTensor;
}

void RigidBody::getInverseInertiaTensor(mat3_t * inverseInertiaTensor) const
{
    *inverseInertiaTensor = RigidBody::inverseInertiaTensor;
}

mat3_t RigidBody::getInverseInertiaTensor(void) const
{
    return inverseInertiaTensor;
}

void RigidBody::getInverseInertiaTensorWorld(mat3_t * inverseInertiaTensor) const
{
    *inverseInertiaTensor = inverseInertiaTensorWorld;
}

mat3_t RigidBody::getInverseInertiaTensorWorld(void) const
{
    return inverseInertiaTensorWorld;
}

void RigidBody::setDamping(const real_t linearDamping,
    const real_t angularDamping)
{
    RigidBody::linearDamping = linearDamping;
    RigidBody::angularDamping = angularDamping;
}

void RigidBody::setLinearDamping(const real_t linearDamping)
{
    RigidBody::linearDamping = linearDamping;
}

real_t RigidBody::getLinearDamping(void) const
{
    return linearDamping;
}

void RigidBody::setAngularDamping(const real_t angularDamping)
{
    RigidBody::angularDamping = angularDamping;
}

real_t RigidBody::getAngularDamping(void) const
{
    return angularDamping;
}

void RigidBody::setPosition(const vec3_t & position)
{
    RigidBody::position = position;
}

void RigidBody::setPosition(const real_t x, const real_t y, const real_t z)
{
    position.x = x;
    position.y = y;
    position.z = z;
}

void RigidBody::getPosition(vec3_t * position) const
{
    *position = RigidBody::position;
}

vec3_t RigidBody::getPosition(void) const
{
    return position;
}

void RigidBody::setOrientation(const quat_t & orientation)
{
    RigidBody::orientation = QuatNormalise( orientation );
}

void RigidBody::setOrientation(const real_t r, const real_t i, const real_t j,
    const real_t k)
{
    orientation.r = r;
    orientation.i = i;
    orientation.j = j;
    orientation.k = k;
    orientation = QuatNormalise( orientation );
}

void RigidBody::getOrientation(quat_t * orientation) const
{
    *orientation = RigidBody::orientation;
}

quat_t RigidBody::getOrientation(void) const
{
    return orientation;
}

void RigidBody::getOrientation(mat3_t * matrix) const
{
    getOrientation(matrix->n);
}

void RigidBody::getOrientation(real_t matrix[9]) const
{
    matrix[0] = transformMatrix.n[ 0];
    matrix[1] = transformMatrix.n[ 1];
    matrix[2] = transformMatrix.n[ 2];

    matrix[3] = transformMatrix.n[ 4];
    matrix[4] = transformMatrix.n[ 5];
    matrix[5] = transformMatrix.n[ 6];

    matrix[6] = transformMatrix.n[ 8];
    matrix[7] = transformMatrix.n[ 9];
    matrix[8] = transformMatrix.n[10];
}

void RigidBody::getTransform(mat4_t * transform) const
{
    memcpy( transform, &transformMatrix.n, sizeof( mat4_t ) );
}

void RigidBody::getTransform(real_t matrix[16]) const
{
    memcpy( matrix, transformMatrix.n, sizeof( real_t ) * 12 );
    matrix[12] = matrix[13] = matrix[14] = 0;
    matrix[15] = 1;
}

void RigidBody::getGLTransform(float matrix[16]) const
{
    matrix[ 0] = (float)transformMatrix.n[ 0];
    matrix[ 1] = (float)transformMatrix.n[ 4];
    matrix[ 2] = (float)transformMatrix.n[ 8];
    matrix[ 3] = 0;

    matrix[ 4] = (float)transformMatrix.n[ 1];
    matrix[ 5] = (float)transformMatrix.n[ 5];
    matrix[ 6] = (float)transformMatrix.n[ 9];
    matrix[ 7] = 0;

    matrix[ 8] = (float)transformMatrix.n[ 2];
    matrix[ 9] = (float)transformMatrix.n[ 6];
    matrix[10] = (float)transformMatrix.n[10];
    matrix[11] = 0;

    matrix[12] = (float)transformMatrix.n[ 3];
    matrix[13] = (float)transformMatrix.n[ 7];
    matrix[14] = (float)transformMatrix.n[11];
    matrix[15] = 1;
}

mat4_t RigidBody::getTransform(void) const
{
    return transformMatrix;
}


vec3_t RigidBody::getPointInLocalSpace(const vec3_t & point) const
{
    return Mat4TransformInverse( point, transformMatrix );
}

vec3_t RigidBody::getPointInWorldSpace(const vec3_t & point) const
{
    return Mat4Transform( point, transformMatrix );
}

vec3_t RigidBody::getDirectionInLocalSpace(const vec3_t & direction) const
{
    return Mat4TransformInverseDirection( direction, transformMatrix );
}

vec3_t RigidBody::getDirectionInWorldSpace(const vec3_t & direction) const
{
    return Mat4TransformDirection( direction, transformMatrix );
}


void RigidBody::setVelocity(const vec3_t & velocity)
{
    RigidBody::velocity = velocity;
}

void RigidBody::setVelocity(const real_t x, const real_t y, const real_t z)
{
    velocity.x = x;
    velocity.y = y;
    velocity.z = z;
}

void RigidBody::getVelocity(vec3_t * velocity) const
{
    *velocity = RigidBody::velocity;
}

vec3_t RigidBody::getVelocity(void) const
{
    return velocity;
}

void RigidBody::addVelocity(const vec3_t & deltaVelocity)
{
    velocity = Vec3Add( velocity, deltaVelocity );
}

void RigidBody::setRotation(const vec3_t & rotation)
{
    RigidBody::rotation = rotation;
}

void RigidBody::setRotation(const real_t x, const real_t y, const real_t z)
{
    rotation.x = x;
    rotation.y = y;
    rotation.z = z;
}

void RigidBody::getRotation(vec3_t * rotation) const
{
    *rotation = RigidBody::rotation;
}

vec3_t RigidBody::getRotation(void) const
{
    return rotation;
}

void RigidBody::addRotation(const vec3_t & deltaRotation)
{
    rotation = Vec3Add( rotation, deltaRotation );
}

///>SetAwake
void RigidBody::setAwake(const bool awake)
{
    if ( awake ) {

        isAwake = true;

        // Add a bit of motion to avoid it falling asleep immediately.
        motion = GetSleepEpsilon() * 2.0f;

    } else {
        
        isAwake = false;
        velocity = Vec3Clear();
        rotation = Vec3Clear();
    }
}
///<SetAwake

void RigidBody::setCanSleep(const bool canSleep)
{
    RigidBody::canSleep = canSleep;

    if ( !canSleep && !isAwake ) {
        setAwake();
    }
}


void RigidBody::getLastFrameAcceleration(vec3_t * acceleration) const
{
    *acceleration = lastFrameAcceleration;
}

vec3_t RigidBody::getLastFrameAcceleration(void) const
{
    return lastFrameAcceleration;
}

///>ClearAccumulators
void RigidBody::clearAccumulators(void)
{
    forceAccum = Vec3Clear();
    torqueAccum = Vec3Clear();
}
///<ClearAccumulators

///>AddForceAtCenter
void RigidBody::addForce(const vec3_t & force)
{
    forceAccum = Vec3Add( forceAccum, force );
///<AddForceAtCenter
    isAwake = true;
///>AddForceAtCenter
}
///<AddForceAtCenter

///>AddForceBody
void RigidBody::addForceAtBodyPoint(const vec3_t & force, const vec3_t & point)
{
    // Convert to coordinates relative to center of mass.
    vec3_t pt = getPointInWorldSpace( point );
    addForceAtPoint( force, pt );
///<AddForceBody

    isAwake = true;
///>AddForceBody
}
///<AddForceBody

void RigidBody::addForceAtPoint(const vec3_t & force, const vec3_t & point)
{
    // Convert to coordinates relative to center of mass.
    vec3_t pt = point;

    pt = Vec3Subtract( pt, position );

    forceAccum = Vec3Add( forceAccum, force );
    torqueAccum = Vec3Add( torqueAccum, Vec3VectorProduct( pt, force ) );

    isAwake = true;
}

void RigidBody::addTorque(const vec3_t & torque)
{
    torqueAccum = Vec3Add( torqueAccum, torque );
    isAwake = true;
}

void RigidBody::setAcceleration(const vec3_t & acceleration)
{
    RigidBody::acceleration = acceleration;
}

void RigidBody::setAcceleration(const real_t x, const real_t y, const real_t z)
{
    acceleration.x = x;
    acceleration.y = y;
    acceleration.z = z;
}

void RigidBody::getAcceleration(vec3_t * acceleration) const
{
    *acceleration = RigidBody::acceleration;
}

vec3_t RigidBody::getAcceleration(void) const
{
    return acceleration;
}
