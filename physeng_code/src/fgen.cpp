/*
 * Implementation file for the rigid body force generators.
 *
 * Part of the Cyclone physics system.
 *
 * Copyright (c) Icosagon 2003. All Rights Reserved.
 *
 * This software is distributed under license. Use of this software
 * implies agreement with all terms and conditions of the accompanying
 * software license.
 */

#include <cyclone/fgen.h>

using namespace cyclone;

void ForceRegistry::updateForces(real_t duration)
{
	Registry::iterator i = registrations.begin();
	for (; i != registrations.end(); i++)
	{
		i->fg->updateForce(i->body, duration);
	}
}

void ForceRegistry::add(RigidBody *body, ForceGenerator *fg)
{
	ForceRegistry::ForceRegistration registration;
	registration.body = body;
	registration.fg = fg;
	registrations.push_back(registration);
}

Buoyancy::Buoyancy(const vec3_t & cOfB, real_t maxDepth, real_t volume,
    real_t waterHeight, real_t liquidDensity /* = 1000.0f */)
{
	centreOfBuoyancy = cOfB;
	Buoyancy::liquidDensity = liquidDensity;
	Buoyancy::maxDepth = maxDepth;
	Buoyancy::volume = volume;
	Buoyancy::waterHeight = waterHeight;
}

void Buoyancy::updateForce(RigidBody * body, real_t duration)
{
	// Calculate the submersion depth
	vec3_t pointInWorld = body->getPointInWorldSpace( centreOfBuoyancy );
	real_t depth = pointInWorld.y;

	// Check if we're out of the water
	if (depth >= waterHeight + maxDepth) return;
	vec3_t force = { 0, 0, 0 };

	// Check if we're at maximum depth
	if ( depth <= waterHeight - maxDepth ) {

		force.y = liquidDensity * volume;
		body->addForceAtBodyPoint( force, centreOfBuoyancy );
		
        return;
	}

	// Otherwise we are partly submerged
	force.y = liquidDensity * volume * ( depth - maxDepth - waterHeight ) /
        ( 2 * maxDepth );
	body->addForceAtBodyPoint( force, centreOfBuoyancy );
}

Gravity::Gravity(const vec3_t & gravity)
: gravity(gravity)
{
}

///>GravityFG
void Gravity::updateForce(RigidBody * body, real_t duration)
{
    // Check that we do not have infinite mass
    if ( !body->hasFiniteMass() ) return;

    // Apply the mass-scaled force to the body
    body->addForce( Vec3Scale( gravity, body->getMass() ) );
}
///<GravityFG

Spring::Spring(const vec3_t & localConnectionPt, RigidBody * other,
    const vec3_t & otherConnectionPt, real_t springConstant, real_t restLength)
: connectionPoint(localConnectionPt),
  otherConnectionPoint(otherConnectionPt),
  other(other),
  springConstant(springConstant),
  restLength(restLength)
{
}

///>SpringFG
void Spring::updateForce(RigidBody * body, real_t duration)
{
    // Calculate the two ends in world space
	vec3_t lws = body->getPointInWorldSpace( connectionPoint );
	vec3_t ows = other->getPointInWorldSpace( otherConnectionPoint );

    // Calculate the vector of the spring
    vec3_t force = Vec3Subtract( lws, ows );

    // Calculate the magnitude of the force
    real_t magnitude = Vec3Magnitude( force );
    magnitude = R_abs( magnitude - restLength );
    magnitude *= springConstant;

    // Calculate the final force and apply it
    force = Vec3Normalise( force );
    force = Vec3Scale( force, -magnitude );

    body->addForceAtPoint( force, lws );
}
///<SpringFG

Aero::Aero(const mat3_t & tensor, const vec3_t & position,
    const vec3_t * windspeed)
{
	Aero::tensor = tensor;
	Aero::position = position;
	Aero::windspeed = windspeed;
}

void Aero::updateForce(RigidBody * body, real_t duration)
{
	Aero::updateForceFromTensor(body, duration, tensor);
}

void Aero::updateForceFromTensor(RigidBody * body, real_t duration,
    const mat3_t & tensor)
{
	// Calculate total velocity (windspeed and body's velocity).
	vec3_t velocity = Vec3Add( body->getVelocity(), *windspeed );

	// Calculate the velocity in body coordinates
    vec3_t bodyVel = Mat4TransformInverseDirection( velocity,
        body->getTransform() );
	
	// Calculate the force in body coordinates
	vec3_t bodyForce = Mat3Transform( bodyVel, tensor );
    vec3_t force = Mat4TransformDirection( bodyForce, body->getTransform() );

	// Apply the force
	body->addForceAtBodyPoint( force, position );
}

AeroControl::AeroControl(const mat3_t & base, const mat3_t & min,
    const mat3_t & max, const vec3_t & position, const vec3_t * windspeed)
:
Aero(base, position, windspeed)
{
	AeroControl::minTensor = min;
	AeroControl::maxTensor = max;
	controlSetting = 0.0f;
}

mat3_t AeroControl::getTensor()
{
	if (controlSetting <= -1.0f) return minTensor;
	else if (controlSetting >= 1.0f) return maxTensor;
	else if (controlSetting < 0)
	{
		return Mat3LinearInterpolate( minTensor, tensor,
            controlSetting + 1.0f );
	}
	else if (controlSetting > 0)
	{
		return Mat3LinearInterpolate( tensor, maxTensor, controlSetting );
	}
	else return tensor;
}

void AeroControl::setControl(real_t value)
{
	controlSetting = value;
}

void AeroControl::updateForce(RigidBody * body, real_t duration)
{
	mat3_t tensor = getTensor();
	Aero::updateForceFromTensor( body, duration, tensor );
}

void Explosion::updateForce(RigidBody * body, real_t duration)
{
    
}

