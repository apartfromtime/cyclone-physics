/*
 * Implementation file for the particle class.
 * 
 * Part of the Cyclone physics system.
 * 
 * Copyright (c) Icosagon 2003. All Rights Reserved.
 *
 * This software is distributed under licence. Use of this software
 * implies agreement with all terms and conditions of the accompanying
 * software licence.
 */

///>ParticleIntegrate
#include <assert.h>
///<ParticleIntegrate
///>ParticleSrc
#include <cyclone/particle.h>

using namespace cyclone;

///<ParticleSrc

/*
 * --------------------------------------------------------------------------
 * FUNCTIONS DECLARED IN HEADER:
 * --------------------------------------------------------------------------
 */

///>ParticleAccum
///>ParticleIntegrate
void Particle::integrate(real_t duration)
{
	// We don't integrate things with zero mass.
	if ( inverseMass <= 0.0f ) {
        return;
    }

    assert(duration > 0.0);

    // Update linear position.
    position = Vec3Add( position, Vec3Scale( velocity, duration ) );

    // Work out the acceleration from the force
    vec3_t resultingAcc = Vec3Add( acceleration,
        Vec3Scale( forceAccum, inverseMass ) );
    
    // Update linear velocity from the acceleration.
    velocity = Vec3Add( velocity, Vec3Scale( resultingAcc, duration ) );

    // Impose drag.
    velocity = Vec3Scale( velocity, R_pow( damping, duration ) );
///<ParticleIntegrate

    // Clear the forces.
    clearAccumulator();
///>ParticleIntegrate
}
///<ParticleIntegrate

///<ParticleAccum


void Particle::setMass(const real_t mass)
{
    assert(mass != 0);
    Particle::inverseMass = ( ( real_t )1.0f ) / mass;
}

real_t Particle::getMass(void) const
{
    if (inverseMass == 0) {
        return REAL_MAX;
    } else {
        return ( ( real_t )1.0f ) / inverseMass;
    }
}

void Particle::setInverseMass(const real_t inverseMass)
{
    Particle::inverseMass = inverseMass;
}

real_t Particle::getInverseMass(void) const
{
    return inverseMass;
}

bool Particle::hasFiniteMass(void) const
{
    return inverseMass >= 0.0f;
}

void Particle::setDamping(const real_t damping)
{
    Particle::damping = damping;
}

real_t Particle::getDamping(void) const
{
    return damping;
}

void Particle::setPosition(const vec3_t & position)
{
    Particle::position = position;
}

void Particle::setPosition(const real_t x, const real_t y, const real_t z)
{
    position.x = x;
    position.y = y;
    position.z = z;
}

void Particle::getPosition(vec3_t * position) const
{
    *position = Particle::position;
}

vec3_t Particle::getPosition(void) const
{
    return position;
}

void Particle::setVelocity(const vec3_t & velocity)
{
    Particle::velocity = velocity;
}

void Particle::setVelocity(const real_t x, const real_t y, const real_t z)
{
    velocity.x = x;
    velocity.y = y;
    velocity.z = z;
}

void Particle::getVelocity(vec3_t * velocity) const
{
    *velocity = Particle::velocity;
}

vec3_t Particle::getVelocity(void) const
{
    return velocity;
}

void Particle::setAcceleration(const vec3_t & acceleration)
{
    Particle::acceleration = acceleration;
}

void Particle::setAcceleration(const real_t x, const real_t y, const real_t z)
{
    acceleration.x = x;
    acceleration.y = y;
    acceleration.z = z;
}

void Particle::getAcceleration(vec3_t * acceleration) const
{
    *acceleration = Particle::acceleration;
}

vec3_t Particle::getAcceleration(void) const
{
    return acceleration;
}

///>ParticleAccum
void Particle::clearAccumulator(void)
{
    forceAccum = Vec3Clear();
}
///<ParticleAccum

///>ParticleAddForce
void Particle::addForce(const vec3_t & force)
{
    forceAccum = Vec3Add( forceAccum, force );
}
///<ParticleAddForce




