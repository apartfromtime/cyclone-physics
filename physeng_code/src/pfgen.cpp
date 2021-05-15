/*
 * Implementation file for the particle force generators.
 * 
 * Part of the Cyclone physics system.
 * 
 * Copyright (c) Icosagon 2003. All Rights Reserved.
 *
 * This software is distributed under licence. Use of this software
 * implies agreement with all terms and conditions of the accompanying
 * software licence.
 */

///>ParticleFGUpdate
#include <cyclone/pfgen.h>

using namespace cyclone;

///<ParticleFGUpdate

///>ParticleFGUpdate
void ParticleForceRegistry::updateForces(real_t duration)
{
    Registry::iterator i = registrations.begin();
    for (; i != registrations.end(); i++)
    {
        i->fg->updateForce(i->particle, duration);
    }
}

void ParticleForceRegistry::add(Particle * particle,
    ParticleForceGenerator * fg)
{
	ParticleForceRegistry::ParticleForceRegistration registration;
	registration.particle = particle;
	registration.fg = fg;
	registrations.push_back(registration);
}
///<ParticleFGUpdate

ParticleGravity::ParticleGravity(const vec3_t & gravity)
: gravity(gravity)
{
}

///>GravityPFG
void ParticleGravity::updateForce(Particle * particle, real_t duration)
{
    // Check that we do not have infinite mass
    if ( !particle->hasFiniteMass() ) {
        return;
    }

    // Apply the mass-scaled force to the particle
    particle->addForce( Vec3Scale( gravity, particle->getMass() ) );
}
///<GravityPFG

ParticleDrag::ParticleDrag(real_t k1, real_t k2)
: k1(k1), k2(k2)
{
}

///>DragPFG
void ParticleDrag::updateForce(Particle * particle, real_t duration)
{
    vec3_t force = particle->getVelocity();

    // Calculate the total drag coefficient
    real_t dragCoeff = Vec3Magnitude( force );
    dragCoeff = k1 * dragCoeff + k2 * dragCoeff * dragCoeff;

    // Calculate the final force and apply it
    force = Vec3Normalise( force );
    force = Vec3Scale( force, -dragCoeff );
    particle->addForce( force );
}
///<DragPFG

ParticleSpring::ParticleSpring(Particle * other, real_t sc, real_t rl)
: other(other), springConstant(sc), restLength(rl)
{
}

///>SpringPFG
void ParticleSpring::updateForce(Particle * particle, real_t duration)
{
    // Calculate the vector of the spring
    vec3_t force = particle->getPosition();
    force = Vec3Subtract( force, other->getPosition() );

    // Calculate the magnitude of the force
    real_t magnitude = Vec3Magnitude( force );
    magnitude = R_abs( magnitude - restLength );
    magnitude *= springConstant;

    // Calculate the final force and apply it
    force = Vec3Normalise( force );
    force = Vec3Scale( force, -magnitude );
    particle->addForce( force );
}
///<SpringPFG

ParticleBuoyancy::ParticleBuoyancy(real_t maxDepth, real_t volume, 
    real_t waterHeight, real_t liquidDensity)
: 
maxDepth(maxDepth), volume(volume), 
waterHeight(waterHeight), liquidDensity(liquidDensity)
{
}

///>BoyancyPFG
void ParticleBuoyancy::updateForce(Particle * particle, real_t duration)
{
    // Calculate the submersion depth
    real_t depth = particle->getPosition().y;
   
    // Check if we're out of the water
    if ( depth >= waterHeight + maxDepth ) {
        return;
    }

    vec3_t force = { 0, 0, 0 };

    // Check if we're at maximum depth
    if ( depth <= waterHeight - maxDepth ) {

        force.y = liquidDensity * volume;
        particle->addForce( force );
        
        return;
    }

    // Otherwise we are partly submerged
    force.y = ( liquidDensity * volume ) * ( depth - maxDepth - waterHeight ) /
        ( 2 * maxDepth );
    particle->addForce( force );
}
///<BoyancyPFG

ParticleBungee::ParticleBungee(Particle * other, real_t sc, real_t rl)
: other(other), springConstant(sc), restLength(rl)
{
}

///>BungeePFG
void ParticleBungee::updateForce(Particle * particle, real_t duration)
{
    // Calculate the vector of the spring
    vec3_t force = particle->getPosition();
    force = Vec3Subtract( force, other->getPosition() );

    // Check if the bungee is compressed
    real_t magnitude = Vec3Magnitude( force );

    if ( magnitude <= restLength ) {
        return;
    }

    // Calculate the magnitude of the force
    magnitude = springConstant * (restLength - magnitude);

    // Calculate the final force and apply it
    force = Vec3Normalise( force );
    force = Vec3Scale( force, -magnitude );
    particle->addForce( force );
}
///<BungeePFG

ParticleFakeSpring::ParticleFakeSpring(vec3_t * anchor, real_t sc, real_t d)
: anchor(anchor), springConstant(sc), damping(d)
{
}

///>FakeSpringPFG
void ParticleFakeSpring::updateForce(Particle * particle, real_t duration)
{
    // Check that we do not have infinite mass
    if ( !particle->hasFiniteMass() ) {
        return;
    }

    // Calculate the relative position of the particle to the anchor
    vec3_t position = particle->getPosition();
    position = Vec3Subtract( position, *anchor );

    // Calculate the constants and check they are in bounds.
    real_t gamma = 0.5f * R_sqrt( 4 * springConstant - damping * damping );

    if ( gamma == 0.0f ) {
        return;
    }

    vec3_t c = Vec3Add( Vec3Scale( position, ( damping / ( 2.0f * gamma ) ) ), 
        Vec3Scale( particle->getVelocity(), ( 1.0f / gamma ) ) );

    // Calculate the target position
    vec3_t target = Vec3Add( Vec3Scale( position, R_cos( gamma * duration ) ), 
        Vec3Scale( c, R_sin( gamma * duration ) ) );

    target = Vec3Scale( target, R_exp( -0.5f * duration * damping ) );

    // Calculate the resulting acceleration and therefore the force
    vec3_t accel = Vec3Subtract( Vec3Scale( Vec3Subtract( target, position ),
        ( 1.0f / duration * duration ) ),  Vec3Scale( particle->getVelocity(),
        duration ) );
    particle->addForce( Vec3Scale( accel, particle->getMass() ) );
}
///<FakeSpringPFG

ParticleAnchoredSpring::ParticleAnchoredSpring(void)
{
}

ParticleAnchoredSpring::ParticleAnchoredSpring(vec3_t * anchor,
    real_t sc, real_t rl)
: anchor(anchor), springConstant(sc), restLength(rl)
{
}

void ParticleAnchoredSpring::init(vec3_t * anchor, real_t springConstant,
    real_t restLength)
{
	ParticleAnchoredSpring::anchor = anchor;
	ParticleAnchoredSpring::springConstant = springConstant;
	ParticleAnchoredSpring::restLength = restLength;
}

///>ASpringPFG
void ParticleAnchoredBungee::updateForce(Particle * particle, real_t duration)
{
	// Calculate the vector of the spring
	vec3_t force = particle->getPosition();
	force = Vec3Subtract( force, *anchor );

	// Calculate the magnitude of the force
	real_t magnitude = Vec3Magnitude( force );

	if ( magnitude < restLength ) {
        return;
    }

	magnitude = magnitude - restLength;
	magnitude *= springConstant;

	// Calculate the final force and apply it
	force = Vec3Normalise( force );
	force = Vec3Scale( force, -magnitude );
	particle->addForce( force );
}

void ParticleAnchoredSpring::updateForce(Particle * particle, real_t duration)
{
    // Calculate the vector of the spring
    vec3_t force = particle->getPosition();
    force = Vec3Subtract( force, *anchor );

    // Calculate the magnitude of the force
    real_t magnitude = Vec3Magnitude( force );
    magnitude = R_abs( magnitude - restLength );
    magnitude *= springConstant;

    // Calculate the final force and apply it
    force = Vec3Normalise( force );
    force = Vec3Scale( force, -magnitude );
    particle->addForce( force );
}


///<ASpringPFG