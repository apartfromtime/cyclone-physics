/*
 * Implementation file for particle links.
 *
 * Part of the Cyclone physics system.
 *
 * Copyright (c) Icosagon 2003. All Rights Reserved.
 *
 * This software is distributed under licence. Use of this software
 * implies agreement with all terms and conditions of the accompanying
 * software licence.
 */

#include <cyclone/plinks.h>

using namespace cyclone;

real_t ParticleLink::currentLength(void) const
{
    vec3_t relativePos = Vec3Subtract( particle[0]->getPosition(),
        particle[1]->getPosition() );
    return Vec3Magnitude( relativePos );
}

unsigned ParticleCable::addContact(ParticleContact * contact,
    unsigned limit) const
{
    // Find the length of the cable
    real_t length = currentLength();

    // Check if we're over-extended
    if ( length < maxLength ) {
        return 0;
    }

    // Otherwise return the contact
    contact->particle[0] = particle[0];
    contact->particle[1] = particle[1];

    // Calculate the normal
    vec3_t normal = Vec3Subtract( particle[1]->getPosition(),
        particle[0]->getPosition() );
    normal = Vec3Normalise( normal );

    contact->contactNormal = normal;
    contact->penetration = length - maxLength;
    contact->restitution = restitution;

    return 1;
}

unsigned ParticleRod::addContact(ParticleContact * contact,
    unsigned limit) const
{
    // Find the length of the rod
    real_t currentLen = currentLength();

    // Check if we're over-extended
    if ( currentLen == length ) {
        return 0;
    }

    // Otherwise return the contact
    contact->particle[0] = particle[0];
    contact->particle[1] = particle[1];

    // Calculate the normal
    vec3_t normal = Vec3Subtract( particle[1]->getPosition(),
        particle[0]->getPosition() );
    normal = Vec3Normalise( normal );

    // The contact normal depends on whether we're extending or compressing
    if ( currentLen > length ) {

        contact->contactNormal = normal;
        contact->penetration = currentLen - length;

    } else {

        contact->contactNormal = Vec3Scale( normal, -1 );
        contact->penetration = length - currentLen;
    }

    // Always use zero restitution (no bounciness)
    contact->restitution = 0;

    return 1;
}

real_t ParticleConstraint::currentLength(void) const
{
	vec3_t relativePos = Vec3Subtract( particle->getPosition(), anchor );
	return Vec3Magnitude( relativePos );
}

unsigned ParticleCableConstraint::addContact(ParticleContact * contact,
    unsigned limit) const
{
	// Find the length of the cable
	real_t length = currentLength();

	// Check if we're over-extended
	if ( length < maxLength ) {
		return 0;
	}

	// Otherwise return the contact
	contact->particle[0] = particle;
	contact->particle[1] = 0;

	// Calculate the normal
	vec3_t normal = Vec3Subtract( anchor, particle->getPosition() );
	normal = Vec3Normalise( normal );

	contact->contactNormal = normal;
	contact->penetration = length - maxLength;
	contact->restitution = restitution;

	return 1;
}

unsigned ParticleRodConstraint::addContact(ParticleContact * contact,
    unsigned limit) const
{
	// Find the length of the rod
	real_t currentLen = currentLength();

	// Check if we're over-extended
	if ( currentLen == length ) {
		return 0;
	}

	// Otherwise return the contact
	contact->particle[0] = particle;
	contact->particle[1] = 0;

	// Calculate the normal
    vec3_t normal = Vec3Subtract( anchor, particle->getPosition() );
    normal = Vec3Normalise( normal );

	// The contact normal depends on whether we're extending or compressing
	if ( currentLen > length ) {

		contact->contactNormal = normal;
		contact->penetration = currentLen - length;

	} else {
		
        contact->contactNormal = Vec3Scale( normal, -1 );
		contact->penetration = length - currentLen;
	}

	// Always use zero restitution (no bounciness)
	contact->restitution = 0;

	return 1;
}
