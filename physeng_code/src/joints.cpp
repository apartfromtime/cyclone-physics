/*
 * Implementation file for joints.
 *
 * Part of the Cyclone physics system.
 *
 * Copyright (c) Icosagon 2003. All Rights Reserved.
 *
 * This software is distributed under licence. Use of this software
 * implies agreement with all terms and conditions of the accompanying
 * software licence.
 */

#include <cyclone/cyclone.h>

using namespace cyclone;

unsigned Joint::addContact(Contact * contact, unsigned limit) const
{
	// Calculate the position of each connection point in world coordinates
	vec3_t a_pos_world = body[0]->getPointInWorldSpace( position[0] );
	vec3_t b_pos_world = body[1]->getPointInWorldSpace( position[1] );

	// Calculate the length of the joint
	vec3_t a_to_b = Vec3Subtract( b_pos_world, a_pos_world );
	vec3_t normal = a_to_b;
	normal = Vec3Normalise( normal );
	real_t length = Vec3Magnitude( a_to_b );

	// Check if it is violated
	if ( R_abs( length ) > error ) {

		contact->body[0] = body[0];
		contact->body[1] = body[1];
		contact->contactNormal = normal;
		contact->contactPoint = Vec3Scale(
			Vec3Add( a_pos_world, b_pos_world ), 0.5f );
		contact->penetration = length-error;
		contact->friction = 1.0f;
		contact->restitution = 0;
		
		return 1;
	}

	return 0;
}

void Joint::set(RigidBody * a, const vec3_t & a_pos, 
				RigidBody * b, const vec3_t & b_pos,
				real_t error)
{
	body[0] = a;
	body[1] = b;

	position[0] = a_pos;
	position[1] = b_pos;

	Joint::error = error;
}