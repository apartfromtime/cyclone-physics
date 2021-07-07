/*
 * Implementation file for the coarse collision detector.
 * 
 * Part of the Cyclone physics system.
 * 
 * Copyright (c) Icosagon 2003. All Rights Reserved.
 *
 * This software is distributed under licence. Use of this software
 * implies agreement with all terms and conditions of the accompanying
 * software licence.
 */


#include <cyclone/collide_coarse.h>

using namespace cyclone;

BoundingSphere cyclone::ConstructBoundingVolumeClass(vec3_t centre,
    real_t radius)
{
    BoundingSphere newSphere;

    newSphere.centre = centre;
    newSphere.radius = radius;
    
    return newSphere;
}

BoundingSphere cyclone::ConstructBoundingVolumeClass(BoundingSphere one,
    BoundingSphere two)
{
    BoundingSphere newSphere = ConstructBoundingVolumeClass( VECTOR3, 0.0f );

    vec3_t centreOffset = Vec3Subtract( two.centre, one.centre );
    real_t distance = Vec3MagnitudeSqr( centreOffset );
    real_t radiusDiff = two.radius - one.radius;

    // Check if the larger sphere encloses the small one
    if ( radiusDiff * radiusDiff >= distance ) {
        if ( one.radius > two.radius ) {

            newSphere.centre = one.centre;
            newSphere.radius = one.radius;
        } else {

            newSphere.centre = two.centre;
            newSphere.radius = two.radius;
        }
    }
    // Otherwise we need to work with partially 
    // overlapping spheres
    else {

        distance = R_sqrt( distance );
        newSphere.radius = ( distance + one.radius + two.radius ) *
            ( ( real_t )0.5f );

        // The new centre is based on one's centre, moved towards
        // two's centre by an ammount proportional to the spheres'
        // radii.
        newSphere.centre = one.centre;

        if ( distance > 0 ) {
            newSphere.centre = Vec3Add( newSphere.centre,
                Vec3Scale( centreOffset,
                ( ( newSphere.radius - one.radius ) / distance ) ) );
        }
    }

    return newSphere;
}

///>SphereBVHOverlap
bool cyclone::Overlaps(BoundingSphere sphere, BoundingSphere other)
{
    real_t distanceSquared = Vec3MagnitudeSqr( Vec3Subtract( sphere.centre,
        other.centre ) );
    return distanceSquared < ( sphere.radius + other.radius ) *
        ( sphere.radius + other.radius );
}
///<SphereBVHOverlap

real_t cyclone::GetGrowth(BoundingSphere sphere, BoundingSphere other)
{
    BoundingSphere newSphere = ConstructBoundingVolumeClass( sphere, other );

    // We return a value proportional to the change in surface
    // area of the sphere.
    return newSphere.radius * newSphere.radius - sphere.radius * sphere.radius;
}

real_t cyclone::GetSize(BoundingSphere sphere)
{
    return ( ( real_t )1.333333 ) * R_PI * sphere.radius * sphere.radius *
        sphere.radius;
}