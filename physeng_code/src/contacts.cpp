/*
 * Implementation file for the contact resolution system.
 * 
 * Part of the Cyclone physics system.
 * 
 * Copyright (c) Icosagon 2003. All Rights Reserved.
 *
 * This software is distributed under licence. Use of this software
 * implies agreement with all terms and conditions of the accompanying
 * software licence.
 */

///>ContactResolver
#include <cyclone/contacts.h>
///<ContactResolver
#include <memory.h>
#include <assert.h>

using namespace cyclone;

// Contact implementation

/* this has to change if the internal structure of the type changes */
const cyclone::Contact cyclone::CONTACT = {
    nullptr, nullptr, 0.0f, 0.0f, Vec3Clear(), Vec3Clear(), 0.0f,
    Mat3Identity(), Vec3Clear(), 0.0f, Vec3Clear()
};

void cyclone::SetBodyData(Contact & contact, RigidBody * one, RigidBody * two,
    real_t friction, real_t restitution)
{
    contact.body[0] = one;
    contact.body[1] = two;
    contact.friction = friction;
    contact.restitution = restitution;
}

///>ContactCalculateInternals
void cyclone::CalculateInternals(Contact & contact, real_t duration)
{
    // Check if the first object is NULL, and swap if it is.
    if ( !contact.body[0] ) {
        SwapBodies( contact );
    }

    assert( contact.body[0] );

    // Calculate an set of axis at the contact point.
    CalculateContactBasis( contact );

    // Store the relative position of the contact relative to each body
    contact.relativeContactPosition[0] = Vec3Subtract( contact.contactPoint,
        contact.body[0]->getPosition() );

    if ( contact.body[1] ) {

        contact.relativeContactPosition[1] = Vec3Subtract( contact.contactPoint,
            contact.body[1]->getPosition() );
    }

    // Find the relative velocity of the bodies at the contact point.
    contact.contactVelocity = CalculateLocalVelocity( contact, 0, duration );

    if ( contact.body[1] ) {
        contact.contactVelocity = Vec3Subtract( contact.contactVelocity,
            CalculateLocalVelocity( contact, 1, duration ) );
    }

    // Calculate the desired change in velocity for resolution
    CalculateDesiredDeltaVelocity( contact, duration );
}
///<ContactCalculateInternals

///>SwapBodies
/**
 * Swaps the bodies in the current contact, so body 0 is at body 1 and vice
 * versa. This also changes the direction of the contact normal, but doesn't
 * update any calculated internal data. If you are calling this method 
 * manually, then call calculateInternals afterwards to make sure the 
 * internal data is up to date.
 */
void cyclone::SwapBodies(Contact & contact)
{
    contact.contactNormal = Vec3Scale( contact.contactNormal, -1 );

    RigidBody * temp = contact.body[0];
    contact.body[0] = contact.body[1];
    contact.body[1] = temp;
}
///<SwapBodies


///>MatchAwakeState
void cyclone::MatchAwakeState(Contact & contact)
{
    // Collisions with the world never cause a body to wake up.
    if ( !contact.body[1] ) return;

    bool body0awake = contact.body[0]->getAwake();
    bool body1awake = contact.body[1]->getAwake();

    // Wake up only the sleeping one
    if ( body0awake ^ body1awake ) {
        if ( body0awake ) {
            contact.body[1]->setAwake();
        } else {
            contact.body[0]->setAwake();
        }
    }
}
///<MatchAwakeState

void cyclone::CalculateDesiredDeltaVelocity(Contact & contact, real_t duration)
{
    const static real_t velocityLimit = ( real_t )0.25f;

///>NewVelocityCalculation
    // Calculate the acceleration induced velocity accumulated this frame
    real_t velocityFromAcc = Vec3ScalarProduct(
        Vec3Scale( contact.body[0]->getLastFrameAcceleration(), duration ),
        contact.contactNormal );

    if ( contact.body[1] ) {
        velocityFromAcc -= Vec3ScalarProduct(
            Vec3Scale( contact.body[1]->getLastFrameAcceleration(), duration ),
            contact.contactNormal );
    }

    // If the velocity is very slow, limit the restitution
    real_t thisRestitution = contact.restitution;
    
    if ( R_abs( contact.contactVelocity.x ) < velocityLimit ) {
        thisRestitution = ( real_t )0.0f;
    }

    // Combine the bounce velocity with the removed
    // acceleration velocity.
    contact.desiredDeltaVelocity =
        -contact.contactVelocity.x
        -thisRestitution * (contact.contactVelocity.x - velocityFromAcc);
///<NewVelocityCalculation            
}

///>CalculateLocalVelocity;CalculateAllLocalVelocity
vec3_t cyclone::CalculateLocalVelocity(Contact & contact, unsigned bodyIndex,
    real_t duration)
{
    RigidBody * thisBody = contact.body[bodyIndex];

    // Work out the velocity of the contact point.
    vec3_t velocity = Vec3VectorProduct( thisBody->getRotation(),
        contact.relativeContactPosition[bodyIndex] );
    velocity = Vec3Add( velocity, thisBody->getVelocity() );

    // Turn the velocity into contact-coordinates.
    vec3_t contactVelocity = Mat3TransformTranspose( velocity,
        contact.contactToWorld );

///<CalculateLocalVelocity
    // Calculate the ammount of velocity that is due to forces without
    // reactions.
    vec3_t accVelocity = Vec3Scale( thisBody->getLastFrameAcceleration(),
        duration );

    // Calculate the velocity in contact-coordinates.
    accVelocity = Mat3TransformTranspose( accVelocity, contact.contactToWorld );

    // We ignore any component of acceleration in the contact normal 
    // direction, we are only interested in planar acceleration
    accVelocity.x = 0;

    // Add the planar velocities - if there's enough friction they will 
    // be removed during velocity resolution
    contactVelocity = Vec3Add( contactVelocity, accVelocity );
///>CalculateLocalVelocity

    // And return it
    return contactVelocity;
}
///<CalculateLocalVelocity;CalculateAllLocalVelocity

///>OrthonormalBasis
/**
 * Constructs an arbitrary orthonormal basis for the contact.
 * This is stored as a 3x3 matrix, where each vector is a column
 * (in other words the matrix transforms contact space into world
 * space). The x direction is generated from the contact normal,
 * and the y and z directionss are set so they are at right angles to
 * it.
 */
///<OrthonormalBasis
///>OrthonormalBasis
void cyclone::CalculateContactBasis(Contact & contact)
{
    vec3_t contactTangent[2];

    // Check whether the Z-axis is nearer to the X or Y axis
    if ( R_abs( contact.contactNormal.x ) > R_abs( contact.contactNormal.y ) ) { 
        
        // Scaling factor to ensure the results are normalised
        const real_t s = ( real_t )1.0f / R_sqrt( contact.contactNormal.z *
            contact.contactNormal.z + contact.contactNormal.x *
            contact.contactNormal.x );

        // The new X-axis is at right angles to the world Y-axis
        contactTangent[0].x = contact.contactNormal.z * s;
        contactTangent[0].y = 0;
        contactTangent[0].z = -contact.contactNormal.x * s;

        // The new Y-axis is at right angles to the new X- and Z- axes
        contactTangent[1].x = contact.contactNormal.y * contactTangent[0].x;
        contactTangent[1].y = contact.contactNormal.z * contactTangent[0].x - 
            contact.contactNormal.x * contactTangent[0].z;
        contactTangent[1].z = -contact.contactNormal.y * contactTangent[0].x;

    } else {
        
        // Scaling factor to ensure the results are normalised
        const real_t s = ( real_t )1.0 / R_sqrt( contact.contactNormal.z *
            contact.contactNormal.z + contact.contactNormal.y *
            contact.contactNormal.y );

        // The new X-axis is at right angles to the world X-axis
        contactTangent[0].x = 0;
        contactTangent[0].y = -contact.contactNormal.z*s;
        contactTangent[0].z = contact.contactNormal.y*s;

        // The new Y-axis is at right angles to the new X- and Z- axes
        contactTangent[1].x = contact.contactNormal.y * contactTangent[0].z - 
            contact.contactNormal.z * contactTangent[0].y;
        contactTangent[1].y = -contact.contactNormal.x * contactTangent[0].z;
        contactTangent[1].z = contact.contactNormal.x * contactTangent[0].y;
    }

    // Make a matrix from the three vectors.
    contact.contactToWorld = Mat3SetComponents( contact.contactNormal,
        contactTangent[0], contactTangent[1] );
}
///<OrthonormalBasis

void cyclone::ApplyVelocityChange(Contact & contact, vec3_t velocityChange[2],
    vec3_t rotationChange[2])
{
    // Get hold of the inverse mass and inverse inertia tensor, both in
    // world coordinates.
    mat3_t inverseInertiaTensor[2];

    contact.body[0]->getInverseInertiaTensorWorld( &inverseInertiaTensor[0] );

    // We will calculate the impulse for each contact axis
    vec3_t impulseContact;

    // Use the short format for frictionless contacts
    if ( contact.friction == ( real_t )0.0f ) {

///>DeltaVelPerImpulse
        // Build a vector that shows the change in velocity in
        // world space for a unit impulse in the direction of the contact 
        // normal.
        vec3_t deltaVelWorld = Vec3VectorProduct(
            contact.relativeContactPosition[0], contact.contactNormal );
        deltaVelWorld = Mat3Transform( deltaVelWorld, inverseInertiaTensor[0] );
        deltaVelWorld = Vec3VectorProduct( deltaVelWorld,
            contact.relativeContactPosition[0] );

        // Work out the change in velocity in contact coordiantes.
        real_t deltaVelocity = Vec3ScalarProduct( deltaVelWorld,
            contact.contactNormal );

        // Add the linear component of velocity change
        deltaVelocity += contact.body[0]->getInverseMass();

        // Check if we need the second body's data
        if ( contact.body[1] ) {

            // Find the inertia tensor for this body
            contact.body[1]->getInverseInertiaTensorWorld(
                &inverseInertiaTensor[1] );

            // Go through the same transformation sequence again
            vec3_t deltaVelWorld = Vec3VectorProduct(
                contact.relativeContactPosition[1], contact.contactNormal );
            deltaVelWorld = Mat3Transform( deltaVelWorld,
                inverseInertiaTensor[1] );
            deltaVelWorld = Vec3VectorProduct( deltaVelWorld,
                contact.relativeContactPosition[1] );

            // Add the change in velocity due to rotation
            deltaVelocity += Vec3ScalarProduct( deltaVelWorld,
                contact.contactNormal );

            // Add the change in velocity due to linear motion
            deltaVelocity += contact.body[1]->getInverseMass();
        }
///<DeltaVelPerImpulse

///>ImpulseContact
        // Calculate the required size of the impulse
        impulseContact.x = contact.desiredDeltaVelocity / deltaVelocity;
        impulseContact.y = 0;
        impulseContact.z = 0;
///<ImpulseContact
    }

    // Otherwise we may have impulses that aren't in the direction of the 
    // contact, so we need to build complete matrices.
    else
    {
///>UpdateVelWithFriction
        real_t inverseMass = contact.body[0]->getInverseMass();

///>AngularVelocityFriction
        // The equivalent of a cross product in matrices is multiplication
        // by a skew symmetric matrix - we build the matrix for converting
        // between linear and angular quantities.
        mat3_t impulseToTorque = Mat3SetSkewSymmetric(
            contact.relativeContactPosition[0] );

        // Build the matrix to convert contact impulse to change in velocity
        // in world coordinates.
        mat3_t deltaVelWorld = impulseToTorque; 
        deltaVelWorld = Mat3Multiply( deltaVelWorld, inverseInertiaTensor[0] );
        deltaVelWorld = Mat3Multiply( deltaVelWorld, impulseToTorque );
        deltaVelWorld = Mat3Scale( deltaVelWorld, -1 );

        // Check if we need to add body two's data
        if ( contact.body[1] ) {

            // Find the inertia tensor for this body
            contact.body[1]->getInverseInertiaTensorWorld(
                &inverseInertiaTensor[1] );

            // Set the cross product matrix
            impulseToTorque =
                Mat3SetSkewSymmetric( contact.relativeContactPosition[1] );

            // Calculate the velocity change matrix
            mat3_t deltaVelWorld2 = impulseToTorque; 
            deltaVelWorld2 = Mat3Multiply( deltaVelWorld2,
                inverseInertiaTensor[1] );
            deltaVelWorld2 = Mat3Multiply( deltaVelWorld2, impulseToTorque );
            deltaVelWorld2 = Mat3Scale( deltaVelWorld2, -1 );

            // Add to the total delta velocity.
            deltaVelWorld = Mat3Add( deltaVelWorld, deltaVelWorld2 );
///<AngularVelocityFriction

            // Add to the inverse mass
            inverseMass += contact.body[1]->getInverseMass();
///>AngularVelocityFriction
        }

        // Do a change of basis to convert into contact coordinates.
        mat3_t deltaVelocity = Mat3Transpose( contact.contactToWorld );
        deltaVelocity = Mat3Multiply( deltaVelocity, deltaVelWorld );
        deltaVelocity = Mat3Multiply( deltaVelocity, contact.contactToWorld );
///<AngularVelocityFriction

        // Add in the linear velocity change
        deltaVelocity.n[0] += inverseMass;
        deltaVelocity.n[4] += inverseMass;
        deltaVelocity.n[8] += inverseMass;

        // Invert to get the impulse needed per unit velocity
        mat3_t impulseMatrix = Mat3Inverse( deltaVelocity );

        // Find the target velocities to kill
        vec3_t velKill = {
            contact.desiredDeltaVelocity, -contact.contactVelocity.y,
            -contact.contactVelocity.z
        };

        // Find the impulse to kill target velocities
        impulseContact = Mat3Transform( velKill, impulseMatrix );

        // Check for exceeding friction
        real_t planarImpulse = R_sqrt( impulseContact.y * impulseContact.y +
            impulseContact.z * impulseContact.z
        );

        if ( planarImpulse > impulseContact.x * contact.friction ) {

            // We need to use dynamic friction
            impulseContact.y /= planarImpulse;
            impulseContact.z /= planarImpulse;

            impulseContact.x = deltaVelocity.n[0] +
                deltaVelocity.n[1] * contact.friction * impulseContact.y +
                deltaVelocity.n[2] * contact.friction * impulseContact.z;
            impulseContact.x = contact.desiredDeltaVelocity / impulseContact.x;
            impulseContact.y *= contact.friction * impulseContact.x;
            impulseContact.z *= contact.friction * impulseContact.x;
        }
///<UpdateVelWithFriction
    }    

///>ImpulseToWorld
    // Convert impulse to world coordinates
    vec3_t impulse = Mat3Transform( impulseContact, contact.contactToWorld );
///<ImpulseToWorld

    // Split the impulse into linear and rotational components
    vec3_t impulsiveTorque = Vec3VectorProduct(
        contact.relativeContactPosition[0], impulse );
    rotationChange[0] = Mat3Transform( impulsiveTorque,
        inverseInertiaTensor[0] );
    velocityChange[0] = Vec3Clear();
    velocityChange[0] = Vec3Add( velocityChange[0],
        Vec3Scale( impulse, contact.body[0]->getInverseMass() ) );

    // Apply the changes
    contact.body[0]->addVelocity( velocityChange[0] );
    contact.body[0]->addRotation( rotationChange[0] );

    if ( contact.body[1] ) {

        // Work out body one's linear and angular changes
        vec3_t impulsiveTorque = Vec3VectorProduct( impulse,
            contact.relativeContactPosition[1] );
        rotationChange[1] = Mat3Transform( impulsiveTorque,
            inverseInertiaTensor[1] );
        velocityChange[1] = Vec3Clear();
        velocityChange[1] = Vec3Add( velocityChange[1],
            Vec3Scale( impulse, -contact.body[1]->getInverseMass() ) );

        // And apply them.
        contact.body[1]->addVelocity( velocityChange[1] );
        contact.body[1]->addRotation( rotationChange[1] );
    }
}

void cyclone::ApplyPositionChange(Contact & contact, vec3_t velocityChange[2],
    vec3_t rotationDirection[2], real_t rotationAmount[2], real_t penetration)
{
    real_t angularLimit = (real_t)1000;//0.1f;
    real_t angularMove[2],linearMove[2];
    int b;

    real_t totalInertia = 0;
    real_t linearInertia[2];
    real_t angularInertia[2];

///>AngularInertia
    // We need to work out the inertia of each object in the direction
    // of the contact normal, due to angular inertia only.
    for (unsigned i = 0; i < 2; i++) {
        if ( contact.body[i] ) {

            mat3_t inverseInertiaTensor;
            contact.body[i]->getInverseInertiaTensorWorld(
                &inverseInertiaTensor );

            // Use the same procedure as for calculating frictionless
            // velocity change to work out the angular inertia.
            vec3_t angularInertiaWorld = Vec3VectorProduct(
                contact.relativeContactPosition[i], contact.contactNormal );
            angularInertiaWorld = Mat3Transform( angularInertiaWorld,
                inverseInertiaTensor );
            angularInertiaWorld = Vec3VectorProduct( angularInertiaWorld,
                contact.relativeContactPosition[i] );
            angularInertia[i] = Vec3ScalarProduct( angularInertiaWorld,
                contact.contactNormal );

            // The linear component is simply the inverse mass
            linearInertia[i] = contact.body[i]->getInverseMass();

            // Keep track of the total inertia from all components
            totalInertia += linearInertia[i] + angularInertia[i];
        }
    }
///<AngularInertia

    real_t inverseMass[2];

    totalInertia = angularInertia[0] + contact.body[0]->getInverseMass();

    if( contact.body[1] ) {

        inverseMass[1] = angularInertia[1] + contact.body[1]->getInverseMass();
        totalInertia += inverseMass[1];

        angularMove[1] = -penetration * angularInertia[1] / totalInertia;
        linearMove[1] = -penetration * contact.body[1]->getInverseMass() /
            totalInertia;

        // To avoid angular projections that are too great (when mass is large
        // but inertia tensor is small) limit the angular move.
        vec3_t projection = contact.relativeContactPosition[1];
        projection = Vec3Add( projection, Vec3Scale( contact.contactNormal,
            -Vec3ScalarProduct( contact.relativeContactPosition[1],
                contact.contactNormal ) ) );
        real_t max = angularLimit * Vec3Magnitude(
            contact.relativeContactPosition[0] );

        if( R_abs( angularMove[1] ) > max ) {

            real_t pp = angularMove[1] + linearMove[1];
            angularMove[1] = angularMove[1] > 0 ? max: -max;
            linearMove[1] = pp - angularMove[1];
        }
    }

    angularMove[0] = penetration * angularInertia[0] / totalInertia;
    linearMove[0] = penetration * contact.body[0]->getInverseMass() /
        totalInertia;

    // To avoid angular projections that are too great (when mass is large
    // but inertia tensor is small) limit the angular move.
    vec3_t projection = contact.relativeContactPosition[0];
    projection = Vec3Add( projection, Vec3Scale( contact.contactNormal,
        -Vec3ScalarProduct( contact.relativeContactPosition[0],
            contact.contactNormal ) ) );
    real_t max = angularLimit * Vec3Magnitude(
        contact.relativeContactPosition[0] );

    if ( R_abs( angularMove[0] ) > max ) {

        real_t pp = angularMove[0] + linearMove[0];
        angularMove[0] = angularMove[0] > 0 ? max: -max;
        linearMove[0] = pp - angularMove[0];
    }

    for(b = 0; b < 2; b++) {
        
        if ( contact.body[b] ) {
            
            if ( angularMove[b] != ( ( real_t )0.0f ) ) {

                vec3_t t = Vec3VectorProduct(
                    contact.relativeContactPosition[b], contact.contactNormal );

                mat3_t inverseInertiaTensor;
                contact.body[b]->getInverseInertiaTensorWorld(
                    &inverseInertiaTensor );

                rotationDirection[b] = Mat3Transform( t, inverseInertiaTensor );
                rotationAmount[b] = angularMove[b] / angularInertia[b];

                assert( rotationAmount[b] != ( ( real_t )0.0f ) );

            } else {
                
                rotationDirection[b] = Vec3Clear();
                rotationAmount[b] = 1;
            }

            velocityChange[b] = contact.contactNormal;
            velocityChange[b] = Vec3Scale( velocityChange[b],
                linearMove[b] / rotationAmount[b] );

            vec3_t pos = contact.body[b]->getPosition();
            pos = Vec3Add( pos, Vec3Scale( contact.contactNormal,
                linearMove[b] ) );
            contact.body[b]->setPosition( pos );

            quat_t q = QuatAddScaledVector( contact.body[b]->getOrientation(),
                rotationDirection[b], rotationAmount[b] * ( ( real_t )0.5f ) );
            contact.body[b]->setOrientation( q );
        }
    }
}


//==============================================================================
//
// Contact resolver implementation
//
//==============================================================================

const cyclone::ContactResolver cyclone::RESOLVER = {
    0,
    0,
    ( real_t )0.1,
    ( real_t )0.1,
    0,
    0,
    false
};

ContactResolver cyclone::ConstructContactResolver(unsigned iterations,
    real_t velocityEpsilon, real_t positionEpsilon)
{
    ContactResolver resolver = RESOLVER;

    resolver.velocityIterations = iterations;
    resolver.positionIterations = iterations;
    resolver.velocityEpsilon = velocityEpsilon;
    resolver.positionEpsilon = positionEpsilon;

    return resolver;
}

ContactResolver cyclone::ConstructContactResolver(unsigned velocityIterations,
    unsigned positionIterations, real_t velocityEpsilon, real_t positionEpsilon)
{
    ContactResolver resolver = RESOLVER;

    resolver.velocityIterations = velocityIterations;
    resolver.positionIterations = positionIterations;
    resolver.velocityEpsilon = velocityEpsilon;
    resolver.positionEpsilon = positionEpsilon;

    return resolver;
}

/** 
 * Returns true if the resolver has valid settings and is ready to go. 
 */
bool cyclone::IsValid(ContactResolver & resolver)
{
    return ( resolver.velocityIterations > 0  ) && 
           ( resolver.positionIterations > 0  ) &&
           ( resolver.positionEpsilon >= 0.0f ) && 
           ( resolver.positionEpsilon >= 0.0f );
}

///>ContactResolverBase
void cyclone::ResolveContacts(ContactResolver & resolver, Contact * contacts,
    unsigned numContacts, real_t duration)
{
    // Make sure we have something to do.
    if ( numContacts == 0 ) {
        return;
    }

///<ContactResolverBase
    if ( !IsValid( resolver ) ) {
        return;
    }
///>ContactResolverBase

    // Prepare the contacts for processing
    PrepareContacts( contacts, numContacts, duration );

    // Resolve the interpenetration problems with the contacts.
    AdjustPositions( resolver, contacts, numContacts, duration );

    // Resolve the velocity problems with the contacts.
    AdjustVelocities( resolver, contacts, numContacts, duration );
}
///<ContactResolverBase

///>PrepareContacts
void cyclone::PrepareContacts(Contact * contacts, unsigned numContacts,
    real_t duration)
{
    // Generate contact velocity and axis information.
    Contact * lastContact = contacts + numContacts;
    
    for(Contact * contact = contacts; contact < lastContact; contact++)
    {
        // Calculate the internal contact data (inertia, basis, etc).
        CalculateInternals( *contact, duration );
    } 
}
///<PrepareContacts

void cyclone::AdjustVelocities(ContactResolver & resolver, Contact * contact,
    unsigned numContacts, real_t duration)
{
    vec3_t velocityChange[2], rotationChange[2];
    vec3_t cp;

    // iteratively handle impacts in order of severity.
    resolver.velocityIterationsUsed = 0;
    
    while ( resolver.velocityIterationsUsed < resolver.velocityIterations ) {

        // Find contact with maximum magnitude of probable velocity change.
        real_t max = resolver.velocityEpsilon;
        unsigned index = numContacts;
        
        for(unsigned i = 0; i < numContacts; i++)
        {
            if ( contact[i].desiredDeltaVelocity > max ) {
                
                max = contact[i].desiredDeltaVelocity;
                index = i;
            }
        }

        if ( index == numContacts ) {
            break;
        }

        // Match the awake state at the contact
        MatchAwakeState( contact[index] );

        // Do the resolution on the contact that came out top.
        ApplyVelocityChange( contact[index], velocityChange, rotationChange );

        // With the change in velocity of the two bodies, the update of 
        // contact velocities means that some of the relative closing 
        // velocities need recomputing.
        for (unsigned i = 0; i < numContacts; i++)
        {
            if ( contact[i].body[0] ) {

                if ( contact[i].body[0] == contact[index].body[0] ) {

                    cp = Vec3VectorProduct( rotationChange[0],
                        contact[i].relativeContactPosition[0] );

                    cp = Vec3Add( cp, velocityChange[0] );

                    contact[i].contactVelocity = Vec3Add(
                        contact[i].contactVelocity, Mat3TransformTranspose( cp,
                            contact[i].contactToWorld ) );

                    CalculateDesiredDeltaVelocity( contact[i], duration );

                } else if ( contact[i].body[0] == contact[index].body[1] ) {

                    cp = Vec3VectorProduct( rotationChange[1],
                        contact[i].relativeContactPosition[0] );

                    cp = Vec3Add( cp, velocityChange[1] );

                    contact[i].contactVelocity = Vec3Add(
                        contact[i].contactVelocity, Mat3TransformTranspose( cp,
                            contact[i].contactToWorld ) );

                    CalculateDesiredDeltaVelocity( contact[i], duration );
                }
            }

            if ( contact[i].body[1] ) {

                if ( contact[i].body[1] == contact[index].body[0] ) {

                    cp = Vec3VectorProduct( rotationChange[0],
                        contact[i].relativeContactPosition[1] );

                    cp = Vec3Add( cp, velocityChange[0] );

                    contact[i].contactVelocity = Vec3Subtract(
                        contact[i].contactVelocity, Mat3TransformTranspose( cp,
                            contact[i].contactToWorld ) );

                    CalculateDesiredDeltaVelocity( contact[i], duration );

                } else if ( contact[i].body[1] == contact[index].body[1] ) {

                    cp = Vec3VectorProduct( rotationChange[1],
                        contact[i].relativeContactPosition[1] );

                    cp = Vec3Add( cp, velocityChange[1] );

                    contact[i].contactVelocity = Vec3Subtract(
                        contact[i].contactVelocity, Mat3TransformTranspose( cp,
                            contact[i].contactToWorld ) );

                    CalculateDesiredDeltaVelocity( contact[i], duration );
                }
            }
        }

        resolver.velocityIterationsUsed++;
    }
}

///>AdjustPositions
void cyclone::AdjustPositions(ContactResolver & resolver, Contact * contact,
    unsigned numContacts, real_t duration)
{
    unsigned i, index;
    vec3_t velocityChange[2], rotationChange[2];
    real_t rotationAmount[2];
    real_t max;
    vec3_t cp;

    // iteratively resolve interpenetration in order of severity.
    resolver.positionIterationsUsed = 0;
    
    while ( resolver.positionIterationsUsed < resolver.positionIterations ) {

        // Find biggest penetration
        max = resolver.positionEpsilon;
        index = numContacts;
        
        for(i = 0; i < numContacts; i++) {

            if( contact[i].penetration > max ) {

                max = contact[i].penetration;
                index = i;
            }
        }

        if ( index == numContacts ) {
            break;
        }

        // Match the awake state at the contact
        MatchAwakeState( contact[index] );

        // Resolve the penetration.
        ApplyPositionChange( contact[index], velocityChange, rotationChange,
            rotationAmount, max );//-positionEpsilon);

        // Again this action may have changed the penetration of other 
        // bodies, so we update contacts.
        for(i = 0; i < numContacts; i++)
        {
            if ( contact[i].body[0] == contact[index].body[0] ) {

                cp = Vec3VectorProduct( rotationChange[0],
                    contact[i].relativeContactPosition[0] );

                cp = Vec3Add( cp, velocityChange[0] );

                contact[i].penetration -= rotationAmount[0] *
                    Vec3ScalarProduct( cp, contact[i].contactNormal );

            } else if( contact[i].body[0] == contact[index].body[1] )
            {
                cp = Vec3VectorProduct( rotationChange[1],
                    contact[i].relativeContactPosition[0] );

                cp = Vec3Add( cp, velocityChange[1] );

                contact[i].penetration -= rotationAmount[1] *
                    Vec3ScalarProduct( cp, contact[i].contactNormal );
            }

            if ( contact[i].body[1] ) {

                if( contact[i].body[1] == contact[index].body[0] ) {

                    cp = Vec3VectorProduct( rotationChange[0],
                        contact[i].relativeContactPosition[1] );

                    cp = Vec3Add( cp, velocityChange[0] );

                    contact[i].penetration += rotationAmount[0] *
                        Vec3ScalarProduct( cp, contact[i].contactNormal );

                } else if( contact[i].body[1] == contact[index].body[1] ) {

                    cp = Vec3VectorProduct( rotationChange[1],
                        contact[i].relativeContactPosition[1] );

                    cp = Vec3Add( cp, velocityChange[1] );

                    contact[i].penetration += rotationAmount[1] *
                        Vec3ScalarProduct( cp, contact[i].contactNormal );
                }
            }
        }

        resolver.positionIterationsUsed++;
    }
}
///<AdjustPositions
