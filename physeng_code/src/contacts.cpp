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

void Contact::setBodyData(RigidBody* one, RigidBody *two,
                          real_t friction, real_t restitution)
{
    Contact::body[0] = one;
    Contact::body[1] = two;
    Contact::friction = friction;
    Contact::restitution = restitution;
}

///>MatchAwakeState
void Contact::matchAwakeState()
{
    // Collisions with the world never cause a body to wake up.
    if (!body[1]) return;

    bool body0awake = body[0]->getAwake();
    bool body1awake = body[1]->getAwake();

    // Wake up only the sleeping one
    if (body0awake ^ body1awake) {
        if (body0awake) body[1]->setAwake();
        else body[0]->setAwake();
    }
}
///<MatchAwakeState

///>SwapBodies
/**
 * Swaps the bodies in the current contact, so body 0 is at body 1 and vice
 * versa. This also changes the direction of the contact normal, but doesn't
 * update any calculated internal data. If you are calling this method 
 * manually, then call calculateInternals afterwards to make sure the 
 * internal data is up to date.
 */
void Contact::swapBodies()
{
    contactNormal = Vec3Scale( contactNormal, -1 );

    RigidBody *temp = body[0];
    body[0] = body[1];
    body[1] = temp;
}
///<SwapBodies

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
inline
///>OrthonormalBasis
void Contact::calculateContactBasis()
{
    vec3_t contactTangent[2];

    // Check whether the Z-axis is nearer to the X or Y axis
    if(R_abs(contactNormal.x) > R_abs(contactNormal.y))
    { 
        // Scaling factor to ensure the results are normalised
        const real_t s = (real_t)1.0f/R_sqrt(contactNormal.z*contactNormal.z +
            contactNormal.x*contactNormal.x);

        // The new X-axis is at right angles to the world Y-axis
        contactTangent[0].x = contactNormal.z*s;
        contactTangent[0].y = 0;
        contactTangent[0].z = -contactNormal.x*s;

        // The new Y-axis is at right angles to the new X- and Z- axes
        contactTangent[1].x = contactNormal.y*contactTangent[0].x;
        contactTangent[1].y = contactNormal.z*contactTangent[0].x - 
            contactNormal.x*contactTangent[0].z;
        contactTangent[1].z = -contactNormal.y*contactTangent[0].x;
    }
    else
    {
        // Scaling factor to ensure the results are normalised
        const real_t s = (real_t)1.0/R_sqrt(contactNormal.z*contactNormal.z + 
            contactNormal.y*contactNormal.y);

        // The new X-axis is at right angles to the world X-axis
        contactTangent[0].x = 0;
        contactTangent[0].y = -contactNormal.z*s;
        contactTangent[0].z = contactNormal.y*s;

        // The new Y-axis is at right angles to the new X- and Z- axes
        contactTangent[1].x = contactNormal.y*contactTangent[0].z - 
            contactNormal.z*contactTangent[0].y;
        contactTangent[1].y = -contactNormal.x*contactTangent[0].z;
        contactTangent[1].z = contactNormal.x*contactTangent[0].y;
    }

    // Make a matrix from the three vectors.
    contactToWorld = Mat3SetComponents( contactNormal, contactTangent[0],
        contactTangent[1] );
}
///<OrthonormalBasis

///>CalculateLocalVelocity;CalculateAllLocalVelocity
vec3_t Contact::calculateLocalVelocity(unsigned bodyIndex, real_t duration)
{
    RigidBody * thisBody = body[bodyIndex];

    // Work out the velocity of the contact point.
    vec3_t velocity = Vec3VectorProduct( thisBody->getRotation(),
        relativeContactPosition[bodyIndex] );
    velocity = Vec3Add( velocity, thisBody->getVelocity() );

    // Turn the velocity into contact-coordinates.
    vec3_t contactVelocity = Mat3TransformTranspose( velocity, contactToWorld );

///<CalculateLocalVelocity
    // Calculate the ammount of velocity that is due to forces without
    // reactions.
    vec3_t accVelocity = Vec3Scale( thisBody->getLastFrameAcceleration(),
        duration );

    // Calculate the velocity in contact-coordinates.
    accVelocity = Mat3TransformTranspose( accVelocity, contactToWorld );

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


void Contact::calculateDesiredDeltaVelocity(real_t duration)
{
    const static real_t velocityLimit = ( real_t )0.25f;

///>NewVelocityCalculation
    // Calculate the acceleration induced velocity accumulated this frame
    real_t velocityFromAcc = Vec3ScalarProduct(
        Vec3Scale( body[0]->getLastFrameAcceleration(), duration ),
        contactNormal );

    if ( body[1] ) {
        velocityFromAcc -= Vec3ScalarProduct(
            Vec3Scale( body[1]->getLastFrameAcceleration(), duration ),
            contactNormal );
    }

    // If the velocity is very slow, limit the restitution
    real_t thisRestitution = restitution;
    
    if ( R_abs( contactVelocity.x ) < velocityLimit ) {
        thisRestitution = ( real_t )0.0f;
    }

    // Combine the bounce velocity with the removed
    // acceleration velocity.
    desiredDeltaVelocity =
        -contactVelocity.x
        -thisRestitution * (contactVelocity.x - velocityFromAcc);
///<NewVelocityCalculation            
}


///>ContactCalculateInternals
void Contact::calculateInternals(real_t duration)
{
    // Check if the first object is NULL, and swap if it is.
    if (!body[0]) swapBodies();
    assert(body[0]);

    // Calculate an set of axis at the contact point.
    calculateContactBasis();

    // Store the relative position of the contact relative to each body
    relativeContactPosition[0] = Vec3Subtract( contactPoint,
        body[0]->getPosition() );

    if ( body[1] ) {

        relativeContactPosition[1] = Vec3Subtract( contactPoint,
            body[1]->getPosition() );
    }

    // Find the relative velocity of the bodies at the contact point.
    contactVelocity = calculateLocalVelocity( 0, duration );

    if ( body[1] ) {
        contactVelocity = Vec3Subtract( contactVelocity,
            calculateLocalVelocity( 1, duration ) );
    }

    // Calculate the desired change in velocity for resolution
    calculateDesiredDeltaVelocity( duration );
}
///<ContactCalculateInternals

void Contact::applyVelocityChange(vec3_t velocityChange[2],
    vec3_t rotationChange[2])
{
    // Get hold of the inverse mass and inverse inertia tensor, both in
    // world coordinates.
    mat3_t inverseInertiaTensor[2];

    body[0]->getInverseInertiaTensorWorld( &inverseInertiaTensor[0] );

    // We will calculate the impulse for each contact axis
    vec3_t impulseContact;

    // Use the short format for frictionless contacts
    if ( friction == ( real_t )0.0f ) {

///>DeltaVelPerImpulse
        // Build a vector that shows the change in velocity in
        // world space for a unit impulse in the direction of the contact 
        // normal.
        vec3_t deltaVelWorld = Vec3VectorProduct( relativeContactPosition[0],
            contactNormal );
        deltaVelWorld = Mat3Transform( deltaVelWorld, inverseInertiaTensor[0] );
        deltaVelWorld = Vec3VectorProduct( deltaVelWorld,
            relativeContactPosition[0] );

        // Work out the change in velocity in contact coordiantes.
        real_t deltaVelocity = Vec3ScalarProduct( deltaVelWorld,
            contactNormal );

        // Add the linear component of velocity change
        deltaVelocity += body[0]->getInverseMass();

        // Check if we need the second body's data
        if ( body[1] ) {

            // Find the inertia tensor for this body
            body[1]->getInverseInertiaTensorWorld( &inverseInertiaTensor[1] );

            // Go through the same transformation sequence again
            vec3_t deltaVelWorld = Vec3VectorProduct(
                relativeContactPosition[1], contactNormal );
            deltaVelWorld = Mat3Transform( deltaVelWorld,
                inverseInertiaTensor[1] );
            deltaVelWorld = Vec3VectorProduct( deltaVelWorld,
                relativeContactPosition[1] );

            // Add the change in velocity due to rotation
            deltaVelocity += Vec3ScalarProduct( deltaVelWorld, contactNormal );

            // Add the change in velocity due to linear motion
            deltaVelocity += body[1]->getInverseMass();
        }
///<DeltaVelPerImpulse

///>ImpulseContact
        // Calculate the required size of the impulse
        impulseContact.x = desiredDeltaVelocity / deltaVelocity;
        impulseContact.y = 0;
        impulseContact.z = 0;
///<ImpulseContact
    }

    // Otherwise we may have impulses that aren't in the direction of the 
    // contact, so we need to build complete matrices.
    else
    {
///>UpdateVelWithFriction
        real_t inverseMass = body[0]->getInverseMass();

///>AngularVelocityFriction
        // The equivalent of a cross product in matrices is multiplication
        // by a skew symmetric matrix - we build the matrix for converting
        // between linear and angular quantities.
        mat3_t impulseToTorque = Mat3SetSkewSymmetric(
            relativeContactPosition[0] );

        // Build the matrix to convert contact impulse to change in velocity
        // in world coordinates.
        mat3_t deltaVelWorld = impulseToTorque; 
        deltaVelWorld = Mat3Multiply( deltaVelWorld, inverseInertiaTensor[0] );
        deltaVelWorld = Mat3Multiply( deltaVelWorld, impulseToTorque );
        deltaVelWorld = Mat3Scale( deltaVelWorld, -1 );

        // Check if we need to add body two's data
        if ( body[1] ) {

            // Find the inertia tensor for this body
            body[1]->getInverseInertiaTensorWorld( &inverseInertiaTensor[1] );

            // Set the cross product matrix
            impulseToTorque =
                Mat3SetSkewSymmetric( relativeContactPosition[1] );

            // Calculate the velocity change matrix
            mat3_t deltaVelWorld2 = impulseToTorque; 
            deltaVelWorld2 = Mat3Multiply( deltaVelWorld2, inverseInertiaTensor[1] );
            deltaVelWorld2 = Mat3Multiply( deltaVelWorld2, impulseToTorque );
            deltaVelWorld2 = Mat3Scale( deltaVelWorld2, -1 );

            // Add to the total delta velocity.
            deltaVelWorld = Mat3Add( deltaVelWorld, deltaVelWorld2 );
///<AngularVelocityFriction

            // Add to the inverse mass
            inverseMass += body[1]->getInverseMass();
///>AngularVelocityFriction
        }

        // Do a change of basis to convert into contact coordinates.
        mat3_t deltaVelocity = Mat3Transpose( contactToWorld );
        deltaVelocity = Mat3Multiply( deltaVelocity, deltaVelWorld );
        deltaVelocity = Mat3Multiply( deltaVelocity, contactToWorld );
///<AngularVelocityFriction

        // Add in the linear velocity change
        deltaVelocity.n[0] += inverseMass;
        deltaVelocity.n[4] += inverseMass;
        deltaVelocity.n[8] += inverseMass;

        // Invert to get the impulse needed per unit velocity
        mat3_t impulseMatrix = Mat3Inverse( deltaVelocity );

        // Find the target velocities to kill
        vec3_t velKill = {
            desiredDeltaVelocity, -contactVelocity.y, -contactVelocity.z
        };

        // Find the impulse to kill target velocities
        impulseContact = Mat3Transform( velKill, impulseMatrix );

        // Check for exceeding friction
        real_t planarImpulse = R_sqrt( impulseContact.y * impulseContact.y +
                                       impulseContact.z * impulseContact.z
        );

        if ( planarImpulse > impulseContact.x * friction ) {

            // We need to use dynamic friction
            impulseContact.y /= planarImpulse;
            impulseContact.z /= planarImpulse;

            impulseContact.x = deltaVelocity.n[0] +
                deltaVelocity.n[1] * friction * impulseContact.y +
                deltaVelocity.n[2] * friction * impulseContact.z;
            impulseContact.x = desiredDeltaVelocity / impulseContact.x;
            impulseContact.y *= friction * impulseContact.x;
            impulseContact.z *= friction * impulseContact.x;
        }
///<UpdateVelWithFriction
    }    

///>ImpulseToWorld
    // Convert impulse to world coordinates
    vec3_t impulse = Mat3Transform( impulseContact, contactToWorld );
///<ImpulseToWorld

    // Split the impulse into linear and rotational components
    vec3_t impulsiveTorque = Vec3VectorProduct( relativeContactPosition[0],
        impulse );
    rotationChange[0] = Mat3Transform( impulsiveTorque,
        inverseInertiaTensor[0] );
    velocityChange[0] = Vec3Clear();
    velocityChange[0] = Vec3Add( velocityChange[0],
        Vec3Scale( impulse, body[0]->getInverseMass() ) );

    // Apply the changes
    body[0]->addVelocity( velocityChange[0] );
    body[0]->addRotation( rotationChange[0] );

    if ( body[1] ) {

        // Work out body one's linear and angular changes
        vec3_t impulsiveTorque = Vec3VectorProduct( impulse,
            relativeContactPosition[1] );
        rotationChange[1] = Mat3Transform( impulsiveTorque,
            inverseInertiaTensor[1] );
        velocityChange[1] = Vec3Clear();
        velocityChange[1] = Vec3Add( velocityChange[1],
            Vec3Scale( impulse, -body[1]->getInverseMass() ) );

        // And apply them.
        body[1]->addVelocity( velocityChange[1] );
        body[1]->addRotation( rotationChange[1] );
    }
}

void Contact::applyPositionChange(vec3_t velocityChange[2],
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
        if ( body[i] ) {

            mat3_t inverseInertiaTensor;
            body[i]->getInverseInertiaTensorWorld( &inverseInertiaTensor );

            // Use the same procedure as for calculating frictionless
            // velocity change to work out the angular inertia.
            vec3_t angularInertiaWorld = Vec3VectorProduct(
                relativeContactPosition[i], contactNormal );
            angularInertiaWorld = Mat3Transform( angularInertiaWorld,
                inverseInertiaTensor );
            angularInertiaWorld = Vec3VectorProduct( angularInertiaWorld,
                relativeContactPosition[i] );
            angularInertia[i] = Vec3ScalarProduct( angularInertiaWorld,
                contactNormal );

            // The linear component is simply the inverse mass
            linearInertia[i] = body[i]->getInverseMass();

            // Keep track of the total inertia from all components
            totalInertia += linearInertia[i] + angularInertia[i];
        }
    }
///<AngularInertia

    real_t inverseMass[2];

    totalInertia = angularInertia[0] + body[0]->getInverseMass();

    if( body[1] ) {

        inverseMass[1] = angularInertia[1] + body[1]->getInverseMass();
        totalInertia += inverseMass[1];

        angularMove[1] = -penetration * angularInertia[1] / totalInertia;
        linearMove[1] = -penetration * body[1]->getInverseMass() / totalInertia;

        // To avoid angular projections that are too great (when mass is large
        // but inertia tensor is small) limit the angular move.
        vec3_t projection = relativeContactPosition[1];
        projection = Vec3Add( projection, Vec3Scale( contactNormal,
            -Vec3ScalarProduct( relativeContactPosition[1], contactNormal ) ) );
        real_t max = angularLimit * Vec3Magnitude( relativeContactPosition[0] );

        if( R_abs( angularMove[1] ) > max ) {

            real_t pp = angularMove[1] + linearMove[1];
            angularMove[1] = angularMove[1] > 0 ? max: -max;
            linearMove[1] = pp - angularMove[1];
        }
    }

    angularMove[0] = penetration * angularInertia[0] / totalInertia;
    linearMove[0] = penetration * body[0]->getInverseMass() / totalInertia;

    // To avoid angular projections that are too great (when mass is large
    // but inertia tensor is small) limit the angular move.
    vec3_t projection = relativeContactPosition[0];
    projection = Vec3Add( projection, Vec3Scale( contactNormal,
        -Vec3ScalarProduct( relativeContactPosition[0], contactNormal ) ) );
    real_t max = angularLimit * Vec3Magnitude( relativeContactPosition[0] );

    if ( R_abs( angularMove[0] ) > max ) {

        real_t pp = angularMove[0] + linearMove[0];
        angularMove[0] = angularMove[0] > 0 ? max: -max;
        linearMove[0] = pp - angularMove[0];
    }

    for(b = 0; b < 2; b++) {
        
        if ( body[b] ) {
            
            if ( angularMove[b] != ( ( real_t )0.0f ) ) {

                vec3_t t = Vec3VectorProduct( relativeContactPosition[b],
                    contactNormal );

                mat3_t inverseInertiaTensor;
                body[b]->getInverseInertiaTensorWorld( &inverseInertiaTensor );

                rotationDirection[b] = Mat3Transform( t, inverseInertiaTensor );
                rotationAmount[b] = angularMove[b] / angularInertia[b];

                assert( rotationAmount[b] != ( ( real_t )0.0f ) );

            } else {
                
                rotationDirection[b] = Vec3Clear();
                rotationAmount[b] = 1;
            }

            velocityChange[b] = contactNormal;
            velocityChange[b] = Vec3Scale( velocityChange[b],
                linearMove[b] / rotationAmount[b] );

            vec3_t pos = body[b]->getPosition();
            pos = Vec3Add( pos, Vec3Scale( contactNormal, linearMove[b] ) );
            body[b]->setPosition( pos );

            quat_t q = QuatAddScaledVector( body[b]->getOrientation(),
                rotationDirection[b], rotationAmount[b] * ( ( real_t )0.5f ) );
            body[b]->setOrientation( q );
        }
    }
}


//==============================================================================
//
// Contact resolver implementation
//
//==============================================================================


ContactResolver::ContactResolver(unsigned iterations, real_t velocityEpsilon,
    real_t positionEpsilon)
{
    setIterations(iterations, iterations);
    setEpsilon(velocityEpsilon, positionEpsilon);
}

ContactResolver::ContactResolver(unsigned velocityIterations,
    unsigned positionIterations, real_t velocityEpsilon, real_t positionEpsilon)
{
    setIterations(velocityIterations);
    setEpsilon(velocityEpsilon, positionEpsilon);
}

void ContactResolver::setIterations(unsigned iterations)
{
    setIterations(iterations, iterations);
}

void ContactResolver::setIterations(unsigned velocityIterations,
    unsigned positionIterations)
{
    ContactResolver::velocityIterations = velocityIterations;
    ContactResolver::positionIterations = positionIterations;
}

void ContactResolver::setEpsilon(real_t velocityEpsilon,
    real_t positionEpsilon)
{
    ContactResolver::velocityEpsilon = velocityEpsilon;
    ContactResolver::positionEpsilon = positionEpsilon;
}

///>ContactResolverBase
void ContactResolver::resolveContacts(Contact * contacts,
    unsigned numContacts, real_t duration)
{
    // Make sure we have something to do.
    if ( numContacts == 0 ) {
        return;
    }

///<ContactResolverBase
    if ( !isValid() ) {
        return;
    }
///>ContactResolverBase

    // Prepare the contacts for processing
    prepareContacts(contacts, numContacts, duration);

    // Resolve the interpenetration problems with the contacts.
    adjustPositions(contacts, numContacts, duration);

    // Resolve the velocity problems with the contacts.
    adjustVelocities(contacts, numContacts, duration);
}
///<ContactResolverBase

///>PrepareContacts
void ContactResolver::prepareContacts(Contact * contacts, unsigned numContacts,
    real_t duration)
{
    // Generate contact velocity and axis information.
    Contact * lastContact = contacts + numContacts;
    
    for(Contact * contact = contacts; contact < lastContact; contact++)
    {
        // Calculate the internal contact data (inertia, basis, etc).
        contact->calculateInternals( duration );
    } 
}
///<PrepareContacts

void ContactResolver::adjustVelocities(Contact * c, unsigned numContacts,
    real_t duration)
{
    vec3_t velocityChange[2], rotationChange[2];
    vec3_t cp;

    // iteratively handle impacts in order of severity.
    velocityIterationsUsed = 0;
    
    while ( velocityIterationsUsed < velocityIterations ) {

        // Find contact with maximum magnitude of probable velocity change.
        real_t max = velocityEpsilon;
        unsigned index = numContacts;
        
        for(unsigned i = 0; i < numContacts; i++)
        {
            if ( c[i].desiredDeltaVelocity > max ) {
                
                max = c[i].desiredDeltaVelocity;
                index = i;
            }
        }

        if ( index == numContacts ) {
            break;
        }

        // Match the awake state at the contact
        c[index].matchAwakeState();

        // Do the resolution on the contact that came out top.
        c[index].applyVelocityChange( velocityChange, rotationChange );

        // With the change in velocity of the two bodies, the update of 
        // contact velocities means that some of the relative closing 
        // velocities need recomputing.
        for (unsigned i = 0; i < numContacts; i++)
        {
            if ( c[i].body[0] ) {

                if ( c[i].body[0] == c[index].body[0] ) {

                    cp = Vec3VectorProduct( rotationChange[0],
                        c[i].relativeContactPosition[0] );

                    cp = Vec3Add( cp, velocityChange[0] );

                    c[i].contactVelocity = Vec3Add( c[i].contactVelocity,
                        Mat3TransformTranspose( cp, c[i].contactToWorld ) );

                    c[i].calculateDesiredDeltaVelocity( duration );

                } else if ( c[i].body[0] == c[index].body[1] ) {

                    cp = Vec3VectorProduct( rotationChange[1],
                        c[i].relativeContactPosition[0] );

                    cp = Vec3Add( cp, velocityChange[1] );

                    c[i].contactVelocity = Vec3Add( c[i].contactVelocity,
                        Mat3TransformTranspose( cp, c[i].contactToWorld ) );

                    c[i].calculateDesiredDeltaVelocity( duration );
                }
            }

            if ( c[i].body[1] ) {

                if ( c[i].body[1] == c[index].body[0] ) {

                    cp = Vec3VectorProduct( rotationChange[0],
                        c[i].relativeContactPosition[1] );

                    cp = Vec3Add( cp, velocityChange[0] );

                    c[i].contactVelocity = Vec3Subtract( c[i].contactVelocity,
                        Mat3TransformTranspose( cp, c[i].contactToWorld ) );

                    c[i].calculateDesiredDeltaVelocity( duration );

                } else if ( c[i].body[1] == c[index].body[1] ) {

                    cp = Vec3VectorProduct( rotationChange[1],
                        c[i].relativeContactPosition[1] );

                    cp = Vec3Add( cp, velocityChange[1] );

                    c[i].contactVelocity = Vec3Subtract( c[i].contactVelocity,
                        Mat3TransformTranspose( cp, c[i].contactToWorld ) );

                    c[i].calculateDesiredDeltaVelocity( duration );
                }
            }
        }

        velocityIterationsUsed++;
    }
}

///>AdjustPositions
void ContactResolver::adjustPositions(Contact * c, unsigned numContacts,
    real_t duration)
{
    unsigned i,index;
    vec3_t velocityChange[2], rotationChange[2];
    real_t rotationAmount[2];
    real_t max;
    vec3_t cp;

    // TODO: Remove this
    for(i=0; i<numContacts; i++) c[i].calculateInternals(duration);

    // iteratively resolve interpenetration in order of severity.
    positionIterationsUsed = 0;
    while(positionIterationsUsed < positionIterations)
    {
        // Find biggest penetration
        max = positionEpsilon;
        index = numContacts;
        
        for(i = 0; i < numContacts; i++) {

            if( c[i].penetration > max ) {

                max = c[i].penetration;
                index = i;
            }
        }

        if ( index == numContacts ) {
            break;
        }

        // Match the awake state at the contact
        c[index].matchAwakeState();

        // Resolve the penetration.
        c[index].applyPositionChange( velocityChange, rotationChange,
            rotationAmount, max );//-positionEpsilon);

        // Again this action may have changed the penetration of other 
        // bodies, so we update contacts.
        for(i = 0; i < numContacts; i++)
        {
            if ( c[i].body[0] == c[index].body[0] ) {

                cp = Vec3VectorProduct( rotationChange[0],
                    c[i].relativeContactPosition[0] );

                cp = Vec3Add( cp, velocityChange[0] );

                c[i].penetration -= rotationAmount[0] *
                    Vec3ScalarProduct( cp, c[i].contactNormal );

            } else if( c[i].body[0] == c[index].body[1] )
            {
                cp = Vec3VectorProduct( rotationChange[1],
                    c[i].relativeContactPosition[0] );

                cp = Vec3Add( cp, velocityChange[1] );

                c[i].penetration -= rotationAmount[1] *
                    Vec3ScalarProduct( cp, c[i].contactNormal );
            }

            if ( c[i].body[1] ) {

                if( c[i].body[1] == c[index].body[0] ) {

                    cp = Vec3VectorProduct( rotationChange[0],
                        c[i].relativeContactPosition[1] );

                    cp = Vec3Add( cp, velocityChange[0] );

                    c[i].penetration += rotationAmount[0] *
                        Vec3ScalarProduct( cp, c[i].contactNormal );

                } else if( c[i].body[1] == c[index].body[1] ) {

                    cp = Vec3VectorProduct( rotationChange[1],
                        c[i].relativeContactPosition[1] );

                    cp = Vec3Add( cp, velocityChange[1] );

                    c[i].penetration += rotationAmount[1] *
                        Vec3ScalarProduct( cp, c[i].contactNormal );
                }
            }
        }

        positionIterationsUsed++;
    }
}
///<AdjustPositions


