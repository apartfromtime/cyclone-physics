/*
 * Implementation file for the fine grained collision detector.
 * 
 * Part of the Cyclone physics system.
 * 
 * Copyright (c) Icosagon 2003. All Rights Reserved.
 *
 * This software is distributed under licence. Use of this software
 * implies agreement with all terms and conditions of the accompanying
 * software licence.
 */


#include <cyclone/collide_fine.h>
#include <memory.h>
#include <assert.h>
#include <cstdlib>

#ifdef _DEBUG
/* For debugging */
#include <gl/glut.h>
#endif /* ifdef _DEBUG */

using namespace cyclone;

const cyclone::CollisionPrimitive cyclone::COLLISION_PRIMITIVE = {
    NULL, Mat4Identity(), Mat4Identity()
};

const cyclone::CollisionSphere cyclone::COLLISION_SPHERE = {
    COLLISION_PRIMITIVE, 0.0f
};

const cyclone::CollisionPlane cyclone::COLLISION_PLANE = {
    COLLISION_PRIMITIVE, VECTOR3, ( real_t )0.0f
};

const cyclone::CollisionBox cyclone::COLLISION_BOX = {
    COLLISION_PRIMITIVE, VECTOR3
};

/**
 * Calculates the internals for the primitive.
 */
void cyclone::CalculateInternals(CollisionPrimitive & primitive)
{
    primitive.transform = Mat4Multiply( primitive.body->getTransform(),
        primitive.offset );
}


/* Intersection Tests */


bool cyclone::SphereAndHalfSpace(const CollisionSphere & sphere,
    const CollisionPlane & plane)
{
    // Find the distance from the origin
    real_t ballDistance = Vec3ScalarProduct( plane.direction,
        Mat4AxisVector( sphere.primitive.transform, 3 ) ) - sphere.radius;

    // Check for the intersection
    return ballDistance <= plane.offset;
}

bool cyclone::SphereAndSphere(const CollisionSphere & one,
    const CollisionSphere & two)
{
    // Find the vector between the objects
    vec3_t midline = Vec3Subtract( Mat4AxisVector( one.primitive.transform, 3 ),
        Mat4AxisVector( two.primitive.transform, 3 ) );

    // See if it is large enough.
    return Vec3MagnitudeSqr( midline ) < ( one.radius + two.radius ) *
        ( one.radius + two.radius );
}

///>BoxBoxSepAxis
static inline real_t TransformToAxis(const CollisionBox & box,
    const vec3_t & axis)
{
    return ( box.halfSize.x * R_abs( Vec3ScalarProduct( axis,
            Mat4AxisVector( box.primitive.transform, 0 ) ) ) +
        box.halfSize.y * R_abs( Vec3ScalarProduct( axis,
            Mat4AxisVector( box.primitive.transform, 1 ) ) ) +
        box.halfSize.z * R_abs( Vec3ScalarProduct( axis,
            Mat4AxisVector( box.primitive.transform, 2 ) ) ) );
}

/**
 * This function checks if the two boxes overlap
 * along the given axis. The final parameter toCentre
 * is used to pass in the vector between the boxes centre
 * points, to avoid having to recalculate it each time.
 */
static inline bool OverlapOnAxis(const CollisionBox & one,
    const CollisionBox & two, const vec3_t & axis, const vec3_t & toCentre)
{
    // Project the half-size of one onto axis
    real_t oneProject = TransformToAxis( one, axis );
    real_t twoProject = TransformToAxis( two, axis );

    // Project this onto the axis
    real_t distance = R_abs( Vec3ScalarProduct( toCentre, axis ) );

    // Check for overlap
    return (distance < oneProject + twoProject);
}
///<BoxBoxSepAxis

// This preprocessor definition is only used as a convenience
// in the boxAndBox intersection  method.
#define TEST_OVERLAP(axis) OverlapOnAxis(one, two, (axis), toCentre)

bool cyclone::BoxAndBox(const CollisionBox & one, const CollisionBox & two)
{
    // Find the vector between the two centres
    vec3_t toCentre = Vec3Subtract( Mat4AxisVector( two.primitive.transform, 3 ),
        Mat4AxisVector( one.primitive.transform, 3 ) );

    return (
        // Check on box one's axes first
        TEST_OVERLAP(Mat4AxisVector( one.primitive.transform, 0 ) ) &&
        TEST_OVERLAP(Mat4AxisVector( one.primitive.transform, 1 ) ) &&
        TEST_OVERLAP(Mat4AxisVector( one.primitive.transform, 2 ) ) &&

        // And on two's
        TEST_OVERLAP(Mat4AxisVector( two.primitive.transform, 0 ) ) &&
        TEST_OVERLAP(Mat4AxisVector( two.primitive.transform, 1 ) ) &&
        TEST_OVERLAP(Mat4AxisVector( two.primitive.transform, 2 ) ) &&

        // Now on the cross products
        TEST_OVERLAP( Vec3VectorProduct(
            Mat4AxisVector( one.primitive.transform, 0 ),
            Mat4AxisVector( two.primitive.transform, 0 ) ) ) &&
        TEST_OVERLAP( Vec3VectorProduct(
            Mat4AxisVector( one.primitive.transform, 0 ),
            Mat4AxisVector( two.primitive.transform, 1 ) ) ) &&
        TEST_OVERLAP( Vec3VectorProduct(
            Mat4AxisVector( one.primitive.transform, 0 ),
            Mat4AxisVector( two.primitive.transform, 2 ) ) ) &&
        TEST_OVERLAP( Vec3VectorProduct(
            Mat4AxisVector( one.primitive.transform, 1 ),
            Mat4AxisVector( two.primitive.transform, 0 ) ) ) &&
        TEST_OVERLAP( Vec3VectorProduct(
            Mat4AxisVector( one.primitive.transform, 1 ),
            Mat4AxisVector( two.primitive.transform, 1 ) ) ) &&
        TEST_OVERLAP( Vec3VectorProduct(
            Mat4AxisVector( one.primitive.transform, 1 ),
            Mat4AxisVector( two.primitive.transform, 2 ) ) ) &&
        TEST_OVERLAP( Vec3VectorProduct(
            Mat4AxisVector( one.primitive.transform, 2 ),
            Mat4AxisVector( two.primitive.transform, 0 ) ) ) &&
        TEST_OVERLAP( Vec3VectorProduct(
            Mat4AxisVector( one.primitive.transform, 2 ),
            Mat4AxisVector( two.primitive.transform, 1 ) ) ) &&
        TEST_OVERLAP( Vec3VectorProduct(
            Mat4AxisVector( one.primitive.transform, 2 ),
            Mat4AxisVector( two.primitive.transform, 2 ) ) )
    );
}
#undef TEST_OVERLAP

bool cyclone::BoxAndHalfSpace(const CollisionBox & box, const CollisionPlane & plane)
{
    // Work out the projected radius of the box onto the plane direction
    real_t projectedRadius = TransformToAxis(box, plane.direction);

    // Work out how far the box is from the origin
    real_t boxDistance = Vec3ScalarProduct( plane.direction,
        Mat4AxisVector( box.primitive.transform, 3 ) ) - projectedRadius;

    // Check for the intersection
    return boxDistance <= plane.offset;
}


/* Collision Detection */


/**
 * Checks if there are more contacts available in the contact
 * data.
 */
bool cyclone::HasContacts(CollisionData & data)
{
    return data.contactsLeft > 0;
}

/**
 * Resets the data so that it has no used contacts recorded.
 */
void cyclone::Reset(CollisionData & data, unsigned int maxContacts)
{
    data.contactsLeft = maxContacts;
    data.contactCount = 0;
    data.contacts = data.contactArray;
}

/**
 * Notifies the data that the given number of contacts have
 * been added.
 */
void cyclone::Add(CollisionData & data, unsigned int count)
{
    // Reduce the number of contacts remaining, add number used and move
    // the array forward
    data.contactsLeft -= count;
    data.contactCount += count;
    data.contacts += count;
}


///>SphereHalfSpace
unsigned cyclone::SphereAndHalfSpace(const CollisionSphere & sphere,
    const CollisionPlane & plane, CollisionData * data)
{
    // Make sure we have contacts
    if ( data->contactsLeft <= 0 ) {
        return 0;
    }

    // Cache the sphere position
    vec3_t position = Mat4AxisVector( sphere.primitive.transform, 3 );
    
    // Find the distance from the plane
    real_t ballDistance = Vec3ScalarProduct( plane.direction, position ) -
        sphere.radius - plane.offset;
    
    if ( ballDistance >= 0 ) {
        return 0;
    }

    // Create the contact - it has a normal in the plane direction.
    Contact * contact = data->contacts;
    contact->contactNormal = plane.direction;
    contact->penetration = -ballDistance;
    contact->contactPoint = Vec3Subtract( position,
        Vec3Scale( plane.direction, ( ballDistance + sphere.radius ) ) );
    SetBodyData( *contact, sphere.primitive.body, NULL,
        data->friction, data->restitution);
    Add( *data, 1 );
    return 1;
}
///<SphereHalfSpace


///>SphereTruePlane
unsigned cyclone::SphereAndTruePlane(const CollisionSphere & sphere,
    const CollisionPlane & plane, CollisionData * data)
{
    // Make sure we have contacts
    if (data->contactsLeft <= 0) {
        return 0;
    }

    // Cache the sphere position
    vec3_t position = Mat4AxisVector( sphere.primitive.transform, 3 );

    // Find the distance from the plane
    real_t centreDistance = Vec3ScalarProduct( plane.direction, position ) -
        plane.offset;

    // Check if we're within radius
    if ( centreDistance * centreDistance > sphere.radius * sphere.radius ) {
        return 0;
    }

    // Check which side of the plane we're on
    vec3_t normal = plane.direction;
    real_t penetration = -centreDistance;

    if (centreDistance < 0) {

        normal = Vec3Scale( normal, -1 );
        penetration = -penetration;
    }

    penetration += sphere.radius;

    // Create the contact - it has a normal in the plane direction.
    Contact* contact = data->contacts;
    contact->contactNormal = normal;
    contact->penetration = penetration;
    contact->contactPoint = Vec3Subtract( position,
        Vec3Scale( plane.direction, centreDistance ) );
    SetBodyData( *contact, sphere.primitive.body, NULL,
        data->friction, data->restitution);

    Add( *data, 1 );
    return 1;
}
///<SphereTruePlane

///>SphereSphere
unsigned cyclone::SphereAndSphere(const CollisionSphere & one,
    const CollisionSphere & two, CollisionData * data)
{
    // Make sure we have contacts
    if ( data->contactsLeft <= 0 ) {
        return 0;
    }

    // Cache the sphere positions
    vec3_t positionOne = Mat4AxisVector( one.primitive.transform, 3 );
    vec3_t positionTwo = Mat4AxisVector( two.primitive.transform, 3 );

    // Find the vector between the objects
    vec3_t midline = Vec3Subtract( positionOne, positionTwo );
    real_t size = Vec3Magnitude( midline );

    // See if it is large enough.
    if ( size <= 0.0f || size >= one.radius + two.radius ) {
        return 0;
    }

    // We manually create the normal, because we have the
    // size to hand.
    vec3_t normal = Vec3Scale( midline, ( ( ( real_t )1.0f) / size ) );

    Contact * contact = data->contacts;
    contact->contactNormal = normal;
    contact->contactPoint = Vec3Add( positionOne,
        Vec3Scale( midline, ( real_t )0.5f ) );
    contact->penetration = ( one.radius + two.radius - size );
    SetBodyData( *contact, one.primitive.body, two.primitive.body, data->friction,
        data->restitution);

    Add( *data, 1 );

    return 1;
}
///<SphereSphere

///>BoxAndHalfSpace
unsigned cyclone::BoxAndHalfSpace(const CollisionBox & box,
    const CollisionPlane & plane, CollisionData * data)
{
    // Make sure we have contacts
    if ( data->contactsLeft <= 0 ) {
        return 0;
    }

    // Check for intersection
    if ( !BoxAndHalfSpace( box, plane ) ) {
        return 0;
    }

    // We have an intersection, so find the intersection points. We can make
    // do with only checking vertices. If the box is resting on a plane
    // or on an edge, it will be reported as four or two contact points.

    // Go through each combination of + and - for each half-size
    static real_t mults[8][3] = {
        { 1, 1, 1 }, {-1, 1, 1 }, { 1,-1, 1 }, {-1,-1, 1 },
        { 1, 1,-1 }, {-1, 1,-1 }, { 1,-1,-1 }, {-1,-1,-1 }
    };

    Contact * contact = data->contacts;
    unsigned contactsUsed = 0;

    for (unsigned i = 0; i < 8; i++) {

        // Calculate the position of each vertex
        vec3_t vertexPos = { mults[i][0], mults[i][1], mults[i][2] };
        vertexPos = Vec3ComponentProduct( vertexPos, box.halfSize );
        vertexPos = Mat4Transform( vertexPos, box.primitive.transform );

///>BoxPlaneTestOne
        // Calculate the distance from the plane
        real_t vertexDistance = Vec3ScalarProduct( vertexPos, plane.direction );

        // Compare this to the plane's distance
        if ( vertexDistance <= plane.offset ) {

            // Create the contact data.

            // The contact point is halfway between the vertex and the
            // plane - we multiply the direction by half the separation 
            // distance and add the vertex location.
            contact->contactPoint = Vec3Add(
                Vec3Scale( plane.direction, ( vertexDistance - plane.offset ) ),
                vertexPos );
            contact->contactNormal = plane.direction;
            contact->penetration = plane.offset - vertexDistance;
///<BoxPlaneTestOne
            
            // Write the appropriate data
            SetBodyData( *contact,  box.primitive.body, NULL,  data->friction,
                data->restitution );

            // Move onto the next contact
            contact++;
            contactsUsed++;
            
            if ( ( int )contactsUsed == data->contactsLeft ) {
                return contactsUsed;
            }

///>BoxPlaneTestOne
        }
///<BoxPlaneTestOne
    }

    Add( *data, contactsUsed );

    return contactsUsed;
}
///<BoxAndHalfSpace


/**
 * This function checks if the two boxes overlap
 * along the given axis, returning the amount of overlap. 
 * The final parameter toCentre
 * is used to pass in the vector between the boxes centre
 * points, to avoid having to recalculate it each time.
 */
static inline real_t PenetrationOnAxis(const CollisionBox & one,
    const CollisionBox & two, const vec3_t & axis, const vec3_t & toCentre)
{
    // Project the half-size of one onto axis
    real_t oneProject = TransformToAxis(one, axis);
    real_t twoProject = TransformToAxis(two, axis);

    // Project this onto the axis
    real_t distance = R_abs( Vec3ScalarProduct( toCentre, axis ) );

    // Return the overlap (i.e. positive indicates
    // overlap, negative indicates separation).
    return oneProject + twoProject - distance;
}

static inline bool TryAxis(const CollisionBox & one, const CollisionBox & two,
    const vec3_t & axis, const vec3_t & toCentre, unsigned index,
    // These values may be updated
    real_t & smallestPenetration, unsigned & smallestCase)
{
    real_t penetration = PenetrationOnAxis(one, two, axis, toCentre);
    if (penetration < 0) return false;
    if (penetration < smallestPenetration) {
        smallestPenetration = penetration;
        smallestCase = index;
    }
    return true;
}

// This preprocessor definition is only used as a convenience
// in the boxAndBox contact generation method.
#define CHECK_OVERLAP(axis, index) \
    if (!TryAxis(one, two, (axis), toCentre, (index), pen, best)) return 0;
   
void FillPointFaceBoxBox(const CollisionBox & one, const CollisionBox & two,
    const vec3_t & toCentre, CollisionData * data, unsigned best, real_t pen)
{
    // This method is called when we know that a vertex from
    // box two is in contact with box one.

    Contact * contact = data->contacts;

    // We know which axis the collision is on (i.e. best), 
    // but we need to work out which of the two faces on 
    // this axis.
    vec3_t normal = Mat4AxisVector( one.primitive.transform, best );

    if ( Vec3ScalarProduct( Mat4AxisVector( one.primitive.transform, best ), toCentre ) > 0) {
        normal = Vec3Scale( normal, -1.0f );
    }

    // Work out which vertex of box two we're colliding with.
    // Using toCentre doesn't work!
    vec3_t vertex = two.halfSize;
    
    if ( Vec3ScalarProduct( Mat4AxisVector( two.primitive.transform, 0 ), normal ) < 0 ) {
        vertex.x = -vertex.x;
    }

    if ( Vec3ScalarProduct( Mat4AxisVector( two.primitive.transform, 1 ), normal ) < 0 ) {
        vertex.y = -vertex.y;
    }

    if ( Vec3ScalarProduct( Mat4AxisVector( two.primitive.transform, 2 ), normal ) < 0 ) {
        vertex.z = -vertex.z;
    }
    
    // Create the contact data
    contact->contactNormal = normal;
    contact->penetration = pen;
    contact->contactPoint = Mat4Transform( vertex,
        two.primitive.transform );
    SetBodyData( *contact,  one.primitive.body, two.primitive.body,
        data->friction, data->restitution );
}

static inline vec3_t ContactPoint(const vec3_t & pOne, const vec3_t & dOne,
    real_t oneSize, const vec3_t & pTwo, const vec3_t & dTwo, real_t twoSize,
    // If this is true, and the contact point is outside
    // the edge (in the case of an edge-face contact) then
    // we use one's midpoint, otherwise we use two's.
    bool useOne)
{
    vec3_t toSt, cOne, cTwo;
    real_t dpStaOne, dpStaTwo, dpOneTwo, smOne, smTwo;
    real_t denom, mua, mub;

    smOne = Vec3MagnitudeSqr( dOne );
    smTwo = Vec3MagnitudeSqr( dTwo );
    dpOneTwo = Vec3ScalarProduct( dTwo, dOne );

    toSt = Vec3Subtract( pOne, pTwo );
    dpStaOne = Vec3ScalarProduct( dOne, toSt );
    dpStaTwo = Vec3ScalarProduct( dTwo, toSt );

    denom = smOne * smTwo - dpOneTwo * dpOneTwo;

    // Zero denominator indicates parrallel lines
    if ( R_abs( denom ) < 0.0001f ) {
        return useOne ? pOne : pTwo;
    }

    mua = ( ( dpOneTwo * dpStaTwo ) - ( smTwo * dpStaOne ) ) / denom;
    mub = ( ( smOne * dpStaTwo ) - ( dpOneTwo * dpStaOne ) ) / denom;

    // If either of the edges has the nearest point out
    // of bounds, then the edges aren't crossed, we have
    // an edge-face contact. Our point is on the edge, which
    // we know from the useOne parameter.
    if ( mua > oneSize || mua < -oneSize || mub > twoSize || mub < -twoSize ) {
        return useOne ? pOne : pTwo;
    } else {

        cOne = Vec3Add( pOne, Vec3Scale( dOne, mua ) );
        cTwo = Vec3Add( pTwo, Vec3Scale( dTwo, mub ) );

        return Vec3Add( Vec3Scale( cOne, 0.5f ), Vec3Scale( cTwo, 0.5f ) );
    }
}

///>BoxAndBox
unsigned cyclone::BoxAndBox(const CollisionBox & one, const CollisionBox & two,
    CollisionData * data)
{
    //if ( !BoxAndBox( one, two ) ) return 0;

    // Find the vector between the two centres
    vec3_t toCentre = Vec3Subtract( Mat4AxisVector( two.primitive.transform, 3 ),
        Mat4AxisVector( one.primitive.transform, 3 ) );

    // We start assuming there is no contact
    real_t pen = REAL_MAX;
    unsigned best = 0xffffff;

    // Now we check each axes, returning if it gives us
    // a separating axis, and keeping track of the axis with
    // the smallest penetration otherwise.
    CHECK_OVERLAP( Mat4AxisVector( one.primitive.transform, 0 ), 0 );
    CHECK_OVERLAP( Mat4AxisVector( one.primitive.transform, 1 ), 1 );
    CHECK_OVERLAP( Mat4AxisVector( one.primitive.transform, 2 ), 2 );

    CHECK_OVERLAP( Mat4AxisVector( two.primitive.transform, 0 ), 3 );
    CHECK_OVERLAP( Mat4AxisVector( two.primitive.transform, 1 ), 4 );
    CHECK_OVERLAP( Mat4AxisVector( two.primitive.transform, 2 ), 5 );

    // Store the best axis-major, in case we run into almost
    // parallel edge collisions later
    unsigned bestSingleAxis = best;

    CHECK_OVERLAP( Vec3VectorProduct( Mat4AxisVector( one.primitive.transform, 0 ),
        Mat4AxisVector( two.primitive.transform, 0 ) ),  6 );
    CHECK_OVERLAP( Vec3VectorProduct( Mat4AxisVector( one.primitive.transform, 0 ),
        Mat4AxisVector( two.primitive.transform, 1 ) ),  7 );
    CHECK_OVERLAP( Vec3VectorProduct( Mat4AxisVector( one.primitive.transform, 0 ),
        Mat4AxisVector( two.primitive.transform, 2 ) ),  8 );
    CHECK_OVERLAP( Vec3VectorProduct( Mat4AxisVector( one.primitive.transform, 1 ),
        Mat4AxisVector( two.primitive.transform, 0 ) ),  9 );
    CHECK_OVERLAP( Vec3VectorProduct( Mat4AxisVector( one.primitive.transform, 1 ),
        Mat4AxisVector( two.primitive.transform, 1 ) ), 10 );
    CHECK_OVERLAP( Vec3VectorProduct( Mat4AxisVector( one.primitive.transform, 1 ),
        Mat4AxisVector( two.primitive.transform, 2 ) ), 11 );
    CHECK_OVERLAP( Vec3VectorProduct( Mat4AxisVector( one.primitive.transform, 2 ),
        Mat4AxisVector( two.primitive.transform, 0 ) ), 12 );
    CHECK_OVERLAP( Vec3VectorProduct( Mat4AxisVector( one.primitive.transform, 2 ),
        Mat4AxisVector( two.primitive.transform, 1 ) ), 13 );
    CHECK_OVERLAP( Vec3VectorProduct( Mat4AxisVector( one.primitive.transform, 2 ),
        Mat4AxisVector( two.primitive.transform, 2 ) ), 14 );
    
    // Make sure we've got a result.
    assert(best != 0xffffff);

    // We now know there's a collision, and we know which
    // of the axes gave the smallest penetration. We now
    // can deal with it in different ways depending on
    // the case.
    if ( best < 3 ) {

        // We've got a vertex of box two on a face of box one.
        FillPointFaceBoxBox( one, two, toCentre, data, best, pen );
        Add( *data, 1 );

        return 1;
    } else if ( best < 6 ) {

        // We've got a vertex of box one on a face of box two.
        // We use the same algorithm as above, but swap around
        // one and two (and therefore also the vector between their 
        // centres).
        FillPointFaceBoxBox( two, one, Vec3Scale( toCentre, -1.0f ), data,
            best - 3, pen );
        Add( *data, 1 );
        
        return 1;

    } else {
        
        // We've got an edge-edge contact. Find out which axes
        best -= 6;
        unsigned oneAxisIndex = best / 3;
        unsigned twoAxisIndex = best % 3;
        vec3_t oneAxis = Mat4AxisVector( one.primitive.transform, oneAxisIndex );
        vec3_t twoAxis = Mat4AxisVector( two.primitive.transform, twoAxisIndex );
        vec3_t axis = Vec3VectorProduct( oneAxis, twoAxis );
        axis = Vec3Normalise( axis );

        // The axis should point from box one to box two.
        if ( Vec3ScalarProduct( axis, toCentre ) > 0 ) {
            axis = Vec3Scale( axis, -1.0f );
        }

        // We have the axes, but not the edges: each axis has 4 edges parallel 
        // to it, we need to find which of the 4 for each object. We do 
        // that by finding the point in the centre of the edge. We know 
        // its component in the direction of the box's collision axis is zero 
        // (its a mid-point) and we determine which of the extremes in each 
        // of the other axes is closest.
        vec3_t ptOnOneEdge = one.halfSize;
        vec3_t ptOnTwoEdge = two.halfSize;
        
        for (unsigned i = 0; i < 3; i++)
        {
            // This doesn't work!
            if ( i == oneAxisIndex ) {
                ptOnOneEdge.n[i] = 0;
            } else if ( Vec3ScalarProduct( Mat4AxisVector( one.primitive.transform, i ),
                axis ) > 0 ) {
                ptOnOneEdge.n[i] = -ptOnOneEdge.n[i];
            }

            if ( i == twoAxisIndex ) {
                ptOnTwoEdge.n[i] = 0;
            } else if ( Vec3ScalarProduct( Mat4AxisVector( two.primitive.transform, i ),
                axis ) < 0 ) {
                ptOnTwoEdge.n[i] = -ptOnTwoEdge.n[i];
            }
        }

        // Move them into world coordinates (they are already oriented
        // correctly, since they have been derived from the axes).
        ptOnOneEdge = Mat4Transform( ptOnOneEdge,
            one.primitive.transform );
        ptOnTwoEdge = Mat4Transform( ptOnTwoEdge,
            two.primitive.transform );

#ifdef _DEBUG

        glPushMatrix();
        glColor3f(1,0,0);
        glTranslatef(ptOnOneEdge.x, ptOnOneEdge.y, ptOnOneEdge.z);
        glutSolidSphere(0.1, 10, 10);
        glPopMatrix();
        glBegin(GL_LINES);
        glVertex3f(ptOnOneEdge.x - oneAxis.x, 
            ptOnOneEdge.y - oneAxis.y,
            ptOnOneEdge.z - oneAxis.z);
        glVertex3f(ptOnOneEdge.x + oneAxis.x, 
            ptOnOneEdge.y + oneAxis.y,
            ptOnOneEdge.z + oneAxis.z);
        glEnd();

        glPushMatrix();
        glColor3f(0,0,1);
        glTranslatef(ptOnTwoEdge.x, ptOnTwoEdge.y, ptOnTwoEdge.z);
        glutSolidSphere(0.1, 10, 10);
        glPopMatrix();
        glBegin(GL_LINES);
        glVertex3f(ptOnTwoEdge.x - twoAxis.x, 
            ptOnTwoEdge.y - twoAxis.y,
            ptOnTwoEdge.z - twoAxis.z);
        glVertex3f(ptOnTwoEdge.x + twoAxis.x, 
            ptOnTwoEdge.y + twoAxis.y,
            ptOnTwoEdge.z + twoAxis.z);
        glEnd();

#endif /* #ifdef _DEBUG */

        // So we have a point and a direction for the colliding edges.
        // We need to find out point of closest approach of the two 
        // line-segments.
        vec3_t vertex = ContactPoint(
            ptOnOneEdge, oneAxis, one.halfSize.n[oneAxisIndex],
            ptOnTwoEdge, twoAxis, two.halfSize.n[twoAxisIndex],
            bestSingleAxis > 2
            );

        // We can fill the contact.
        Contact * contact = data->contacts;

        contact->penetration = pen;
        contact->contactNormal = axis;
        contact->contactPoint = vertex;
        SetBodyData( *contact,  one.primitive.body, two.primitive.body,
            data->friction, data->restitution );
        Add( *data, 1 );
        
        return 1;
    }

    /*return 0;*/
}
///<BoxAndBox

///>BoxAndPoint
unsigned cyclone::BoxAndPoint(const CollisionBox & box, const vec3_t & point,
    CollisionData * data)
{
    // Transform the point into box coordinates
    vec3_t relPt = Mat4TransformInverse( point, box.primitive.transform );

    vec3_t normal = Vec3Clear();

    // Check each axis, looking for the axis on which the
    // penetration is least deep.
    real_t min_depth = box.halfSize.x - R_abs( relPt.x );
    
    if ( min_depth < 0 ) {
        return 0;
    }

    normal = Vec3Scale( Mat4AxisVector( box.primitive.transform, 0 ),
        ( real_t )( ( relPt.x < 0 ) ? -1 : 1 ) );

    real_t depth = box.halfSize.y - R_abs( relPt.y );

    if ( depth < 0 ) {
        return 0;
    } else if ( depth < min_depth ) {

        min_depth = depth;
        normal = Vec3Scale( Mat4AxisVector( box.primitive.transform, 1 ),
            ( real_t )( ( relPt.y < 0 ) ? -1 : 1 ) );
    }

    depth = box.halfSize.z - R_abs( relPt.z );

    if ( depth < 0 ) {
        return 0;
    } else if ( depth < min_depth ) {

        min_depth = depth;
        normal = Vec3Scale( Mat4AxisVector( box.primitive.transform, 2 ),
            ( real_t )( ( relPt.z < 0 ) ? -1 : 1 ) );
    }

    // Compile the contact
    Contact * contact = data->contacts;
    contact->contactNormal = normal;
    contact->contactPoint = point;
    contact->penetration = min_depth;

    // Note that we don't know what rigid body the point
    // belongs to, so we just use NULL. Where this is called
    // this value can be left, or filled in.
    SetBodyData( *contact,  box.primitive.body, NULL, data->friction,
        data->restitution );
    Add( *data, 1 );

    return 1;
}
///<BoxAndPoint


///>BoxAndSphere
unsigned cyclone::BoxAndSphere(const CollisionBox & box,
    const CollisionSphere & sphere, CollisionData * data)
{
///>SphereBoxCentreToBox
    // Transform the centre of the sphere into box coordinates
    vec3_t centre = Mat4AxisVector( sphere.primitive.transform, 3 );
    vec3_t relCentre = Mat4TransformInverse( Mat4AxisVector( sphere.primitive.transform, 3 ),
        box.primitive.transform );
///<SphereBoxCentreToBox

///>SphereBoxEarlyOut
    // Early out check to see if we can exclude the contact
    if ( R_abs( relCentre.x ) - sphere.radius > box.halfSize.x ||
         R_abs( relCentre.y ) - sphere.radius > box.halfSize.y ||
         R_abs( relCentre.z ) - sphere.radius > box.halfSize.z ) {
        return 0;
    }
///<SphereBoxEarlyOut
    
///>SphereBoxClampAndTest
    vec3_t closestPt = { 0, 0, 0 };
    real_t dist;

    // Clamp each coordinate to the box.
    dist = relCentre.x;
    
    if ( dist > box.halfSize.x ) {
        dist = box.halfSize.x;
    }

    if ( dist < -box.halfSize.x ) {
        dist = -box.halfSize.x;
    }

    closestPt.x = dist;

    dist = relCentre.y;
    
    if ( dist > box.halfSize.y ) {
        dist = box.halfSize.y;
    }

    if ( dist < -box.halfSize.y ) {
        dist = -box.halfSize.y;
    }

    closestPt.y = dist;

    dist = relCentre.z;
    
    if ( dist > box.halfSize.z ) {
        dist = box.halfSize.z;
    }

    if ( dist < -box.halfSize.z ) {
        dist = -box.halfSize.z;
    }

    closestPt.z = dist;

    // Check we're in contact
    dist = Vec3MagnitudeSqr( Vec3Subtract( closestPt, relCentre ) );
    
    if ( dist > sphere.radius * sphere.radius ) {
        return 0;
    }

///<SphereBoxClampAndTest

///>SphereBoxUntransformClosestPoint
    // Compile the contact
    vec3_t closestPtWorld = Mat4Transform( closestPt,
        box.primitive.transform );
///<SphereBoxUntransformClosestPoint

    Contact * contact = data->contacts;
    contact->contactNormal = Vec3Normalise(
        Vec3Subtract( closestPtWorld, centre ) );
    contact->contactPoint = closestPtWorld;
    contact->penetration = sphere.radius - R_sqrt(dist);
    SetBodyData( *contact,  box.primitive.body, sphere.primitive.body,
        data->friction, data->restitution );

    Add( *data, 1 );

    return 1;
}
///<BoxAndSphere
