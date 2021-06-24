/*
 * The fracture demo.
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
#include <gl/glut.h>
#include "../app.h"
#include "../timing.h"

#include <stdio.h>

#define MAX_BLOCKS 9

cyclone::Random global_random;

class Block : public cyclone::CollisionBox
{
public:
    bool exists;

    Block()
    :
    exists(false)
    {
        primitive = cyclone::COLLISION_PRIMITIVE;
        primitive.body = new cyclone::RigidBody();
    }

    ~Block()
    {
        delete primitive.body;
    }

    /** Draws the block. */
    void render()
    {
        // Get the OpenGL transformation
        GLfloat mat[16];
        primitive.body->getGLTransform(mat);

        if (primitive.body->getAwake()) glColor3f(1.0,0.7,0.7);
        else glColor3f(0.7,0.7,1.0);

        glPushMatrix();
        glMultMatrixf(mat);
        glScalef(halfSize.x*2, halfSize.y*2, halfSize.z*2);
        glutSolidCube(1.0f);
        glPopMatrix();
    }

    /** Sets the block to a specific location. */
    void setState(const cyclone::vec3_t & position, 
                  const cyclone::quat_t & orientation, 
                  const cyclone::vec3_t & extents,
                  const cyclone::vec3_t & velocity)
    {
        primitive.body->setPosition(position);
        primitive.body->setOrientation(orientation);
        primitive.body->setVelocity(velocity);
        primitive.body->setRotation( cyclone::Vec3Clear() );
        halfSize = extents;

        cyclone::real_t mass = halfSize.x * halfSize.y * halfSize.z * 8.0f;
        primitive.body->setMass(mass);

        cyclone::mat3_t tensor = Mat3SetBlockInertiaTensor( halfSize, mass );
        primitive.body->setInertiaTensor(tensor);

        primitive.body->setLinearDamping(0.95f);
        primitive.body->setAngularDamping(0.8f);
        primitive.body->clearAccumulators();
        primitive.body->setAcceleration(0,-10.0f,0);

        //primitive.body->setCanSleep(false);
        primitive.body->setAwake();

        primitive.body->calculateDerivedData();
    }

///>CalculateBlockTensor
    /**
     * Calculates and sets the mass and inertia tensor of this block,
     * assuming it has the given constant density.
     */
    void calculateMassProperties(cyclone::real_t invDensity)
    {
        // Check for infinite mass
        if (invDensity <= 0)
        {
            // Just set zeros for both mass and inertia tensor
            primitive.body->setInverseMass(0);
            primitive.body->setInverseInertiaTensor( cyclone::Mat3Identity() );
        }
        else
        {
            // Otherwise we need to calculate the mass
            cyclone::real_t volume = cyclone::Vec3Magnitude( halfSize ) * 2.0;
            cyclone::real_t mass = volume / invDensity;

            primitive.body->setMass(mass);

            // And calculate the inertia tensor from the mass and size
            mass *= 0.333f;
            cyclone::mat3_t tensor = cyclone::Mat3SetInertiaTensorCoeffs(
                mass * halfSize.y*halfSize.y + halfSize.z*halfSize.z,
                mass * halfSize.y*halfSize.x + halfSize.z*halfSize.z,
                mass * halfSize.y*halfSize.x + halfSize.z*halfSize.y
                );
            primitive.body->setInertiaTensor(tensor);
        }

    }
///<CalculateBlockTensor

///>DivideBlock
    /** 
     * Performs the division of the given block into four, writing the 
     * eight new blocks into the given blocks array. The blocks array can be
     * a pointer to the same location as the target pointer: since the 
     * original block is always deleted, this effectively reuses its storage.
     * The algorithm is structured to allow this reuse.
     */
    void divideBlock(const cyclone::Contact& contact, 
        Block* target, Block* blocks)
    {
        // Find out if we're block one or two in the contact structure, and
        // therefore what the contact normal is.
        cyclone::vec3_t normal = contact.contactNormal;
        cyclone::RigidBody *body = contact.body[0];

        if (body != target->primitive.body) 
        {
            normal = Vec3Invert( normal );
            body = contact.body[1];
        }

        // Work out where on the body (in body coordinates) the contact is
        // and its direction.
        cyclone::vec3_t point = primitive.body->getPointInLocalSpace(contact.contactPoint);
        normal = primitive.body->getDirectionInLocalSpace(normal);

        // Work out the centre of the split: this is the point coordinates
        // for each of the axes perpendicular to the normal, and 0 for the
        // axis along the normal.
        point = cyclone::Vec3Subtract( point, cyclone::Vec3Scale( normal,
            cyclone::Vec3ScalarProduct( point, normal ) ) );

        // Take a copy of the half size, so we can create the new blocks.
        cyclone::vec3_t size = target->halfSize;

        // Take a copy also of the body's other data.
        cyclone::RigidBody tempBody;
        tempBody.setPosition(primitive.body->getPosition());
        tempBody.setOrientation(primitive.body->getOrientation());
        tempBody.setVelocity(primitive.body->getVelocity());
        tempBody.setRotation(primitive.body->getRotation());
        tempBody.setLinearDamping(primitive.body->getLinearDamping());
        tempBody.setAngularDamping(primitive.body->getAngularDamping());
        tempBody.setInverseInertiaTensor(primitive.body->getInverseInertiaTensor());
        tempBody.calculateDerivedData();

        // Remove the old block
        target->exists = false;

        // Work out the inverse density of the old block
        cyclone::real_t invDensity = cyclone::Vec3Magnitude( halfSize ) * 8 *
            primitive.body->getInverseMass();

        // Now split the block into eight.
        for (unsigned i = 0; i < 8; i++)
        {
            // Find the minimum and maximum extents of the new block
            // in old-block coordinates
            cyclone::vec3_t min, max;
            if ((i & 1) == 0) {
                min.x = -size.x;
                max.x = point.x;
            } else {
                min.x = point.x;
                max.x = size.x;
            }
            if ((i & 2) == 0) {
                min.y = -size.y;
                max.y = point.y;
            } else {
                min.y = point.y;
                max.y = size.y;
            }
            if ((i & 4) == 0) {
                min.z = -size.z;
                max.z = point.z;
            } else {
                min.z = point.z;
                max.z = size.z;
            }

            // Get the origin and half size of the block, in old-body 
            // local coordinates.
            cyclone::vec3_t halfSize = cyclone::Vec3Scale(
                cyclone::Vec3Subtract( max, min ), 0.5f );
            cyclone::vec3_t newPos = cyclone::Vec3Add( halfSize, min );

            // Convert the origin to world coordinates.
            newPos = tempBody.getPointInWorldSpace(newPos);

			// Work out the direction to the impact.
			cyclone::vec3_t direction = cyclone::Vec3Subtract( newPos,
                contact.contactPoint );
			direction = cyclone::Vec3Normalise( direction );

            // Set the body's properties (we assume the block has a body
            // already that we're going to overwrite).
            blocks[i].primitive.body->setPosition(newPos);
            blocks[i].primitive.body->setVelocity( cyclone::Vec3Add(
                tempBody.getVelocity(),
                cyclone::Vec3Scale( direction, 10.0f ) ) );
            blocks[i].primitive.body->setOrientation(tempBody.getOrientation());
            blocks[i].primitive.body->setRotation(tempBody.getRotation());
            blocks[i].primitive.body->setLinearDamping(tempBody.getLinearDamping());
            blocks[i].primitive.body->setAngularDamping(tempBody.getAngularDamping());
			blocks[i].primitive.body->setAwake(true);
			blocks[i].primitive.body->setAcceleration(cyclone::GRAVITY);
			blocks[i].primitive.body->clearAccumulators();
			blocks[i].primitive.body->calculateDerivedData();
            blocks[i].primitive.offset = cyclone::Mat4Identity();
            blocks[i].exists = true;
			blocks[i].halfSize = halfSize;

            // Finally calculate the mass and inertia tensor of the new block
            blocks[i].calculateMassProperties(invDensity);
        }
    }
///<DivideBlock
};

/**
 * The main demo class definition.
 */
class FractureDemo : public RigidBodyApplication
{
	/** Tracks if a block has been hit. */
	bool hit;
	bool ball_active;
	unsigned fracture_contact;

	/** Handle random numbers. */
	cyclone::Random random;

    /** Holds the bodies. */
    Block blocks[MAX_BLOCKS];

    /** Holds the projectile. */
    cyclone::CollisionSphere ball;

    /** Processes the contact generation code. */
	virtual void generateContacts();

    /** Processes the objects in the simulation forward in time. */
    virtual void updateObjects(cyclone::real_t duration);

    /** Resets the position of all the blocks. */
    virtual void reset();

	/** Processes the physics. */
	virtual void update();

public:
    /** Creates a new demo object. */
    FractureDemo();

    /** Returns the window title for the demo. */
    virtual const char* getTitle();

    /** Display the particle positions. */
    virtual void display();
};

// Method definitions
FractureDemo::FractureDemo()
    :
    RigidBodyApplication()
{
	// Create the ball.
    ball = cyclone::COLLISION_SPHERE;
	ball.primitive.body = new cyclone::RigidBody();
	ball.radius = 0.25f;
	ball.primitive.body->setMass(5.0f);
	ball.primitive.body->setDamping(0.9f, 0.9f);
	cyclone::mat3_t it = cyclone::Mat3SetDiagonal( 5.0f, 5.0f, 5.0f );
	ball.primitive.body->setInertiaTensor(it);
	ball.primitive.body->setAcceleration(cyclone::GRAVITY);

	ball.primitive.body->setCanSleep(false);
	ball.primitive.body->setAwake(true);

    // Set up the initial block
    reset();
}

const char* FractureDemo::getTitle()
{
    return "Cyclone > Fracture Demo";
}

void FractureDemo::generateContacts()
{
	hit = false;

	// Create the ground plane data
	cyclone::CollisionPlane plane = cyclone::COLLISION_PLANE;
    cyclone::vec3_t direction = { 0, 1, 0 };
	plane.direction = direction;
	plane.offset = 0;

	// Set up the collision data structure
	cData.reset(maxContacts);
	cData.friction = (cyclone::real_t)0.9;
	cData.restitution = (cyclone::real_t)0.2;
	cData.tolerance = (cyclone::real_t)0.1;

	// Perform collision detection
	for (Block *block = blocks; block < blocks+MAX_BLOCKS; block++)
	{
		if (!block->exists) continue;

		// Check for collisions with the ground plane
		if (!cData.hasMoreContacts()) return;

        /* collision detection */
		cyclone::BoxAndHalfSpace( *block, plane, &cData );

		if (ball_active)
		{
			// And with the sphere
			if (!cData.hasMoreContacts()) return;

            /* collision detection */
			if ( cyclone::BoxAndSphere( *block, ball, &cData ) ) {

				hit = true;
				fracture_contact = cData.contactCount-1;
			}
		}

		// Check for collisions with each other box
		for (Block *other = block+1; other < blocks+MAX_BLOCKS; other++)
		{
			if (!other->exists) continue;

			if (!cData.hasMoreContacts()) return;

            /* collision detection */
			cyclone::BoxAndBox( *block, *other, &cData );
		}
	}

	// Check for sphere ground collisions
	if (ball_active)
	{
		if (!cData.hasMoreContacts()) return;

        /* collision detection */
		cyclone::SphereAndHalfSpace( ball, plane, &cData );
	}
}

void FractureDemo::reset()
{
	// Only the first block exists
	blocks[0].exists = true;
	for (Block *block = blocks+1; block < blocks+MAX_BLOCKS; block++)
	{
		block->exists = false;
	}

	// Set the first block
    cyclone::vec3_t halfSize = { 4, 4, 4 };
	blocks[0].halfSize = halfSize;
	blocks[0].primitive.body->setPosition(0, 7, 0);
	blocks[0].primitive.body->setOrientation(1,0,0,0);
	blocks[0].primitive.body->setVelocity(0,0,0);
	blocks[0].primitive.body->setRotation(0,0,0);
	blocks[0].primitive.body->setMass(100.0f);
	cyclone::mat3_t it = Mat3SetBlockInertiaTensor( blocks[0].halfSize,
        100.0f );
	blocks[0].primitive.body->setInertiaTensor(it);
	blocks[0].primitive.body->setDamping(0.9f, 0.9f);
	blocks[0].primitive.body->calculateDerivedData();
	CalculateInternals( blocks[0].primitive );

	blocks[0].primitive.body->setAcceleration(cyclone::GRAVITY);
	blocks[0].primitive.body->setAwake(true);
	blocks[0].primitive.body->setCanSleep(true);


	ball_active = true;

	// Set up the ball
	ball.primitive.body->setPosition(0,5.0f,20.0f);
	ball.primitive.body->setOrientation(1,0,0,0);
	ball.primitive.body->setVelocity(
		random.randomBinomial(4.0f),
		random.randomReal(1.0f, 6.0f),
		-20.0f
		);
	ball.primitive.body->setRotation(0,0,0);
	ball.primitive.body->calculateDerivedData();
	ball.primitive.body->setAwake(true);
	CalculateInternals( ball.primitive );

	hit = false;

	// Reset the contacts
	cData.contactCount = 0;
}

void FractureDemo::update()
{
	RigidBodyApplication::update();

	// Handle fractures.
	if (hit)
	{
		blocks[0].divideBlock(
			cData.contactArray[fracture_contact], 
			blocks, 
			blocks+1
			);
		ball_active = false;
	}
}

void FractureDemo::updateObjects(cyclone::real_t duration)
{
    for (Block *block = blocks; block < blocks+MAX_BLOCKS; block++)
    {
        if (block->exists)
        {
            block->primitive.body->integrate(duration);
            CalculateInternals( block->primitive );
        }
    }

	if (ball_active)
	{
		ball.primitive.body->integrate(duration);
		CalculateInternals( ball.primitive );
	}
}

void FractureDemo::display()
{
	const static GLfloat lightPosition[] = {0.7f,1,0.4f,0};

	RigidBodyApplication::display();

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
	glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	glEnable(GL_NORMALIZE);
	for (Block *block = blocks; block < blocks+MAX_BLOCKS; block++)
	{
		if (block->exists) block->render();
	}
	glDisable(GL_NORMALIZE);

	if (ball_active)
	{
		glColor3f(0.4f, 0.7f, 0.4f);
		glPushMatrix();
		cyclone::vec3_t pos = ball.primitive.body->getPosition();
		glTranslatef(pos.x, pos.y, pos.z);
		glutSolidSphere(0.25f, 16, 8);
		glPopMatrix();
	}

	glDisable(GL_LIGHTING);
	glDisable(GL_COLOR_MATERIAL);

	// Draw some scale circles
	glColor3f(0.75, 0.75, 0.75);
	for (unsigned i = 1; i < 20; i++)
	{
		glBegin(GL_LINE_LOOP);
		for (unsigned j = 0; j < 32; j++)
		{
			float theta = 3.1415926f * j / 16.0f;
			glVertex3f(i*cosf(theta),0.0f,i*sinf(theta));
		}
		glEnd();
	}
	glBegin(GL_LINES);
	glVertex3f(-20,0,0);
	glVertex3f(20,0,0);
	glVertex3f(0,0,-20);
	glVertex3f(0,0,20);
	glEnd();

	RigidBodyApplication::drawDebug();
}

/**
 * Called by the common demo framework to create an application
 * object (with new) and return a pointer.
 */
Application* getApplication()
{
    return new FractureDemo();
}

