/*
 * The explosion demo.
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

#define OBJECTS 5

// Holds a transform matrix for rendering objects
// reflected in the floor.
GLfloat floorMirror[16] = 
{
    1, 0, 0, 0,
    0, -1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
};

class Ball : public cyclone::CollisionSphere
{
public:
    Ball()
    {
        primitive = cyclone::COLLISION_PRIMITIVE;
        primitive.body = new cyclone::RigidBody;
    }

    ~Ball()
    {
        delete primitive.body;
    }

    /** Draws the box, excluding its shadow. */
    void render()
    {
        // Get the OpenGL transformation
        GLfloat mat[16];
        primitive.body->getGLTransform(mat);

        if (primitive.body->getAwake()) glColor3f(1.0,0.7,0.7);
        else glColor3f(0.7,0.7,1.0);

        glPushMatrix();
        glMultMatrixf(mat);
        glutSolidSphere(radius, 20, 20);
        glPopMatrix();
    }

    /** Draws the ground plane shadow for the box. */
    void renderShadow()
    {
        // Get the OpenGL transformation
        GLfloat mat[16];
        primitive.body->getGLTransform(mat);

        glPushMatrix();
        glScalef(1.0f, 0, 1.0f);
        glMultMatrixf(mat);
        glutSolidSphere(radius, 20, 20);
        glPopMatrix();
    }

    /** Sets the box to a specific location. */
    void setState(cyclone::vec3_t position, 
                  cyclone::quat_t orientation, 
                  cyclone::real_t radius,
                  cyclone::vec3_t velocity)
    {
        cyclone::vec3_t rotation = { 0, 0, 0 };

        primitive.body->setPosition(position);
        primitive.body->setOrientation(orientation);
        primitive.body->setVelocity(velocity);
        primitive.body->setRotation(rotation);
        Ball::radius = radius;

        cyclone::real_t mass = 4.0f*0.3333f*3.1415f * radius*radius*radius;
        primitive.body->setMass(mass);

        cyclone::mat3_t tensor = cyclone::Mat3Identity();
        cyclone::real_t coeff = 0.4f*mass*radius*radius;
        tensor = cyclone::Mat3SetInertiaTensorCoeffs( coeff, coeff, coeff );
        primitive.body->setInertiaTensor(tensor);

        primitive.body->setLinearDamping(0.95f);
        primitive.body->setAngularDamping(0.8f);
        primitive.body->clearAccumulators();
        primitive.body->setAcceleration(0,-10.0f,0);

        //primitive.body->setCanSleep(false);
        primitive.body->setAwake();

        primitive.body->calculateDerivedData();
    }

    /** Positions the box at a random location. */
    void random(cyclone::Random *random)
    {
        const static cyclone::vec3_t minPos = {-5,  5,-5 };
        const static cyclone::vec3_t maxPos = { 5, 10, 5 };

        setState(
            random->randomVector(minPos, maxPos),
            random->randomQuaternion(),
            random->randomReal(0.5f, 1.5f),
            cyclone::Vec3Clear()
        );
    }
};

class Box : public cyclone::CollisionBox
{
public:
    bool isOverlapping;

    Box()
    {
        primitive = cyclone::COLLISION_PRIMITIVE;
        primitive.body = new cyclone::RigidBody;
    }

    ~Box()
    {
        delete primitive.body;
    }

    /** Draws the box, excluding its shadow. */
    void render()
    {
        // Get the OpenGL transformation
        GLfloat mat[16];
        primitive.body->getGLTransform(mat);

        if (isOverlapping) glColor3f(0.7,1.0,0.7);
        else if (primitive.body->getAwake()) glColor3f(1.0,0.7,0.7);
        else glColor3f(0.7,0.7,1.0);

        glPushMatrix();
        glMultMatrixf(mat);
        glScalef(halfSize.x*2, halfSize.y*2, halfSize.z*2);
        glutSolidCube(1.0f);
        glPopMatrix();
    }

    /** Draws the ground plane shadow for the box. */
    void renderShadow()
    {
        // Get the OpenGL transformation
        GLfloat mat[16];
        primitive.body->getGLTransform(mat);

        glPushMatrix();
        glScalef(1.0f, 0, 1.0f);
        glMultMatrixf(mat);
        glScalef(halfSize.x*2, halfSize.y*2, halfSize.z*2);
        glutSolidCube(1.0f);
        glPopMatrix();
    }

    /** Sets the box to a specific location. */
    void setState(const cyclone::vec3_t &position, 
                  const cyclone::quat_t &orientation, 
                  const cyclone::vec3_t &extents,
                  const cyclone::vec3_t &velocity)
    {
        cyclone::vec3_t rotation = { 0, 0, 0 };

        primitive.body->setPosition(position);
        primitive.body->setOrientation(orientation);
        primitive.body->setVelocity(velocity);
        primitive.body->setRotation(rotation);
        halfSize = extents;

        cyclone::real_t mass = halfSize.x * halfSize.y * halfSize.z * 8.0f;
        primitive.body->setMass(mass);

        cyclone::vec3_t hs = { halfSize.x, halfSize.y, halfSize.z };
        cyclone::mat3_t tensor = cyclone::Mat3Identity();
        tensor = cyclone::Mat3SetBlockInertiaTensor(hs, mass);
        primitive.body->setInertiaTensor(tensor);

        primitive.body->setLinearDamping(0.95f);
        primitive.body->setAngularDamping(0.8f);
        primitive.body->clearAccumulators();
        primitive.body->setAcceleration(0,-10.0f,0);

        //primitive.body->setCanSleep(false);
        primitive.body->setAwake();

        primitive.body->calculateDerivedData();
    }

    /** Positions the box at a random location. */
    void random(cyclone::Random *random)
    {
        const static cyclone::vec3_t minPos = {-5, 5, -5 };
        const static cyclone::vec3_t maxPos = { 5, 10, 5 };
        const static cyclone::vec3_t minSize = { 0.5f, 0.5f, 0.5f };
        const static cyclone::vec3_t maxSize = { 4.5f, 1.5f, 1.5f };

        setState(
            random->randomVector(minPos, maxPos),
            random->randomQuaternion(),
            random->randomVector(minSize, maxSize),
            cyclone::Vec3Clear()
        );
    }
};

/**
 * The main demo class definition.
 */
class ExplosionDemo : public RigidBodyApplication
{
    bool editMode, upMode;

    /** 
     * Holds the number of boxes in the simulation.
     */
    const static unsigned boxes = OBJECTS;

    /** Holds the box data. */
    Box boxData[boxes];

    /** 
     * Holds the number of balls in the simulation.
     */
    const static unsigned balls = 1;//OBJECTS;

    /** Holds the ball data. */
    Ball ballData[balls];


    /** Detonates the explosion. */ 
    void fire();

    /** Resets the position of all the boxes and primes the explosion. */
    virtual void reset();

	/** Processes the contact generation code. */
	virtual void generateContacts();

    /** Processes the objects in the simulation forward in time. */
    virtual void updateObjects(cyclone::real_t duration);

public:
    /** Creates a new demo object. */
    ExplosionDemo();

    /** Sets up the rendering. */
    virtual void initGraphics();

    /** Returns the window title for the demo. */
    virtual const char* getTitle();

    /** Display the particle positions. */
    virtual void display();

    /** Handles a key press. */
    virtual void key(unsigned char key);

    /** Handle a mouse drag */
    virtual void mouseDrag(int x, int y);
};

// Method definitions
ExplosionDemo::ExplosionDemo() 
    : 
    RigidBodyApplication(),
    editMode(false),
    upMode(false)
{   
    // Reset the position of the boxes
    reset();
}

const char* ExplosionDemo::getTitle()
{ 
    return "Cyclone > Explosion Demo"; 
}

void ExplosionDemo::fire()
{
    cyclone::vec3_t position = ballData[0].primitive.body->getPosition();
    position = Vec3Normalise( position );
    ballData[0].primitive.body->addForce( Vec3Scale( position, -1000.0f ) );
}

void ExplosionDemo::reset()
{
    Box *box = boxData;

    cyclone::vec3_t pos0 = { 0, 3, 0 };
    cyclone::quat_t ori0 = { 1.0f, 0.0f, 0.0f, 0.0f };
    cyclone::vec3_t ext0 = { 4, 1, 1 };
    cyclone::vec3_t vel0 = { 0, 1, 0 };

    box++->setState( pos0, ori0, ext0, vel0 );

    cyclone::vec3_t pos1 = { 0, 4.75, 2 };
    cyclone::quat_t ori1 = { 1.0f, 0.1f, 0.05f, 0.01f };
    cyclone::vec3_t ext1 = { 1, 1, 4 };
    cyclone::vec3_t vel1 = { 0, 0, 0 };

    box++->setState( pos1, ori1, ext1, vel1 );

    // Create the random objects
    cyclone::Random random;
    for (; box < boxData+boxes; box++) 
    {
        box->random(&random);
    }

    for (Ball *ball = ballData; ball < ballData+balls; ball++) 
    {
        ball->random(&random);
    }

	// Reset the contacts
	cData.contactCount = 0;
}

void ExplosionDemo::generateContacts()
{
	// Note that this method makes a lot of use of early returns to avoid 
	// processing lots of potential contacts that it hasn't got room to 
	// store.

    // Create the ground plane data
    cyclone::CollisionPlane plane = cyclone::COLLISION_PLANE;
    plane.direction.x = 0;
    plane.direction.y = 1;
    plane.direction.z = 0;
    plane.offset = 0;

    // Set up the collision data structure
    cData.reset(maxContacts);
    cData.friction = (cyclone::real_t)0.9;
    cData.restitution = (cyclone::real_t)0.6;
    cData.tolerance = (cyclone::real_t)0.1;

    // Perform exhaustive collision detection
    for (Box *box = boxData; box < boxData+boxes; box++)
    {
        // Check for collisions with the ground plane
		if (!cData.hasMoreContacts()) return;

        /* collision detection */
        cyclone::BoxAndHalfSpace( *box, plane, &cData );

        // Check for collisions with each other box
        for (Box *other = box+1; other < boxData+boxes; other++)
        {
    		if (!cData.hasMoreContacts()) return;

            /* collision detection */
			cyclone::BoxAndBox( *box, *other, &cData );

            /* intersection test */
            if ( cyclone::BoxAndBox( *box, *other ) ) {
                box->isOverlapping = other->isOverlapping = true;
            }
		}

		// Check for collisions with each ball
        for (Ball *other = ballData; other < ballData+balls; other++)
        { 
    		if (!cData.hasMoreContacts()) return;

            /* collision detection */
            cyclone::BoxAndSphere( *box, *other, &cData );
		}
    }

    for (Ball *ball = ballData; ball < ballData+balls; ball++)
    {
        // Check for collisions with the ground plane
        if (!cData.hasMoreContacts()) return;

        /* collision detection */
        cyclone::SphereAndHalfSpace( *ball, plane, &cData );

        for (Ball *other = ball+1; other < ballData+balls; other++)
        {
            // Check for collisions with the ground plane
    		if (!cData.hasMoreContacts()) return;

            cyclone::SphereAndSphere( *ball, *other, &cData );
        }
    }
}

void ExplosionDemo::updateObjects(cyclone::real_t duration)
{
///>ExplosionUpdate
    // Update the physics of each box in turn
    for (Box *box = boxData; box < boxData+boxes; box++)
    {
        // Run the physics
        box->primitive.body->integrate(duration);
        CalculateInternals( box->primitive );
        box->isOverlapping = false;
    }

    // Update the physics of each ball in turn
    for (Ball *ball = ballData; ball < ballData+balls; ball++)
    {
        // Run the physics
        ball->primitive.body->integrate(duration);
        CalculateInternals( ball->primitive );
    }
///<ExplosionUpdate
}


void ExplosionDemo::initGraphics()
{
    GLfloat lightAmbient[] = {0.8,0.8,0.8,1.0};
    GLfloat lightDiffuse[] = {0.9,0.95,1.0,1.0};

    glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);

    glEnable(GL_LIGHT0);

    Application::initGraphics();
}

void ExplosionDemo::display()
{
	const static GLfloat lightPosition[] = {1,-1,0,0};
	const static GLfloat lightPositionMirror[] = {1,1,0,0};

    // Update the physics of each box in turn
    for (Box *box = boxData; box < boxData+boxes; box++)
    {
        // Run the physics
        CalculateInternals( box->primitive );
        box->isOverlapping = false;
    }

    // Update the physics of each ball in turn
    for (Ball *ball = ballData; ball < ballData+balls; ball++)
    {
        // Run the physics
        CalculateInternals( ball->primitive );
    }

    // Clear the viewport and set the camera direction
	RigidBodyApplication::display();

    // Render each element in turn as a shadow
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
	glPushMatrix();
    glMultMatrixf(floorMirror);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glColor3f(1,0,0);
    for (Box *box = boxData; box < boxData+boxes; box++)
    {
        box->render();
    }
    for (Ball *ball = ballData; ball < ballData+balls; ball++)
    {
        ball->render();
    }
    glPopMatrix();
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

    // Render each shadow in turn
	glEnable(GL_BLEND);
    glColor4f(0,0,0,0.1f);
    glDisable(GL_DEPTH_TEST);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
    for (Box *box = boxData; box < boxData+boxes; box++)
    {
        box->renderShadow();
    }
    for (Ball *ball = ballData; ball < ballData+balls; ball++)
    {
        ball->renderShadow();
    }
	glDisable(GL_BLEND);

    // Render the boxes themselves
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glLightfv(GL_LIGHT0, GL_POSITION, lightPositionMirror);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    for (Box *box = boxData; box < boxData+boxes; box++)
    {
        box->render();
    }
    for (Ball *ball = ballData; ball < ballData+balls; ball++)
    {
        ball->render();
    }
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    // Finish the frame, rendering any additional information
	drawDebug();
}

void ExplosionDemo::mouseDrag(int x, int y)
{
    if (editMode)
    {
        cyclone::vec3_t mPos = {
            (x-last_x) * 0.125f,
            0,
            (y-last_y) * 0.125f
        };

        boxData[0].primitive.body->setPosition( Vec3Add( boxData[0].primitive.body->getPosition(),
            mPos ) );
        boxData[0].primitive.body->calculateDerivedData();
    }
    else if (upMode)
    {
        cyclone::vec3_t mPos = {
            0,
            (y-last_y) * 0.125f,
            0
        };

        boxData[0].primitive.body->setPosition( Vec3Add( boxData[0].primitive.body->getPosition(),
            mPos ) );
        boxData[0].primitive.body->calculateDerivedData();
    }
    else
    {
        RigidBodyApplication::mouseDrag(x, y);
    }

    // Remember the position
    last_x = x;
    last_y = y;
}


void ExplosionDemo::key(unsigned char key)
{
    switch(key)
    {
    case 'e': case 'E':
        editMode = !editMode;
        upMode = false;
        return;

    case 't': case 'T':
        upMode = !upMode;
        editMode = false;
        return;

    case 'w': case 'W':
        for (Box *box = boxData; box < boxData+boxes; box++) 
            box->primitive.body->setAwake();
        for (Ball *ball = ballData; ball < ballData+balls; ball++)
            ball->primitive.body->setAwake();
        return;
    }

    RigidBodyApplication::key(key);
}

/**
 * Called by the common demo framework to create an application
 * object (with new) and return a pointer.
 */
Application* getApplication()
{
    return new ExplosionDemo();
} 