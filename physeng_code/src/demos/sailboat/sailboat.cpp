/*
 * The sailboat demo.
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
#include <cassert>


/**
 * The main demo class definition.
 */
class SailboatDemo : public Application
{
	cyclone::Buoyancy buoyancy;

	cyclone::Aero sail;
	cyclone::RigidBody sailboat;
	cyclone::ForceRegistry registry;

	cyclone::Random r;
	cyclone::vec3_t windspeed;

	float sail_control;

public:
    /** Creates a new demo object. */
    SailboatDemo();
    virtual ~SailboatDemo();

    /** Returns the window title for the demo. */
    virtual const char* getTitle();

	/** Display the particles. */
	virtual void display();

	/** Update the particle positions. */
	virtual void update();
	 
	/** Handle a key press. */
    virtual void key(unsigned char key);
};

// Method definitions
SailboatDemo::SailboatDemo()
:
Application()
{

	sail = cyclone::Aero( cyclone::Mat3Construct(0,0,0, 0,0,0, 0,0,-1.0f),
		cyclone::Vec3Construct(2.0f, 0, 0), &windspeed );

	buoyancy = cyclone::Buoyancy( cyclone::Vec3Construct(0.0f, 0.5f, 0.0f),
		1.0f, 3.0f, 1.6f);
	sail_control = 0;
	windspeed = cyclone::Vec3Construct(0,0,0);

	// Set up the boat's rigid body.
	sailboat.setPosition(0, 1.6f, 0);
	sailboat.setOrientation(1,0,0,0);

	sailboat.setVelocity(0,0,0);
	sailboat.setRotation(0,0,0);

	sailboat.setMass(200.0f);
	cyclone::mat3_t it = cyclone::Mat3SetBlockInertiaTensor(
		cyclone::Vec3Construct(2,1,1), 100.0f);
	sailboat.setInertiaTensor(it);

	sailboat.setDamping(0.8f, 0.8f);

	sailboat.setAcceleration(cyclone::GRAVITY);
	sailboat.calculateDerivedData();

	sailboat.setAwake();
	sailboat.setCanSleep(false);

	registry.add(&sailboat, &sail);
	registry.add(&sailboat, &buoyancy);
}

SailboatDemo::~SailboatDemo()
{
}

static void drawBoat()
{
	// Left Hull
	glPushMatrix();
	glTranslatef(0, 0, -1.0f);
	glScalef(2.0f, 0.4f, 0.4f);
	glutSolidCube(1.0f);
	glPopMatrix();

	// Right Hull
	glPushMatrix();
	glTranslatef(0, 0, 1.0f);
	glScalef(2.0f, 0.4f, 0.4f);
	glutSolidCube(1.0f);
	glPopMatrix();

	// Deck
	glPushMatrix();
	glTranslatef(0, 0.3f, 0);
	glScalef(1.0f, 0.1f, 2.0f);
	glutSolidCube(1.0f);
	glPopMatrix();

	// Mast
	glPushMatrix();
	glTranslatef(0, 1.8f, 0);
	glScalef(0.1f, 3.0f, 0.1f);
	glutSolidCube(1.0f);
	glPopMatrix();

}

void SailboatDemo::display()
{
	// Clear the view port and set the camera direction
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	cyclone::vec3_t pos = sailboat.getPosition();
	cyclone::vec3_t offset = cyclone::Vec3Construct(4.0f, 0, 0);
	offset = cyclone::Mat4TransformDirection( offset, sailboat.getTransform() );
	gluLookAt(pos.x+offset.x, pos.y+5.0f, pos.z+offset.z,  
		      pos.x, pos.y, pos.z,  
			  0.0, 1.0, 0.0);

	glColor3f(0.6f,0.6f,0.6f);
	int bx = int(pos.x);
	int bz = int(pos.z);
	glBegin(GL_QUADS);
	for (int x = -20; x <= 20; x++) for (int z = -20; z <= 20; z++)
	{
		glVertex3f(bx+x-0.1f, 0, bz+z-0.1f);
		glVertex3f(bx+x-0.1f, 0, bz+z+0.1f);
		glVertex3f(bx+x+0.1f, 0, bz+z+0.1f);
		glVertex3f(bx+x+0.1f, 0, bz+z-0.1f);
	}
	glEnd();

	// Set the transform matrix for the aircraft
	cyclone::mat4_t transform = sailboat.getTransform();
	cyclone::mat4_t gl_transform;
	gl_transform = cyclone::Mat4FillGLArray( transform );
	glPushMatrix();
	glMultMatrixf(gl_transform.n);

	// Draw the boat
	glColor3f(0,0,0);
	drawBoat();
	glPopMatrix();

	char buffer[256];
	sprintf(
		buffer, 
		"Speed %.1f", 
		cyclone::Vec3Magnitude( sailboat.getVelocity() )
		);
	glColor3f(0,0,0);
	renderText(10.0f, 24.0f, buffer);

	sprintf(
		buffer, 
		"Sail Control: %.1f",
		sail_control
		);
	renderText(10.0f, 10.0f, buffer);
}

void SailboatDemo::update()
{
	// Find the duration of the last frame in seconds
	float duration = (float)TimingData::get().lastFrameDuration * 0.001f;
	if (duration <= 0.0f) return;

	// Start with no forces or acceleration.
	sailboat.clearAccumulators();

	// Add the forces acting on the boat.
	registry.updateForces(duration);

	// Update the boat's physics.
	sailboat.integrate(duration);

	// Change the wind speed.
	windspeed = cyclone::Vec3Add( cyclone::Vec3Scale( windspeed, 0.9f ),
		r.randomXZVector( 1.0f ) );

	Application::update();
}

const char* SailboatDemo::getTitle()
{
	return "Cyclone > Sail Boat Demo";
}

void SailboatDemo::key(unsigned char key)
{
	switch(key)
	{
	case 'q': case 'Q':
		sail_control -= 0.1f;
		break;

	case 'e': case 'E':
		sail_control += 0.1f;
		break;

	case 'w': case 'W':
		sail_control = 0.0f;
		break;

	default:
		Application::key(key);
	}

	// Make sure the controls are in range
	if (sail_control < -1.0f) sail_control = -1.0f;
	else if (sail_control > 1.0f) sail_control = 1.0f;

	// Update the control surfaces
	//sail.setControl(sail_control);
}

/**
 * Called by the common demo framework to create an application
 * object (with new) and return a pointer.
 */
Application* getApplication()
{
    return new SailboatDemo();
}
