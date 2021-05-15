/*
 * Interface file for the force generators.
 *
 * Part of the Cyclone physics system.
 *
 * Copyright (c) Icosagon 2003. All Rights Reserved.
 *
 * This software is distributed under licence. Use of this software
 * implies agreement with all terms and conditions of the accompanying
 * software licence.
 */

/**
 * @file
 *
 * This file contains the interface and sample force generators.
 */
#ifndef CYCLONE_PFGEN_H
#define CYCLONE_PFGEN_H

#include "core.h"
#include "particle.h"
///>ParticleFGRegistry
#include <vector>

namespace cyclone {
///<ParticleFGRegistry

///>ParticleFGInterface
    /**     
     * A force generator can be asked to add a force to one or more
     * particles.
     */
    class ParticleForceGenerator
    {
    public:

        /**
         * Overload this in implementations of the interface to calculate
         * and update the force applied to the given particle.
         */
        virtual void updateForce(Particle * particle, real_t duration) = 0;
    };
///<ParticleFGInterface

///>GravityPFG
    /**
     * A force generator that applies a gravitational force. One instance
     * can be used for multiple particles.
     */
    class ParticleGravity : public ParticleForceGenerator
    {
        /** Holds the acceleration due to gravity. */
        vec3_t gravity;

    public:

        /** Creates the generator with the given acceleration. */
        ParticleGravity(const vec3_t & gravity);

        /** Applies the gravitational force to the given particle. */
        virtual void updateForce(Particle * particle, real_t duration);
    };
///<GravityPFG

///>DragPFG
    /**
     * A force generator that applies a drag force. One instance
     * can be used for multiple particles.
     */
    class ParticleDrag : public ParticleForceGenerator
    {
        /** Holds the velocity drag coeffificent. */
        real_t k1;

        /** Holds the velocity squared drag coeffificent. */
        real_t k2;

    public:

        /** Creates the generator with the given coefficients. */
        ParticleDrag(real_t k1, real_t k2);

        /** Applies the drag force to the given particle. */
        virtual void updateForce(Particle * particle, real_t duration);
    };
///<DragPFG

///>ASpringPFG
	/**
     * A force generator that applies a Spring force, where
     * one end is attached to a fixed point in space.
     */
    class ParticleAnchoredSpring : public ParticleForceGenerator
    {
	protected:
        /** The location of the anchored end of the spring. */
        vec3_t * anchor;

        /** Holds the sprint constant. */
        real_t springConstant;

        /** Holds the rest length of the spring. */
        real_t restLength;

    public:
		ParticleAnchoredSpring();

        /** Creates a new spring with the given parameters. */
        ParticleAnchoredSpring(vec3_t * anchor, real_t springConstant, 
			real_t restLength);

		/** Retrieve the anchor point. */
		const vec3_t * getAnchor(void) const { return anchor; }

		/** Set the spring's properties. */
		void init(vec3_t * anchor, real_t springConstant,
            real_t restLength);

        /** Applies the spring force to the given particle. */
        virtual void updateForce(Particle * particle, real_t duration);
    };
///<ASpringPFG

	/**
	* A force generator that applies a bungee force, where
	* one end is attached to a fixed point in space.
	*/
	class ParticleAnchoredBungee : public ParticleAnchoredSpring
	{
	public:
		/** Applies the spring force to the given particle. */
		virtual void updateForce(Particle * particle, real_t duration);
	};

///>FakeSpringPFG
    /**
     * A force generator that fakes a stiff spring force, and where
     * one end is attached to a fixed point in space.
     */
    class ParticleFakeSpring : public ParticleForceGenerator
    {
        /** The location of the anchored end of the spring. */
        vec3_t * anchor;

        /** Holds the sprint constant. */
        real_t springConstant;

        /** Holds the damping on the oscillation of the spring. */
        real_t damping;

    public:

        /** Creates a new spring with the given parameters. */
        ParticleFakeSpring(vec3_t * anchor, real_t springConstant,
            real_t damping);

        /** Applies the spring force to the given particle. */
        virtual void updateForce(Particle * particle, real_t duration);
    };
///<FakeSpringPFG

///>SpringPFG
    /**
     * A force generator that applies a Spring force.
     */
    class ParticleSpring : public ParticleForceGenerator
    {
        /** The particle at the other end of the spring. */
        Particle * other;

        /** Holds the sprint constant. */
        real_t springConstant;

        /** Holds the rest length of the spring. */
        real_t restLength;

    public:

        /** Creates a new spring with the given parameters. */
        ParticleSpring(Particle * other, real_t springConstant,
            real_t restLength);

        /** Applies the spring force to the given particle. */
        virtual void updateForce(Particle * particle, real_t duration);
    };
///<SpringPFG

///>BungeePFG
    /**
     * A force generator that applies a spring force only
     * when extended.
     */
    class ParticleBungee : public ParticleForceGenerator
    {
        /** The particle at the other end of the spring. */
        Particle * other;

        /** Holds the sprint constant. */
        real_t springConstant;

        /**
         * Holds the length of the bungee at the point it begins to
         * generate a force.
         */
        real_t restLength;

    public:

        /** Creates a new bungee with the given parameters. */
        ParticleBungee(Particle * other, real_t springConstant,
            real_t restLength);

        /** Applies the spring force to the given particle. */
        virtual void updateForce(Particle * particle, real_t duration);
    };
///<BungeePFG

///>BuoyancyPFG
    /**
     * A force generator that applies a buoyancy force for a plane of
     * liquid parrallel to XZ plane.
     */
    class ParticleBuoyancy : public ParticleForceGenerator
    {
        /**
         * The maximum submersion depth of the object before
         * it generates its maximum boyancy force.
         */
        real_t maxDepth;

        /**
         * The volume of the object.
         */
        real_t volume;

        /**
         * The height of the water plane above y=0. The plane will be
         * parrallel to the XZ plane.
         */
        real_t waterHeight;

        /**
         * The density of the liquid. Pure water has a density of
         * 1000kg per cubic meter.
         */
        real_t liquidDensity;

    public:

        /** Creates a new buoyancy force with the given parameters. */
        ParticleBuoyancy(real_t maxDepth, real_t volume, real_t waterHeight,
            real_t liquidDensity = 1000.0f);

        /** Applies the buoyancy force to the given particle. */
        virtual void updateForce(Particle * particle, real_t duration);
    };
///<BuoyancyPFG

///>ParticleFGRegistry
    /**
     * Holds all the force generators and the particles they apply to.
     */
    class ParticleForceRegistry
    {
    protected:

        /**
         * Keeps track of one force generator and the particle it
         * applies to.
         */
        struct ParticleForceRegistration
        {
            Particle * particle;
            ParticleForceGenerator * fg;
        };

        /**
         * Holds the list of registrations.
         */
        typedef std::vector<ParticleForceRegistration> Registry;
        Registry registrations;

    public:
        /**
         * Registers the given force generator to apply to the
         * given particle.
         */
        void add(Particle * particle, ParticleForceGenerator * fg);

        /**
         * Removes the given registered pair from the registry.
         * If the pair is not registered, this method will have
         * no effect.
         */
        void remove(Particle * particle, ParticleForceGenerator * fg);

        /**
         * Clears all registrations from the registry. This will
         * not delete the particles or the force generators
         * themselves, just the records of their connection.
         */
        void clear();

        /**
         * Calls all the force generators to update the forces of
         * their corresponding particles.
         */
        void updateForces(real_t duration);
    };
}
///<ParticleFGRegistry

#endif // CYCLONE_PFGEN_H
