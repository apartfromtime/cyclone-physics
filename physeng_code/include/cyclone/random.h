/*
 * Interface file for the random number generator.
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
 * This file contains the definitions for a random number generator.
 */
#ifndef CYCLONE_RANDOM_H
#define CYCLONE_RANDOM_H

#include "core.h"

namespace cyclone {


	/**
	 * Keeps track of one random stream: i.e. a seed and its output.
	 * This is used to get random numbers. Rather than a function, this
	 * allows there to be several streams of repeatable random numbers
	 * at the same time. Uses the RandRotB algorithm.
	 */
	class Random
	{
	public:

		/** 
		 * Creates a new random number stream with a seed based on
		 * timing data.
		 */
		Random(void);

		/** 
		 * Creates a new random stream with the given seed.
		 */
		Random(unsigned seed);

		/**
		 * Sets the seed value for the random stream.
		 */
		void seed(unsigned seed);

		/**
		 * Returns the next random bitstring from the stream. This is
		 * the fastest method.
		 */
		unsigned randomBits(void);

		/**
		 * Returns a random floating point number between 0 and 1.
		 */
		real_t randomReal(void);

        /**
         * Returns a random floating point number between 0 and scale.
         */
        real_t randomReal(real_t scale);

        /**
         * Returns a random floating point number between min and max.
         */
        real_t randomReal(real_t min, real_t max);

		/**
		 * Returns a random integer less than the given value.
		 */
		unsigned randomInt(unsigned max);

		/**
		 * Returns a random binomially distributed number between -scale 
		 * and +scale.
		 */
		real_t randomBinomial(real_t scale);

		/**
		 * Returns a random vector where each component is binomially
		 * distributed in the range (-scale to scale) [mean = 0.0f].
		 */
		vec3_t randomVector(real_t scale);

		/**
		 * Returns a random vector where each component is binomially
		 * distributed in the range (-scale to scale) [mean = 0.0f],
		 * where scale is the corresponding component of the given
		 * vector.
		 */
		vec3_t randomVector(const vec3_t & scale);

        /**
         * Returns a random vector in the cube defined by the given
         * minimum and maximum vectors. The probability is uniformly
         * distributed in this region.
         */
        vec3_t randomVector(const vec3_t & min, const vec3_t & max);

		/**
		 * Returns a random vector where each component is binomially
		 * distributed in the range (-scale to scale) [mean = 0.0f],
		 * except the y coordinate which is zero.
		 */
		vec3_t randomXZVector(real_t scale);

        /**
         * Returns a random orientation (i.e. normalized) quaternion.
         */
        quat_t randomQuaternion(void);

	private:
		// Internal mechanics
		int p1, p2;
		unsigned buffer[17];
	};

} // namespace cyclone

#endif // CYCLONE_BODY_H
