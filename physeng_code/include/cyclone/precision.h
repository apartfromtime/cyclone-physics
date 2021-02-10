/*
 * Interface file for code that changes when the core's precision is
 * altered.
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
 * Because Cyclone is designed to work at either single or double
 * precision, mathematical functions such as sqrt cannot be used
 * in the source code or headers. This file provides defines for
 * the real number type and mathematical formulae that work on it.
 *
 * @note All the contents of this file need to be changed to compile
 * Cyclone at a different precision.
 */

#ifndef CYCLONE_PRECISION_H
#define CYCLONE_PRECISION_H

#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>

///>VectorIntro
namespace cyclone {

#ifdef __cplusplus
extern "C" {
#endif          /* defined __cplusplus */

/** 
 * Defines we're in single precision mode, for any code
 * that needs to be conditionally compiled.
 */

/* default to single precision mode */
#if ( defined SINGLE_PRECISION ) || ( !defined DOUBLE_PRECISION )

///<SinglePrecision
typedef float real32_t;
typedef real32_t real_t;
typedef real32_t real;

/** Defines the highest value for the real number. */
#define REAL_MAX        FLT_MAX
/** Defines the precision of the square root operator. */
#define R_sqrt          sqrtf
/** Defines the precision of the absolute magnitude operator. */
#define R_abs           fabsf
/** Defines the precision of the sine operator. */
#define R_sin           sinf
/** Defines the precision of the cosine operator. */
#define R_cos           cosf
/** Defines the precision of the exponent operator. */
#define R_exp           expf
/** Defines the precision of the power operator. */
#define R_pow           powf
/** Defines the precision of the floating point modulo operator. */
#define R_mod           fmodf
/** Defines the precision of PI. */
#define R_PI            ( real32_t )M_PI     /* 3.1415927f */
///<SinglePrecision

#else

/** 
 * Defines we're in double precision mode, for any code
 * that needs to be conditionally compiled.
 */

///>DoublePrecision
typedef double real64_t;
typedef real64_t real_t;
typedef real64_t real;

#define REAL_MAX        DBL_MAX
#define R_sqrt          sqrt
#define R_abs           fabs
#define R_sin           sin
#define R_cos           cos
#define R_exp           exp
#define R_pow           pow
#define R_mod           fmod
#define R_PI            M_PI            /* 3.14159265358979 */
///<DoublePrecision

#endif /* #if ( defined SINGLE_PRECISION ) || ( !defined DOUBLE_PRECISION ) */

#ifdef __cplusplus
};
#endif /* defined __cplusplus */

} /* namespace cyclone */
///<VectorIntro

#endif /* CYCLONE_PRECISION_H */
