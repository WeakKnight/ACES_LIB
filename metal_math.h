#pragma once

#ifndef __METAL_STDLIB

#include <cmath>

#ifndef M_PI_F
#define M_PI_F 3.14159265358979323846264338327950288f   /* pi */
#endif

#ifndef HALF_MIN
#define HALF_MIN 6.103515625e-5
#endif

#ifndef HALF_MAX
#define HALF_MAX 65504.0
#endif

#define xyz get_swizzle(0, 1, 2)
#define xzy get_swizzle(0, 2, 1)
#define yxz get_swizzle(1, 0, 2)
#define yzx get_swizzle(1, 2, 0)
#define zxy get_swizzle(2, 0, 1)
#define zyx get_swizzle(2, 1, 0)

#define constant static
#define max fmax
#define min fmin

class float2;
class float3;
class float3x3;

class float2
{
public:
    float2(float x, float y);
    
    float& operator[](int index);
    
    float x, y;
};

class float3
{
public:
    float3();
    float3(float x, float y, float z);

    float& operator[](int index);
    float3 get_swizzle(int m0, int m1, int m2);

    float x,y,z;
};

class float2x2
{
public:
    float2x2(float2 c0, float2 c1);
    float2& operator[](int index);
    
    float2 column0;
    float2 column1;
};

class float3x3
{
public:
    float3x3();
    float3x3(float3 c0, float3 c1, float3 c2);
    float3& operator[](int index);
    float3x3 operator* (float3x3& m);
    
    float3 column0;
    float3 column1;
    float3 column2;
};

float3 operator*(float3 v, float3x3 m);
float3 operator*(float3x3 m, float3 v);

float3 operator*(float scaler, float3 v);
float3x3 operator*(float scaler, float3x3 m);

float dot(const float3& a, const float3& b);
float determinant(float3x3 m);
float3x3 transpose(float3x3 m);
float3x3 adjugate(float3x3 m);
float3x3 inverse(float3x3 m);
float clamp(float val, float minVal, float maxVal);
float3 clamp(float3 val, float minVal, float maxVal);
float sign( float x);

#else

float3x3 adjugate(float3x3 m)
{
    float a0 = +(m[1][1] * m[2][2] - m[2][1] * m[1][2]);
    float b0 = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]);
    float c0 = +(m[0][1] * m[1][2] - m[0][2] * m[1][1]);
    
    float a1 = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]);
    float b1 = +(m[0][0] * m[2][2] - m[0][2] * m[2][0]);
    float c1 = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]);
    
    float a2 = +(m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    float b2 = -(m[0][0] * m[2][1] - m[2][0] * m[0][1]);
    float c2 = +(m[0][0] * m[1][1] - m[1][0] * m[0][1]);
    
    return {{a0, b0, c0},{a1, b1, c1},{a2, b2, c2}};
}

float3x3 inverse(float3x3 m)
{
    float det_rcp = 1.0f / determinant(m);
    return det_rcp * adjugate(m);
}

#endif


