//
//  metal_math.cpp
//  HDRViewer-macOS
//
//  Created by 李天宇 on 6/29/20.
//  Copyright © 2020 Apple. All rights reserved.
//

#include "metal_math.h"

//
// float2
//

float2::float2(float x, float y):
x(x),
y(y)
{
}

float& float2::operator[](int index)
{
    return reinterpret_cast<float*>(this)[index];
}

//
// float3
//

float3::float3()
:
x(0.0f),
y(0.0f),
z(0.0f)
{
}

float3::float3(float x, float y, float z):
x(x),
y(y),
z(z)
{
}

float& float3::operator[](int index)
{
    return reinterpret_cast<float*>(this)[index];
}

float3 float3::get_swizzle(int m0, int m1, int m2)
{
    return
    {
        reinterpret_cast<float*>(this)[m0],
        reinterpret_cast<float*>(this)[m1],
        reinterpret_cast<float*>(this)[m2]
    };
}

//
// float2x2
//
float2x2::float2x2(float2 c0, float2 c1):
column0(c0),
column1(c1)
{
}

float2& float2x2::operator[](int index)
{
    return reinterpret_cast<float2*>(this)[index];
}

//
// float3x3
//

float3x3::float3x3():
column0({}),
column1({}),
column2({})
{
}

float3x3::float3x3(float3 c0, float3 c1, float3 c2):
column0(c0),
column1(c1),
column2(c2)
{
}

float3& float3x3::operator[](int index)
{
    return reinterpret_cast<float3*>(this)[index];
}

float3x3 float3x3::operator* (float3x3& m)
{
    return {(*this) * m[0], (*this) * m[1], (*this) * m[2]};
}

//
// Arith
//

float3 operator*(float3 v, float3x3 m)
{
    return {dot(v, m[0]), dot(v, m[1]), dot(v, m[2])};
}

float3 operator*(float3x3 m, float3 v)
{
    return
    {
        m[0].x * v.x + m[1].x * v.y + m[2].x * v.z,
        m[0].y * v.x + m[1].y * v.y + m[2].y * v.z,
        m[0].z * v.x + m[1].z * v.y + m[2].z * v.z
    };
}

float3 operator*(float3 v1, float3 v2)
{
    return {v1.x * v2.x, v1.y * v2.y, v1.z * v2.z};
}

float3 operator*(float scaler, float3 v)
{
    return {scaler * v.x, scaler * v.y, scaler * v.z};
}

float3x3 operator*(float scaler, float3x3 m)
{
    return {
        {scaler * m.column0.x, scaler * m.column0.y, scaler * m.column0.z},
        {scaler * m.column1.x, scaler * m.column1.y, scaler * m.column1.z},
        {scaler * m.column2.x, scaler * m.column2.y, scaler * m.column2.z}};
}

//
// math function
//

float dot(const float3& a, const float3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

float determinant(float3x3 m)
{
    return m.column0[0] * m.column1[1] * m.column2[2]
    + m.column1[0] * m.column2[1] * m.column0[2]
    + m.column2[0] * m.column0[1] * m.column1[2]
    - m.column2[0] * m.column1[1] * m.column0[2]
    - m.column1[0] * m.column0[1] * m.column2[2]
    - m.column0[0] * m.column2[1] * m.column1[2];
}

/*
0,0 1,0 2,0      0,0 0,1 0,2
0,1 1,1 2,1      1,0 1,1 1,2
0,2 1,2 2,2      2,0 2,1 2,2
*/
float3x3 transpose(float3x3 m)
{
    return
    {
        {m.column0[0], m.column1[0], m.column2[0]},
        {m.column0[1], m.column1[1], m.column2[1]},
        {m.column0[2], m.column1[2], m.column2[2]}
    };
}

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

float clamp(float val, float minVal, float maxVal)
{
    if(val < minVal)
    {
        return minVal;
    }
    
    if(val > maxVal)
    {
        return maxVal;
    }
    
    return val;
}

float3 clamp(float3 val, float minVal, float maxVal)
{
    return {clamp(val.x, minVal, maxVal), clamp(val.y, minVal, maxVal), clamp(val.z, minVal, maxVal)};
}

float sign( float x)
{
    float y;
    if (x < 0.0)
    {
        y = -1.0;
    }
    else if (x > 0.0)
    {
        y = 1.0;
    }
    else
    {
        y = 0.0;
    }
    return y;
}

