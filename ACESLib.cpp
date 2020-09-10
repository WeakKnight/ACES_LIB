#include "ACESLib.h"
#include <iostream>
#include <fstream>
#include <cstdint>
#include <stdio.h>
#include "OutputTransform.h"

static float* dataPtr = nullptr; 

void print_f3(float3 v)
{
    printf("X: %f Y: %f Z: %f\n", v.x, v.y, v.z);
}

float remapping_uv_to_linear(float v)
{
    // 56 [0, 0.875] [0, 0.2]
    // 8  [0.875, 1.0] [0.2, 1.0]
    return v;
}

float remapping_linear_to_uv(float v)
{
    return v;
}

extern "C"
{
    void Test_Val()
    {
        {
            float3 test = {50.0f, 60.0f, 70.0f};
            float3 res1 = ForwardWrappingPrecise(test);
            float3 res2 = ACESHDR10Standard(test, 15.0f);
            float3 backward = BackwardWrappingPrecise(res1);
            print_f3(res1);
            print_f3(res2);
            print_f3(backward);
        }

        {
            float3 test = {10.0f, 16.0f, 17.0f};
            float3 res1 = ForwardWrappingPrecise(test);
            float3 res2 = ACESHDR10Standard(test, 15.0f);
            float3 backward = BackwardWrappingPrecise(res1);
            print_f3(res1);
            print_f3(res2);
            print_f3(backward);
        }

        {
            float3 test = {0.2f, 1.1f, 0.3f};
            float3 res1 = ForwardWrappingPrecise(test);
            float3 res2 = ACESHDR10Standard(test, 15.0f);
            float3 backward = BackwardWrappingPrecise(res1);
            print_f3(res1);
            print_f3(res2);
            print_f3(backward);
        }

        {
            float3 test = {0.0f, 0.1f, 0.03f};
            float3 res1 = ForwardWrappingPrecise(test);
            float3 res2 = ACESHDR10Standard(test, 15.0f);
            float3 backward = BackwardWrappingPrecise(res1);
            print_f3(res1);
            print_f3(res2);
            print_f3(backward);
        }
    }

    void* GenLut(float mid, int lutSize)
    {
        if(dataPtr != nullptr)
        {
            delete dataPtr;
        }

        dataPtr = new float[lutSize * lutSize * lutSize * 4];

        size_t lutLength = lutSize * lutSize * lutSize;
        for(int x = 0; x < lutLength; x++)
        {
            // from 1d to 3d
            int i = int(x / (lutSize * lutSize ));
            int j = int((x - (i * lutSize * lutSize)) / lutSize);
            int k = x - (i * lutSize * lutSize + j * lutSize);
            
            float denominator = (float)lutSize - 1.0f;

            float rIn = float(k)/denominator;
            float gIn = float(j)/denominator;
            float bIn = float(i)/denominator;
            
            float3 cIn = float3(rIn, gIn, bIn);
            float3 lutValue = ComputeLutValue(cIn, mid);
            
            int index = (int)(i * lutSize * lutSize + j * lutSize + k);
            dataPtr[index * 4 + 0] = lutValue.x;
            dataPtr[index * 4 + 1] = lutValue.y;
            dataPtr[index * 4 + 2] = lutValue.z;
            dataPtr[index * 4 + 3] = 0.0f;
        }

        return (void*)dataPtr;
    }

    void* GenRRTLut(int lutSize)
    {
        if(dataPtr != nullptr)
        {
            delete dataPtr;
        }

        dataPtr = new float[lutSize * lutSize * lutSize * 4];

        size_t lutLength = lutSize * lutSize * lutSize;
        for(int x = 0; x < lutLength; x++)
        {
            // from 1d to 3d
            int i = int(x / (lutSize * lutSize ));
            int j = int((x - (i * lutSize * lutSize)) / lutSize);
            int k = x - (i * lutSize * lutSize + j * lutSize);
            
            float denominator = (float)lutSize - 1.0f;

            float rIn = float(k)/denominator;
            float gIn = float(j)/denominator;
            float bIn = float(i)/denominator;
            
            float3 cIn = float3(remapping_uv_to_linear(rIn), remapping_uv_to_linear(gIn), remapping_uv_to_linear(bIn));
            float3 lutValue = ComputeRRTLutValue(cIn);
            
            int index = (int)(i * lutSize * lutSize + j * lutSize + k);
            dataPtr[index * 4 + 0] = lutValue.x;
            dataPtr[index * 4 + 1] = lutValue.y;
            dataPtr[index * 4 + 2] = lutValue.z;
            dataPtr[index * 4 + 3] = 0.0f;
        }

        return (void*)dataPtr;
    }

    ACES_LIB_API void* GenLinearTestLut(int lutSize)
    {
        if(dataPtr != nullptr)
        {
            delete dataPtr;
        }

        dataPtr = new float[lutSize * lutSize * lutSize * 4];

        size_t lutLength = lutSize * lutSize * lutSize;
        for(int x = 0; x < lutLength; x++)
        {
            // from 1d to 3d
            int i = int(x / (lutSize * lutSize ));
            int j = int((x - (i * lutSize * lutSize)) / lutSize);
            int k = x - (i * lutSize * lutSize + j * lutSize);
            
            float denominator = (float)lutSize - 1.0f;

            float rIn = float(k)/denominator;
            float gIn = float(j)/denominator;
            float bIn = float(i)/denominator;
            
            float3 cIn = float3(rIn, gIn, bIn);
            float3 lutValue = cIn;
            
            int index = (int)(i * lutSize * lutSize + j * lutSize + k);
            dataPtr[index * 4 + 0] = lutValue.x;
            dataPtr[index * 4 + 1] = lutValue.y;
            dataPtr[index * 4 + 2] = lutValue.z;
            dataPtr[index * 4 + 3] = 0.0f;
        }

        return (void*)dataPtr;
    }
}
