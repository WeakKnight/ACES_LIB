#include "ACESLib.h"
#include <iostream>
#include <fstream>
#include <cstdint>
#include "OutputTransform.h"

static float* dataPtr = nullptr; 

extern "C"
{
    void* GenLut(float mid)
    {
        if(dataPtr != nullptr)
        {
            delete dataPtr;
        }

        dataPtr = new float[64 * 64 * 64 * 4];

        size_t lutLength = 64 * 64 * 64;
        for(int x = 0; x < lutLength; x++)
        {
            // from 1d to 3d
            int i = int(x / (64 * 64 ));
            int j = int((x - (i * 64 * 64)) / 64);
            int k = x - (i * 64 * 64 + j * 64);
            
            float rIn = float(k)/63.0f;
            float gIn = float(j)/63.0f;
            float bIn = float(i)/63.0f;
            
            float3 cIn = float3(rIn, gIn, bIn);
            float3 lutValue = ComputeLutValue(cIn, mid);
            
            int index = (int)(i * 64 * 64 + j * 64 + k);
            dataPtr[index * 4 + 0] = lutValue.x;
            dataPtr[index * 4 + 1] = lutValue.y;
            dataPtr[index * 4 + 2] = lutValue.z;
            dataPtr[index * 4 + 3] = 0.0f;
        }

        return (void*)dataPtr;
    }
}
