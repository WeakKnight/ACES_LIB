#include <iostream>
#include <fstream>
#include <cstdint>
#include "OutputTransform.h"

extern "C"
{
    void GenLut(const char* path, float mid)
    {
        float data[64 * 64 * 64 * 4];
        float *dataPtr = &data[0];

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

        std::fstream file;
        file.open(path, std::ios::app | std::ios::binary);
        file.write(reinterpret_cast<char*>(dataPtr), sizeof(float) * 64 * 64 * 64 * 4); // ideally, you should memcpy it to a char buffer.
        file.close();
    }
}
