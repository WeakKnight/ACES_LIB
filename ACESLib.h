#pragma once
#define ACES_LIB_API __declspec(dllexport) 

extern "C" 
{
    ACES_LIB_API void* GenLut(float mid, int lutSize);
    ACES_LIB_API void* GenRRTLut(int lutSize);
    ACES_LIB_API void* GenLinearTestLut(int lutSize);
    ACES_LIB_API void Test_Val();
}

