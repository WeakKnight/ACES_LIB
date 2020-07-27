#pragma once
#define ACES_LIB_API __declspec(dllexport) 

extern "C" 
{
    ACES_LIB_API void* GenLut(float mid);
}