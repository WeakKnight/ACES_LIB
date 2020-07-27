#include <iostream>
#include "OutputTransform.h"

extern "C"
{
    typedef struct vec3
    {
        float x;
        float y;
        float z;
    } vec3;

    vec3 ComputeLutValue(float rpq, float gpq, float bpq, float mid)
    {
        float3 result;
        float3 pq = float3(rpq, gpq, bpq);
        float3 lin = ST2084_2_Y_f3(pq);
        
        const float Y_MIN = 0.0001;                     // black luminance (cd/m^2)
        const float Y_MAX = 1000.0;                     // peak white luminance (cd/m^2)
        
        const Chromaticities DISPLAY_PRI = REC2020_PRI; // encoding primaries (device setup)
        const Chromaticities LIMITING_PRI = REC2020_PRI;// limiting primaries
        
        result = OutputTransform(lin, Y_MIN, mid, Y_MAX, DISPLAY_PRI, LIMITING_PRI);
        vec3 res;
        res.x = result.x;
        res.y = result.y;
        res.z = result.z;
        
        return res;
    }
}
// void say_hello(){
//     std::cout << "Hello, from ACESLib!\n";
// }
