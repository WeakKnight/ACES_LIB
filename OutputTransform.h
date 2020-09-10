//
//  OutputTransform.h
//  HDRViewer-macOS
//
//  Created by 李天宇 on 6/29/20.
//  Copyright © 2020 Apple. All rights reserved.
//

#ifndef OutputTransform_h
#define OutputTransform_h

#include "metal_math.h"

/* ---- Chromaticities of some common primary sets ---- */

constant float3x3 M1 =
{
  {  0.5, -1.0, 0.5 },
  { -1.0,  1.0, 0.5 },
  {  0.5,  0.0, 0.0 }
};

struct Chromaticities
{
    float2 R,G,B,White;
};

constant Chromaticities AP0 = // ACES Primaries from SMPTE ST2065-1
{
  { 0.73470,  0.26530},
  { 0.00000,  1.00000},
  { 0.00010, -0.07700},
  { 0.32168,  0.33767}
};

constant Chromaticities AP1 = // Working space and rendering primaries for ACES 1.0
{
  { 0.713,    0.293},
  { 0.165,    0.830},
  { 0.128,    0.044},
  { 0.32168,  0.33767}
};

constant Chromaticities REC2020_PRI =
{
  { 0.70800,  0.29200},
  { 0.17000,  0.79700},
  { 0.13100,  0.04600},
  { 0.31270,  0.32900}
};

constant float3x3 AP0_2_AP1_MAT = //mul( AP0_2_XYZ_MAT, XYZ_2_AP1_MAT );
{
    {1.4514393161, -0.2365107469, -0.2149285693},
    {-0.0765537734,  1.1762296998, -0.0996759264},
    {0.0083161484, -0.0060324498,  0.9977163014}
};

constant float3x3 AP1_2_XYZ_MAT =
{
    {0.6624541811, 0.1340042065, 0.1561876870},
    {0.2722287168, 0.6740817658, 0.0536895174},
    {-0.0055746495, 0.0040607335, 1.0103391003},
};

constant float3 AP1_RGB2Y =
{
    AP1_2_XYZ_MAT[1][0],
    AP1_2_XYZ_MAT[1][1],
    AP1_2_XYZ_MAT[1][2]
};

constant float3x3 D60_2_D65_CAT =
{
    {0.987224,   -0.00611327, 0.0159533},
    {-0.00759836,  1.00186,    0.00533002},
    {0.00307257, -0.00509595, 1.08168}
};

float3 XYZ_2_xyY(float3 XYZ)
{
    float3 xyY;
    float divisor = (XYZ[0] + XYZ[1] + XYZ[2]);
    if (divisor == 0.) divisor = 1e-10;
    xyY[0] = XYZ[0] / divisor;
    xyY[1] = XYZ[1] / divisor;
    xyY[2] = XYZ[1];

    return xyY;
}

float3 xyY_2_XYZ( float3 xyY )
{
    float3 XYZ;
    XYZ[0] = xyY[0] * xyY[2] / max( xyY[1], 1e-10);
    XYZ[1] = xyY[2];
    XYZ[2] = (1.0 - xyY[0] - xyY[1]) * xyY[2] / max( xyY[1], 1e-10);

    return XYZ;
}

float3x3 RGBtoXYZ(Chromaticities N)
{
    float3x3 M = float3x3(
                      float3(N.R.x, N.R.y, 1.0 - (N.R.x + N.R.y)),
                      float3(N.G.x, N.G.y, 1.0 - (N.G.x + N.G.y)),
                      float3(N.B.x, N.B.y, 1.0 - (N.B.x + N.B.y))
                      );
    
    float3 wh = float3(N.White.x / N.White.y, 1.0, (1.0 - (N.White.x + N.White.y)) / N.White.y);
    float3x3 invM = inverse(M);
    wh = invM * wh;
    float3x3 WH = float3x3(
                           float3(wh.x, 0.0, 0.0),
                           float3(0.0, wh.y, 0.0),
                           float3(0.0, 0.0, wh.z)
                           );
    M = M * WH;
    return M;
}

float3x3 XYZtoRGB( Chromaticities N)
{
    float3x3 M = inverse(RGBtoXYZ(N));
    return M;
}

struct TsPoint
{
    float x;        // ACES
    float y;        // luminance
    float slope;    //
};

struct float5
{
    float a,b,c,d,e;
};

struct float6
{
    float a,b,c,d,e,f;
};

struct TsParams
{
    TsPoint Min;
    TsPoint Mid;
    TsPoint Max;
    float6 coefsLow;
    float6 coefsHigh;
};

constant float MIN_STOP_SDR = -6.5;
constant float MAX_STOP_SDR = 6.5;

constant float MIN_STOP_RRT = -15.;
constant float MAX_STOP_RRT = 18.;

constant float MIN_LUM_SDR = 0.02;
constant float MAX_LUM_SDR = 48.0;

constant float MIN_LUM_RRT = 0.0001;
constant float MAX_LUM_RRT = 10000.0;

float interpolate1D( float2x2 table, float p)
{
    if( p <= table[0][0] )
    {
        return table[0][1];
    }
    
    if( p >= table[1][0] )
    {
        return table[1][1];
    }
    
    if( table[0][0] <= p && p < table[1][0])
    {
        float s = (p - table[0][0]) / (table[1][0] - table[0][0]);
        return table[0][1] * ( 1.0 - s ) + table[1][1] * s;
    }
    
    return 0.0;
}

float lookup_ACESmin( float minLum )
{
    const float2x2 minTable =
    {
        { log10(MIN_LUM_RRT), MIN_STOP_RRT },
        { log10(MIN_LUM_SDR), MIN_STOP_SDR }
    };

    return 0.18*pow( 2., interpolate1D( minTable, log10( minLum)));
}

float lookup_ACESmax( float maxLum )
{
    const float2x2 maxTable =
    {
        { log10(MAX_LUM_SDR), MAX_STOP_SDR },
        { log10(MAX_LUM_RRT), MAX_STOP_RRT }
    };

    return 0.18*pow( 2., interpolate1D( maxTable, log10( maxLum)));
}

float5 init_coefsLow(
    TsPoint TsPointLow,
    TsPoint TsPointMid
)
{
    float5 coefsLow;

    float knotIncLow = (log10(TsPointMid.x) - log10(TsPointLow.x)) / 3.;
    // float halfKnotInc = (log10(TsPointMid.x) - log10(TsPointLow.x)) / 6.;

    // Determine two lowest coefficients (straddling minPt)
    coefsLow.a = (TsPointLow.slope * (log10(TsPointLow.x)-0.5*knotIncLow)) + ( log10(TsPointLow.y) - TsPointLow.slope * log10(TsPointLow.x));
    coefsLow.b = (TsPointLow.slope * (log10(TsPointLow.x)+0.5*knotIncLow)) + ( log10(TsPointLow.y) - TsPointLow.slope * log10(TsPointLow.x));
    // NOTE: if slope=0, then the above becomes just
        // coefsLow[0] = log10(TsPointLow.y);
        // coefsLow[1] = log10(TsPointLow.y);
    // leaving it as a variable for now in case we decide we need non-zero slope extensions

    // Determine two highest coefficients (straddling midPt)
    coefsLow.d = (TsPointMid.slope * (log10(TsPointMid.x)-0.5*knotIncLow)) + ( log10(TsPointMid.y) - TsPointMid.slope * log10(TsPointMid.x));
    coefsLow.e = (TsPointMid.slope * (log10(TsPointMid.x)+0.5*knotIncLow)) + ( log10(TsPointMid.y) - TsPointMid.slope * log10(TsPointMid.x));
    
    // Middle coefficient (which defines the "sharpness of the bend") is linearly interpolated
    float2x2 bendsLow =
    {
        {MIN_STOP_RRT, 0.18},
        {MIN_STOP_SDR, 0.35}
    };
    
    float pctLow = interpolate1D( bendsLow, log2(TsPointLow.x/0.18));
    coefsLow.c = log10(TsPointLow.y) + pctLow*(log10(TsPointMid.y)-log10(TsPointLow.y));

    return coefsLow;
}

float5 init_coefsHigh(
    TsPoint TsPointMid,
    TsPoint TsPointMax
)
{
    float5 coefsHigh;

    float knotIncHigh = (log10(TsPointMax.x) - log10(TsPointMid.x)) / 3.;
    // float halfKnotInc = (log10(TsPointMax.x) - log10(TsPointMid.x)) / 6.;

    // Determine two lowest coefficients (straddling midPt)
    coefsHigh.a = (TsPointMid.slope * (log10(TsPointMid.x)-0.5*knotIncHigh)) + ( log10(TsPointMid.y) - TsPointMid.slope * log10(TsPointMid.x));
    coefsHigh.b = (TsPointMid.slope * (log10(TsPointMid.x)+0.5*knotIncHigh)) + ( log10(TsPointMid.y) - TsPointMid.slope * log10(TsPointMid.x));

    // Determine two highest coefficients (straddling maxPt)
    coefsHigh.d = (TsPointMax.slope * (log10(TsPointMax.x)-0.5*knotIncHigh)) + ( log10(TsPointMax.y) - TsPointMax.slope * log10(TsPointMax.x));
    coefsHigh.e = (TsPointMax.slope * (log10(TsPointMax.x)+0.5*knotIncHigh)) + ( log10(TsPointMax.y) - TsPointMax.slope * log10(TsPointMax.x));
    // NOTE: if slope=0, then the above becomes just
        // coefsHigh[0] = log10(TsPointHigh.y);
        // coefsHigh[1] = log10(TsPointHigh.y);
    // leaving it as a variable for now in case we decide we need non-zero slope extensions
    
    // Middle coefficient (which defines the "sharpness of the bend") is linearly interpolated
    float2x2 bendsHigh =
    {
        {MAX_STOP_SDR, 0.89},
        {MAX_STOP_RRT, 0.90}
    };
    
    float pctHigh = interpolate1D( bendsHigh, log2(TsPointMax.x/0.18));
    coefsHigh.c = log10(TsPointMid.y) + pctHigh*(log10(TsPointMax.y)-log10(TsPointMid.y));
    
    return coefsHigh;
}

float shift( float input, float expShift)
{
    return pow(2.,(log2(input)-expShift));
}

float pow10( float x)
{
    return pow(10.0, x);
}

float inv_ssts
(
    float y,
    TsParams C
)
{
    const int N_KNOTS_LOW = 4;
    const int N_KNOTS_HIGH = 4;

    const float KNOT_INC_LOW = (log10(C.Mid.x) - log10(C.Min.x)) / (N_KNOTS_LOW - 1.);
    const float KNOT_INC_HIGH = (log10(C.Max.x) - log10(C.Mid.x)) / (N_KNOTS_HIGH - 1.);

    // KNOT_Y is luminance of the spline at each knot
    float KNOT_Y_LOW[ N_KNOTS_LOW];
    KNOT_Y_LOW[0] = ( C.coefsLow.a + C.coefsLow.b) / 2.;
    KNOT_Y_LOW[1] = ( C.coefsLow.b + C.coefsLow.c) / 2.;
    KNOT_Y_LOW[2] = ( C.coefsLow.c + C.coefsLow.d) / 2.;
    KNOT_Y_LOW[3] = ( C.coefsLow.d + C.coefsLow.e) / 2.;

    float KNOT_Y_HIGH[ N_KNOTS_HIGH];
    KNOT_Y_HIGH[0] = ( C.coefsHigh.a + C.coefsHigh.b) / 2.;
    KNOT_Y_HIGH[1] = ( C.coefsHigh.b + C.coefsHigh.c) / 2.;
    KNOT_Y_HIGH[2] = ( C.coefsHigh.c + C.coefsHigh.d) / 2.;
    KNOT_Y_HIGH[3] = ( C.coefsHigh.d + C.coefsHigh.e) / 2.;

    float logy = log10( max(y,1e-10));

    float logx;
    if (logy <= log10(C.Min.y))
    {
        logx = log10(C.Min.x);
    }
    else if ( (logy > log10(C.Min.y)) && (logy <= log10(C.Mid.y)) )
    {
        unsigned int j = 0;
        float3 cf = {0.0f, 0.0f, 0.0f};
        if ( logy > KNOT_Y_LOW[ 0] && logy <= KNOT_Y_LOW[ 1])
        {
            cf[0] = C.coefsLow.a;
            cf[ 1] = C.coefsLow.b;
            cf[ 2] = C.coefsLow.c;
            j = 0;
        }
        else if ( logy > KNOT_Y_LOW[ 1] && logy <= KNOT_Y_LOW[ 2])
        {
            cf[0] = C.coefsLow.b;
            cf[1] = C.coefsLow.c;
            cf[2] = C.coefsLow.d;
            j = 1;
        }
        else if ( logy > KNOT_Y_LOW[ 2] && logy <= KNOT_Y_LOW[ 3])
        {
            cf[0] = C.coefsLow.c;
            cf[1] = C.coefsLow.d;
            cf[2] = C.coefsLow.e;
            j = 2;
        }

        float3 tmp = M1 * cf;

        float a = tmp[ 0];
        float b = tmp[ 1];
        float c = tmp[ 2];
        c = c - logy;

        float d = sqrt( b * b - 4. * a * c);

        float t = ( 2. * c) / ( -d - b);

        logx = log10(C.Min.x) + ( t + j) * KNOT_INC_LOW;

    }
    else if ( (logy > log10(C.Mid.y)) && (logy < log10(C.Max.y)) )
    {
        unsigned int j = 0;
        float3 cf = {0.0f, 0.0f, 0.0f};
        if ( logy >= KNOT_Y_HIGH[ 0] && logy <= KNOT_Y_HIGH[ 1])
        {
            cf[ 0] = C.coefsHigh.a;
            cf[ 1] = C.coefsHigh.b;
            cf[ 2] = C.coefsHigh.c;
            j = 0;
        }
        else if ( logy > KNOT_Y_HIGH[ 1] && logy <= KNOT_Y_HIGH[ 2])
        {
            cf[ 0] = C.coefsHigh.b;
            cf[ 1] = C.coefsHigh.c;
            cf[ 2] = C.coefsHigh.d;
            j = 1;
        }
        else if ( logy > KNOT_Y_HIGH[ 2] && logy <= KNOT_Y_HIGH[ 3])
        {
            cf[ 0] = C.coefsHigh.c;
            cf[ 1] = C.coefsHigh.d;
            cf[ 2] = C.coefsHigh.e;
            j = 2;
        }

        float3 tmp = M1 * cf;

        float a = tmp[ 0];
        float b = tmp[ 1];
        float c = tmp[ 2];
        c = c - logy;

        float d = sqrt( b * b - 4. * a * c);

        float t = ( 2. * c) / ( -d - b);

        logx = log10(C.Mid.x) + ( t + j) * KNOT_INC_HIGH;

    }
    else
    { //if ( logy >= log10(C.Max.y) ) {
        logx = log10(C.Max.x);
    }

    return pow10( logx);
}

TsParams init_TsParams(
    float minLum,
    float maxLum,
    float expShift = 0
)
{
    TsPoint MIN_PT = { lookup_ACESmin(minLum), minLum, 0.0};
    TsPoint MID_PT = { 0.18, 4.8, 1.55};
    TsPoint MAX_PT = { lookup_ACESmax(maxLum), maxLum, 0.0};
    
    float5 cLow = init_coefsLow( MIN_PT, MID_PT);
    float5 cHigh = init_coefsHigh( MID_PT, MAX_PT);
    
    MIN_PT.x = shift(lookup_ACESmin(minLum),expShift);
    MID_PT.x = shift(0.18,expShift);
    MAX_PT.x = shift(lookup_ACESmax(maxLum),expShift);
    
    TsParams P = {
        {MIN_PT.x, MIN_PT.y, MIN_PT.slope},
        {MID_PT.x, MID_PT.y, MID_PT.slope},
        {MAX_PT.x, MAX_PT.y, MAX_PT.slope},
        {cLow.a, cLow.b, cLow.c, cLow.d, cLow.e, cLow.e},
        {cHigh.a, cHigh.b, cHigh.c, cHigh.d, cHigh.e, cHigh.e}
    };
         
    return P;
}

float rgb_2_saturation( float3 rgb )
{
    float minrgb = min( min(rgb.x, rgb.y ), rgb.z );
    float maxrgb = max( max(rgb.x, rgb.y ), rgb.z );
    return ( max( maxrgb, 1e-10 ) - max( minrgb, 1e-10 ) ) / max( maxrgb, 1e-2 );
}

float rgb_2_yc( float3 rgb, float ycRadiusWeight = 1.75)
{
    // Converts RGB to a luminance proxy, here called YC
    // YC is ~ Y + K * Chroma
    // Constant YC is a cone-shaped surface in RGB space, with the tip on the
    // neutral axis, towards white.
    // YC is normalized: RGB 1 1 1 maps to YC = 1
    //
    // ycRadiusWeight defaults to 1.75, although can be overridden in function
    // call to rgb_2_yc
    // ycRadiusWeight = 1 -> YC for pure cyan, magenta, yellow == YC for neutral
    // of same value
    // ycRadiusWeight = 2 -> YC for pure red, green, blue  == YC for  neutral of
    // same value.

    float r = rgb[0];
    float g = rgb[1];
    float b = rgb[2];
  
    float chroma = sqrt(b*(b-g)+g*(g-r)+r*(r-b));

    return ( b + g + r + ycRadiusWeight * chroma) / 3.;
}

float sigmoid_shaper( float x)
{
    // Sigmoid function in the range 0 to 1 spanning -2 to +2.

    float t = max( 1.0 - abs( 0.5 * x ), 0.0 );
    float y = 1.0 + sign(x) * (1.0 - t*t);
    return 0.5 * y;
}

// ------- Glow module functions
float glow_fwd( float ycIn, float glowGainIn, float glowMid)
{
   float glowGainOut;

   if (ycIn <= 2./3. * glowMid)
   {
       glowGainOut = glowGainIn;
   }
   else if ( ycIn >= 2. * glowMid)
   {
       glowGainOut = 0.;
   }
   else
   {
       glowGainOut = glowGainIn * (glowMid / ycIn - 1./2.);
   }

   return glowGainOut;
}

// Transformations from RGB to other color representations
float rgb_2_hue( float3 rgb )
{
    // Returns a geometric hue angle in degrees (0-360) based on RGB values.
    // For neutral colors, hue is undefined and the function will return a quiet NaN value.
    float hue;
    if (rgb[0] == rgb[1] && rgb[1] == rgb[2])
    {
        //hue = FLT_NAN; // RGB triplets where RGB are equal have an undefined hue
        hue = 0;
    }
    else
    {
        hue = (180. / M_PI_F) * atan2(sqrt(3.0)*(rgb[1] - rgb[2]), 2 * rgb[0] - rgb[1] - rgb[2]);
    }

    if (hue < 0.)
        hue = hue + 360;

    return clamp( hue, 0., 360. );
}

float center_hue( float hue, float centerH)
{
    float hueCentered = hue - centerH;
    if (hueCentered < -180.)
        hueCentered += 360;
    else if (hueCentered > 180.)
        hueCentered -= 360;
    return hueCentered;
}

// ------- Red modifier functions
float cubic_basis_shaper
(
  float x,
  float w   // full base width of the shaper function (in degrees)
)
{
    //return Square( smoothstep( 0, 1, 1 - abs( 2 * x/w ) ) );

    float M[4][4] =
    {
        { -1./6,  3./6, -3./6,  1./6 },
        {  3./6, -6./6,  3./6,  0./6 },
        { -3./6,  0./6,  3./6,  0./6 },
        {  1./6,  4./6,  1./6,  0./6 }
    };
  
    float knots[5] = { -0.5f * w, -0.25f * w, 0.0f, 0.25f * w, 0.5f * w };
  
    float y = 0;
    if ((x > knots[0]) && (x < knots[4]))
    {
        float knot_coord = (x - knots[0]) * 4.0 / w;
        int j = knot_coord;
        float t = knot_coord - j;
      
        float monomials[4] = { t*t*t, t*t, t, 1.0 };

        // (if/else structure required for compatibility with CTL < v1.5.)
        if ( j == 3) {
            y = monomials[0] * M[0][0] + monomials[1] * M[1][0] +
                monomials[2] * M[2][0] + monomials[3] * M[3][0];
        } else if ( j == 2) {
            y = monomials[0] * M[0][1] + monomials[1] * M[1][1] +
                monomials[2] * M[2][1] + monomials[3] * M[3][1];
        } else if ( j == 1) {
            y = monomials[0] * M[0][2] + monomials[1] * M[1][2] +
                monomials[2] * M[2][2] + monomials[3] * M[3][2];
        } else if ( j == 0) {
            y = monomials[0] * M[0][3] + monomials[1] * M[1][3] +
                monomials[2] * M[2][3] + monomials[3] * M[3][3];
        } else {
            y = 0.0;
        }
    }
  
    return y * 1.5;
}


float3x3 calc_sat_adjust_matrix(float sat, float3 rgb2Y)
{
  //
  // This function determines the terms for a 3x3 saturation matrix that is
  // based on the luminance of the input.
  //
  float3x3 M;
  M[0][0] = (1.0 - sat) * rgb2Y[0] + sat;
  M[1][0] = (1.0 - sat) * rgb2Y[0];
  M[2][0] = (1.0 - sat) * rgb2Y[0];
  
  M[0][1] = (1.0 - sat) * rgb2Y[1];
  M[1][1] = (1.0 - sat) * rgb2Y[1] + sat;
  M[2][1] = (1.0 - sat) * rgb2Y[1];
  
  M[0][2] = (1.0 - sat) * rgb2Y[2];
  M[1][2] = (1.0 - sat) * rgb2Y[2];
  M[2][2] = (1.0 - sat) * rgb2Y[2] + sat;

  M = transpose(M);
  return M;
}


float3 rrt_sweeteners(float3 aces)
{
    const float RRT_GLOW_GAIN = 0.05;
    const float RRT_GLOW_MID = 0.08;
    
    // Red modifier constants
    const float RRT_RED_SCALE = 0.82;
    const float RRT_RED_PIVOT = 0.03;
    const float RRT_RED_HUE = 0.;
    const float RRT_RED_WIDTH = 135.;
    
    const float RRT_SAT_FACTOR = 0.96;
    float3x3 RRT_SAT_MAT = calc_sat_adjust_matrix( RRT_SAT_FACTOR, AP1_RGB2Y);
    
    // --- Glow module --- //
    float saturation = rgb_2_saturation( aces);
    float ycIn = rgb_2_yc( aces);
    float s = sigmoid_shaper( (saturation - 0.4) / 0.2);
    float addedGlow = 1. + glow_fwd( ycIn, RRT_GLOW_GAIN * s, RRT_GLOW_MID);

    aces = addedGlow * aces;

    // --- Red modifier --- //
    float hue = rgb_2_hue( aces);
    float centeredHue = center_hue( hue, RRT_RED_HUE);
    float hueWeight = cubic_basis_shaper( centeredHue, RRT_RED_WIDTH);

    aces[0] = aces[0] + hueWeight * saturation * (RRT_RED_PIVOT - aces[0]) * (1. - RRT_RED_SCALE);

    // --- ACES to RGB rendering space --- //
    aces = clamp( aces, 0., 65535.);
    float3 rgbPre = aces * AP0_2_AP1_MAT;
    rgbPre = clamp( rgbPre, 0., 65504.);
    
    // --- Global desaturation --- //
    rgbPre = RRT_SAT_MAT * rgbPre;
    return rgbPre;
}

float ssts
(
    float x,
    TsParams C
)
{
    const int N_KNOTS_LOW = 4;
    const int N_KNOTS_HIGH = 4;

    // Check for negatives or zero before taking the log. If negative or zero,
    // set to HALF_MIN.
    float logx = log10( max(x, 6.103515625e-5 ));

    float logy;

    if ( logx <= log10(C.Min.x) )
    {
        logy = logx * C.Min.slope + ( log10(C.Min.y) - C.Min.slope * log10(C.Min.x) );
    }
    else if (( logx > log10(C.Min.x) ) && ( logx < log10(C.Mid.x) ))
    {
        float knot_coord = (N_KNOTS_LOW-1) * (logx-log10(C.Min.x))/(log10(C.Mid.x)-log10(C.Min.x));
        int j = knot_coord;
        float t = knot_coord - j;

        float coefs[6] = {C.coefsLow.a, C.coefsLow.b, C.coefsLow.c, C.coefsLow.d, C.coefsLow.e, C.coefsLow.f};
        float3 cf = { coefs[ j], coefs[ j + 1], coefs[ j + 2]};

        float3 monomials = { t * t, t, 1. };
        logy = dot( monomials,  M1 * cf);

    }
    else if (( logx >= log10(C.Mid.x) ) && ( logx < log10(C.Max.x) ))
    {

        float knot_coord = (N_KNOTS_HIGH-1) * (logx-log10(C.Mid.x))/(log10(C.Max.x)-log10(C.Mid.x));
        int j = knot_coord;
        float t = knot_coord - j;
        
        float coefs[6] = {C.coefsHigh.a, C.coefsHigh.b, C.coefsHigh.c, C.coefsHigh.d, C.coefsHigh.e, C.coefsHigh.f};
        float3 cf = { coefs[ j], coefs[ j + 1], coefs[ j + 2]};

        float3 monomials = { t * t, t, 1. };
        logy = dot(monomials,  M1 * cf);
    }
    else
    {
        //if ( logIn >= log10(C.Max.x) ) {
        logy = logx * C.Max.slope + ( log10(C.Max.y) - C.Max.slope * log10(C.Max.x) );
    }

    return pow10(logy);

}

float3 ssts_f3
(
    float3 x,
    TsParams C
)
{
    float3 res;
    res[0] = ssts( x[0], C);
    res[1] = ssts( x[1], C);
    res[2] = ssts( x[2], C);

    return res;
}

float Y_2_linCV(float Y, float Ymax, float Ymin)
{
    return (Y - Ymin) / (Ymax - Ymin);
}

float3 Y_2_linCV_f3( float3 Y, float Ymax, float Ymin)
{
    float3 linCV;
    linCV[0] = Y_2_linCV( Y[0], Ymax, Ymin);
    linCV[1] = Y_2_linCV( Y[1], Ymax, Ymin);
    linCV[2] = Y_2_linCV( Y[2], Ymax, Ymin);
    return linCV;
}

float3 limit_to_primaries
(
    float3 XYZ,
    Chromaticities LIMITING_PRI
)
{
    float3x3 XYZ_2_LIMITING_PRI_MAT = XYZtoRGB( LIMITING_PRI);
    float3x3 LIMITING_PRI_2_XYZ_MAT = RGBtoXYZ( LIMITING_PRI);

    // XYZ to limiting primaries
    float3 rgb =  XYZ_2_LIMITING_PRI_MAT * XYZ;

    // Clip any values outside the limiting primaries
    float3 limitedRgb = clamp( rgb, 0., 1.);
    
    // Convert limited RGB to XYZ
    return LIMITING_PRI_2_XYZ_MAT * limitedRgb;
}

float linCV_2_Y( float linCV, float Ymax, float Ymin)
{
    return linCV * (Ymax - Ymin) + Ymin;
}

float3 linCV_2_Y_f3( float3 linCV, float Ymax, float Ymin)
{
    float3 Y;
    Y[0] = linCV_2_Y( linCV[0], Ymax, Ymin);
    Y[1] = linCV_2_Y( linCV[1], Ymax, Ymin);
    Y[2] = linCV_2_Y( linCV[2], Ymax, Ymin);
    return Y;
}

// Converts from linear cd/m^2 to the non-linear perceptually quantized space
// Note that this is in float, and assumes normalization from 0 - 1
// (0 - pq_C for linear) and does not handle the integer coding in the Annex
// sections of SMPTE ST 2084-2014
float Y_2_ST2084( float C )
//pq_r
{
    // Constants from SMPTE ST 2084-2014
    const float pq_m1 = 0.1593017578125; // ( 2610.0 / 4096.0 ) / 4.0;
    const float pq_m2 = 78.84375; // ( 2523.0 / 4096.0 ) * 128.0;
    const float pq_c1 = 0.8359375; // 3424.0 / 4096.0 or pq_c3 - pq_c2 + 1.0;
    const float pq_c2 = 18.8515625; // ( 2413.0 / 4096.0 ) * 32.0;
    const float pq_c3 = 18.6875; // ( 2392.0 / 4096.0 ) * 32.0;

    const float pq_C = 10000.0;
    // Note that this does NOT handle any of the signal range
    // considerations from 2084 - this returns full range (0 - 1)
    float L = C / pq_C;
    float Lm = pow( L, pq_m1 );
    float N = ( pq_c1 + pq_c2 * Lm ) / ( 1.0 + pq_c3 * Lm );
    N = pow( N, pq_m2 );
    return N;
}

float3 Y_2_ST2084_f3( float3 input)
{
    // converts from linear cd/m^2 to PQ code values
    float3 result;
    result[0] = Y_2_ST2084( input[0]);
    result[1] = Y_2_ST2084( input[1]);
    result[2] = Y_2_ST2084( input[2]);

    return result;
}

float ST2084_2_Y(float N)
{
    // Constants from SMPTE ST 2084-2014
    const float pq_m1 = 0.1593017578125; // ( 2610.0 / 4096.0 ) / 4.0;
    const float pq_m2 = 78.84375; // ( 2523.0 / 4096.0 ) * 128.0;
    const float pq_c1 = 0.8359375; // 3424.0 / 4096.0 or pq_c3 - pq_c2 + 1.0;
    const float pq_c2 = 18.8515625; // ( 2413.0 / 4096.0 ) * 32.0;
    const float pq_c3 = 18.6875; // ( 2392.0 / 4096.0 ) * 32.0;

    const float pq_C = 10000.0;
    
    float Np = pow( N, 1.0 / pq_m2 );
    float L = Np - pq_c1;
    if(L < 0.0)
    {
        L = 0.0;
    }
    L = L / ( pq_c2 - pq_c3 * Np );
    L = pow( L, 1.0 / pq_m1 );
    return L * pq_C;
}

float3 ST2084_2_Y_f3(float3 val)
{
    float3 res;
    res[0] = ST2084_2_Y(val[0]);
    res[1] = ST2084_2_Y(val[1]);
    res[2] = ST2084_2_Y(val[2]);
    return res;
}

float3 OutputTransform(float3 input,
                     float Y_MIN,
                     float Y_MID,
                     float Y_MAX,
                     Chromaticities DISPLAY_PRI,
                     Chromaticities LIMITING_PRI)
{
    float3x3 XYZ_2_DISPLAY_PRI_MAT = XYZtoRGB( DISPLAY_PRI);
    TsParams PARAMS_DEFAULT = init_TsParams( Y_MIN, Y_MAX);
    float expShift = log2(inv_ssts(Y_MID, PARAMS_DEFAULT))-log2(0.18);
    TsParams PARAMS = init_TsParams( Y_MIN, Y_MAX, expShift);
    
    // RRT sweeteners
    float3 rgbPre = rrt_sweeteners(input);
    float3 rgbPost = ssts_f3( rgbPre, PARAMS);
    
    // At this point data encoded AP1, scaled absolute luminance (cd/m^2)

    /*  Scale absolute luminance to linear code value  */
    float3 linearCV = Y_2_linCV_f3( rgbPost, Y_MAX, Y_MIN);
    // Rendering primaries to XYZ
    float3 XYZ = linearCV * AP1_2_XYZ_MAT;
    XYZ = XYZ * D60_2_D65_CAT;
    // CIE XYZ to display encoding primaries
    linearCV = XYZ_2_DISPLAY_PRI_MAT * XYZ;
    linearCV = clamp(linearCV, 0., 65535.);
    float3 outputCV = Y_2_ST2084_f3( clamp( linCV_2_Y_f3(linearCV, Y_MAX, 0.0), 0.0, 65535.) );
    
    return outputCV;
}

float3 RRT(float3 input)
{
    float3 rgbPre = rrt_sweeteners(input);
    return rgbPre;
}

float3 ComputeRRTLutValue(float3 lin)
{
    return RRT(lin);
}

float3 ColorExpanding(float3 ap0)
{
    
    return ap0;
}

float3 ACESHDR10Standard(float3 lin, float mid)
{
    const float Y_MIN = 0.0001f;                     // black luminance (cd/m^2)
    const float Y_MAX = 1000.0f;                     // peak white luminance (cd/m^2)
    
    const Chromaticities DISPLAY_PRI = REC2020_PRI; // encoding primaries (device setup)
    const Chromaticities LIMITING_PRI = REC2020_PRI;// limiting primaries

    // Color Expanding


    float3 result = OutputTransform(1.5f * lin, Y_MIN, mid, Y_MAX, DISPLAY_PRI, LIMITING_PRI);

    return result;
}

float3 ComputeLutValue(float3 pq, float mid)
{
    float3 lin = ST2084_2_Y_f3(pq);
    
    return ACESHDR10Standard(lin, mid);
}

constant float3x3 ACESInputMatrix = //AP0_2_AP1_MAT * CONST_RRT_SAT_MAT;
{
    {0.865855216, 0.262915015, -0.128670529},
    {0.0409506075,1.03673851,  -0.0777018442},
    {0.0309127532, 0.11267148, 0.856339156}
};

constant float3x3 ACESOutputMatrix = //AP1_2_XYZ_MAT * D60_2_D65_CAT * XYZ_2_Rec2020_MAT;
{
    {1.02579927, -0.0200525038, -0.00577139854},
    {-0.0022350112, 1.00458252, -0.00235230662},
    {-0.00501400419, -0.0252933856, 1.03044021}
};

float ssts_fitting(float x)
{
    static const float contrast = 3.18016763;
    static const float shoulder = 1.15319486;
    static const float b = 0.30002518;
    static const float c = 132.42972136;

    float logx = log10(x);

    float z = pow(logx + 6.0f, contrast);
    float logy = z / (pow(z, shoulder) * b + c);

    float y = pow(10.0, logy * 7.0f - 4.0f);
    y = clamp(y, 0.001f, 65535.0f);
    return y;
}

float3 ssts_fitting_f3(float3 val)
{
    float3 res;
    res[0] = ssts_fitting(val[0]);
    res[1] = ssts_fitting(val[1]);
    res[2] = ssts_fitting(val[2]);
    return res; 
}

float3 inv_ssts_f3(float3 val, TsParams C)
{
    float3 res;
    res[0] = inv_ssts(val[0], C);
    res[1] = inv_ssts(val[1], C);
    res[2] = inv_ssts(val[2], C);
    return res; 
}

float3 rrt_sweeteners_fitting(float3 aces)
{    
    const float RRT_SAT_FACTOR = 0.96;
    float3x3 RRT_SAT_MAT = calc_sat_adjust_matrix( RRT_SAT_FACTOR, AP1_RGB2Y);
    
    // --- ACES to RGB rendering space --- //
    aces = clamp( aces, 0., 65535.);
    float3 rgbPre = aces * AP0_2_AP1_MAT;
    rgbPre = clamp( rgbPre, 0., 65504.);
    
    // --- Global desaturation --- //
    rgbPre = RRT_SAT_MAT * rgbPre;
    return rgbPre;
}


float3 rrt_sweeteners_fitting_backward(float3 rgbPre)
{    
    const float RRT_SAT_FACTOR = 0.96;
    float3x3 RRT_SAT_MAT = calc_sat_adjust_matrix( RRT_SAT_FACTOR, AP1_RGB2Y);
    
    rgbPre = inverse(RRT_SAT_MAT) * rgbPre;
    float3 aces = rgbPre * inverse(AP0_2_AP1_MAT);

    return aces;
}

float3 ForwardWrappingPrecise(float3 lin)
{
    const float Y_MIN = 0.0001;                     // black luminance (cd/m^2)
    const float Y_MAX = 1000.0;                     // peak white luminance (cd/m^2)
    const float mid   = 15.0;
    float3x3 XYZ_2_DISPLAY_PRI_MAT = XYZtoRGB(REC2020_PRI);
    TsParams PARAMS_DEFAULT = init_TsParams( Y_MIN, Y_MAX);
    float expShift = log2(inv_ssts(mid, PARAMS_DEFAULT))-log2(0.18);
    TsParams PARAMS = init_TsParams( Y_MIN, Y_MAX, expShift);

    float3 rgbPre = rrt_sweeteners_fitting(lin);
    rgbPre = clamp(rgbPre, 0.00001f, 65535.0f);
    float3 rgbPost = ssts_f3(rgbPre, PARAMS);

    rgbPost = rgbPost * AP1_2_XYZ_MAT;
    rgbPost = rgbPost * D60_2_D65_CAT;

    float3 linearCV = XYZtoRGB(REC2020_PRI) * rgbPost;
    linearCV = clamp(linearCV, 0.00000001f, 65535.0f);
    float3 outputCV = Y_2_ST2084_f3(linearCV);
    return outputCV;
}

float3 BackwardWrappingPrecise(float3 outputCV)
{
    const float Y_MIN = 0.0001;                     // black luminance (cd/m^2)
    const float Y_MAX = 1000.0;                     // peak white luminance (cd/m^2)
    const float mid   = 15.0;
    float3x3 XYZ_2_DISPLAY_PRI_MAT = XYZtoRGB(REC2020_PRI);
    TsParams PARAMS_DEFAULT = init_TsParams( Y_MIN, Y_MAX);
    float expShift = log2(inv_ssts(mid, PARAMS_DEFAULT))-log2(0.18);
    TsParams PARAMS = init_TsParams( Y_MIN, Y_MAX, expShift);

    float3 linearCV = ST2084_2_Y_f3(outputCV);
    float3 rgbPost = RGBtoXYZ(REC2020_PRI) * linearCV;
    rgbPost = rgbPost * inverse(D60_2_D65_CAT);
    rgbPost = rgbPost * inverse(AP1_2_XYZ_MAT);

    float3 rgbPre = inv_ssts_f3(rgbPost, PARAMS);
    float3 lin = rrt_sweeteners_fitting_backward(rgbPre);

    return lin;
}

float3 ComputeLutValuePrecise(float3 pq)
{
    float3 result;
    float3 lin = ST2084_2_Y_f3(pq);
    
    const float Y_MIN = 0.0001;                     // black luminance (cd/m^2)
    const float Y_MAX = 1000.0;                     // peak white luminance (cd/m^2)
    const float mid   = 15.0;
    
    const Chromaticities DISPLAY_PRI = REC2020_PRI; // encoding primaries (device setup)
    const Chromaticities LIMITING_PRI = REC2020_PRI;// limiting primaries
    
    result = OutputTransform(lin, Y_MIN, mid, Y_MAX, DISPLAY_PRI, LIMITING_PRI);
    
    return result;
}


#endif /* OutputTransform_h */
