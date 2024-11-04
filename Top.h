#ifndef TOP_H
#define TOP_H
#include <iostream>
#include "FixedPoint.h"    // fixed_t
#include <hls_math.h>      // HLS math functions
#include <hls_stream.h>    // HLS stream definitions
#include <ap_fixed.h>      // Fixed-point types

// Define constants
#define In_2D_Array_Index 102
#define Out_1D_Array_Index 306
#define Angle_Index 30
#define Max_Rot_Index 1000
#define Max_U_V_Index 1000
#define Max_UV_Cross 1000
#define Max_Cross_Mag 500
#define Max_Cross_Norm 500
#define Max_Area_Index 500
#define P_Score_Index 500
#define Pos_Score_Ind 500
#define Rot_Points_Index 10000
#define Element_Index 1300

extern "C" {

    void Top(fixed_t Ordered_X[], fixed_t Ordered_Y[], fixed_t Ordered_Z[],
             int Num_Rows, int Angle_Steps, fixed_t *Best_Score, fixed_t *Best_X, fixed_t *Best_Y);
}

#endif // TOP_H
