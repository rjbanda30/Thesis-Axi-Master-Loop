// Top_sw.cpp
#include "Top_sw.h"
#include "FixedPoint.h"
#include <cmath>
#include <cstring>  // For memset
#include <iostream>

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

// Software reference implementation of the Top function
void Top_sw(fixed_t Ordered_X[], fixed_t Ordered_Y[], fixed_t Ordered_Z[],
            int Num_Rows, int Angle_Steps, fixed_t *Best_Score, fixed_t *Best_X, fixed_t *Best_Y) {
#ifndef __SYNTHESIS__
    std::cout << "Entering TopSW Function" << std::endl<< std::endl<< std::endl<< std::endl<< std::endl<< std::endl<< std::endl;
#endif

    // Initialize local arrays
    fixed_t Out_1D_Array_X[Out_1D_Array_Index];
    fixed_t Out_1D_Array_Y[Out_1D_Array_Index];
    fixed_t Out_1D_Array_Z[Out_1D_Array_Index];

    // Zero initialize the output arrays
    memset(Out_1D_Array_X, 0, sizeof(Out_1D_Array_X));
    memset(Out_1D_Array_Y, 0, sizeof(Out_1D_Array_Y));
    memset(Out_1D_Array_Z, 0, sizeof(Out_1D_Array_Z));

    // Copy input arrays to local arrays
    for (int i = 0; i < Num_Rows; i++) {
        if (i < Out_1D_Array_Index) {
            Out_1D_Array_X[i] = Ordered_X[i];
            Out_1D_Array_Y[i] = Ordered_Y[i];
            Out_1D_Array_Z[i] = Ordered_Z[i];
        }
    }

    // Compute Rotation Angles
    fixed_t X_Rotation_Angle[Angle_Index];
    fixed_t Y_Rotation_Angle[Angle_Index];

    double Two_Pi = 6.283185307179;
    double Factor = Two_Pi / Angle_Steps;
    double Neg_Two_Pi = -Two_Pi;

    for (int i = 0; i < Angle_Index; i++) {
        double Factor_Mul = Factor * i;
        double Angle_0 = Neg_Two_Pi + Factor_Mul;
        X_Rotation_Angle[i] = double_to_fixed(Angle_0);
        Y_Rotation_Angle[i] = double_to_fixed(Angle_0);
    }

    // Compute Rotation Matrices
    fixed_t Rot_Matrix_X[Max_Rot_Index];
    fixed_t Rot_Matrix_Y[Max_Rot_Index];
    fixed_t Rot_Matrix_Z[Max_Rot_Index];

    memset(Rot_Matrix_X, 0, sizeof(Rot_Matrix_X));
    memset(Rot_Matrix_Y, 0, sizeof(Rot_Matrix_Y));
    memset(Rot_Matrix_Z, 0, sizeof(Rot_Matrix_Z));

    fixed_t Element_X[Element_Index];
    fixed_t Element_Y[Element_Index];
    memset(Element_X, 0, sizeof(Element_X));
    memset(Element_Y, 0, sizeof(Element_Y));

    int X_Y_Index = 0;

    for (int idx = 0; idx < Angle_Steps * Angle_Steps && idx < Element_Index; idx++) {
        int j = idx / Angle_Steps;
        int i = idx % Angle_Steps;

        Element_X[X_Y_Index] = X_Rotation_Angle[j];
        Element_Y[X_Y_Index] = Y_Rotation_Angle[i];
        X_Y_Index++;

        double cos_X = std::cos(fixed_to_double(X_Rotation_Angle[j]));
        double cos_Y = std::cos(fixed_to_double(Y_Rotation_Angle[i]));
        double sin_X = std::sin(fixed_to_double(X_Rotation_Angle[j]));
        double sin_Y = std::sin(fixed_to_double(Y_Rotation_Angle[i]));

        int d = 3 * idx;
        if (d + 2 < Max_Rot_Index) {
            Rot_Matrix_X[d] = double_to_fixed(cos_Y);
            Rot_Matrix_Y[d] = double_to_fixed(sin_X * sin_Y);
            Rot_Matrix_Z[d] = double_to_fixed(cos_X * sin_Y);

            Rot_Matrix_X[d + 1] = double_to_fixed(0.0); // Fixed-point 0
            Rot_Matrix_Y[d + 1] = double_to_fixed(cos_X);
            Rot_Matrix_Z[d + 1] = double_to_fixed(-sin_X);

            Rot_Matrix_X[d + 2] = double_to_fixed(-sin_Y);
            Rot_Matrix_Y[d + 2] = double_to_fixed(cos_Y * sin_X);
            Rot_Matrix_Z[d + 2] = double_to_fixed(cos_X * cos_Y);
        }
    }

    // Initialize variables for scoring
    *Best_Score = double_to_fixed(10000.0); // Initialize to a high value
    double Temp_Best_Score = 10000.0;
    int Best_Score_Element_Out = 10000;
    // Define Num_Faces and Score_Index
    int Num_Faces = Num_Rows / 3;
    int Score_Index = 0;

    const int End_Rot_Loop = (Angle_Steps + 1) * (Angle_Steps + 1);

    // Rotated Points Arrays
    fixed_t Rotated_X_Points_0[Rot_Points_Index] = {0};
    fixed_t Rotated_Y_Points_0[Rot_Points_Index] = {0};
    fixed_t Rotated_Z_Points_0[Rot_Points_Index] = {0};

    fixed_t Rotated_X_Points_1[Rot_Points_Index] = {0};
    fixed_t Rotated_Y_Points_1[Rot_Points_Index] = {0};
    fixed_t Rotated_Z_Points_1[Rot_Points_Index] = {0};

    fixed_t Rotated_X_Points_2[Rot_Points_Index] = {0};
    fixed_t Rotated_Y_Points_2[Rot_Points_Index] = {0};
    fixed_t Rotated_Z_Points_2[Rot_Points_Index] = {0};

    // Compute Scores
    for (int Rot_M_Indx = 0; Rot_M_Indx < End_Rot_Loop && Rot_M_Indx < Max_Rot_Index / 3; Rot_M_Indx++) {
        // Produces Rotated Points
        for (int i = 0; i < (Num_Rows / 3) && i < Rot_Points_Index; i++) {
            // Compute Rotated Points for Set 0
            Rotated_X_Points_0[i] = fixed_mul(Rot_Matrix_X[Rot_M_Indx * 3], Out_1D_Array_X[3 * i]) +
                                    fixed_mul(Rot_Matrix_Y[Rot_M_Indx * 3], Out_1D_Array_Y[3 * i]) +
                                    fixed_mul(Rot_Matrix_Z[Rot_M_Indx * 3], Out_1D_Array_Z[3 * i]);

            Rotated_Y_Points_0[i] = fixed_mul(Rot_Matrix_X[Rot_M_Indx * 3 + 1], Out_1D_Array_X[3 * i]) +
                                    fixed_mul(Rot_Matrix_Y[Rot_M_Indx * 3 + 1], Out_1D_Array_Y[3 * i]) +
                                    fixed_mul(Rot_Matrix_Z[Rot_M_Indx * 3 + 1], Out_1D_Array_Z[3 * i]);

            Rotated_Z_Points_0[i] = fixed_mul(Rot_Matrix_X[Rot_M_Indx * 3 + 2], Out_1D_Array_X[3 * i]) +
                                    fixed_mul(Rot_Matrix_Y[Rot_M_Indx * 3 + 2], Out_1D_Array_Y[3 * i]) +
                                    fixed_mul(Rot_Matrix_Z[Rot_M_Indx * 3 + 2], Out_1D_Array_Z[3 * i]);

            // Compute Rotated Points for Set 1
            Rotated_X_Points_1[i] = fixed_mul(Rot_Matrix_X[Rot_M_Indx * 3], Out_1D_Array_X[(3 * i) + 1]) +
                                    fixed_mul(Rot_Matrix_Y[Rot_M_Indx * 3], Out_1D_Array_Y[(3 * i) + 1]) +
                                    fixed_mul(Rot_Matrix_Z[Rot_M_Indx * 3], Out_1D_Array_Z[(3 * i) + 1]);

            Rotated_Y_Points_1[i] = fixed_mul(Rot_Matrix_X[Rot_M_Indx * 3 + 1], Out_1D_Array_X[(3 * i) + 1]) +
                                    fixed_mul(Rot_Matrix_Y[Rot_M_Indx * 3 + 1], Out_1D_Array_Y[(3 * i) + 1]) +
                                    fixed_mul(Rot_Matrix_Z[Rot_M_Indx * 3 + 1], Out_1D_Array_Z[(3 * i) + 1]);

            Rotated_Z_Points_1[i] = fixed_mul(Rot_Matrix_X[Rot_M_Indx * 3 + 2], Out_1D_Array_X[(3 * i) + 1]) +
                                    fixed_mul(Rot_Matrix_Y[Rot_M_Indx * 3 + 2], Out_1D_Array_Y[(3 * i) + 1]) +
                                    fixed_mul(Rot_Matrix_Z[Rot_M_Indx * 3 + 2], Out_1D_Array_Z[(3 * i) + 1]);

            // Compute Rotated Points for Set 2
            Rotated_X_Points_2[i] = fixed_mul(Rot_Matrix_X[Rot_M_Indx * 3], Out_1D_Array_X[(3 * i) + 2]) +
                                    fixed_mul(Rot_Matrix_Y[Rot_M_Indx * 3], Out_1D_Array_Y[(3 * i) + 2]) +
                                    fixed_mul(Rot_Matrix_Z[Rot_M_Indx * 3], Out_1D_Array_Z[(3 * i) + 2]);

            Rotated_Y_Points_2[i] = fixed_mul(Rot_Matrix_X[Rot_M_Indx * 3 + 1], Out_1D_Array_X[(3 * i) + 2]) +
                                    fixed_mul(Rot_Matrix_Y[Rot_M_Indx * 3 + 1], Out_1D_Array_Y[(3 * i) + 2]) +
                                    fixed_mul(Rot_Matrix_Z[Rot_M_Indx * 3 + 1], Out_1D_Array_Z[(3 * i) + 2]);

            Rotated_Z_Points_2[i] = fixed_mul(Rot_Matrix_X[Rot_M_Indx * 3 + 2], Out_1D_Array_X[(3 * i) + 2]) +
                                    fixed_mul(Rot_Matrix_Y[Rot_M_Indx * 3 + 2], Out_1D_Array_Y[(3 * i) + 2]) +
                                    fixed_mul(Rot_Matrix_Z[Rot_M_Indx * 3 + 2], Out_1D_Array_Z[(3 * i) + 2]);
        }

        // Compute U and V vectors
        double Comp_U_0[Max_U_V_Index];
        double Comp_U_1[Max_U_V_Index];
        double Comp_U_2[Max_U_V_Index];

        double Comp_V_0[Max_U_V_Index];
        double Comp_V_1[Max_U_V_Index];
        double Comp_V_2[Max_U_V_Index];
        const int Count_6 = Num_Faces;

        for (int i = 0; i < Count_6 && i < Max_U_V_Index; i++) {
            // U = Point1 - Point0
            Comp_U_0[i] = fixed_to_double(Rotated_X_Points_1[i] - Rotated_X_Points_0[i]);
            Comp_U_1[i] = fixed_to_double(Rotated_Y_Points_1[i] - Rotated_Y_Points_0[i]);
            Comp_U_2[i] = fixed_to_double(Rotated_Z_Points_1[i] - Rotated_Z_Points_0[i]);
        }

        for (int j = 0; j < Count_6 && j < Max_U_V_Index; j++) {
            // V = Point2 - Point0
            Comp_V_0[j] = fixed_to_double(Rotated_X_Points_2[j] - Rotated_X_Points_0[j]);
            Comp_V_1[j] = fixed_to_double(Rotated_Y_Points_2[j] - Rotated_Y_Points_0[j]);
            Comp_V_2[j] = fixed_to_double(Rotated_Z_Points_2[j] - Rotated_Z_Points_0[j]);
        }

        // Compute U Cross V
        double U_Cross_V_X[Max_UV_Cross];
        double U_Cross_V_Y[Max_UV_Cross];
        double U_Cross_V_Z[Max_UV_Cross];
        double Temp_X = 0;
        double Temp_Y = 0;
        double Temp_Z = 0;

        const int Count_7 = Num_Faces;

        for (int i = 0; i < Count_7 && i < Max_UV_Cross; i++) {
            // Computes cross product
            Temp_X = (Comp_U_1[i] * Comp_V_2[i]) - (Comp_U_2[i] * Comp_V_1[i]);
            Temp_Y = (Comp_U_0[i] * Comp_V_2[i]) - (Comp_U_2[i] * Comp_V_0[i]);
            Temp_Z = (Comp_U_0[i] * Comp_V_1[i]) - (Comp_U_1[i] * Comp_V_0[i]);

            U_Cross_V_X[i] = Temp_X;
            U_Cross_V_Y[i] = Temp_Y;
            U_Cross_V_Z[i] = Temp_Z;
        }

        // Compute Cross Magnitude
        double Cross_Magnitude_Array_0[Max_Cross_Mag];
        const int Count_8 = Num_Faces;

        for (int i = 0; i < Count_8 && i < Max_Cross_Mag; i++) {
            Cross_Magnitude_Array_0[i] = std::sqrt(std::pow(U_Cross_V_X[i], 2) +
                                                  std::pow(U_Cross_V_Y[i], 2) +
                                                  std::pow(U_Cross_V_Z[i], 2));
        }

        // Compute U Cross V Normal
        double U_V_Cross_Norm_Array_0[Max_Cross_Norm];
        double U_V_Cross_Norm_Array_1[Max_Cross_Norm];
        double U_V_Cross_Norm_Array_2[Max_Cross_Norm];

        const int Count_9 = Num_Faces;

        for (int i = 0; i < Count_9 && i < Max_Cross_Norm; i++) {
            if (Cross_Magnitude_Array_0[i] != 0) {
                U_V_Cross_Norm_Array_0[i] = U_Cross_V_X[i] / Cross_Magnitude_Array_0[i];
                U_V_Cross_Norm_Array_1[i] = U_Cross_V_Y[i] / Cross_Magnitude_Array_0[i];
                U_V_Cross_Norm_Array_2[i] = U_Cross_V_Z[i] / Cross_Magnitude_Array_0[i];
            }
            else {

                U_V_Cross_Norm_Array_0[i] = 0;
                U_V_Cross_Norm_Array_1[i] = 0;
                U_V_Cross_Norm_Array_2[i] = 0;

                // std::cout << "Zero magnitude at index " << i << ". Points are colinear or duplicate." << std::endl;
            }
        }

        // Compute Cross Area
        double Cross_Area_Out_0[Max_Area_Index];
        const int Count_10 = Num_Faces;

        for (int i = 0; i < Count_10 && i < Max_Area_Index; i++) {
            Cross_Area_Out_0[i] = 0.5 * Cross_Magnitude_Array_0[i];
        }

        // Compute Point Scores
        double Point_Score_Out_0[P_Score_Index] = {0};
        int Unit_X[3] = {1, 0, 0};
        int Unit_Y[3] = {0, 1, 0};
        int Unit_Z[3] = {0, 0, 1};
        const int Count_11 = Num_Faces;

        for(int i = 0; i < Count_11 && i < P_Score_Index; i++) {
            double dot_X = U_V_Cross_Norm_Array_0[i] * Unit_X[0] +
                           U_V_Cross_Norm_Array_1[i] * Unit_X[1] +
                           U_V_Cross_Norm_Array_2[i] * Unit_X[2];
            double dot_Y = U_V_Cross_Norm_Array_0[i] * Unit_Y[0] +
                           U_V_Cross_Norm_Array_1[i] * Unit_Y[1] +
                           U_V_Cross_Norm_Array_2[i] * Unit_Y[2];
            double dot_Z = U_V_Cross_Norm_Array_0[i] * Unit_Z[0] +
                           U_V_Cross_Norm_Array_1[i] * Unit_Z[1] +
                           U_V_Cross_Norm_Array_2[i] * Unit_Z[2];

            double abs_sum = std::abs(dot_X) + std::abs(dot_Y) + std::abs(dot_Z);
            Point_Score_Out_0[i] = Cross_Area_Out_0[i] * abs_sum;
            std::cout << "Point_Score_Out_0: " << Point_Score_Out_0[i] << std::endl;
        }

        // Compute Total Scores and Find Best Score
        double Position_Score_Out[Pos_Score_Ind] = {0};
        const int Count_12 = Count_6;

        double Temp_Score = 0.0;
        double Temp_Avg_Score = 0.0;

        for(int rot_idx = 0; rot_idx < End_Rot_Loop && rot_idx < Pos_Score_Ind; rot_idx++) {
            for(int i = 0; i < Count_6 && i < P_Score_Index; i++) {
                Temp_Score += Point_Score_Out_0[i];
            }

            Temp_Avg_Score = (Num_Faces > 0) ? (Temp_Score / Num_Faces) : 0.0;
            Position_Score_Out[Score_Index] = Temp_Avg_Score;

            std::cout << "Score_Index: " << Score_Index << std::endl;
            std::cout << "Position_Score_Out: " << Position_Score_Out[Score_Index]  << std::endl;

            if(Position_Score_Out[Score_Index] < Temp_Best_Score) {
                Best_Score_Element_Out = Score_Index;
                *Best_Score = double_to_fixed(Position_Score_Out[Score_Index]);
                Temp_Best_Score = Position_Score_Out[Score_Index];
                // Debugging Best Score
                std::cout << "SW_Best Score: " << fixed_to_double(*Best_Score) << std::endl;
            }

            Score_Index++;
            Temp_Score = 0.0;
            Temp_Avg_Score = 0.0;
        }


    }
    // Assign Best_X and Best_Y based on Best_Score_Element_Out
            if(Best_Score_Element_Out < X_Y_Index) {
                *Best_X = Element_X[Best_Score_Element_Out];
                *Best_Y = Element_Y[Best_Score_Element_Out];
            }
}
