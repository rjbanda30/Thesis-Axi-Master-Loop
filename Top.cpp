// Top.cpp
#include "Top.h"
#include <cfloat>
#include <hls_math.h>      // For HLS math functions
#include <ap_fixed.h>      // For fixed-point data types
#include <iostream>        // For debugging



extern "C" {
void Top(fixed_t Ordered_X[], fixed_t Ordered_Y[], fixed_t Ordered_Z[],
         int Num_Rows, int Angle_Steps, fixed_t *Best_Score, fixed_t *Best_X, fixed_t *Best_Y) {
	#pragma HLS ARRAY_PARTITION dim=1 type=complete variable=Ordered_Z
	#pragma HLS ARRAY_PARTITION dim=1 type=complete variable=Ordered_Y
	#pragma HLS ARRAY_PARTITION dim=1 type=complete variable=Ordered_X





    #pragma HLS INTERFACE mode=m_axi depth=1 port=Ordered_X
    #pragma HLS INTERFACE mode=m_axi depth=1 port=Ordered_Y
    #pragma HLS INTERFACE mode=m_axi depth=1 port=Ordered_Z

    #pragma HLS INTERFACE mode=s_axilite bundle=control port=Num_Rows
    #pragma HLS INTERFACE mode=s_axilite bundle=control port=Angle_Steps
    #pragma HLS INTERFACE mode=s_axilite bundle=control port=Best_Score
    #pragma HLS INTERFACE mode=s_axilite bundle=control port=Best_X
    #pragma HLS INTERFACE mode=s_axilite bundle=control port=Best_Y
    #pragma HLS INTERFACE s_axilite port=return bundle=control

    //***************************************
    // Converts 2D array to 1D
    //***************************************

    // Temp is used to reduce array partition factor
    fixed_t Out_1D_Array_X[Out_1D_Array_Index] = {0};
    fixed_t Out_1D_Array_Y[Out_1D_Array_Index] = {0};
    fixed_t Out_1D_Array_Z[Out_1D_Array_Index] = {0};
    fixed_t Temp_1D_X = 0;
    fixed_t Temp_1D_Y = 0;
    fixed_t Temp_1D_Z = 0;
    const int Count_0 = Num_Rows;

    //***************************************
    // Converts 2D array to 1D
    //***************************************

    Input_Loop:
    for (int i = 0; i < Count_0; i++) {
        #pragma HLS PIPELINE II=4

        // Assists in reducing the number of times In_2D_Array needs to be accessed
        Temp_1D_X = Ordered_X[i];
        Temp_1D_Y = Ordered_Y[i];
        Temp_1D_Z = Ordered_Z[i];

        // Write the row to the output stream
        Out_1D_Array_X[i] = Temp_1D_X;
        Out_1D_Array_Y[i] = Temp_1D_Y;
        Out_1D_Array_Z[i] = Temp_1D_Z;
    }

    //***************************************
    // Compute_Rotation_Angle
    //***************************************

    const int Count_1 = Angle_Steps;
    const fixed_t Two_Pi = fixed_t(6.283185307179); // 2 * PI as fixed_t
    const fixed_t Factor = Two_Pi / Angle_Steps;
    const fixed_t Neg_Two_Pi = fixed_t(-Two_Pi);
    fixed_t Angle_0 = 0;
    fixed_t Angle_1 = 0;

    fixed_t Factor_Mul = 0;
    fixed_t X_Rotation_Angle[Angle_Index];
    fixed_t Y_Rotation_Angle[Angle_Index];

    //***************************************
    // Compute_Rotation_Angle
    //***************************************

    Compute_Angles_Array_Loop:
    for(int i = 0; i <= Count_1; i++){
        #pragma HLS LOOP_TRIPCOUNT min=1 max=50
        #pragma HLS PIPELINE II=1
        #pragma HLS UNROLL factor=3

        Factor_Mul = Factor * i;
        Angle_0 = Neg_Two_Pi + Factor_Mul;
        X_Rotation_Angle[i] = Angle_0;
        Y_Rotation_Angle[i] = Angle_0;
    }

    //***************************************
    // Compute Rotation Matrices
    //***************************************

    fixed_t Rot_Matrix_X[Max_Rot_Index];
    fixed_t Rot_Matrix_Y[Max_Rot_Index];
    fixed_t Rot_Matrix_Z[Max_Rot_Index];

    const int Count_2 = Angle_Steps;
    const int Count_3 = Angle_Steps;
    int X_Y_Index=0;
    fixed_t Element_X[Element_Index] = {0};
    fixed_t Element_Y[Element_Index] = {0};

    int c = 0;
    int d = 0;

    Rot_Matrix_Man_Flattened_Loop:
    for (int idx = 0; idx < Count_2 * Count_3; idx++) {
        #pragma HLS PIPELINE II=2
        #pragma HLS LOOP_TRIPCOUNT min=1 max=2500

        int j = idx / Count_3;  // Original outer loop index
        int i = idx % Count_3;  // Original inner loop index

        // Perform local computations
        Element_X[X_Y_Index] = X_Rotation_Angle[j];
        Element_Y[X_Y_Index] = Y_Rotation_Angle[i];

        X_Y_Index++;

        fixed_t cos_X = hls::cos(X_Rotation_Angle[j]);
        fixed_t cos_Y = hls::cos(Y_Rotation_Angle[i]);
        fixed_t sin_X = hls::sin(X_Rotation_Angle[j]);
        fixed_t sin_Y = hls::sin(Y_Rotation_Angle[i]);
        fixed_t neg_sin_Y = fixed_t(-sin_Y);
        fixed_t neg_sin_X = fixed_t(-sin_X);

        fixed_t sin_X_sin_Y = sin_X * sin_Y;
        fixed_t cos_X_sin_Y = cos_X * sin_Y;
        fixed_t cos_Y_sin_X = cos_Y * sin_X;
        fixed_t cos_X_cos_Y = cos_X * cos_Y;

        // Assign values to Rot_Matrix
        int d_idx = 3 * idx;  // Adjusted calculation for `d` based on the flattened loop index

        // Ensure indices do not exceed array bounds
        if (d_idx + 2 >= Max_Rot_Index) {

            #ifndef __SYNTHESIS__
            std::cout << "Rotation matrix index out of bounds: " << d_idx + 2 << std::endl;
            #endif
            continue; // Skip if out of bounds
        }

        Rot_Matrix_X[d_idx] = cos_Y;
        Rot_Matrix_Y[d_idx] = sin_X_sin_Y;
        Rot_Matrix_Z[d_idx] = cos_X_sin_Y;

        Rot_Matrix_X[d_idx + 1] = 0;
        Rot_Matrix_Y[d_idx + 1] = cos_X;
        Rot_Matrix_Z[d_idx + 1] = neg_sin_X;

        Rot_Matrix_X[d_idx + 2] = neg_sin_Y;
        Rot_Matrix_Y[d_idx + 2] = cos_Y_sin_X;
        Rot_Matrix_Z[d_idx + 2] = cos_X_cos_Y;
    }

    //***************************************
    // Initialize Best_Score Before Score_Loop
    //***************************************

    *Best_Score = 10000;
    fixed_t Temp_Best_Score = 10000; // Initialize to large value
    int Best_Score_Element_Out = 10000; // Initialize to large value

    //***************************************
    // For Loop to Control Flow of Scores
    //***************************************

    const int End_Rot_Loop = (Angle_Steps + 1) * (Angle_Steps + 1);
    const int Best_Elem_Loop = (Angle_Steps + 1) * (Angle_Steps + 1);

    // Initialize rotated points arrays
    fixed_t Rotated_X_Points_0[Rot_Points_Index] = {0};
    fixed_t Rotated_Y_Points_0[Rot_Points_Index] = {0};
    fixed_t Rotated_Z_Points_0[Rot_Points_Index] = {0};

    fixed_t Rotated_X_Points_1[Rot_Points_Index] = {0};
    fixed_t Rotated_Y_Points_1[Rot_Points_Index] = {0};
    fixed_t Rotated_Z_Points_1[Rot_Points_Index] = {0};

    fixed_t Rotated_X_Points_2[Rot_Points_Index] = {0};
    fixed_t Rotated_Y_Points_2[Rot_Points_Index] = {0};
    fixed_t Rotated_Z_Points_2[Rot_Points_Index] = {0};

    fixed_t Rotated_X_0 = 0;
    fixed_t Rotated_Y_0 = 0;
    fixed_t Rotated_Z_0 = 0;

    fixed_t Rotated_X_1 = 0;
    fixed_t Rotated_Y_1 = 0;
    fixed_t Rotated_Z_1 = 0;

    fixed_t Rotated_X_2 = 0;
    fixed_t Rotated_Y_2 = 0;
    fixed_t Rotated_Z_2 = 0;

    const int Num_Faces = Num_Rows / 3;
    const int In_1D_Limit = Num_Faces;
    int Score_Index = 0;

    Score_Loop:
    for(int Rot_M_Indx = 0; Rot_M_Indx < End_Rot_Loop; Rot_M_Indx++){
        #pragma HLS LOOP_FLATTEN
        #pragma HLS LOOP_TRIPCOUNT min=121 max=2601 avg=961

        //***************************************
        // Produces Rotated Points
        //***************************************

        Rotate_Points_Loop:
        for(int i = 0; i < In_1D_Limit; i++){
            #pragma HLS PIPELINE II=3

            Rotated_X_0 = Rot_Matrix_X[Rot_M_Indx] * Out_1D_Array_X[(3*i)]   + Rot_Matrix_Y[Rot_M_Indx] * Out_1D_Array_Y[(3*i)]   + Rot_Matrix_Z[Rot_M_Indx] * Out_1D_Array_Z[(3*i)];
            Rotated_Y_0 = Rot_Matrix_X[Rot_M_Indx + 1] * Out_1D_Array_X[3*i] + Rot_Matrix_Y[Rot_M_Indx + 1] * Out_1D_Array_Y[(3*i)] + Rot_Matrix_Z[Rot_M_Indx + 1] * Out_1D_Array_Z[(3*i)];
            Rotated_Z_0 = Rot_Matrix_X[Rot_M_Indx + 2] * Out_1D_Array_X[3*i] + Rot_Matrix_Y[Rot_M_Indx + 2] * Out_1D_Array_Y[(3*i)] + Rot_Matrix_Z[Rot_M_Indx + 2] * Out_1D_Array_Z[(3*i)];

            Rotated_X_1 = Rot_Matrix_X[Rot_M_Indx] * Out_1D_Array_X[(3*i)+1] + Rot_Matrix_Y[Rot_M_Indx] * Out_1D_Array_Y[(3*i)+1] + Rot_Matrix_Z[Rot_M_Indx] * Out_1D_Array_Z[(3*i)+1];
            Rotated_Y_1 = Rot_Matrix_X[Rot_M_Indx + 1] * Out_1D_Array_X[(3*i)+1] + Rot_Matrix_Y[Rot_M_Indx + 1] * Out_1D_Array_Y[(3*i)+1] + Rot_Matrix_Z[Rot_M_Indx + 1] * Out_1D_Array_Z[(3*i)+1];
            Rotated_Z_1 = Rot_Matrix_X[Rot_M_Indx + 2] * Out_1D_Array_X[(3*i)+1] + Rot_Matrix_Y[Rot_M_Indx + 2] * Out_1D_Array_Y[(3*i)+1] + Rot_Matrix_Z[Rot_M_Indx + 2] * Out_1D_Array_Z[(3*i)+1];

            Rotated_X_2 = Rot_Matrix_X[Rot_M_Indx] * Out_1D_Array_X[(3*i)+2] + Rot_Matrix_Y[Rot_M_Indx] * Out_1D_Array_Y[(3*i)+2] + Rot_Matrix_Z[Rot_M_Indx] * Out_1D_Array_Z[(3*i)+2];
            Rotated_Y_2 = Rot_Matrix_X[Rot_M_Indx + 1] * Out_1D_Array_X[(3*i)+2] + Rot_Matrix_Y[Rot_M_Indx + 1] * Out_1D_Array_Y[(3*i)+2] + Rot_Matrix_Z[Rot_M_Indx + 1] * Out_1D_Array_Z[(3*i)+2];
            Rotated_Z_2 = Rot_Matrix_X[Rot_M_Indx + 2] * Out_1D_Array_X[(3*i)+2] + Rot_Matrix_Y[Rot_M_Indx + 2] * Out_1D_Array_Y[(3*i)+2] + Rot_Matrix_Z[Rot_M_Indx + 2] * Out_1D_Array_Z[(3*i)+2];

            Rotated_X_Points_0[i] = Rotated_X_0;
            Rotated_Y_Points_0[i] = Rotated_Y_0;
            Rotated_Z_Points_0[i] = Rotated_Z_0;

            Rotated_X_Points_1[i] = Rotated_X_1;
            Rotated_Y_Points_1[i] = Rotated_Y_1;
            Rotated_Z_Points_1[i] = Rotated_Z_1;

            Rotated_X_Points_2[i] = Rotated_X_2;
            Rotated_Y_Points_2[i] = Rotated_Y_2;
            Rotated_Z_Points_2[i] = Rotated_Z_2;
        }

        //***************************************
        // Computes U and V that will be used in the cross product
        //***************************************

        fixed_t Comp_U_0[Max_U_V_Index];
        fixed_t Comp_U_1[Max_U_V_Index];
        fixed_t Comp_U_2[Max_U_V_Index];

        fixed_t Comp_V_0[Max_U_V_Index];
        fixed_t Comp_V_1[Max_U_V_Index];
        fixed_t Comp_V_2[Max_U_V_Index];
        const int Count_6 = Num_Faces;

        Compute_U_Loop:
        for(int i = 0;i<Count_6; i++){
            #pragma HLS PIPELINE II=1

            // This is U = Point 1 - Point 0
            Comp_U_0[i] = Rotated_X_Points_1[i] - Rotated_X_Points_0[i];
            Comp_U_1[i] = Rotated_Y_Points_1[i] - Rotated_Y_Points_0[i];
            Comp_U_2[i] = Rotated_Z_Points_1[i] - Rotated_Z_Points_0[i];
        }

        Compute_V_Loop:
        for(int j = 0;j<Count_6; j++){
            #pragma HLS PIPELINE II=1

            // This is V = Point 2 - Point 0
            Comp_V_0[j] = Rotated_X_Points_2[j] - Rotated_X_Points_0[j];
            Comp_V_1[j] = Rotated_Y_Points_2[j] - Rotated_Y_Points_0[j];
            Comp_V_2[j] = Rotated_Z_Points_2[j] - Rotated_Z_Points_0[j];
        }

        //***************************************
        // Computes U Cross V
        //***************************************

        fixed_t U_Cross_V_X[Max_UV_Cross];
        fixed_t U_Cross_V_Y[Max_UV_Cross];
        fixed_t U_Cross_V_Z[Max_UV_Cross];
        fixed_t Temp_X = 0;
        fixed_t Temp_Y = 0;
        fixed_t Temp_Z = 0;

        const int Count_7 = Num_Faces;

        U_Cross_V_Loop:
        for(int i = 0; i < Count_7 ; i++){
            #pragma HLS PIPELINE II=1

            // Computes cross product components
            Temp_X = (Comp_U_1[i] * Comp_V_2[i]) - (Comp_U_2[i] * Comp_V_1[i]);
            Temp_Y = (Comp_U_0[i] * Comp_V_2[i]) - (Comp_U_2[i] * Comp_V_0[i]);
            Temp_Z = (Comp_U_0[i] * Comp_V_1[i]) - (Comp_U_1[i] * Comp_V_0[i]);

            U_Cross_V_X[i] = Temp_X;
            U_Cross_V_Y[i] = Temp_Y;
            U_Cross_V_Z[i] = Temp_Z;
        }

        //***************************************
        // Computes U Cross V Magnitude
        //***************************************

        fixed_t Cross_Magnitude_Array_0[Max_Cross_Mag];

        const int Count_8 = Num_Faces;

        Cross_Magnitude_Loop:
        for(int i = 0;i<Count_8; i++)
        {
            #pragma HLS PIPELINE II=1


            Cross_Magnitude_Array_0[i] = hls::sqrt( (U_Cross_V_X[i] * U_Cross_V_X[i]) +
                                                  (U_Cross_V_Y[i] * U_Cross_V_Y[i]) +
                                                  (U_Cross_V_Z[i] * U_Cross_V_Z[i]) );
        }

        //***************************************
        // Computes U Cross V Normal
        //***************************************

        fixed_t U_V_Cross_Norm_Array_0[Max_Cross_Norm];
        fixed_t U_V_Cross_Norm_Array_1[Max_Cross_Norm];
        fixed_t U_V_Cross_Norm_Array_2[Max_Cross_Norm];

        int h = 0;
        int Count_9 = Num_Faces;

        U_V_Cross_Norm_Loop:

        for(int i = 0;i<Count_9; i++)
        {
            #pragma HLS PIPELINE II=1
             if (Cross_Magnitude_Array_0[i] != 0) {

                U_V_Cross_Norm_Array_0[i] = U_Cross_V_X[i] / Cross_Magnitude_Array_0[i];
                U_V_Cross_Norm_Array_1[i] = U_Cross_V_Y[i] / Cross_Magnitude_Array_0[i];
                U_V_Cross_Norm_Array_2[i] = U_Cross_V_Z[i] / Cross_Magnitude_Array_0[i];

             }
             else {

                 // Debuging Rotated Points
#ifndef __SYNTHESIS__

                 std::cout << "Zero magnitude at index " << i << ". Points are colinear or duplicate." << std::endl;
                 std::cout << "Points: (" << Rotated_X_Points_0[i] << ", " << Rotated_Y_Points_0[i] << ", " << Rotated_Z_Points_0[i] << "), "
                           << "(" << Rotated_X_Points_1[i] << ", " << Rotated_Y_Points_1[i] << ", " << Rotated_Z_Points_1[i] << "), "
                           << "(" << Rotated_X_Points_2[i] << ", " << Rotated_Y_Points_2[i] << ", " << Rotated_Z_Points_2[i] << ")" << std::endl;
#endif


                 U_V_Cross_Norm_Array_0[i] = 0;
                 U_V_Cross_Norm_Array_1[i] = 0;
                 U_V_Cross_Norm_Array_2[i] = 0;
             }
        }

        //***************************************
        // Computes U Cross V Area
        //***************************************

        fixed_t Cross_Area_Out_0[Max_Area_Index];

        // Initializing Local Variables
        const int Count_10 = Num_Faces;

        Cross_Area_Loop:

        for(int i = 0;i<Count_10; i++){
            #pragma HLS PIPELINE II=1


            Cross_Area_Out_0[i] = fixed_t(0.5) * Cross_Magnitude_Array_0[i];
        }

        //***************************************
        // Computes the scores of each face at current print position
        //***************************************

        fixed_t Point_Score_Out_0[P_Score_Index] = {0};

        int Unit_X[3] = {1, 0, 0};
        int Unit_Y[3] = {0, 1, 0};
        int Unit_Z[3] = {0, 0, 1};
        const int Count_11 = Num_Faces;

        Point_Score_Loop:

        for(int i = 0;i<Count_11; i++)
        {
            #pragma HLS PIPELINE II=1

            // Compute dot products with unit vectors
            fixed_t dot_X = (U_V_Cross_Norm_Array_0[i] * Unit_X[0]) +
                            (U_V_Cross_Norm_Array_1[i] * Unit_X[1]) +
                            (U_V_Cross_Norm_Array_2[i] * Unit_X[2]);
            fixed_t dot_Y = (U_V_Cross_Norm_Array_0[i] * Unit_Y[0]) +
                            (U_V_Cross_Norm_Array_1[i] * Unit_Y[1]) +
                            (U_V_Cross_Norm_Array_2[i] * Unit_Y[2]);
            fixed_t dot_Z = (U_V_Cross_Norm_Array_0[i] * Unit_Z[0]) +
                            (U_V_Cross_Norm_Array_1[i] * Unit_Z[1]) +
                            (U_V_Cross_Norm_Array_2[i] * Unit_Z[2]);

            // Compute absolute sum of dot products using hls::abs
            fixed_t abs_X = hls::abs(dot_X);
            fixed_t abs_Y = hls::abs(dot_Y);
            fixed_t abs_Z = hls::abs(dot_Z);

            fixed_t abs_sum = abs_X + abs_Y + abs_Z;

            // Compute point score
            Point_Score_Out_0[i] = Cross_Area_Out_0[i] * abs_sum;

            // Debugging Point Score Out
#ifndef __SYNTHESIS__
            std::cout << "Point_Score_Out_0: " << Point_Score_Out_0[i] << " for i: " << i << std::endl;
#endif

        }

        //***************************************
        // Computes the TOTAL SCORE of each print position
        //***************************************

        fixed_t Position_Score_Out[Pos_Score_Ind] = {0};
        const int Count_12 = Num_Faces;

        // Used to add all scores for current position
        fixed_t Temp_Score = 0;
        fixed_t Temp_Avg_Score = 0;

        // K is used to iterate through the current position face scores;
        // Each input row is a 3D point; Three points is a face; Adding all faces of a print position will give current print position score
        // All position scores will be compared to find the best
        int k = 0;

        Position_Score_Loop:
        for(int i = 0; i < Count_12; i++, k++)
        {
            #pragma HLS PIPELINE II=7

            Temp_Score = Temp_Score + Point_Score_Out_0[k];
        }

        Temp_Avg_Score = Temp_Score / Num_Faces;
        Position_Score_Out[Score_Index] = Temp_Avg_Score;

        // Debugging Score index and Position Score Out
        #ifndef __SYNTHESIS__

        std::cout << "Score_Index: " << Score_Index << std::endl;
        std::cout << "Position_Score_Out: " << Position_Score_Out[Score_Index]  << std::endl;
		#endif

        If_Position:
        if(Position_Score_Out[Score_Index] < Temp_Best_Score){

        	Best_Score_Element_Out = Score_Index;
            std::cout << "Score_Index In Position: " << Score_Index << std::endl;

                *Best_Score = Position_Score_Out[Score_Index];
                Temp_Best_Score = Position_Score_Out[Score_Index];


                //Debugging Best Score
                #ifndef __SYNTHESIS__
                std::cout << "Best Score: " << *Best_Score << std::endl;
				#endif
        }

        Score_Index++;
        Temp_Score = 0;
        Temp_Avg_Score = 0;

    }

//    const int Best_Elem_Loop = (Angle_Steps + 1) * (Angle_Steps + 1);

    Best_Score_Loop:
    for (int i = 0; i < Best_Elem_Loop; i++)
    {
        if (i == Best_Score_Element_Out)
        {
            *Best_X = Element_X[i];
            *Best_Y = Element_Y[i];
        }

    }

    	// Debugging Best X and Best Y
		#ifndef __SYNTHESIS__

    	std::cout << "Best X: " << *Best_X << " Best Y: " << *Best_Y<< std::endl;


		#endif

}
}
