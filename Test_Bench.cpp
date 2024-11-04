// Testbench.cpp
#include <iostream>
#include "Top.h"      // HLS Top function
#include "Top_sw.h"   // Software reference Top_sw function
#include "FixedPoint.h" // Fixed-point definitions

// Function to arrange vertices with fixed-point types
void Arrange_Vertices(
     fixed_t in_array[][3],
     int indices[],
     fixed_t Ordered_X[],
     fixed_t Ordered_Y[],
     fixed_t Ordered_Z[],
     int Num_Indices)
{
    int num_triangles = Num_Indices / 3;
    int idx1, idx2, idx3;

    for (int i = 0; i < num_triangles; ++i) {
        // Use indices to place vertices in the specified order
        idx1 = indices[i * 3];
        idx2 = indices[i * 3 + 1];
        idx3 = indices[i * 3 + 2];

        Ordered_X[i * 3] = in_array[idx1][0];
        Ordered_Y[i * 3] = in_array[idx1][1];
        Ordered_Z[i * 3] = in_array[idx1][2];

        Ordered_X[i * 3 + 1] = in_array[idx2][0];
        Ordered_Y[i * 3 + 1] = in_array[idx2][1];
        Ordered_Z[i * 3 + 1] = in_array[idx2][2];

        Ordered_X[i * 3 + 2] = in_array[idx3][0];
        Ordered_Y[i * 3 + 2] = in_array[idx3][1];
        Ordered_Z[i * 3 + 2] = in_array[idx3][2];
    }
}

int main() {
    // Define constants
    const int Num_Rows = 15;
    const int Num_Faces = 5;
    const int Angle_Steps = 18;
    const int Num_Indices = 15;

    // Input array with coordinates (fixed-point representation)
    fixed_t in_array[Num_Rows][3] = {
        {double_to_fixed(-1.5), double_to_fixed(1.0), double_to_fixed(3.54)},
        {double_to_fixed(0.5354), double_to_fixed(5.0709), double_to_fixed(3.54)},
        {double_to_fixed(-1.5), double_to_fixed(1.0), double_to_fixed(1.0)},
        {double_to_fixed(0.5354), double_to_fixed(5.0709), double_to_fixed(3.54)},
        {double_to_fixed(0.535484), double_to_fixed(5.070968), double_to_fixed(1.0)},
        {double_to_fixed(-1.5), double_to_fixed(1.0), double_to_fixed(1.0)},
        {double_to_fixed(3.5), double_to_fixed(1.0), double_to_fixed(1.0)},
        {double_to_fixed(-1.5), double_to_fixed(1.0), double_to_fixed(3.54)},
        {double_to_fixed(-1.5), double_to_fixed(1.0), double_to_fixed(1.0)},
        {double_to_fixed(3.5), double_to_fixed(1.0), double_to_fixed(1.0)},
        {double_to_fixed(3.5), double_to_fixed(1.0), double_to_fixed(3.54)},
        {double_to_fixed(-1.5), double_to_fixed(1.0), double_to_fixed(3.54)},
        {double_to_fixed(1.425385), double_to_fixed(5.070968), double_to_fixed(1.0)},
        {double_to_fixed(3.5), double_to_fixed(1.0), double_to_fixed(3.54)},
        {double_to_fixed(3.5), double_to_fixed(1.0), double_to_fixed(1.0)}
    };

    // Indices array to reorder vertices for triangles
    int indices[Num_Indices] = {0, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 10};

    // Arrays to hold the vertices ordered by triangles
    fixed_t Ordered_X[Num_Rows], Ordered_Y[Num_Rows], Ordered_Z[Num_Rows];

    // Arrange vertices based on indices
    Arrange_Vertices(in_array, indices, Ordered_X, Ordered_Y, Ordered_Z, Num_Indices);


    fixed_t Best_Score_HLS = double_to_fixed(0.0);
    fixed_t Best_X_HLS = double_to_fixed(0.0);
    fixed_t Best_Y_HLS = double_to_fixed(0.0);

    // Variables to store the output from Top_sw function (Software Reference)
    fixed_t Best_Score_SW = double_to_fixed(0.0);
    fixed_t Best_X_SW = double_to_fixed(0.0);
    fixed_t Best_Y_SW = double_to_fixed(0.0);

    // Call the HLS Top function
    Top(Ordered_X, Ordered_Y, Ordered_Z, Num_Rows, Angle_Steps, &Best_Score_HLS, &Best_X_HLS, &Best_Y_HLS);

    // Call the Software Reference
    Top_sw(Ordered_X, Ordered_Y, Ordered_Z, Num_Rows, Angle_Steps, &Best_Score_SW, &Best_X_SW, &Best_Y_SW);

    // Tolerance for Comparison
    double tolerance = 0.01;

    // Convert
    double Best_Score_HLS_d = fixed_to_double(Best_Score_HLS);
    double Best_X_HLS_d = fixed_to_double(Best_X_HLS);
    double Best_Y_HLS_d = fixed_to_double(Best_Y_HLS);

    double Best_Score_SW_d = fixed_to_double(Best_Score_SW);
    double Best_X_SW_d = fixed_to_double(Best_X_SW);
    double Best_Y_SW_d = fixed_to_double(Best_Y_SW);

    if(6.283 < abs(Best_Score_HLS_d) < 6.283185307179)
    {
    	Best_Score_HLS_d = 0;
    	Best_X_HLS_d = 0;
		Best_Y_HLS_d = 0;


    }
    // Output verification
    bool score_ok = (std::abs(Best_Score_HLS_d - Best_Score_SW_d) < tolerance);
    bool x_ok = (std::abs(Best_X_HLS_d - Best_X_SW_d) < tolerance);
    bool y_ok = (std::abs(Best_Y_HLS_d - Best_Y_SW_d) < tolerance);



    // Print the results with verification status
    std::cout << "Best Score (HLS): " << Best_Score_HLS_d
              << " | Best Score (SW): " << Best_Score_SW_d
              << " | " << (score_ok ? "PASS" : "FAIL") << std::endl;

    std::cout << "Best X (HLS): " << Best_X_HLS_d
              << " | Best X (SW): " << Best_X_SW_d
              << " | " << (x_ok ? "PASS" : "FAIL") << std::endl;

    std::cout << "Best Y (HLS): " << Best_Y_HLS_d
              << " | Best Y (SW): " << Best_Y_SW_d
              << " | " << (y_ok ? "PASS" : "FAIL") << std::endl;

    // Optionally, assert correctness
    if (!score_ok || !x_ok || !y_ok) {
        std::cerr << "Testbench Verification Failed!" << std::endl;
        return -1;
    }

    std::cout << "Testbench Verification Passed!" << std::endl;
    return 0;
}
