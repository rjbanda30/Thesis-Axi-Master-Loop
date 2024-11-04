// Top_sw.h
#ifndef TOP_SW_H
#define TOP_SW_H

#include "FixedPoint.h"

// Software reference implementation
void Top_sw(fixed_t Ordered_X[], fixed_t Ordered_Y[], fixed_t Ordered_Z[],
            int Num_Rows, int Angle_Steps, fixed_t *Best_Score, fixed_t *Best_X, fixed_t *Best_Y);

#endif // TOP_SW_H
