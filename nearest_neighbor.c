/* The MIT License

   Copyright (c) 2018 Giulio Genovese

   Author: Giulio Genovese <giulio.genovese@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

#include <stdlib.h>

int elementsInBin[12];
int *binData[12];
int elementsInShiftedBin[11];
int *binDataShifted[11];

int findClosestSitesToPointsAlongAxis(
    int n_raw,
    float *raw_x,
    float *raw_y,
    int n_axis,
    float *axis_x,
    float *axis_y,
    int *ret)
{
    int i;
    float *raw_a = NULL;
    float *raw_b = NULL;
    float *axis_a = NULL;
    float axis_max_val;
    float bin_width;
    int bin_idx;
    float quotient;
    float reminder;
    int *curr_bin_data;
    int curr_bin_size;
    float curr_axis_x;
    float curr_axis_y;
    float x_dist;
    float y_dist;
    double best_val;
    int best_idx;
    int j;
    int curr_idx;
    double sq_dist;
    double axis_max_dist;
    int use_y = 1;
    int use_x = 1;

    for (i=0; i<n_axis; i++)
    {
        if ( axis_x[i] > 0.0001 )
        {
            use_y = 0;
            break;
        }
    }

    for (i=0; i<n_axis; i++)
    {
        if ( axis_y[i] > 0.0001 )
        {
            use_x = 0;
            break;
        }
    }

    if ( use_y )
    {
        raw_a = raw_y;
        raw_b = raw_x;
        axis_a = axis_y;
    }
    else if ( use_x )
    {
        raw_a = raw_x;
        raw_b = raw_y;
        axis_a = axis_x;
    }
    else
    {
        return -1;
    }

    axis_max_val = axis_a[n_axis - 1];
    bin_width = axis_max_val / 12.0f;
    axis_max_dist = (double)bin_width;

    for (i=0; i<n_raw; i++)
    {
      if ( (double)raw_b[i] > axis_max_dist ) continue;
        bin_idx = (int)(raw_a[i] / bin_width);
        if ( bin_idx < 0 ) bin_idx = 0;
        if ( bin_idx > 11 ) bin_idx = 11;
        elementsInBin[bin_idx]++;
        bin_idx = (int)(raw_a[i] / bin_width - 0.5f);
        if ( bin_idx < 0 ) bin_idx = 0;
        if ( bin_idx > 10 ) bin_idx = 10;
        elementsInShiftedBin[bin_idx]++;
     }

    for (i=0; i<=11; i++)
    {
        binData[i] = (int *)malloc( (size_t)elementsInBin[i] * sizeof(int) );
        elementsInBin[i] = 0;
        if ( i == 11 ) continue;
        binDataShifted[i] = (int *)malloc( (size_t)elementsInShiftedBin[i] * sizeof(int) );
        elementsInShiftedBin[i] = 0;
    }

    for (i=0; i<n_raw; i++)
    {
      if ( (double)raw_b[i] > axis_max_dist ) continue;
        bin_idx = (int)(raw_a[i] / bin_width);
        if ( bin_idx < 0 ) bin_idx = 0;
        if ( bin_idx > 11 ) bin_idx = 11;
        binData[bin_idx][elementsInBin[bin_idx]] = i;
        elementsInBin[bin_idx]++;
        bin_idx = (int)(raw_a[i] / bin_width - 0.5f);
        if ( bin_idx < 0 ) bin_idx = 0;
        if ( bin_idx > 10 ) bin_idx = 10;
        binDataShifted[bin_idx][elementsInShiftedBin[bin_idx]] = i;
        elementsInShiftedBin[bin_idx]++;
    }

    for (i=0; i<n_axis; i++)
    {
        quotient = axis_a[i] / bin_width;
        bin_idx = (int)quotient;
        reminder = quotient - (float)bin_idx;
        curr_bin_data = NULL;
        curr_bin_size = 0;
        if ( bin_idx < 0 ) bin_idx = 0;
        if ( bin_idx > 11 ) bin_idx = 11;

        if ( 0.25f <= reminder && reminder <= 0.75f )
        {
            curr_bin_data = binData[bin_idx];
            curr_bin_size = elementsInBin[bin_idx];
        }
        else
        {
            if ( reminder < 0.25f )
            {
                if ( bin_idx == 0 )
                {
                    curr_bin_data = binData[bin_idx];
                    curr_bin_size = elementsInBin[bin_idx];
                }
                else
                {
                    curr_bin_data = binDataShifted[bin_idx-1];
                    curr_bin_size = elementsInShiftedBin[bin_idx-1];
                }
            }
            else if ( bin_idx == 11 )
            {
                curr_bin_data = binData[bin_idx];
                curr_bin_size = elementsInBin[bin_idx];
            }
            else
            {
                curr_bin_data = binDataShifted[bin_idx];
                curr_bin_size = elementsInShiftedBin[bin_idx];
            }
        }

        curr_axis_x = axis_x[i];
        curr_axis_y = axis_y[i];
        best_val = 1e20;
        best_idx = -1;

        for (j=0; j<curr_bin_size; j++)
        {
            curr_idx = curr_bin_data[j];
            x_dist = raw_x[curr_idx] - curr_axis_x;
            y_dist = raw_y[curr_idx] - curr_axis_y;
            sq_dist = (double)(x_dist * x_dist + y_dist * y_dist);
            if ( sq_dist < best_val )
            {
                best_val = sq_dist;
                best_idx = curr_idx;
            }
        }

        ret[i] = best_idx;
    }

    for (i=0; i<=11; i++)
    {
        free( (void *)binData[i] );
        elementsInBin[i] = 0;
        if ( i > 10 ) continue;
        free( (void *)binDataShifted[i] );
        elementsInShiftedBin[i] = 0;
    }

    return 0;
}
