#include <stdio.h>
#include <math.h>

#include "pan_tompkins.h"
#include "kernels.def"
#include "data.h"

// The signal array is where the most recent samples are kept. 
    // The other arrays are the outputs of each filtering module
dataType signal_data[BUFFER_SIZE], lowpass_data[BUFFER_SIZE], highpass_data[BUFFER_SIZE], derivative_data[BUFFER_SIZE], squaring_data[BUFFER_SIZE], integrated_data[BUFFER_SIZE];
int R_loc[N];

dataType peaki = 0, spki = 0, npki = 0, threshold1 = 0, threshold2 = 0;
__uint32_t previous_peak = 0;
__uint32_t peak_counter = 0;

int searchback_end = 0;

dataType single_convolution(dataType *x, int n, dataType *h, int nc);
dataType max(dataType *x, int n);

int ecg_processing (dataType *input_data, int sample, int *output_data)
{
    int i = 0, j = 0;

    // Signal status variables
    int current = 0;
    
    // Detection parameters
    int searchback_start = 0;

    bool_enum is_qrs = false;

    float min_rr_width = 0.2*FS;
    int max_rr_width = BUFFER_SIZE;

    int dist_previous_peak = 0;
    int dist_searchback_end = 0;

        // Perform a shifting inside the data buffers after BUFFER_SIZE elements
        if (sample >= BUFFER_SIZE)
        {
            for (i = 0; i < BUFFER_SIZE - 1; i++)
            {
                signal_data[i] = signal_data[i+1];
				lowpass_data[i] = lowpass_data[i+1];
				highpass_data[i] = highpass_data[i+1];
				derivative_data[i] = derivative_data[i+1];
				squaring_data[i] = squaring_data[i+1];
				integrated_data[i] = integrated_data[i+1];
            }
            current = BUFFER_SIZE - 1;
        }
        else
        {
            current = sample;
        }


        // Check validity
        if (isfinite(input_data[current]) == 0)
        {
            input_data[current] = input_data[current-1];
        }

        // Cancellation of DC component and normalization
        if (current >= 1)
			signal_data[current] = input_data[current] - input_data[current-1] + 0.995*signal_data[current-1];
		else
			signal_data[current] = 0;

        signal_data[current] /= MAX_VALUE;

         
        // Low Pass Filtering
        if (current > NC_Lo){
            lowpass_data[current] = single_convolution(signal_data + (current - NC_Lo), NC_Lo, h_Lo, NC_Lo);
        }

                
        // High Pass Filtering
        if (current > NC_Hi){
            highpass_data[current] = single_convolution(lowpass_data + (current - NC_Hi), NC_Hi, h_Hi, NC_Hi);
        }


        // Derivative
        if (current > NC_Der){
            derivative_data[current] = single_convolution(highpass_data + (current - NC_Der), NC_Der, h_Der, NC_Der);
        }
        
        // Squaring
        squaring_data[current] = derivative_data[current]*derivative_data[current];
        
        // Moving window integration
        if (current > NC_Int){
            integrated_data[current] = single_convolution(squaring_data + (current - NC_Int), NC_Int, h_Int, NC_Int);
        }
        //printf("integrated_data[%d]:%f\n", current, integrated_data[current]);


        
        // Find R points
        if (current < 2)
            return 0;

        dist_previous_peak = sample - previous_peak;
        dist_searchback_end = sample - searchback_end;

        if( dist_previous_peak < 0 )
        {
            dist_previous_peak += BUFFER_SIZE;   
        }
        
        if( dist_searchback_end < 0)
        {
            dist_searchback_end += BUFFER_SIZE;
        }
        
        // Check if a seachback is required
        if (dist_previous_peak > max_rr_width && dist_searchback_end > max_rr_width) // we use threshold2
        {
            searchback_start = current - dist_previous_peak;
            searchback_end = current;
            for (i = searchback_start; i <= searchback_end; i++)
            {
                // update peaki with current integrated value
                peaki = integrated_data[i]; 
                // Look for a QRS
                is_qrs = false;
                // Found a peak in searchback (threshold2)
                if (peaki > threshold2)
                {
                    spki = 0.750 * spki + 0.250 * peaki;
                    is_qrs = true;
                }
                // Update data structures
                if (is_qrs)
                {
                    if (peak_counter == 0 || dist_previous_peak >= min_rr_width)
                    {
                        // Add a new peak
                        peak_counter = peak_counter + 1;
                        #ifdef DEBUG
                        R_loc[peak_counter-1] = sample;
                        #else
                            R_loc[(peak_counter-1)%N_AVG] = sample;
                        #endif
                    }
                    else if (integrated_data[searchback_start] < peaki)
                    {
                        // Replace previous peak
                        #ifdef DEBUG
                        R_loc[peak_counter-1] = sample;
                        #else
                        R_loc[(peak_counter-1)%N_AVG] = sample;
                        #endif
                    }
                    previous_peak = sample;
                }
                else
                {
                    // Update npki
                    npki = 0.875 * npki + 0.125 * peaki;
                }
                // Adjust thresholds
                threshold1 = npki + 0.25 * (spki - npki);
                threshold2 = 0.5 * threshold1;
            }
        }
        else // we use threshold1
        {
            // update peaki with current integrated value
            peaki = integrated_data[current]; 
            // Look for a QRS
            is_qrs = false;
            // Found a peak (threshold1)
            if (peaki > threshold1)
            {
                spki = 0.875 * spki + 0.125 * peaki;
                is_qrs = true;
            }
            // Update data structures
            if (is_qrs)
            {
                if (peak_counter == 0 || dist_previous_peak >= min_rr_width)
                {
                    // Add a new peak
                    peak_counter = peak_counter + 1;
                    #ifdef DEBUG
                    R_loc[peak_counter-1] = sample;
                    #else
                    R_loc[(peak_counter-1)%N_AVG] = sample;
                    #endif
                }
                else if (integrated_data[current - dist_previous_peak] < peaki)
                {
                    // Replace previous peak
                    #ifdef DEBUG
                    R_loc[peak_counter-1] = sample;
                    #else
                    R_loc[(peak_counter-1)%N_AVG] = sample;
                    #endif
                }
                previous_peak = sample;
            }
            else
            {
                // Update npki
                npki = 0.875 * npki + 0.125 * peaki;
            }
            // Adjust thresholds
            threshold1 = npki + 0.25 * (spki - npki);
            threshold2 = 0.5 * threshold1;
        }
        
        #ifndef DEBUG
        float hr = 0;
        for (int i=1; i<N_AVG; i++)
        {
            int idx = (((peak_counter-1)%N_AVG) + 1 + i) % N_AVG;
            dataType temp = R_loc[idx] - R_loc[(idx+N_AVG-1)% N_AVG];
            if(temp < 0)
            {
                temp += BUFFER_SIZE;    
            }
            hr += temp;
        }
        hr /= (N_AVG - 1);
        //printf("HR @ %d = %f\n", sample, 60.0f/(hr/FS));

        #endif

    return hr;

}


int main()
{
    float hr;
    int sample = 0;  
    
    while (sample < N){
        hr = ecg_processing (sample < BUFFER_SIZE ? input_data : input_data+sample-BUFFER_SIZE+1, sample, R_loc);
        sample++;
        printf("HR @ %d = %f\n", sample, 60.0f/(hr/FS));
    }

    // #ifdef DEBUG
    // while (sample < N){
    //     c = ecg_processing (input_data, sample, R_loc);
    //     sample++;
    // for(i=0; i<c; i++)
    //     printf("R_loc[%d]: %d\n", i, R_loc[i]-DELAY);
    // #endif

    return 0;
}
