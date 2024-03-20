#ifndef PAN_TOMPKINS_H
#define PAN_TOMPKINS_H

typedef float dataType; //int or float
typedef enum {false, true} bool_enum;

//#define DEBUG 1  //ifndef DEBUG we compute the HR avg 
                    //ifdef DEBUG we see the R_locs

#define N 5528         // ECG length
#define NC_Lo 13		// Number of low-pass filter components
#define NC_Hi 33        // Number of high-pass filter components
#define NC_Der 5        // Number of components in the derivative transfer function
#define NC_Int 31       // Number of components in the window integration
#define FS 128			// Sampling frequency
#define DELAY 52
#define MAX_VALUE 2400
#define N_AVG 8

#define MIN(a, b) (((a)<(b))? (a): (b))
#define MAX(a, b) (((a)>(b))? (a): (b))

#define BUFFER_SIZE 205

#endif