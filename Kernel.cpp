#include "SPH.hpp"
#include <float.h>
# define M_PI           3.14159265358979323846  /* pi */


// Smoothing function
double Wab(std::vector<double> pos, int partA, int partB, double h, size_t choice)
{   
    // Constant normalisation
    double alphaD;

    // Support domain
    double k; 

    // Compute the distance btween two particle A and B
    double r =  sqrt(distance(pos, partA, partB));

    switch (choice){

    case 1 : // Gausian Kernel
        alphaD = 1.0/(pow(M_PI,3/2)*pow(h,3));
        k = DBL_MAX;
        return alphaD*exp(-pow(r/h,2));
    break;

    case 2 : // Bell-shaped Kernel
        alphaD = 105.0/(16.0*M_PI*pow(h,3));
        k = 1;
        if(r/h <= 1.0)
            return alphaD*((1.0+3.0*(r/h)) * pow((1.0-(r/h)),3));
        else
            return 0.0;
    break;

    case 3 : // Cubic spline Kernel
        alphaD = 3.0/(2.0*M_PI*pow(h,3)); 
        k = 2;
        if(0.0 <= r/h && r/h < 1.0)
            return alphaD*((2.0/3.0) - pow((r/h),2) + (1.0/2.0)*pow((r/h),3));
        else if (1.0 <= r/h && r/h < 2.0)
            return alphaD*((1.0/6.0) * pow((1.0-(r/h)),3));
        else
            return 0.0;
    break;

    case 4 : // Quadratic Kernel
        alphaD = 5.0/(4.0*M_PI*pow(h,3)); 
        k = 2;
        if(0.0 <= r/h && r/h <= 2.0)
            return alphaD*((3.0/16.0)*pow((r/h),2) - (3.0/4.0)*(r/h) + (3.0/4.0));
        else
            return 0.0;
    break;

    case 5 : // Quintic Kernel
        alphaD = 21.0/(16.0*M_PI*pow(h,3)); 
        k = 2;
        if(0.0 <= r/h && r/h <= 2.0)
            return alphaD*(pow((1.0-(1.0/2.0)*(r/h)),4) * (2.0*(r/h) + 1.0));
        else
            return 0.0;
    break;

    case 6 : // Quintic spline Kernel
        alphaD = 3.0/(359.0*M_PI*pow(h,3)); 
        k = 3;
        if(0.0 <= r/h && r/h < 1.0)
            return alphaD*(pow((3.0-(r/h)),5) - 6.0*pow((2.0-(r/h)),5) + 15.0*pow((1.0-(r/h)),5));
        else if (1.0 <= r/h && r/h < 2.0)
            return alphaD*(pow((3.0-(r/h)),5) - 6.0*pow((2.0-(r/h)),5));
        else if (2.0 <= r/h && r/h < 3.0)
            return alphaD*(pow((3.0-(r/h)),5));
        else
            return 0.0;
    break;

    default: return 0.0;

    }

}



// Darivative of the smoothing function
double grad_Wab(std::vector<double> pos, int partA, int partB, double h, size_t choice)
{   
    // Constant normalisation
    double alphaD;

    // Support domain
    double k; 

    // Compute the distance btween two particle A and B
    double r =  sqrt(distance(pos, partA, partB));

    switch (choice){

    case 1 : // Gausian Kernel
        alphaD = 1.0/(pow(M_PI,3/2)*pow(h,3));
        k = DBL_MAX;
        return (alphaD/h)*(-2.0*(r/h))*exp(-pow(r/h,2));
    break;

    case 2 : // Bell-shaped Kernel
        alphaD = 105.0/(16.0*M_PI*pow(h,3));
        k = 1;
        if(r/h <= 1.0)
            return (alphaD/h)*3*(pow((1.0-(r/h)),3) - ((1.0+3.0*(r/h)) * pow((1.0-(r/h)),2)));
        else
            return 0.0;
    break;

    case 3 : // Cubic spline Kernel
        alphaD = 3.0/(2.0*M_PI*pow(h,3)); 
        k = 2;
        if(0.0 <= r/h && r/h < 1.0)
            return (alphaD/h)*((2.0/3.0)*pow((r/h),2) - 2.0*(r/h));
        else if (1.0 <= r/h && r/h < 2.0)
            return (alphaD/h)*((-1.0/2.0) * pow((2.0-(r/h)),2));
        else
            return 0.0;
    break;

    case 4 : // Quadratic Kernel
        alphaD = 5.0/(4.0*M_PI*pow(h,3)); 
        k = 2;
        if(0.0 <= r/h && r/h <= 2.0)
            return (alphaD/h)*((3.0/8.0)*(r/h) - (3.0/4.0));
        else
            return 0.0;
    break;

    case 5 : // Quintic Kernel
        alphaD = 21.0/(16.0*M_PI*pow(h,3)); 
        k = 2;
        if(0.0 <= r/h && r/h <= 2.0)
            return (alphaD/h)*((-5.0*(r/h))*pow((1.0-(1.0/2.0)*(r/h)),3));
        else
            return 0.0;
    break;

    case 6 : // Quintic spline Kernel
        alphaD = 3.0/(359.0*M_PI*pow(h,3)); 
        k = 3;
        if(0.0 <= r/h && r/h < 1.0)
            return (alphaD/h)*((-5.0)*pow((3.0-(r/h)),4) + 30.0*pow((2.0-(r/h)),4) - 75.0*pow((1.0-(r/h)),4));
        else if (1.0 <= r/h && r/h < 2.0)
            return (alphaD/h)*((-5.0)*pow((3.0-(r/h)),4) + 30.0*pow((2.0-(r/h)),4));
        else if (2.0 <= r/h && r/h < 3.0)
            return (alphaD/h)*((-5.0)*pow((3.0-(r/h)),4));
        else
            return 0.0;
    break;

    default: return 0.0;

    }

}

