#ifndef SHAPE_H
#define SHAPE_H 

#include <cmath>

float get_eta(float zeta, int geo) {

    float eta = 0;      // dimensionless radius

    // parabola
    if (geo == 0) {     

        eta = 1.0 - (zeta * zeta);

    }

    // giromill
    else if (geo == 1) { 

        eta = 1.0;
    }

    return eta;

}

float get_del(float zeta, float B, int geo) {

    float del = 0;      // blade's angle

    // parabola
    if (geo == 0) {     

        del = atan(2.0 * B * zeta);

    }

    // giromill
    else if (geo == 1) {

        del = 0.0;
    }

    return del;

}

float get_ar(float B, float H, float C, int geo) {

    float AR = 0;       // blade's aspect ratio

    // parabola
    if (geo == 0) {

        float t1 = std::sqrt(1.0 + 4 * B * B);
        float t2 = 1.0 / (2.0 * B);
        float t3 = std::log(2 * B + t1);
        float len = H * (t1 + (t2 * t3));

        AR = len / C;

    }

    // giromill
    else if (geo == 1) {

        AR = 2 * H / C;
    }

    return AR;
}

#endif