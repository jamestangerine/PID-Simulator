#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <emscripten/emscripten.h>

// Signal buffer format: First values are the initial measurements in the order of a, v, and x (if the PID controller controls position); all following values correspond to indices in the indices buffer
// Indices buffer format: First value is the total length of the simulation; all following values correspond to values in the indices buffer
EMSCRIPTEN_KEEPALIVE
void pid_compute(double* signals, uint16_t* indices, double* output, double* canvasOutput, double kP, double kI, double kD, double dt, double maxOutput, double kVoltage, uint8_t kMeasurement, uint8_t intClamp, double intClampMin, double intClampMax, uint8_t intZone, double intZoneThresh, double viewMax, double vRes) {
    const double kdProd = kD/dt;
    const double kiProd = kI*dt;
    double prevError = 0;
    double integral = 0;
    double a;
    if (kMeasurement) {
        double v = *signals;
        uint16_t length = *(indices++);
        for (uint16_t current = 0; current < length; current++) {
            if (current == *indices) {
                signals++;
                indices++;
            }
            double error = *signals - v;
            if (intZone && ((error > 0 && error > intZoneThresh) || (error < 0 && -error > intZoneThresh))) {
                integral = 0;
                a = (kP*error + kdProd*(error-prevError));
            } else {
                integral += error;
                if (intClamp) {
                    double computedITerm = kiProd*integral;
                    a = kP*error + (computedITerm > intClampMax ? intClampMax : (computedITerm < intClampMin ? intClampMin : computedITerm)) + kdProd*(error-prevError);
                } else
                    a = kP*error + kiProd*integral + kdProd*(error-prevError);
            }
            if (a > maxOutput)
                a = maxOutput;
            else if (a < -maxOutput)
                a = -maxOutput;
            a *= kVoltage*12;
            v += a*dt;
            *(output++) = v;
            *(canvasOutput++) = (viewMax-v)*vRes;
            prevError = error;
        }
    } else {
        double v = *(signals++);
        double x = *signals;
        uint16_t length = *(indices++);
        for (uint16_t current = 0; current < length; current++) {
            if (current == *indices) {
                signals++;
                indices++;
            }
            double error = *signals - x;
            if (intZone && ((error > 0 && error > intZoneThresh) || (error < 0 && -error > intZoneThresh))) {
                integral = 0;
                a = (kP*error + kdProd*(error-prevError));
            } else {
                integral += error;
                if (intClamp) {
                    double computedITerm = kiProd*integral;
                    a = kP*error + (computedITerm > intClampMax ? intClampMax : (computedITerm < intClampMin ? intClampMin : computedITerm)) + kdProd*(error-prevError);
                } else
                    a = kP*error + kiProd*integral + kdProd*(error-prevError);
            }
            if (a > maxOutput)
                a = maxOutput;
            else if (a < -maxOutput)
                a = -maxOutput;
            a *= kVoltage*12;
            v += a*dt;
            x += v*dt;
            *(output++) = x;
            *(canvasOutput++) = (viewMax-x)*vRes;
            prevError = error;
        }
    }
}