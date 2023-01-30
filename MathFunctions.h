#ifndef MATHFUNCTIONS_H_
#define MATHFUNCTIONS_H_

#include <cfloat>

#include "Job.h"

class MathFunctions {
public:
    static void copyVector(const float* _a, float* _b, const int _dims) {
        for(int i=0; i<_dims; i++) {
            _b[i] = _a[i];
        }
    }
    
    static void printVectors(std::vector<Job>* jobs) {
        int dimensions = jobs->at(0).getDimensions();
     //   Add a comment to this line
        for(int j = 0; j < jobs->size(); j++) {
            printf("Job %.3d id: %d vec: (", j, jobs->at(j).getID());
            for(int k = 0; k < dimensions; k++) {
                printf("%f ", jobs->at(j).getKPIVec()[k]);
            }
            printf(")\n");
        }
    }
            
            static void VectorAddition(float* _a, const float* _b, const int _dims) {
                for(int i=0; i<_dims; i++) {
                    _a[i] += _b[i];
                }
            }
            
            static void VectorSubtraction(float* _a, const float* _b, const int _dims) {
                for(int i=0; i<_dims; i++) {
                    _a[i] -= _b[i];
                }
            }
            
            static void VecDiv(float* _a, const float _b, const int _dims) {
                if(_b != 0) {
                    for(unsigned i=0; i<_dims; i++) {
                        _a[i] /= _b;
                    }
                }
            }
            
            static bool areVectorsEqual(const float* _a, const float* _b, const int _dims) {
                for(int i=0; i<_dims; i++) {
                    if(fabs(_a[i] - _b[i]) > 0.0000001) {
                        return false;
                    }
                }
                return true;
            }
            
            static bool areFloatsEqual(const float _a, const float _b) {
                return fabs(_a - _b) < 0.0000001;
            }
            
            static float DotProd(float* _a, float* _b, float* _c, int _dims) {
                float res = 0.0;
                for(unsigned i=0; i<_dims; i++) {
                    res += _a[i] * _b[i] * _c[i];
                }
                return res;
            }
            
            static float DotProd(float* _a, float* _b, int _dims) {
                float res = 0.0;
                for(unsigned i=0; i<_dims; i++) {
                    res += _a[i] * _b[i];
                }
                return res;
            }
            
            static float* VecCoordProd(float* _res, float* _a, const float* _b, int _dims) {
                for(unsigned i=0; i<_dims; i++) {
                    _res[i] = _a[i] * _b[i];
                }
                return _res;
            };
            
            
            /**
             * Vector b is subtracted from vector a. If any component is smaller than 0
             * (or rather -0.00001 to tolerate rounding effects), this function returns
             * false, otherwise it returns true.
             * @param a first vector (usually: remaining capacity of machine).
             * @param b second vector (usually: job requirements).
             * @param _dims (dimension of vectors).
             * @return true if and only if job fits into machine.
             */
            static bool isVectorDifferencePositive(float* a, const float* b, int _dims) {
                float t;
                for(int i=0; i<_dims; i++) {
                    t = a[i] - b[i];
                    if(t < -0.00001) { // not 0 because of rounding effects
                        return false;
                    }
                }
                return true;
            };
            
            static float Linf(float* vec, int _dims) {
                float max = -FLT_MAX;
                for (unsigned i = 0; i < _dims; ++i) {
                    if (max < vec[i]) {
                        max = vec[i];
                    }
                }
                return max;
            };
            
            static bool LInfDecSort(Job &a, Job &b) {
                int dims = a.getDimensions();
                return (Linf(a.getKPIVec(), dims) > Linf(b.getKPIVec(), dims));
            };
            
            static float IntervalRounding(float _left, float _num, bool _leftOpen, float _epsilon) {
                if (_leftOpen) {
                    if ((_left < _num) && (_num <= ((1.0 + _epsilon) * _left))) {
                        return (_left);
                    } else {
                        return IntervalRounding(((1.0 + _epsilon) * _left), _num, true, _epsilon);
                    }
                } else {
                    if ((_left <= _num) && (_num < ((1.0 + _epsilon) * _left))) {
                        return (_left);
                    } else {
                        return IntervalRounding(((1.0 + _epsilon) * _left), _num, false, _epsilon);
                    }
                }
            };
            
            /**
             * This function returns true if at least one coordinate of the given vector is negative.
             * @param _vec vector that is checked.
             * @param _dimensions dimensions of vector.
             * @return true if at least one coordinate is negative, false otherwise.
             */
            static bool isVectorNegative(const float* _vec, int _dimensions) {
                for(unsigned i = 0; i < _dimensions; ++i) {
                    if(_vec[i] < 0.00001) {
                        return true;
                    }
                }
                return false;
            }
            
            /**
             * This function sums up the components of the given vector and returns the sum.
             * @param v vector of floats.
             * @param length length of vector.
             * @return sum of vector components.
             */
            static float sum(float* v, int length) {
                float res = 0.;
                for(int i=0; i<length; i++) {
                    res += v[i];
                }
                return res;
            }
            static float max(float* v, int dimensions) {
                float max = 0.;
                for(int i=0; i<dimensions; i++) {
                    if(v[i]>max)
                        max=v[i];
                }
                return max;
            }
        };
#endif
