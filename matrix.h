#ifndef MATRIX_H
#define MATRIX_H

#include "vector.h"

typedef struct {
  float m[4][4];
} mat4_t;

mat4_t mat4_identity(void);

mat4_t mat4_make_scale(float, float, float);

mat4_t mat4_make_translate(float, float, float);

mat4_t mat4_make_rotation_x(float);
mat4_t mat4_make_rotation_y(float);
mat4_t mat4_make_rotation_z(float);

vec4_t mat4_mul_vec4(mat4_t m, vec4_t v);
mat4_t mat4_mul_mat4(mat4_t, mat4_t);  

#endif
