#ifndef MESH_H
#define MESH_H

#include "vector.h"
#include "triangle.h"

#define N_MESH_VERTICES 5
#define N_MESH_FACES 12
#define S_COMET_VERTICES 8
// in mesh.h we declare our array of vertices we use
// extern keyword becuase we will initialize and define
// them in mesh.c
extern vec3_t mesh_vertices[N_MESH_VERTICES];

//declare an array of faces
extern face_t mesh_faces[N_MESH_FACES];

#endif