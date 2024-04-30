#include <stdio.h>
#include <SDL2/SDL.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include "vector.h"
#include "mesh.h"
#include "triangle.h"
#include "matrix.h"

// Define constants
#define N_MESH_FACES 12
#define FPS 60
#define FRAME_TARGET_TIME (1000 / FPS)
// Global variables
SDL_Renderer *renderer = NULL;
SDL_Window *window = NULL;
SDL_Texture *texture = NULL;
uint32_t *color_buffer = NULL;
bool is_running = false;
int window_width;
int window_height;
int scaling_factor = 700;
int previous_frame_time = 0;

// Define functions
bool initialize_windowing_system();
void clean_up_windowing_system();
void process_keyboard_input();
void setup_memory_buffer();
void clear_color_buffer(uint32_t color);
void draw_pixel(int x, int y, uint32_t color);
void draw_line(int x0, int y0, int x1, int y1, uint32_t color);
void draw_triangle(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t color);
void draw_rectangle(int x, int y, int width, int height, uint32_t color);
void project_model();
void update_state();
vec2_t perspective_projection_point(vec3_t point_3d);
vec3_t vec3_rotate_z(vec3_t v, float angle);
vec3_t vec3_rotate_y(vec3_t v, float angle);
vec3_t vec3_rotate_x(vec3_t v, float angle);
float vec3_length(vec3_t v);
vec3_t vec3_add(vec3_t a, vec3_t b);
vec3_t vec3_sub(vec3_t a, vec3_t b);
vec3_t vec3_mul(vec3_t v, float scalar);
vec3_t vec3_div(vec3_t v, float scalar);
vec3_t vec3_cross(vec3_t a, vec3_t b);
float vec3_dot(vec3_t a, vec3_t b);
float vec2_length(vec2_t v);
vec2_t vec2_add(vec2_t a, vec2_t b);
vec2_t vec2_sub(vec2_t a, vec2_t b);
vec2_t vec2_mul(vec2_t v, float scalar);
vec2_t vec2_div(vec2_t v, float scalar);
float vec2_dot(vec2_t a, vec2_t b);
mat4_t mat4_identity();
mat4_t mat4_make_scale(float sx, float sy, float sz);
mat4_t mat4_make_translate(float tx, float ty, float tz);
mat4_t mat4_make_rotation_x(float angle);
mat4_t mat4_make_rotation_y(float angle);
mat4_t mat4_make_rotation_z(float angle);
vec4_t mat4_mul_vec4(mat4_t M, vec4_t v);
vec4_t vec4_from_vec3(vec3_t v);
vec3_t vec3_from_vec4(vec4_t v);

// Global variables for the cube
vec3_t camera_position = {.x = 0, .y = 0, .z = 0};
int t_cnt = 0;
triangle_t triangles_to_render[1000];

vec3_t cube_scale = {.x = 1, .y = 1, .z = 1};
vec3_t cube_translate = {.x = 1, .y = 1, .z = 7};
vec3_t cube_rotation = {.x = 0, .y = 0, .z = 0};

vec3_t mesh_vertices[N_MESH_VERTICES] = {
    {.x = -1, .y = -1, .z = -1}, // 1
    {.x = -1, .y = -1, .z = 1},  // 2
    {.x = 1, .y = -1, .z = -1},   // 3
    {.x = 1, .y = -1, .z = 1},    // 4
    {.x = 0, .y = 1, .z = 0},     // 5
};

face_t mesh_faces[N_MESH_FACES] = {
    // front
    {.a = 1, .b = 2, .c = 3}, // bottom-left, top-left, top-right
    {.a = 1, .b = 3, .c = 4}, // bottom-left, top-right, bottom-right

    // right
    {.a = 4, .b = 3, .c = 5}, // bottom-right, top-right, top-left
    {.a = 4, .b = 5, .c = 6}, // bottom-right, top-left, bottom-left

    // back
    {.a = 6, .b = 5, .c = 7}, // bottom-left, top-left, top-right
    {.a = 6, .b = 7, .c = 8}, // bottom-left, top-right, bottom-right

    // left
    {.a = 8, .b = 7, .c = 2}, // bottom-right, top-right, top-left
    {.a = 8, .b = 2, .c = 1}, // bottom-right, top-left, bottom-left

    // top
    {.a = 2, .b = 7, .c = 5}, // bottom-left, top-right, top-left
    {.a = 2, .b = 5, .c = 3}, // bottom-left, bottom-right, top-right

    // bottom
    {.a = 6, .b = 8, .c = 1}, // bottom-right, top-right, bottom-left
    {.a = 6, .b = 1, .c = 4}  // bottom-right, bottom-left, top-left
};

bool initialize_windowing_system()
{
    if (SDL_Init(SDL_INIT_EVERYTHING) != 0)
    {
        fprintf(stderr, "SDL_Init failed\n");
        return false;
    }

    // Query SDL for display resolution
    SDL_DisplayMode display_mode;
    SDL_GetCurrentDisplayMode(0, &display_mode);
    window_width = display_mode.w;
    window_height = display_mode.h;

    window = SDL_CreateWindow(NULL, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, window_width, window_height, SDL_WINDOW_BORDERLESS);

    if (!window)
    {
        fprintf(stderr, "SDL_CreateWindow() Failed\n");
        return false;
    }

    renderer = SDL_CreateRenderer(window, -1, 0);
    if (!renderer)
    {
        fprintf(stderr, "SDL_CreateRenderer() Failed\n");
        return false;
    }

    // SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN);
    return true;
}

void clean_up_windowing_system()
{
    free(color_buffer);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

void process_keyboard_input(void)
{
    SDL_Event event;       // create the union
    SDL_PollEvent(&event); // fill the union with the current state
    // check the state of the union event
    switch (event.type)
    {
    case SDL_QUIT:
        is_running = false;
        break;
    case SDL_KEYDOWN:
        if (event.key.keysym.sym == SDLK_ESCAPE)
            is_running = false;
        break;
    }
}

void setup_memory_buffer(void)
{
    // Allocate the required memory in bytes to hold the color buffer
    color_buffer = (uint32_t *)malloc(window_width * window_height * sizeof(uint32_t));

    // Creating a SDL texture that will be use to display the color buffer
    texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STREAMING, window_width, window_height);
}

void clear_color_buffer(uint32_t color)
{
    for (int i = 0; i < window_width * window_height; i++)
        color_buffer[i] = color;
}

void draw_triangle(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t color)
{
    draw_line(x0, y0, x1, y1, color);
    draw_line(x1, y1, x2, y2, color);
    draw_line(x2, y2, x0, y0, color);
}

void draw_line(int x0, int y0, int x1, int y1, uint32_t color)
{
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy, e2; /* error value e_xy */

    for (;;)
    { /* loop */
        draw_pixel(x0, y0, 0xFFFFFFFF);
        if (x0 == x1 && y0 == y1)
        {
            break;
        }
        e2 = 2 * err;
        if (e2 >= dy)
        {
            err += dy;
            x0 += sx;
        } /* e_xy+e_x > 0 */
        if (e2 <= dx)
        {
            err += dx;
            y0 += sy;
        } /* e_xy+e_y < 0 */
    }
}

void draw_rectangle(int x, int y, int width, int height, uint32_t color)
{
    for (int i = y; i < height + y; i++)
    {
        for (int j = x; j < width + x; j++)
        {
            draw_pixel(j, i, color);
        }
    }
}

void draw_pixel(int x, int y, uint32_t color)
{
    // confirm pixel is in the visible window space
    if (x >= 0 && x < window_width && y >= 0 && y < window_height)
        color_buffer[(y * window_width) + x] = color;
}

vec2_t perspective_projection_point(vec3_t point_3d)
{
    vec2_t projected_point = {.x = (scaling_factor * point_3d.x) / point_3d.z,
                              .y = (scaling_factor * point_3d.y) / point_3d.z};
    return projected_point;
}

void project_model(void)
{
    for (int i = 0; i < N_MESH_FACES; i++)
    {
        face_t mesh_face = mesh_faces[i];

        vec3_t face_vertices[3];
        face_vertices[0] = mesh_vertices[mesh_face.a - 1];
        face_vertices[1] = mesh_vertices[mesh_face.b - 1];
        face_vertices[2] = mesh_vertices[mesh_face.c - 1];

        mat4_t scale_matrix = mat4_make_scale(cube_scale.x, cube_scale.y, cube_scale.z);
        mat4_t rotation_matrix_x = mat4_make_rotation_x(cube_rotation.x); // pass the angle as float
        mat4_t rotation_matrix_y = mat4_make_rotation_y(cube_rotation.y);
        mat4_t rotation_matrix_z = mat4_make_rotation_z(cube_rotation.z);
        mat4_t translate_matrix = mat4_make_translate(cube_translate.x, cube_translate.y, cube_translate.z);

        vec3_t transformed_vertices[3];

        mat4_t world_matrix = mat4_identity();
        world_matrix = mat4_mul_mat4(scale_matrix, world_matrix);
        world_matrix = mat4_mul_mat4(rotation_matrix_x, world_matrix);
        world_matrix = mat4_mul_mat4(rotation_matrix_y, world_matrix);
        world_matrix = mat4_mul_mat4(rotation_matrix_z, world_matrix);
        world_matrix = mat4_mul_mat4(translate_matrix, world_matrix);

        /* TRANSFORMAION */
        // Sub Loop, apply transformations to each vertex of the current face
        for (int j = 0; j < 3; j++)
        {

            vec4_t transformed_vertex = vec4_from_vec3(face_vertices[j]);

            // all we need is this singlue transofrmation
            transformed_vertex = mat4_mul_vec4(world_matrix, transformed_vertex);

            transformed_vertices[j] = vec3_from_vec4(transformed_vertex);
        }

        /* BACKFACE CULLING  */
        vec3_t vertex_a = transformed_vertices[0]; /*    A    */
        vec3_t vertex_b = transformed_vertices[1]; /*   / \   */
        vec3_t vertex_c = transformed_vertices[2]; /*  C---B  */

        // Get the vector subtractiono of B-A and C-A
        vec3_t vector_ab = vec3_sub(vertex_b, vertex_a);
        vec3_t vector_ac = vec3_sub(vertex_c, vertex_a);

        // Compute the face normal (using corss product to find perpendiculiar vector)
        vec3_t normal = vec3_cross(vector_ab, vector_ac);

        // Find the vector between a point in the triangle and the camera origin.
        vec3_t camera_ray = vec3_sub(camera_position, vertex_a);

        // Calculate how aligned the camera ray is with the face normal (using the dot product)
        float dot_normal_camera = vec3_dot(camera_ray, normal);

        // Bypass triangles that are looking away from the camera by continuing to next face
        // in main loop
        if (dot_normal_camera < 0)
        {
            continue;
        }

        /* PROJECTION */

        triangle_t projected_triangle;

        // Sub Loop, project the vertices of curent face
        for (int j = 0; j < 3; j++)
        {
            vec2_t projected_point = perspective_projection_point(transformed_vertices[j]);

            /* translate projected vertex to center of screen */
            projected_point.x += (window_width / 2);
            projected_point.y += (window_height / 2);

            projected_triangle.points[j] = projected_point;
        }
        // save the projected triangle
        triangles_to_render[t_cnt++] = projected_triangle;
    }
     for (int i = 0; i < t_cnt; i++)
    {
        triangle_t triangle = triangles_to_render[i];
        draw_rectangle(triangle.points[0].x, triangle.points[0].y, 2, 2, 0xFF00FF00);
        draw_rectangle(triangle.points[1].x, triangle.points[1].y, 2, 2, 0xFF00FF00);
        draw_rectangle(triangle.points[2].x, triangle.points[2].y, 2, 2, 0xFF00FF00);
        draw_triangle(
            triangle.points[0].x,
            triangle.points[0].y,
            triangle.points[1].x,
            triangle.points[1].y,
            triangle.points[2].x,
            triangle.points[2].y,
            0xFF00FF00);
    }
    t_cnt = 0;
}

vec3_t vec3_rotate_z(vec3_t v, float angle)
{
    vec3_t rotated_vector = {
        .x = v.x * cos(angle) - v.y * sin(angle),
        .y = v.x * sin(angle) + v.y * cos(angle),
        .z = v.z};
    return rotated_vector;
}

vec3_t vec3_rotate_y(vec3_t v, float angle)
{
    vec3_t rotated_vector = {
        .x = v.x * cos(angle) - v.z * sin(angle),
        .y = v.y,
        .z = v.x * sin(angle) + v.z * cos(angle),
    };
    return rotated_vector;
}

vec3_t vec3_rotate_x(vec3_t v, float angle)
{
    vec3_t rotated_vector = {
        .x = v.x,
        .y = v.y * cos(angle) - v.z * sin(angle),
        .z = v.y * sin(angle) + v.z * cos(angle),
    };
    return rotated_vector;
}

float vec3_length(vec3_t v)
{
    return sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
}

vec3_t vec3_add(vec3_t a, vec3_t b)
{
    vec3_t result = {
        .x = a.x + b.x,
        .y = a.y + b.y,
        .z = a.z + b.z};
    return result;
}

vec3_t vec3_sub(vec3_t a, vec3_t b)
{
    vec3_t result = {
        .x = a.x - b.x,
        .y = a.y - b.y,
        .z = a.z - b.z};
    return result;
}

vec3_t vec3_mul(vec3_t v, float scalar)
{
    vec3_t result = {
        .x = v.x * scalar,
        .y = v.y * scalar,
        .z = v.z * scalar};
    return result;
}

vec3_t vec3_div(vec3_t v, float scalar)
{
    vec3_t result = {
        .x = v.x / scalar,
        .y = v.y / scalar,
        .z = v.z / scalar};
    return result;
}

vec3_t vec3_cross(vec3_t a, vec3_t b)
{
    vec3_t result = {
        .x = a.y * b.z - a.z * b.y,
        .y = a.z * b.x - a.x * b.z,
        .z = a.x * b.y - a.y * b.x};
    return result;
}

// vec.c implementation
float vec3_dot(vec3_t a, vec3_t b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);

    // we should do these for vec2_t as well
}

// Vec2 methods
float vec2_length(vec2_t v)
{
    return sqrt((v.x * v.x) + (v.y * v.y));
}

vec2_t vec2_add(vec2_t a, vec2_t b)
{
    vec2_t result = {
        .x = a.x + b.x,
        .y = a.y + b.y};
    return result;
}

vec2_t vec2_sub(vec2_t a, vec2_t b)
{
    vec2_t result = {
        .x = a.x - b.x,
        .y = a.y - b.y};
    return result;
}

vec2_t vec2_mul(vec2_t v, float scalar)
{
    vec2_t result = {
        .x = v.x * scalar,
        .y = v.y * scalar};
    return result;
}

vec2_t vec2_div(vec2_t v, float scalar)
{
    vec2_t result = {
        .x = v.x / scalar,
        .y = v.y / scalar};
    return result;
}

float vec2_dot(vec2_t a, vec2_t b)
{
    return (a.x * b.x) + (a.y * b.y);
}

mat4_t mat4_identity(void)
{

    mat4_t m = {{{1, 0, 0, 0},
                 {0, 1, 0, 0},
                 {0, 0, 1, 0},
                 {0, 0, 0, 1}}};

    return m;
}

mat4_t mat4_make_scale(float sx, float sy, float sz)
{
    // | sx  0  0  0 |
    // |  0 sy  0  0 |
    // |  0  0 sz  0 |
    // |  0  0  0  1 |
    mat4_t m = mat4_identity();
    m.m[0][0] = sx;
    m.m[1][1] = sy;
    m.m[2][2] = sz;

    return m;
}

mat4_t mat4_make_translate(float tx, float ty, float tz)
{
    // | 1  0  0  tx |
    // |  0 1  0  ty |
    // |  0  0 1  tz |
    // |  0  0  0  1 |

    mat4_t M = mat4_identity();

    M.m[0][3] = tx;
    M.m[1][3] = ty;
    M.m[2][3] = tz;

    return M;
}

mat4_t mat4_make_rotation_x(float angle)
{
    float c = cos(angle);
    float s = sin(angle);

    // | 1  0   0  0 |
    // | 0  c  -s  0 |
    // | 0  s   c  0 |
    // | 0  0   0  1 |

    mat4_t M = mat4_identity();
    M.m[1][1] = c;
    M.m[1][2] = -s;
    M.m[2][1] = s;
    M.m[2][2] = c;

    return M;
}

mat4_t mat4_make_rotation_y(float angle)
{
    float c = cos(angle);
    float s = sin(angle);

    // |  c  0   s  0 |
    // |  0  1   0  0 |
    // | -s  s   c  0 |
    // |  0  0   0  1 |

    mat4_t M = mat4_identity();
    M.m[0][0] = c;
    M.m[0][2] = -s;
    M.m[2][0] = s;
    M.m[2][2] = c;

    return M;
}

mat4_t mat4_make_rotation_z(float angle)
{
    float c = cos(angle);
    float s = sin(angle);

    // |  c  -s   0  0 |
    // |  s   c   0  0 |
    // |  0   0   0  0 |
    // |  0   0   0  1 |

    mat4_t M = mat4_identity();
    M.m[0][0] = c;
    M.m[0][2] = -s;
    M.m[2][0] = s;
    M.m[2][2] = c;

    return M;
}

vec4_t mat4_mul_vec4(mat4_t M, vec4_t v)
{
    vec4_t result;
    result.x = M.m[0][0] * v.x + M.m[0][1] * v.y + M.m[0][2] * v.z + M.m[0][3] * v.w;
    result.y = M.m[1][0] * v.x + M.m[1][1] * v.y + M.m[1][2] * v.z + M.m[1][3] * v.w;
    result.z = M.m[2][0] * v.x + M.m[2][1] * v.y + M.m[2][2] * v.z + M.m[2][3] * v.w;
    result.w = M.m[3][0] * v.x + M.m[3][1] * v.y + M.m[3][2] * v.z + M.m[3][3] * v.w;

    return result;
}

vec4_t vec4_from_vec3(vec3_t v)
{
    vec4_t result = {v.x, v.y, v.z, 1.0};

    return result;
}

vec3_t vec3_from_vec4(vec4_t v)
{
    vec3_t result = {v.x, v.y, v.z};

    return result;
}

mat4_t mat4_mul_mat4(mat4_t A, mat4_t B)
{
    mat4_t M;
    // apply dot product betwwen rows and colums for every position of the matrix.
    // loop all rows and columns
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            M.m[i][j] = (A.m[i][0] * B.m[0][j]) +
                        (A.m[i][1] * B.m[1][j]) +
                        (A.m[i][2] * B.m[2][j]) +
                        (A.m[i][3] * B.m[3][j]);
    return M;
}

void update_state()
{
    // rotation
    cube_rotation.x += 0.01;
    cube_rotation.y += 0.01;
    cube_rotation.z += 0.01;

    // transforming
    cube_translate.y += .009;
    cube_translate.x += .03;
    cube_scale.x -= .01;
    cube_scale.y -= .01;
    cube_scale.z += .01;

    int time_to_wait = FRAME_TARGET_TIME - (SDL_GetTicks() - previous_frame_time);

    if (time_to_wait > 0 && time_to_wait <= FRAME_TARGET_TIME)
    {
        SDL_Delay((Uint32)time_to_wait);
    }

    // Update previous_frame_time after the delay
    previous_frame_time = SDL_GetTicks();

    clear_color_buffer(0xFF000000);

    project_model();
}

void run_render_pipeline()
{
    SDL_UpdateTexture(texture, NULL, color_buffer, (int)(window_width * sizeof(uint32_t)));
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_RenderPresent(renderer);
}

int main(void)
{
    is_running = initialize_windowing_system();
    setup_memory_buffer();

    while (is_running)
    {
        process_keyboard_input();
        update_state();
        run_render_pipeline();
    }

    clean_up_windowing_system();
    return 0;
}
