#include <stdio.h>
#include <SDL2/SDL.h>
#include <stdbool.h>
#include <stdint.h>
#include "vector.h"
#include "mesh.h"
#include "triangle.h"

SDL_Renderer *renderer = NULL;
SDL_Window *window = NULL;
SDL_Texture *texture = NULL;
uint32_t *color_buffer = NULL;
bool is_running = false;
int window_width;
int window_height;
void draw_pixel(int x, int y, uint32_t color);
void draw_line(int x0, int y0, int x1, int y1, uint32_t color);
void draw_triangle(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t color);
void draw_rectangle(int x, int y, int width, int height, uint32_t color);
bool initialize_windowing_system();
void clean_up_windowing_system();
void project_model(void);

int scaling_factor = 700;
int previous_frame_time = 0;

#define FPS 60
#define FRAME_TARGET_TIME (1000 / FPS)

vec3_t camera_position = {.x = 0, .y = 0, .z = 0};
int t_cnt = 0;
vec3_t cube_rotation = {.x = 0, .y = 0, .z = 0};
triangle_t triangles_to_render[1000];

vec3_t mesh_vertices[8] = {
    {.x = -1, .y = -1, .z = -1}, // 1
    {.x = -1, .y = 1, .z = -1},  // 2
    {.x = 1, .y = 1, .z = -1},   // 3
    {.x = 1, .y = -1, .z = -1},  // 4
    {.x = 1, .y = 1, .z = 1},    // 5
    {.x = 1, .y = -1, .z = 1},   // 6
    {.x = -1, .y = 1, .z = 1},   // 7
    {.x = -1, .y = -1, .z = 1}   // 8
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
    int delta_x = x1 - x0;
    int delta_y = y1 - y0;

    int steps = abs(delta_x) >= abs(delta_y) ? abs(delta_x) : abs(delta_y);

    float x_inc = delta_x / (float)steps;
    float y_inc = delta_y / (float)steps;

    float cur_x = x0;
    float cur_y = y0;
    for (int i = 0; i < steps; i++)
    {
        draw_pixel(round(cur_x), round(cur_y), color);
        cur_x += x_inc;
        cur_y += y_inc;
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

        // New array to store transformed vertices
        vec3_t transformed_vertices[3];

        /* TRANSFORMAION */
        // Sub Loop, apply transformations to each vertex of the current face
        for (int j = 0; j < 3; j++)
        {
            vec3_t transformed_vertex = face_vertices[j];

            transformed_vertex = vec3_rotate_x(transformed_vertex, cube_rotation.x);
            transformed_vertex = vec3_rotate_y(transformed_vertex, cube_rotation.y);
            transformed_vertex = vec3_rotate_z(transformed_vertex, cube_rotation.z);

            // translate the vertex away from the camera
            transformed_vertex.z += 5;

            // save transformed vertex in an array of transformed vertices
            transformed_vertices[j] = transformed_vertex;
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
}

void update_state()
{
    cube_rotation.x += 0.01;
    cube_rotation.y += 0.01;
    cube_rotation.z += 0.01;
    int time_to_wait = FRAME_TARGET_TIME - (SDL_GetTicks() - previous_frame_time);

    if (time_to_wait > 0 && time_to_wait <= FRAME_TARGET_TIME)
    {
        SDL_Delay((Uint32)time_to_wait);
    }

    // Update previous_frame_time after the delay
    previous_frame_time = SDL_GetTicks();

    clear_color_buffer(0xFF000000);

    project_model();

    for (int i = 0; i < t_cnt; i++)
    {
        triangle_t triangle = triangles_to_render[i];
        draw_rectangle(triangle.points[0].x, triangle.points[0].y, 2,2, 0xFF00FF00);
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
