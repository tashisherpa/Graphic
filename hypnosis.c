#include <stdio.h>
#include <SDL2/SDL.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include "vector.h"
#include "mesh.h"
#include "triangle.h"
#include "matrix.h"

SDL_Renderer *renderer = NULL;
SDL_Window *window = NULL;
SDL_Texture *texture = NULL;
uint32_t *color_buffer = NULL;
bool is_running = false;

int window_width;
int window_height;

int scaling_factor = 700;
int previous_frame_time = 0;
bool isEnable = true;

#define FPS 60
#define FRAME_TARGET_TIME (1000 / FPS)

bool initialize_windowing_system();
void clean_up_windowing_system();
void process_keyboard_input();
void setup_memory_buffer();
void clear_color_buffer(uint32_t color, bool isEnable);
void draw_line(int x0, int y0, int x1, int y1, uint32_t color);
void draw_rectangle(int x, int y, int width, int height, uint32_t color);
void draw_triangle(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t color);
void draw_hollow_rectangle(int x, int y, int width, int height, uint32_t color);
void draw_pixel(int x, int y, uint32_t color);
void update_state();
void run_render_pipeline();
void project_pyramid(uint32_t color); 

// Array of colors
uint32_t colors[] = {
    0xFF0000FF, // Red
    0xFFFF0000, // Blue
    0xFF00FF00, // Green
    0xFFFFFF00, // Yellow
    0xFFFF00FF, // Purple
    0xFF00FFFF, // Cyan
    0xFF800080, // Purple
    0xFFFFA500, // Orange
};

uint32_t string_color = 0xFF0000FF; // Red color for red string animation 

int num_colors = sizeof(colors) / sizeof(colors[0]);

// Global variables for the cube
vec3_t camera_position = {.x = 0, .y = 0, .z = 0};
int t_cnt = 0;
triangle_t triangles_to_render[1000];

vec3_t scale = {.x = 1, .y = 1, .z = 1};
vec3_t translate = {.x = 0, .y = 0, .z = 2}; // Translate the pyramid along the z-axis to move it away from the camera
vec3_t rotation = {.x = 0, .y = 0, .z = 0};

vec3_t mesh_vertices[N_MESH_VERTICES] = {
    {.x = -0.1, .y = -0.1, .z = -0.1}, // Bottom-left-front vertex
    {.x = 0.1, .y = -0.1, .z = -0.1},  // Bottom-right-front vertex
    {.x = 0.1, .y = -0.1, .z = 0.1},   // Bottom-right-back vertex
    {.x = -0.1, .y = -0.1, .z = 0.1},  // Bottom-left-back vertex
    {.x = 0, .y = 0.1, .z = 0},        // Top vertex
};

face_t mesh_faces[N_MESH_FACES] = {
    {.a = 5, .b = 2, .c = 1}, // front (clockwise)
    {.a = 5, .b = 3, .c = 2}, // right (clockwise)
    {.a = 5, .b = 4, .c = 3}, // back (clockwise)
    {.a = 5, .b = 1, .c = 4}, // left (clockwise)
    {.a = 3, .b = 4, .c = 1}, // bottom (clockwise)
    {.a = 2, .b = 3, .c = 1}, // bottom (clockwise)
};

vec3_t scale_big_comet = {.x = 2, .y = 2, .z = 2};
vec3_t translate_big_comet = {.x = 1.5, .y = -0.5, .z = 2}; // Translate the pyramid along the z-axis to move it away from the camera
vec3_t rotation_big_comet = {.x = 0, .y = 0, .z = 0};

vec3_t big_comet_vertices[N_MESH_VERTICES] = {
    {.x = -0.1, .y = -0.1, .z = -0.1}, // Bottom-left-front vertex
    {.x = 0.1, .y = -0.1, .z = -0.1},  // Bottom-right-front vertex
    {.x = 0.1, .y = -0.1, .z = 0.1},   // Bottom-right-back vertex
    {.x = -0.1, .y = -0.1, .z = 0.1},  // Bottom-left-back vertex
    {.x = 0, .y = 0.1, .z = 0},        // Top vertex
};

face_t big_comet_faces[N_MESH_FACES] = {
    {.a = 5, .b = 2, .c = 1}, // front (clockwise)
    {.a = 5, .b = 3, .c = 2}, // right (clockwise)
    {.a = 5, .b = 4, .c = 3}, // back (clockwise)
    {.a = 5, .b = 1, .c = 4}, // left (clockwise)
    {.a = 3, .b = 4, .c = 1}, // bottom (clockwise)
    {.a = 2, .b = 3, .c = 1}, // bottom (clockwise)
};

vec3_t scale_small_comet = {.x = 0.5, .y = 0.5, .z = 0.5};
vec3_t translate_small_comet = {.x = 1.5, .y = -0.5, .z = 2}; // Translate the pyramid along the z-axis to move it away from the camera
vec3_t rotation_small_comet = {.x = 0, .y = 0, .z = 0};
vec3_t small_comet_vertices[S_COMET_VERTICES] = {
    {.x = -0.1, .y = -0.1, .z = -0.1}, // Bottom-left-front vertex
    {.x = 0.1, .y = -0.1, .z = -0.1},  // Bottom-right-front vertex
    {.x = 0.1, .y = -0.1, .z = 0.1},   // Bottom-right-back vertex
    {.x = -0.1, .y = -0.1, .z = 0.1},  // Bottom-left-back vertex
    {.x = 0, .y = 0.1, .z = 0},        // Top vertex
};

face_t small_comet_faces[N_MESH_FACES] = {
    {.a = 5, .b = 2, .c = 1}, // front (clockwise)
    {.a = 5, .b = 3, .c = 2}, // right (clockwise)
    {.a = 5, .b = 4, .c = 3}, // back (clockwise)
    {.a = 5, .b = 1, .c = 4}, // left (clockwise)
    {.a = 3, .b = 4, .c = 1}, // bottom (clockwise)
    {.a = 2, .b = 3, .c = 1}, // bottom (clockwise)
};

// Draw a hollow rectangle with specified dimensions and color
void draw_hollow_rectangle(int x, int y, int width, int height, uint32_t color)
{
    // Draw top and bottom edges
    for (int i = x; i < x + width; i++)
    {
        draw_pixel(i, y, color); // Top edge
        draw_pixel(i, y + height - 1, color); // Bottom edge
    }

    // Draw left and right edges
    for (int i = y; i < y + height; i++)
    {
        draw_pixel(x, i, color); // Left edge
        draw_pixel(x + width - 1, i, color); // Right edge
    }
}

// Draw a single pixel at the specified coordinates with the given color
void draw_pixel(int x, int y, uint32_t color)
{
    // Ensure the pixel is within the visible window space
    if (x >= 0 && x < window_width && y >= 0 && y < window_height)
        color_buffer[(y * window_width) + x] = color;
}

// Draw a line between two points with the specified color
void draw_line(int x0, int y0, int x1, int y1, uint32_t color)
{
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy, e2; /* error value e_xy */

    for (;;)
    { /* loop */
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
        draw_pixel(x0, y0, color);
    }
}

// Draw a filled rectangle with specified dimensions and color
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

// Draw a triangle between three points with the specified color
void draw_triangle(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t color)
{
    draw_line(x0, y0, x1, y1, color); // Draw line from point 0 to point 1
    draw_line(x1, y1, x2, y2, color); // Draw line from point 1 to point 2
    draw_line(x2, y2, x0, y0, color); // Draw line from point 2 to point 0
}

// Draw a single pixel at the specified coordinates with the given color for sparkle animation
void draw_sparkle_pixel(int x, int y, uint32_t color)
{
    if (x >= 0 && x < window_width && y >= 0 && y < window_height)
        color_buffer[(y * window_width) + x] = color;
}

//Makes the bursts animation 
void sparkle_bursts_animation()
{
    const int num_bursts = 60;          // Number of sparkle bursts
    const float explosion_speed = 5.0f; // Speed multiplier for explosion effect

    int animation_duration_ms = 10000; // Animation duration in milliseconds
    int start_time = SDL_GetTicks();  // Get the current time in milliseconds

    while ((SDL_GetTicks() - start_time <= animation_duration_ms))
    {
        int burst_x = rand() % window_width;
        int burst_y = rand() % window_height;
        uint32_t color = colors[rand() % num_colors]; 

        for (float progress = 0.0f; progress < 1.0f; progress += 0.01f * explosion_speed)
        {
            int sparkle_radius = (int)(progress * 50); // Adjust the multiplier for the burst size

            for (int angle_deg = 0; angle_deg < 360; angle_deg += 10) // Adjust the angle step for sparkles
            {
                float angle_rad = angle_deg * (M_PI / 180.0f);
                int sparkle_x = burst_x + (int)(sparkle_radius * cos(angle_rad));
                int sparkle_y = burst_y + (int)(sparkle_radius * sin(angle_rad));

                // Calculate fade-out based on distance from burst center and progress
                float distance_factor = 1.0f - (progress / 1.0f);
                float fade_factor = distance_factor * distance_factor; // Square the distance_factor for a stronger fade-out effect
                uint32_t faded_color = (color & 0x00FFFFFF) | ((uint32_t)(fade_factor * 255.0f) << 24);

                draw_sparkle_pixel(sparkle_x, sparkle_y, faded_color);
            }

            SDL_UpdateTexture(texture, NULL, color_buffer, window_width * sizeof(uint32_t));
            SDL_RenderClear(renderer);
            SDL_RenderCopy(renderer, texture, NULL, NULL);
            SDL_RenderPresent(renderer);

            SDL_Delay(5); // Adjust the delay between frames for the explosion effect
        }

        // Clear the burst by redrawing the background color after the explosion ends
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); 
        SDL_RenderClear(renderer);
        SDL_RenderPresent(renderer);

        SDL_Delay(100); // Delay before starting the next burst
    }
}

// Function to add two vectors together
vec3_t vec3_add(vec3_t a, vec3_t b)
{
    vec3_t result = {
        .x = a.x + b.x,
        .y = a.y + b.y,
        .z = a.z + b.z};
    return result;
}

// Function to subtract one vector from another
vec3_t vec3_sub(vec3_t a, vec3_t b)
{
    vec3_t result = {
        .x = a.x - b.x,
        .y = a.y - b.y,
        .z = a.z - b.z};
    return result;
}

// Function to multiply a vector by a scalar value
vec3_t vec3_mul(vec3_t v, float scalar)
{
    vec3_t result = {
        .x = v.x * scalar,
        .y = v.y * scalar,
        .z = v.z * scalar};
    return result;
}

// Function to divide a vector by a scalar value
vec3_t vec3_div(vec3_t v, float scalar)
{
    vec3_t result = {
        .x = v.x / scalar,
        .y = v.y / scalar,
        .z = v.z / scalar};
    return result;
}

// Function to calculate the cross product of two vectors
vec3_t vec3_cross(vec3_t a, vec3_t b)
{
    vec3_t result = {
        .x = a.y * b.z - a.z * b.y,
        .y = a.z * b.x - a.x * b.z,
        .z = a.x * b.y - a.y * b.x};
    return result;
}

// Function to calculate the dot product of two vectors
float vec3_dot(vec3_t a, vec3_t b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);

    // we should do these for vec2_t as well
}

// Function to calculate the length of a 2D vector
float vec2_length(vec2_t v)
{
    return sqrt((v.x * v.x) + (v.y * v.y));
}

// Function to create an identity matrix
mat4_t mat4_identity(void)
{

    mat4_t m = {{{1, 0, 0, 0},
                 {0, 1, 0, 0},
                 {0, 0, 1, 0},
                 {0, 0, 0, 1}}};

    return m;
}

// Function to create a scaling matrix
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

// Function to create a translation matrix
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

// Function to create a rotation matrix around the X-axis
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

// Function to create a rotation matrix around the Y-axis
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

// Function to create a rotation matrix around the Z-axis
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

// Function to multiply a 4x4 matrix by a 4D vector
vec4_t mat4_mul_vec4(mat4_t M, vec4_t v)
{
    vec4_t result;
    result.x = M.m[0][0] * v.x + M.m[0][1] * v.y + M.m[0][2] * v.z + M.m[0][3] * v.w;
    result.y = M.m[1][0] * v.x + M.m[1][1] * v.y + M.m[1][2] * v.z + M.m[1][3] * v.w;
    result.z = M.m[2][0] * v.x + M.m[2][1] * v.y + M.m[2][2] * v.z + M.m[2][3] * v.w;
    result.w = M.m[3][0] * v.x + M.m[3][1] * v.y + M.m[3][2] * v.z + M.m[3][3] * v.w;

    return result;
}

// Convert vec3 to vec4 by adding 1 as the fourth component
vec4_t vec4_from_vec3(vec3_t v)
{
    vec4_t result = {v.x, v.y, v.z, 1.0};

    return result;
}

// Convert vec4 to vec3 by ignoring the fourth component
vec3_t vec3_from_vec4(vec4_t v)
{
    vec3_t result = {v.x, v.y, v.z};

    return result;
}

mat4_t mat4_mul_mat4(mat4_t A, mat4_t B)
{
    mat4_t M;
    // Apply matrix multiplication between A and B, storing the result in M
    // Loop through each row and column to calculate each element of the resulting matrix
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            M.m[i][j] = (A.m[i][0] * B.m[0][j]) +
                        (A.m[i][1] * B.m[1][j]) +
                        (A.m[i][2] * B.m[2][j]) +
                        (A.m[i][3] * B.m[3][j]);
    return M;
}

vec2_t perspective_projection_point(vec3_t point_3d)
{
    vec2_t projected_point = {.x = (scaling_factor * point_3d.x) / point_3d.z,
                              .y = (scaling_factor * point_3d.y) / point_3d.z};
    return projected_point;
}

//Function to initialized SDL, window creation and properties, and renderers 
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

    SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN);

    return true;
}

// Free memory and destroy SDL window and renderer
void clean_up_windowing_system()
{
    free(color_buffer);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

void process_keyboard_input(void)
{
    // Poll SDL events for keyboard input
    SDL_Event event;
    SDL_PollEvent(&event);
    // Set is_running to false upon receiving SDL_QUIT or pressing ESC key
    //This is to stop the program by pressing ESC on keyboard 
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

// Allocate memory for color buffer and create SDL texture
void setup_memory_buffer(void)
{
    color_buffer = (uint32_t *)malloc(window_width * window_height * sizeof(uint32_t));
    texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STREAMING, window_width, window_height);
}

// Clear color buffer with specified color if isEnable is true
void clear_color_buffer(uint32_t color, bool isEnable)
{
    if (isEnable)
    {
        for (int i = 0; i < window_width * window_height; i++)
            color_buffer[i] = color;
    }
}

// Clear color buffer and run render pipeline
void update_state()
{
    clear_color_buffer(0xFF000000, isEnable);
    run_render_pipeline();
}

// Update SDL texture with color buffer data and render to window
void run_render_pipeline()
{
    SDL_UpdateTexture(texture, NULL, color_buffer, (int)(window_width * sizeof(uint32_t)));
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_RenderPresent(renderer);
}

// Project and render a pyramid object with backface culling and perspective projection
void project_pyramid(uint32_t color)
{
    for (int i = 0; i < N_MESH_FACES; i++)
    {
        face_t mesh_face = mesh_faces[i];

        vec3_t face_vertices[3];
        face_vertices[0] = mesh_vertices[mesh_face.a - 1];
        face_vertices[1] = mesh_vertices[mesh_face.b - 1];
        face_vertices[2] = mesh_vertices[mesh_face.c - 1];

        mat4_t scale_matrix = mat4_make_scale(scale.x, scale.y, scale.z);
        mat4_t rotation_matrix_x = mat4_make_rotation_x(rotation.x); // pass the angle as float
        mat4_t rotation_matrix_y = mat4_make_rotation_y(rotation.y);
        mat4_t rotation_matrix_z = mat4_make_rotation_z(rotation.z);
        mat4_t translate_matrix = mat4_make_translate(translate.x, translate.y, translate.z);

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
        // draw_rectangle(triangle.points[0].x, triangle.points[0].y, 2, 2, 0xFF00FF00);
        // draw_rectangle(triangle.points[1].x, triangle.points[1].y, 2, 2, 0xFF00FF00);
        // draw_rectangle(triangle.points[2].x, triangle.points[2].y, 2, 2, 0xFF00FF00);
        draw_triangle(
            triangle.points[0].x,
            triangle.points[0].y,
            triangle.points[1].x,
            triangle.points[1].y,
            triangle.points[2].x,
            triangle.points[2].y,
            color);
    }
    t_cnt = 0;
}

// Animate spinning flowers with rotation and drawing triangles
void spinning_flower(float flower_index, u_int32_t color)
{
    float rotation_angle = 0.0f; // Initial rotation angle
    int num_triangles = 10;      // Number of triangles in the flower
    int current_frame_time = SDL_GetTicks();
    float delta_time = (current_frame_time - previous_frame_time) / 2000.0f; // Convert milliseconds to seconds

    // Update the rotation angle based on elapsed time to make it rotate
    float rotation_speed = 0.5f; // Adjust rotation speed as needed
    rotation_angle += rotation_speed * delta_time;

    // Center of the screen
    int x_center = (flower_index + 1) * (window_width / 4);
    int y_center = window_height / 2;

    // Radius of the flower
    int radius = 100 + flower_index * 50;

    for (int i = 0; i < num_triangles; i++)
    {
        // Calculate the angle for this triangle
        float angle = (2 * M_PI / num_triangles) * i + rotation_angle;

        // Calculate the coordinates of the vertices for the triangle
        int x0 = x_center;
        int y0 = y_center;
        int x1 = x_center + radius * cos(angle);
        int y1 = y_center + radius * sin(angle);
        int x2 = x_center + radius * cos(angle + 0.5 * M_PI);
        int y2 = y_center + radius * sin(angle + 0.5 * M_PI);

        // Draw the triangle
        draw_triangle(x0, y0, x1, y1, x2, y2, color);
    }
}

// Animate the drawing of a red string of fate and heart shape using pixel drawing
void red_string_of_fate_animation()
{
    const int line_growth_speed = 5; // Adjust line growth speed as needed
    const int heart_speed = 1;       // Adjust heart drawing speed as needed

    int start_x = -100;                    // Starting x-coordinate of the line (left side of screen off screen)
    int mid_x = window_width / 2;          // middle x-coordinate of the line (middle of the screen)
    int end_x = window_width + 100;        // Ending x-coordinate of the line (right side of screen off screen)
    int line_y = window_height / 2;        // Y-coordinate of the line (middle of the screen)
    int heart_x = mid_x;                   // X-coordinate of the heart starts from the mid_x
    int heart_y = window_height / 2 - 170; // Y-coordinate of the heart

    // Grow the line from left to the middle of the screen
    for (int x = start_x; x <= mid_x; x += line_growth_speed)
    {
        draw_line(start_x, line_y, x, line_y, string_color);
        SDL_UpdateTexture(texture, NULL, color_buffer, (int)window_width * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
        SDL_Delay(10); // Adjust delay for line growth speed
    }

    // Draw the heart shape starting from the bottom right and ending at the bottom left
    for (int i = 180; i >= 0; i -= heart_speed)
    {
        float angle = i * 3.14159f / 180.0f;
        int x = heart_x + 160 * (sin(angle) * sin(angle) * sin(angle));                                       // Increase scale for x-coordinate
        int y = heart_y - 130 * cos(angle) + 50 * cos(2 * angle) + 20 * cos(3 * angle) + 10 * cos(4 * angle); // Increase scale for y-coordinate
        draw_pixel(x - 1, y + 1, string_color);
        draw_pixel(x - 1, y - 1, string_color);
        draw_pixel(x - 1, y, string_color); // Draw the pixel to form the heart shape
        draw_pixel(x, y, string_color);     // Draw the pixel to form the heart shape
        draw_pixel(x + 1, y, string_color); // Draw the pixel to form the heart shape
        draw_pixel(x + 1, y - 1, string_color);
        draw_pixel(x + 1, y + 1, string_color);
        SDL_UpdateTexture(texture, NULL, color_buffer, (int)window_width * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
        SDL_Delay(10); // Adjust delay for heart drawing speed
    }

    for (int i = 360; i >= 180; i -= heart_speed)
    {
        float angle = i * 3.14159f / 180.0f;
        int x = heart_x + 160 * (sin(angle) * sin(angle) * sin(angle));                                       // Increase scale for x-coordinate
        int y = heart_y - 130 * cos(angle) + 50 * cos(2 * angle) + 20 * cos(3 * angle) + 10 * cos(4 * angle); // Increase scale for y-coordinate

        draw_pixel(x - 1, y + 1, string_color);
        draw_pixel(x - 1, y - 1, string_color);
        draw_pixel(x - 1, y, string_color);
        draw_pixel(x, y, string_color);
        draw_pixel(x + 1, y, string_color);
        draw_pixel(x + 1, y - 1, string_color);
        draw_pixel(x + 1, y + 1, string_color);

        SDL_UpdateTexture(texture, NULL, color_buffer, (int)window_width * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
        SDL_Delay(10); // Adjust delay for heart drawing speed
    }

    // Move the line to the left
    for (int x = mid_x; x <= end_x; x += line_growth_speed)
    {
        draw_line(mid_x, line_y, x, line_y, string_color);
        SDL_UpdateTexture(texture, NULL, color_buffer, (int)window_width * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
        SDL_Delay(10); // Adjust delay for line growth speed
    }
}

// Draw animated stairs using rectangles to create an ascending and descending stairs crossing paths effect
void draw_stairs(int step_width, int step_height)
{
    int x = 80;                             // Initial x-coordinate of the stairs
    int y_up = window_height - step_height; // Initial y-coordinate of the stairs going up
    int y_down = step_height * 3;           // Initial y-coordinate of the stairs going down

    // Define the delay between drawing each step
    int delay_ms = 500; // Adjust this value as needed

    while (x < window_width + step_width)
    {
        draw_rectangle(x, y_up, step_width, step_height, 0xFFFFFFFF); // Draw each step using draw_rectangle
        x += step_width;                                              // Move to the next step position
        y_up -= step_height;                                          // Move upwards for the next step

        // Update the screen after drawing each step
        SDL_UpdateTexture(texture, NULL, color_buffer, window_width * sizeof(uint32_t));
        SDL_RenderClear(renderer);
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);

        draw_rectangle(x - step_width - 80, y_down, step_width, step_height, 0xFFFFFFFF); // Draw each step using draw_rectangle, manipulate x to match with the stairs going up 
        y_down += step_height;                                                            // Move downwards for the next step

        // Update the screen after drawing each step
        SDL_UpdateTexture(texture, NULL, color_buffer, window_width * sizeof(uint32_t));
        SDL_RenderClear(renderer);
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);

        // Delay before drawing the next step
        SDL_Delay(400);
    }
}

// Project and render a pyramid shape representing a comet object with transformations and backface culling
void project_comit(vec3_t comet_vertices[N_MESH_VERTICES], face_t comet_faces[N_MESH_FACES], vec3_t scale_comet, vec3_t translate_comet, vec3_t rotation_comet)
{
    for (int i = 0; i < N_MESH_FACES; i++)
    {
        face_t commet_face = comet_faces[i];

        vec3_t face_vertices[3];
        face_vertices[0] = comet_vertices[commet_face.a - 1];
        face_vertices[1] = comet_vertices[commet_face.b - 1];
        face_vertices[2] = comet_vertices[commet_face.c - 1];

        mat4_t scale_matrix = mat4_make_scale(scale_comet.x, scale_comet.y, scale_comet.z);
        mat4_t rotation_matrix_x = mat4_make_rotation_x(rotation_comet.x); // pass the angle as float
        mat4_t rotation_matrix_y = mat4_make_rotation_y(rotation_comet.y);
        mat4_t rotation_matrix_z = mat4_make_rotation_z(rotation_comet.z);
        mat4_t translate_matrix = mat4_make_translate(translate_comet.x, translate_comet.y, translate_comet.z);

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
        // draw_rectangle(triangle.points[0].x, triangle.points[0].y, 2, 2, 0xFF00FF00);
        // draw_rectangle(triangle.points[1].x, triangle.points[1].y, 2, 2, 0xFF00FF00);
        // draw_rectangle(triangle.points[2].x, triangle.points[2].y, 2, 2, 0xFF00FF00);
        draw_triangle(
            triangle.points[0].x,
            triangle.points[0].y,
            triangle.points[1].x,
            triangle.points[1].y,
            triangle.points[2].x,
            triangle.points[2].y,
            0xFFFFFFFF);
    }
    t_cnt = 0;
}

int main()
{
    is_running = initialize_windowing_system();

    if (!is_running)
    {
        return -1;
    }

    float zoom_factor = 1.0f; // Initial zoom factor for hypnosis effect for Act 1 Scene 1
    float zoom_speed = 0.09f; // Zoom speed for Act 1 Scene 1
    int animation_duration_ms = 6000; // Animation duration in milliseconds for Act 1 Scene 1
    int start_time = SDL_GetTicks();  // Get the current time in milliseconds to keep track of the scenes and animations 

    float zoomin_factor = 1.0f; // Initial zoom factor for zooming in rectangles effect for Act 1 Scene 2
    float zoomin_speed = 0.005f; // Zoom speed for Act 1 Scene 2

    setup_memory_buffer();

    while (is_running)
    {
        process_keyboard_input();

        isEnable = true;
        clear_color_buffer(0xFF000000, isEnable); // Clear color buffer with black color

        // Act 1: Scene 1
        int numRectangles = 19; // Number of rectangles to draw
        for (int i = 0; i < numRectangles; ++i)
        {
            // Calculate rectangle properties based on zoom factor and index
            int rectWidth = 10 * zoom_factor + (i * 80);
            int rectHeight = 10 * zoom_factor + (i * 80);
            int x = (window_width - rectWidth) / 2 - (numRectangles * 2);   // Adjust position to center vanishing point
            int y = (window_height - rectHeight) / 2 - (numRectangles * 2); // Adjust position to center vanishing point

            // Draw hollow rectangle with calculated properties
            draw_hollow_rectangle(x, y, rectWidth, rectHeight, 0xFF00FF00);
            zoom_factor += zoom_speed; // Increase zoom factor
            SDL_Delay(5); // Delay for smooth animation

            // Reset zoom and speed for looping effect
            if (rectWidth > (window_width + (i * 5)) && rectHeight > (window_height + (i * 5))) 
            {
                if ((SDL_GetTicks() - start_time <= animation_duration_ms))
                {
                    zoom_factor = 1.0f;
                }
                else if ((SDL_GetTicks() - start_time <= animation_duration_ms))
                {
                    zoom_speed = 1.0f;
                }
            }
        }

        // Act 1: scene 2 and 3
        if ((SDL_GetTicks() - start_time >= 17000) && (SDL_GetTicks() - start_time <= 24000))
        {
            const int numRects = 25; // Number of rectangles for scene 2 and 3
            for (int i = 0; i < numRects; ++i)
            {
                // Randomly select a color from the array
                uint32_t color = colors[rand() % num_colors];
                int rectWidth = window_width * zoomin_factor + (i * 40);
                int rectHeight = 530 * zoomin_factor + (i * 40);
                int x = (window_width - rectWidth) / 2 - (numRects * 2);   // Adjust position to center vanishing point
                int y = (window_height - rectHeight) / 2 - (numRects * 2); // Adjust position to center vanishing point

                // Draw hollow rectangle with calculated properties
                draw_hollow_rectangle(x, y, rectWidth, rectHeight, color);
                zoomin_factor -= zoomin_speed; // Decrease zoom-in factor for zooming in effect 
            }

            if (zoomin_factor >= 2.0f) // Reset zoom for looping effect
            {
                zoomin_factor = 1.0f;
            }

            SDL_Delay(150); // Delay for smooth animation
        }

        isEnable = false; // Disable rendering for next frame
        clear_color_buffer(0xFF000000, isEnable); // Clear color buffer with black color

        // Act 2: Scene 1
        if ((SDL_GetTicks() - start_time >= 25000) && (SDL_GetTicks() - start_time <= 40000))
        {
            // Update color and transformation based on time
            uint32_t  color = 0xFFFFFFFF;
            if (SDL_GetTicks() - start_time <= 28000)
            {
                rotation.y += 0.15;
            }

            //Update scaling, rotation, and translation factors to control the animation as the scene plays 
            if ((SDL_GetTicks() - start_time >= 28000) && (SDL_GetTicks() - start_time <= 32000))
            {
                scale.x += 0.1;
                scale.y += 0.1;
                scale.z += 0.1;
                rotation.x += 0.1;
                rotation.z += 0.1;
                rotation.y += 0.05;
            }

            if ((SDL_GetTicks() - start_time >= 32000) && (SDL_GetTicks() - start_time <= 37000))
            {
                rotation.y += 0.1;
                translate.y += 0.01;
            }

            if ((SDL_GetTicks() - start_time >= 37000) && (SDL_GetTicks() - start_time <= 40000))
            {
                color = colors[rand() % num_colors];
                rotation.y += 0.3;
                translate.y -= 0.09;
            }
            project_pyramid(color); // Project pyramid with updated properties
        }

        // Act 2: Scene 2: Spinning triangle flowers 
        if ((SDL_GetTicks() - start_time >= 41000) && (SDL_GetTicks() - start_time <= 50000))
        {
            // Array of flower colors with different opacities
            uint32_t flower_colors[] = {
                0xFFFFC0CB, // Pink with 100% opacity
                0xE6FFC0CB, // Pink with 90% opacity
                0xCCFFC0CB, // Pink with 80% opacity
                0xB3FFC0CB, // Pink with 70% opacity
                0x99FFC0CB, // Pink with 60% opacity
                0x80FFC0CB, // Pink with 50% opacity
                0x66FFC0CB, // Pink with 40% opacity
                0x4DFFC0CB, // Pink with 30% opacity
                0x33FFC0CB, // Pink with 20% opacity
                0x1AFFC0CB  // Pink with 10% opacity
            };

            int num_flower_colors = sizeof(colors) / sizeof(colors[0]);

            // Act 2: Scene 2: Spinning triangle flowers transitions 
            // Draw each flower one by one
            u_int32_t color = 0xFFFFFFF;
            if (SDL_GetTicks() - start_time > 48000)
            {
                color = flower_colors[rand() % num_flower_colors];
            }

            if (SDL_GetTicks() - start_time >= 42000)
            {
                spinning_flower(0, color); // Draw first spinning flower
            }
            if (SDL_GetTicks() - start_time >= 45000)
            {
                spinning_flower(0.85, color); // Draw second spinning flower
            }
            if (SDL_GetTicks() - start_time >= 48000)
            {
                spinning_flower(2, color); // Draw third spinning flower
            }
        }

        // Act 2: Scene 3: bursting sparkles 
        if (SDL_GetTicks() - start_time >= 52000 && (SDL_GetTicks() - start_time <= 62000))
        {
            sparkle_bursts_animation(); //plays the bursting sparkles animation 
        }

        // Act 3: Scene 1: stairs animation 
        if (SDL_GetTicks() - start_time >= 63000 && (SDL_GetTicks() - start_time <= 68000))
        {
            int step_width = 80;
            int step_height = 28;
            draw_stairs(step_width, step_height); // Draw animated stairs
        }

        // Act 3: Scene 2: comet animation 
        if (SDL_GetTicks() - start_time >= 68000 && (SDL_GetTicks() - start_time <= 82000))
        {
            // Update properties for comet animation
            rotation_big_comet.z += 1.0;
            translate_big_comet.x-=0.05;
            translate_big_comet.y+=0.01;

            // Additional animations for smaller comet
            if (SDL_GetTicks() - start_time >= 68000 && (SDL_GetTicks() - start_time <= 74000)){
                rotation_small_comet.z += 1.0;
                translate_small_comet.x-=0.05;
                translate_small_comet.y+=0.01;
            }
            if(SDL_GetTicks() - start_time >= 74000 && (SDL_GetTicks() - start_time <= 82000)){
                rotation_small_comet.z += 1;
                translate_small_comet.x -=0.03;
                translate_small_comet.y+=0.05;
                translate_small_comet.z+=0.01;
            }
            // Project big comet
            project_comit(big_comet_vertices, big_comet_faces, scale_big_comet, translate_big_comet, rotation_big_comet); 
            // Project small comet
            project_comit(small_comet_vertices, small_comet_faces, scale_small_comet, translate_small_comet, rotation_small_comet);
        }

        // Act 3 Scene 3: Red string of fate going across the screen and forming heart 
        //Final scene based on time
        if (SDL_GetTicks() - start_time >= 82000 && (SDL_GetTicks() - start_time <= 90000))
        {
            red_string_of_fate_animation(); // Play red string of fate scene 
        }

        SDL_UpdateTexture(texture, NULL, color_buffer, (int)window_width * sizeof(uint32_t)); // Update SDL texture with color buffer
        SDL_RenderCopy(renderer, texture, NULL, NULL); // Render texture to window
        SDL_RenderPresent(renderer); // Present rendered frame
    }

    clean_up_windowing_system(); 
    return 0;
}