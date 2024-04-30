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
void project_pyramid(); // Declaration of the function

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

void draw_hollow_rectangle(int x, int y, int width, int height, uint32_t color)
{
    // Draw top and bottom edges
    for (int i = x; i < x + width; i++)
    {
        draw_pixel(i, y, color);
        draw_pixel(i, y + height - 1, color);
    }

    // Draw left and right edges
    for (int i = y; i < y + height; i++)
    {
        draw_pixel(x, i, color);
        draw_pixel(x + width - 1, i, color);
    }
}

void draw_pixel(int x, int y, uint32_t color)
{
    // Confirm pixel is in the visible window space
    if (x >= 0 && x < window_width && y >= 0 && y < window_height)
        color_buffer[(y * window_width) + x] = color;
}

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
        draw_pixel(x0, y0, 0xFFFFFFFF);
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
void draw_triangle(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t color)
{
    draw_line(x0, y0, x1, y1, color);
    draw_line(x1, y1, x2, y2, color);
    draw_line(x2, y2, x0, y0, color);
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

vec2_t perspective_projection_point(vec3_t point_3d)
{
    vec2_t projected_point = {.x = (scaling_factor * point_3d.x) / point_3d.z,
                              .y = (scaling_factor * point_3d.y) / point_3d.z};
    return projected_point;
}

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
    SDL_Event event;
    SDL_PollEvent(&event);
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
    color_buffer = (uint32_t *)malloc(window_width * window_height * sizeof(uint32_t));
    texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STREAMING, window_width, window_height);
}

void clear_color_buffer(uint32_t color, bool isEnable)
{
    if(isEnable){
        for (int i = 0; i < window_width * window_height; i++)
            color_buffer[i] = color;
    }
}

void update_state()
{
    clear_color_buffer(0xFF000000, true);
    run_render_pipeline();
}

void run_render_pipeline()
{
    SDL_UpdateTexture(texture, NULL, color_buffer, (int)(window_width * sizeof(uint32_t)));
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_RenderPresent(renderer);
}

void project_pyramid(void)
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
            0xFF800080);
    }
    t_cnt = 0;
}

void update_pyramid_state(int start_time)
{
    // Calculate the elapsed time since the start of the animation
    int current_time = SDL_GetTicks();
    int elapsed_time = current_time - start_time;

    // Determine the scaling factor based on elapsed time
    float scale_factor = 1.0f + elapsed_time / 1000.0f; // Adjust the scaling speed as needed

    // Update the scale of the pyramid only if the elapsed time is less than a certain duration
    if (elapsed_time < 4000)
    { // Change 5000 to the desired duration in milliseconds
        scale.x = scale.y = scale.z = scale_factor;
    }
    // Rotate the pyramid
    rotation.z += 0.1;

    clear_color_buffer(0xFF000000, true);
    // Project and render the pyramid
    project_pyramid();
}

void spinning_triangle()
{
    float rotation_angle = 0.0f; // Initial rotation angle
    int current_frame_time = SDL_GetTicks();
    float delta_time = (current_frame_time - previous_frame_time) / 2000.0f; // Convert milliseconds to seconds

    // Update the rotation angle based on elapsed time to make it rotate
    float rotation_speed = 0.5f; // Adjust rotation speed as needed
    rotation_angle += rotation_speed * delta_time;

    // Calculate the coordinates of the vertices for an equilateral triangle
    int triangle_size = 100;                                       // Size of the triangle
    int x_center = window_width / 2;                               // x-coordinate of the center of the screen
    int y_center = window_height / 2;                              // y-coordinate of the center of the screen
    int half_triangle_height = (int)(triangle_size * sqrt(3) / 2); // Half of the height of the equilateral triangle

    // Calculate the coordinates of the vertices
    int x0 = x_center;                            // x-coordinate of the top vertex
    int y0 = y_center - half_triangle_height;     // y-coordinate of the top vertex
    int x1 = x_center + triangle_size / 2;        // x-coordinate of the bottom-right vertex
    int y1 = y_center + half_triangle_height / 2; // y-coordinate of the bottom-right vertex
    int x2 = x_center - triangle_size / 2;        // x-coordinate of the bottom-left vertex
    int y2 = y_center + half_triangle_height / 2; // y-coordinate of the bottom-left vertex

    // Apply rotation transformation to the triangle vertices
    float cos_angle = cos(rotation_angle);
    float sin_angle = sin(rotation_angle);

    // Rotate each vertex of the triangle around the center of the screen
    int rotated_x0 = cos_angle * (x0 - x_center) - sin_angle * (y0 - y_center) + x_center;
    int rotated_y0 = sin_angle * (x0 - x_center) + cos_angle * (y0 - y_center) + y_center;

    int rotated_x1 = cos_angle * (x1 - x_center) - sin_angle * (y1 - y_center) + x_center;
    int rotated_y1 = sin_angle * (x1 - x_center) + cos_angle * (y1 - y_center) + y_center;

    int rotated_x2 = cos_angle * (x2 - x_center) - sin_angle * (y2 - y_center) + x_center;
    int rotated_y2 = sin_angle * (x2 - x_center) + cos_angle * (y2 - y_center) + y_center;

    // Draw the equilateral triangle with the rotated vertices
    draw_triangle(rotated_x0, rotated_y0, rotated_x1, rotated_y1, rotated_x2, rotated_y2, 0xFF800080);
}

int main()
{
    is_running = initialize_windowing_system();

    if (!is_running)
    {
        return -1;
    }

    float zoom_factor = 1.0f;
    float zoom_speed = 0.09f;
    int animation_duration_ms = 6000; // Animation duration in milliseconds
    int start_time = SDL_GetTicks();  // Get the current time in milliseconds

    float zoomin_factor = 1.0f;
    float zoomin_speed = 0.005f;

    // Array of colors
    uint32_t colors[] = {
        0xFF0000FF, // Blue
        0xFFFF0000, // Red
        0xFF00FF00, // Green
        0xFFFFFF00, // Yellow
        0xFFFF00FF, // Purple
        0xFF00FFFF, // Cyan
        0xFF800080, // Purple
        0xFFFFA500, // Orange
    };
    int num_colors = sizeof(colors) / sizeof(colors[0]);

    setup_memory_buffer();

    while (is_running)
    {
        process_keyboard_input();

        clear_color_buffer(0xFF000000, false);

        // Act 1: Scene 1
        // int numRectangles = 19;
        // for (int i = 0; i < numRectangles; ++i)
        // {
        //     int rectWidth = 10 * zoom_factor + (i * 80);
        //     int rectHeight = 10 * zoom_factor + (i * 80);
        //     int x = (window_width - rectWidth) / 2 - (numRectangles * 2);   // Adjust position to center vanishing point
        //     int y = (window_height - rectHeight) / 2 - (numRectangles * 2); // Adjust position to center vanishing point
        //     draw_hollow_rectangle(x, y, rectWidth, rectHeight, 0xFF00FF00);
        //     zoom_factor += zoom_speed;
        //     SDL_Delay(5);
        //     if (rectWidth > (window_width + (i * 5)) && rectHeight > (window_height + (i * 5))) // Reset zoom for looping effect
        //     {
        //         if ((SDL_GetTicks() - start_time <= animation_duration_ms))
        //         {
        //             zoom_factor = 1.0f;
        //         }
        //         else if ((SDL_GetTicks() - start_time <= animation_duration_ms))
        //         {
        //             zoom_speed = 1.0f;
        //         }
        //     }
        // }
        // // Act 1: scene 2 and 3
        // if ((SDL_GetTicks() - start_time >= 17000))
        // {
        //     const int numRects = 20;
        //     for (int i = 0; i < numRects; ++i)
        //     {
        //         // Randomly select a color from the array
        //         uint32_t color = colors[rand() % num_colors];
        //         int rectWidth = window_width * zoomin_factor + (i * 40);
        //         int rectHeight = 530 * zoomin_factor + (i * 40);
        //         int x = (window_width - rectWidth) / 2 - (numRects * 2);   // Adjust position to center vanishing point
        //         int y = (window_height - rectHeight) / 2 - (numRects * 2); // Adjust position to center vanishing point
        //         draw_hollow_rectangle(x, y, rectWidth, rectHeight, color);
        //         zoomin_factor -= zoomin_speed;
        //     }

        //     if (zoomin_factor >= 2.0f) // Reset zoom for looping effect
        //     {
        //         zoomin_factor = 1.0f;
        //     }

        //     SDL_Delay(100);
        // }

        // clear_color_buffer(0xFF000000, false);

        // Act 2: Scene 1
        // if ((SDL_GetTicks() - start_time >= 24000))
        // {
            // update_pyramid_state(start_time);

            // Act 2: Scene 2: Spinning triangle
            spinning_triangle();
        // }

        SDL_UpdateTexture(texture, NULL, color_buffer, (int)window_width * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
    }

    clean_up_windowing_system();
    return 0;
}
