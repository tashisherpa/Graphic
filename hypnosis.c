#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <SDL2/SDL.h>
#include <math.h>

SDL_Renderer *renderer = NULL;
SDL_Window *window = NULL;
SDL_Texture *texture = NULL;
uint32_t *color_buffer = NULL;
bool is_running = false;

int window_width;
int window_height;

bool initialize_windowing_system();
void clean_up_windowing_system();
void process_keyboard_input();
void setup_memory_buffer();
void clear_color_buffer(uint32_t color);
void draw_hollow_rectangle(int x, int y, int width, int height, uint32_t color);
void draw_pixel(int x, int y, uint32_t color);
void update_state();
void run_render_pipeline();

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

void clear_color_buffer(uint32_t color)
{
    for (int i = 0; i < window_width * window_height; i++)
        color_buffer[i] = color;
}

void update_state()
{
    clear_color_buffer(0xFF000000);
    run_render_pipeline();
}

void run_render_pipeline()
{
    SDL_UpdateTexture(texture, NULL, color_buffer, (int)(window_width * sizeof(uint32_t)));
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_RenderPresent(renderer);
}

int main()
{
    is_running = initialize_windowing_system();

    if (!is_running)
    {
        return -1;
    }

    float zoom_factor = 1.0f;
    float zoom_speed = 0.1f;
    int animation_duration_ms = 6000; // Animation duration in milliseconds
    int start_time = SDL_GetTicks(); // Get the current time in milliseconds

    float zoomin_factor = 1.0f;
    float zoomin_speed = 0.005f;

    setup_memory_buffer();

    while (is_running)
    {
        process_keyboard_input();
        // update_state();

        clear_color_buffer(0xFF000000);

        int numRectangles = 30;
        for (int i = 0; i < numRectangles; ++i)
        {
            uint32_t color = colors[rand() % num_colors];
            int rectWidth = 10 * zoom_factor + (i*50);
            int rectHeight = 10 * zoom_factor + (i*50);
            int x = (window_width - rectWidth) / 2 ; // Adjust position to center vanishing point
            int y = (window_height - rectHeight) / 2 ; // Adjust position to center vanishing point
            draw_hollow_rectangle(x, y, rectWidth, rectHeight, color);
            zoom_factor += zoom_speed;
            SDL_Delay(5);
            if (rectWidth>(window_width+(i*3)) && rectHeight>(window_height+(i*3))) // Reset zoom for looping effect
            {
                if((SDL_GetTicks() - start_time <= animation_duration_ms)){
                    zoom_factor = 1.0f;
                }
            }
        }

        if((SDL_GetTicks() - start_time >= 15000)){
            const int numRects = 20;
            for (int i = 0; i < numRects; ++i)
            {
                // Randomly select a color from the array
                uint32_t color = colors[rand() % num_colors];
                int rectWidth = window_width * zoomin_factor + (i * 40);
                int rectHeight = 500 * zoomin_factor + (i * 40);
                int x = (window_width - rectWidth) / 2 - (numRects * 2);   // Adjust position to center vanishing point
                int y = (window_height - rectHeight) / 2 - (numRects * 2); // Adjust position to center vanishing point
                draw_hollow_rectangle(x, y, rectWidth, rectHeight, color);
                zoomin_factor -= zoomin_speed;
            }

            if (zoomin_factor >= 2.0f) // Reset zoom for looping effect
            {
                zoomin_factor = 1.0f;
            }

            SDL_Delay(100);
        }

        SDL_UpdateTexture(texture, NULL, color_buffer, (int)window_width * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
    }

    clean_up_windowing_system();
    return 0;
}