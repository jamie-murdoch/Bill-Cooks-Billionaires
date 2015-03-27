#ifndef GRAPHICS_H
#define GRAPHICS_H
#include <SDL2/SDL.h>
#include <vector>

#include "algebra.h"
#include "graph.h"

using namespace std;

SDL_Renderer* renderer;
int width = 800, height = 800;
int minX, minY;
inline void setup_sdl(const Graph &graph) {
    SDL_Window *window;                    // Declare a pointer
    SDL_Init(SDL_INIT_VIDEO);              // Initialize SDL2
    // Create an application window with the following settings:
    window = SDL_CreateWindow(
        "TSP Edge Elimination",                  // window title
        SDL_WINDOWPOS_UNDEFINED,           // initial x position
        SDL_WINDOWPOS_UNDEFINED,           // initial y position
        width,                               // width, in pixels
        height,                               // height, in pixels
        SDL_WINDOW_OPENGL                  // flags - see below
    );
    // Check that the window was successfully made
    if (window == NULL) {
        // In the event that the window could not be made...
        printf("Could not create window: %s\n", SDL_GetError());
    }

    renderer = SDL_CreateRenderer(window, -1, 0);

    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderClear(renderer);
    SDL_RenderPresent(renderer);

    double minX, minY, maxX, maxY;
    graph.get_bounding_box(minX, minY, maxX, maxY);

}

inline SDL_Point to_sdl_point(const Point2D &point) {
	SDL_Point sdl_p = {nint(point.x()), nint(point.y())};
	return sdl_p; 
}

inline void draw_points(const vector<Point2D> &points, int r, int g, int b) {
	SDL_SetRenderDrawColor(renderer, r, g, b, 255);
	vector<SDL_Point> sdl_points(points.size());
	for(int i = 0; i < (int)points.size(); i++) {
		sdl_points[i] = to_sdl_point(points[i]);
		cout << points[i] << endl;
	}

	SDL_RenderDrawPoints(renderer, &sdl_points[0], sdl_points.size());
	SDL_RenderPresent(renderer);
}

#endif