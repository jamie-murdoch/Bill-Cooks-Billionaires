#ifndef GRAPHICS_H
#define GRAPHICS_H
#include <SDL2/SDL.h>
#include <SDL2/SDL_test_font.h>
#include <vector>
#include <stdlib.h> 

#include "algebra.h"
#include "graph.h"

using namespace std;

SDL_Renderer* renderer;
int width = 700, height = 700;
double minOffset, scale;

inline void draw_circle(int n_cx, int n_cy, int radius)
{
    // if the first pixel in the screen is represented by (0,0) (which is in sdl)
    // remember that the beginning of the circle is not in the middle of the pixel
    // but to the left-top from it:
 
    double error = (double)-radius;
    double x = (double)radius -0.5;
    double y = (double)0.5;
    double cx = n_cx - 0.5;
    double cy = n_cy - 0.5;
 
    while (x >= y)
    {
        SDL_RenderDrawPoint(renderer, (int)(cx + x), (int)(cy + y));
        SDL_RenderDrawPoint(renderer, (int)(cx + y), (int)(cy + x));
 
        if (x != 0)
        {
            SDL_RenderDrawPoint(renderer, (int)(cx - x), (int)(cy + y));
            SDL_RenderDrawPoint(renderer, (int)(cx + y), (int)(cy - x));
        }
 
        if (y != 0)
        {
            SDL_RenderDrawPoint(renderer, (int)(cx + x), (int)(cy - y));
            SDL_RenderDrawPoint(renderer, (int)(cx - y), (int)(cy + x));
        }
 
        if (x != 0 && y != 0)
        {
            SDL_RenderDrawPoint(renderer, (int)(cx - x), (int)(cy - y));
            SDL_RenderDrawPoint(renderer, (int)(cx - y), (int)(cy - x));
        }
 
        error += y;
        ++y;
        error += y;
 
        if (error >= 0)
        {
            --x;
            error -= x;
            error -= x;
        }
    }
}

inline void draw_fill_circle(int cx, int cy, int radius)
{
    // Note that there is more to altering the bitrate of this 
    // method than just changing this value.  See how pixels are
    // altered at the following web page for tips:
    //   http://www.libsdl.org/intro.en/usingvideo.html
    //static const int BPP = 4;
 
    double r = (double)radius;
 
    for (double dy = -1.0; dy <= r; dy += 1.0)
    {
        // This loop is unrolled a bit, only iterating through half of the
        // height of the circle.  The result is used to draw a scan line and
        // its mirror image below it.
 
        // The following formula has been simplified from our original.  We
        // are using half of the width of the circle because we are provided
        // with a center and we need left/right coordinates.
 
        double dx = nint(sqrt((2.0 * r * dy) - (dy * dy)));
        int x = cx - dx;
 
        // Grab a pointer to the left-most pixel for each half of the circle
        // Uint8 *target_pixel_a = (Uint8 *)surface->pixels + ((int)(cy + r - dy)) * surface->pitch + x * BPP;
        // Uint8 *target_pixel_b = (Uint8 *)surface->pixels + ((int)(cy - r + dy)) * surface->pitch + x * BPP;
 
        for (; x <= cx + dx; x++)
        {
        	SDL_RenderDrawPoint(renderer, x, (int)(cy + r - dy));
        	SDL_RenderDrawPoint(renderer, x, (int)(cy - r + dy));
            // *(Uint32 *)target_pixel_a = pixel;
            // *(Uint32 *)target_pixel_b = pixel;
            // target_pixel_a += BPP;
            // target_pixel_b += BPP;
        }
    }
}

inline void setup_sdl(const Graph &graph) {
    SDL_Window *window;                    // Declare a pointer
    SDL_Init(SDL_INIT_VIDEO);              // Initialize SDL2
    // Create an application window with the following settings:
    window = SDL_CreateWindow(
        "TSP Edge Elimination",                  // window title
        0,           // initial x position
        0,           // initial y position
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
    minOffset = min(minX, minY);
    scale = max(maxY, maxX) - minOffset;
}

inline int to_screen(double v) {
	return nint((v - minOffset)*width/scale);
}
inline SDL_Point to_sdl_point(const Point2D &point) {
	SDL_Point sdl_p = {to_screen(point.x()), to_screen(point.y())};
	return sdl_p; 
}


inline void draw_points(const vector<Point2D> &points, int r, int g, int b, bool use_labels = true) {
	int radius = 3;
	char label[100];
	SDL_SetRenderDrawColor(renderer, r, g, b, 255);
	for(int i = 0; i < (int)points.size(); i++) {
		SDL_Point p = to_sdl_point(points[i]);
		draw_circle(p.x, p.y, radius);

	}

	if(use_labels) {
		for(int i = 0; i < (int)points.size(); i++) {
			SDL_Point p = to_sdl_point(points[i]);
			sprintf(label,"%d",i);
			SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
			SDLTest_DrawString(renderer, p.x + radius, p.y - radius*4, label);
		}
	}
}

inline void draw_tree(KdNode *node, int minX = 0, int minY = 0, int maxX = width, int maxY = height) {
	if(!node->loson) return;

 	//x draw vertical
 	if(node->cutdim == 0) {
 		int x = to_screen(node->cutval);

 		SDL_RenderDrawLine(renderer, x, minY, x, maxY);

 		draw_tree(node->loson, minX, minY, x, maxY);
 		draw_tree(node->hison, x, minY, maxX, maxY);
 	}

 	//y horiz
 	if(node->cutdim == 1) {
 		int y = to_screen(node->cutval);

 		SDL_RenderDrawLine(renderer, minX, y, maxX, y);

 		draw_tree(node->loson, minX, minY, maxX, y);
 		draw_tree(node->hison, minX, y, maxX, maxY);
 	}
}

#endif





