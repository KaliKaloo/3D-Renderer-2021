#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

#define WIDTH 320
#define HEIGHT 240

void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

//Task 2: Single Element Numerical Interpolation
std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
    float interval = (to - from)/(numberOfValues-1);
    std::vector<float> result(numberOfValues);
    for(size_t i=0 ; i< result.size(); i++){
        result[i] = from+(interval*i);
    }
    return result;
}

//Task 4: Three Element Numerical Interpolation
std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
    //creates vector of size 3. For each corresponding value from[i] and to[i], calculate the intervel. In total there should be three.
    std::vector<float> interval(3); 
    for (size_t i=0; i<interval.size(); i++){
        interval[i] = (to[i] - from[i])/(numberOfValues-1);
    }

    std::vector<glm::vec3> result(numberOfValues);
    for(size_t i=0; i<result.size(); i++){
        result[i] = glm::vec3(from[0]+(interval[0]*i), from[1]+(interval[1]*i),from[2]+(interval[2]*i));
    }

    return result;
}

//Task 3: Single Dimension Greyscale Interpolation
void gray_interpolate(DrawingWindow &window) {
	window.clearPixels();
    std::vector<float> result = interpolateSingleFloats(0, 255, window.width);

	for (size_t x = 0; x < window.width; x++) {
        float grayPixel = result[x];
		uint32_t colour = (255 << 24) + (int(grayPixel) << 16) + (int(grayPixel) << 8) + int(grayPixel);
		for (size_t y = 0; y < window.height; y++) {
			window.setPixelColour(x, y, colour); //width x height
		}
	}
}

//Task 5: Two Dimensional Colour Interpolation
void rgb_interpolate(DrawingWindow &window) {
	window.clearPixels();
    glm::vec3 topLeft(255, 0, 0);        // red 
    glm::vec3 topRight(0, 0, 255);       // blue 
    glm::vec3 bottomRight(0, 255, 0);    // green 
    glm::vec3 bottomLeft(255, 255, 0);   // yellow

    std::vector<glm::vec3> leftValues = interpolateThreeElementValues(topLeft, bottomLeft, window.width);
    std::vector<glm::vec3> rightValues = interpolateThreeElementValues(topRight, bottomRight, window.width);

	for (size_t x = 0; x < window.height; x++) {
        //interpolate the pixel values from one end to the other as you go down the window
        std::vector<glm::vec3> rowValues = interpolateThreeElementValues(leftValues[x], rightValues[x], window.width);
        
		for (size_t y = 0; y < window.width; y++) {
            uint32_t colour = (255 << 24) + (int(rowValues[y].x) << 16) + (int(rowValues[y].y) << 8) + int(rowValues[y].z);
            //.x .y .z is a away of accessing vectors by array indexing
			window.setPixelColour(y, x, colour); //width x height
		}
	}
}

int main(int argc, char *argv[]) {
     //TASK 1
    // std::vector<float> result;
    // result = interpolateSingleFloats(2.2, 8.5, 7);
    // for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
    // std::cout << std::endl;

    //TASK 4
    // std::vector<glm::vec3> result;
    // glm::vec3 from(1, 4, 9.2);
    // glm::vec3 to(4, 1, 9.8);
    // result = interpolateThreeElementValues(from, to, 4);
    // for(size_t i=0; i<result.size(); i++) {
    //     std::cout << glm::to_string(result[i]) << " ";}
    // std::cout << std::endl;

	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);

		// draw(window);
        // gray_interpolate(window); //TASK 2
        rgb_interpolate(window); //TASk 5
        
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
