#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <CanvasPoint.h>
#include <Colour.h>
#include <TextureMap.h>

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

void drawLine(CanvasPoint from, CanvasPoint to, Colour colour, DrawingWindow &window){
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float numberOfSteps= std::max(abs(xDiff),abs(yDiff));

    float xStepSize= xDiff/numberOfSteps;
    float yStepSize = yDiff/numberOfSteps;

    //get colour
    uint32_t c = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);

    for (float i=0.0; i<numberOfSteps; i++){
        float x = from.x + (xStepSize*i);
        float y = from.y + (yStepSize*i);
        window.setPixelColour(round(x), round(y), c);
    }
}

void strokedTriangle(CanvasTriangle t, Colour colour, DrawingWindow &window){
    CanvasPoint p_1 = t[0];
    CanvasPoint p_2 = t[1];
    CanvasPoint p_3 = t[2];

    drawLine(p_1,p_2,colour, window);
    drawLine(p_2,p_3,colour, window);
    drawLine(p_3,p_1,colour, window);
    
}
void fill_half_triangle(CanvasTriangle t, Colour colour, DrawingWindow &window){
    CanvasPoint top = t[0];
    CanvasPoint mid = t[1];
    CanvasPoint bot = t[2];

	float x1_diff = mid.x - top.x;
	float y1_diff = mid.y - top.y;
	float x2_diff = bot.x - top.x;
	float y2_diff = bot.y - top.y;

	float steps1 = std::max(abs(x1_diff),abs(y1_diff));
	float steps2 = std::max(abs(x2_diff),abs(y2_diff));

	float x1_step = x1_diff/steps1;
	float x2_step = x2_diff/steps2;
	float y1_step = y1_diff/steps1;
	float y2_step = y2_diff/steps2;

	for(float i = 0.0; i < steps1; i++) {
		for(float j = 0.0; j < steps2; j++) {
			float x1 = top.x + (x1_step * i);
			float y1 = top.y + (y1_step * i);
			float x2 = top.x + (x2_step * j);
			float y2 = top.y + (y2_step * j);
       
			drawLine(CanvasPoint(round(x1),round(y1)),CanvasPoint(round(x2),round(y2)), colour, window);
		}
	}
}

void filledTriangle(CanvasTriangle t, Colour colour, DrawingWindow &window){
    fill_half_triangle(t, colour, window);
    strokedTriangle(t, Colour(255,255,255), window);
}

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
    float interval = (to - from)/(numberOfValues-1);
    std::vector<float> result(numberOfValues);
    for(size_t i=0 ; i< result.size(); i++){
        result[i] = from+(interval*i);
    }
    return result;
}

template <typename T>
std::vector<CanvasPoint> interpolateRoundFloats(T from, T to, int numberOfValues){
    std::vector<CanvasPoint> points;
    std::vector<float> xs = interpolateSingleFloats(from.x, to.x, numberOfValues);
    std::vector<float> ys = interpolateSingleFloats(from.y, to.y, numberOfValues);
    for(int i=0; i<numberOfValues; i++){
        points.push_back(CanvasPoint(round(xs[i]), round(ys[i])));
    }
    return points;
}

void sortTriangleVertices(CanvasTriangle &triangle) {
	if (triangle[0].y > triangle[1].y) {
		std::swap(triangle[0], triangle[1]);
	}
	if (triangle[1].y > triangle[2].y) {
		std::swap(triangle[1], triangle[2]);
		if (triangle[0].y > triangle[1].y) {
			std::swap(triangle[0], triangle[1]);
		}	
	}
}

void textureTriangle(CanvasTriangle t, TextureMap texture, DrawingWindow &window){
    sortTriangleVertices(t);

    std::vector<CanvasPoint> xStart = interpolateRoundFloats(t[0], t[1], t[1].y - t[0].y+1);
    xStart.pop_back();
    std::vector<CanvasPoint> xStart2 = interpolateRoundFloats(t[1], t[2], t[2].y - t[1].y+1);
    xStart.insert(xStart.end(), xStart2.begin(), xStart2.end());
    //begin() function is used to return an iterator pointing to the first element of the vector
    std::vector<CanvasPoint> xEnd = interpolateRoundFloats(t[0], t[2], t[2].y - t[0].y+1);


    std::vector<CanvasPoint> textureStart = interpolateRoundFloats(t[0].texturePoint, t[1].texturePoint, t[1].y - t[0].y+1);
    xStart.pop_back();
    std::vector<CanvasPoint> textureStart2 = interpolateRoundFloats(t[1].texturePoint, t[2].texturePoint, t[2].y - t[1].y+1);
    textureStart.insert(textureStart.end(), textureStart2.begin(), textureStart2.end());
    std::vector<CanvasPoint> textureEnd = interpolateRoundFloats(t[0].texturePoint, t[2].texturePoint, t[2].y - t[0].y+1);

    for(int i=0; i<=t[2].y - t[0].y; i++){
        float numberOfSteps = abs(xStart[i].x -xEnd[i].x);
        std::vector<CanvasPoint> texpoint = interpolateRoundFloats(textureStart[i], textureEnd[i],numberOfSteps+1);
        std::vector<CanvasPoint> points = interpolateRoundFloats(xStart[i], xEnd[i],numberOfSteps+1);

        for(float i=0.0; i<=numberOfSteps; i++){
            uint32_t col = texture.pixels[(round(texpoint[i].y)*texture.width) + round(texpoint[i].x)];
            window.setPixelColour(points[i].x, points[i].y, col);
        }  
    }
	strokedTriangle(t, Colour(255,255,255), window);

}

void random_triangle(DrawingWindow &window, int event){
    //task 2: make random traingle
    CanvasPoint t_1(rand()%(window.width-1), rand()%(window.height-1));
    CanvasPoint t_2(rand()%(window.width-1), rand()%(window.height-1));
    CanvasPoint t_3(rand()%(window.width-1), rand()%(window.height-1));
    CanvasTriangle t(t_1,t_2,t_3);
    Colour colour = Colour(rand()%256, rand()%256, rand()%256);

    if(event == SDLK_u) {
		strokedTriangle(CanvasTriangle(t_1,t_2,t_3), colour, window);
	}
	if(event == SDLK_f) {
		filledTriangle(CanvasTriangle(t_1,t_2,t_3), colour, window);
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
        else if (event.key.keysym.sym == SDLK_u || event.key.keysym.sym == SDLK_f) random_triangle(window, event.key.keysym.sym);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

    TextureMap texture("src/texture.ppm");
    CanvasPoint p_1(160,10);
    CanvasPoint p_2(300,230);
    CanvasPoint p_3(10,150);
    p_1.texturePoint = TexturePoint(195,5);
    p_2.texturePoint = TexturePoint(395,380);
    p_3.texturePoint = TexturePoint(65,330);
    CanvasTriangle t(p_1,p_2,p_3);


	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);

        // drawLine(p_1, p_2, Colour(255,255,255), window);
        textureTriangle(t, texture, window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
