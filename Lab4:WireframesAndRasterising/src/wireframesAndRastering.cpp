
#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <CanvasPoint.h>
#include <Colour.h>
#include <TextureMap.h>
#include <ModelTriangle.h>
#include <unordered_map>
#include <string>

using namespace std;
using namespace glm;

#define WIDTH 640
#define HEIGHT 400

vec3 camera(0.0,0.0,4.0);

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
    std::vector<float> zs = interpolateSingleFloats(from.depth, to.depth, numberOfValues);
    for(int i=0; i<numberOfValues; i++){
        points.push_back(CanvasPoint(round(xs[i]), round(ys[i]), (zs[i]))); //don't round z
    }
    return points;
}

void drawLine(CanvasPoint from, CanvasPoint to, Colour colour, vector<vector<float>> &depthBuffer, DrawingWindow &window){ //depth buffer ampersand (same for others)
    float numberOfSteps= fmax(fmax(abs(to.x - from.x), abs(to.y - from.y)), 1);
    vector<CanvasPoint> points = interpolateRoundFloats(from, to, numberOfSteps+1);

    for (float i=0.0; i<=numberOfSteps; i++){
        float pointDepth = 1/-points[i].depth;
        if(pointDepth > depthBuffer[round(points[i].y)][round(points[i].x)]){
            depthBuffer[round(points[i].y)][round(points[i].x)] = pointDepth;
            uint32_t c = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
            window.setPixelColour(points[i].x, points[i].y, c);
        }
        // else{cout << to <<endl;}
    }
}


void strokedTriangle(CanvasTriangle t, Colour colour, vector<vector<float>> &depthBuffer, DrawingWindow &window){
    CanvasPoint p_1 = t[0];
    CanvasPoint p_2 = t[1];
    CanvasPoint p_3 = t[2];

    drawLine(p_1,p_2,colour, depthBuffer, window);
    drawLine(p_2,p_3,colour, depthBuffer ,window);
    drawLine(p_3,p_1,colour, depthBuffer, window);
    
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

void filledTriangle(CanvasTriangle t, Colour colour,vector<vector<float>> &depthBuffer, DrawingWindow &window){
    sortTriangleVertices(t);

    std::vector<CanvasPoint> xStart = interpolateRoundFloats(t[0], t[1], t[1].y - t[0].y+1);
    if(t[2].y - t[1].y+1 >1){
        xStart.pop_back();
        std:: vector<CanvasPoint> start2 = interpolateRoundFloats(t[1], t[2], t[2].y - t[1].y+1); //t[1] not t[2]
        xStart.insert(xStart.end(), start2.begin(), start2.end());
    }
    std::vector<CanvasPoint> xEnd = interpolateRoundFloats(t[0], t[2], t[2].y - t[0].y+1);

    for(int i=0; i<=t[2].y - t[0].y; i++){
        drawLine(xStart[i],xEnd[i], colour, depthBuffer, window); 
    }
	strokedTriangle(t, colour, depthBuffer, window);
}


// void textureTriangle(CanvasTriangle t, TextureMap texture, DrawingWindow &window){
//     sortTriangleVertices(t);

//     std::vector<CanvasPoint> xStart = interpolateRoundFloats(t[0], t[1], t[1].y - t[0].y+1);
//     xStart.pop_back();
//     std::vector<CanvasPoint> xStart2 = interpolateRoundFloats(t[1], t[2], t[2].y - t[1].y+1);
//     xStart.insert(xStart.end(), xStart2.begin(), xStart2.end());
//     //begin() function is used to return an iterator pointing to the first element of the vector
//     std::vector<CanvasPoint> xEnd = interpolateRoundFloats(t[0], t[2], t[2].y - t[0].y+1);

//     std::vector<CanvasPoint> textureStart = interpolateRoundFloats(t[0].texturePoint, t[1].texturePoint, t[1].y - t[0].y+1);
//     xStart.pop_back();
//     std::vector<CanvasPoint> textureStart2 = interpolateRoundFloats(t[1].texturePoint, t[2].texturePoint, t[2].y - t[1].y+1);
//     textureStart.insert(textureStart.end(), textureStart2.begin(), textureStart2.end());
//     std::vector<CanvasPoint> textureEnd = interpolateRoundFloats(t[0].texturePoint, t[2].texturePoint, t[2].y - t[0].y+1);

//     for(int i=0; i<=t[2].y - t[0].y; i++){
//         float numberOfSteps = abs(xStart[i].x -xEnd[i].x);
//         std::vector<CanvasPoint> texpoint = interpolateRoundFloats(textureStart[i], textureEnd[i],numberOfSteps+1);
//         std::vector<CanvasPoint> points = interpolateRoundFloats(xStart[i], xEnd[i],numberOfSteps+1);

//         for(float i=0.0; i<=numberOfSteps; i++){
//             uint32_t col = texture.pixels[(round(texpoint[i].y)*texture.width) + round(texpoint[i].x)];
//             window.setPixelColour(points[i].x, points[i].y, col);
//         }  
//     }
// 	strokedTriangle(t, Colour(255,255,255), window);
// }

void random_triangle(DrawingWindow &window, int event){
    //task 2: make random traingle
    CanvasPoint t_1(rand()%(window.width-1), rand()%(window.height-1));
    CanvasPoint t_2(rand()%(window.width-1), rand()%(window.height-1));
    CanvasPoint t_3(rand()%(window.width-1), rand()%(window.height-1));
    CanvasTriangle t(t_1,t_2,t_3);
    Colour colour = Colour(rand()%256, rand()%256, rand()%256);
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 vertex, DrawingWindow &window){
    int planemultiplyer = 450;
    float focal = 2.0;

	// glm::vec3 vertex = vertexPosition - camera;
    float u = -round(planemultiplyer*focal * (vertex.x-camera.x)/(vertex.z-camera.z)) + (window.width/2);
    float v = round(planemultiplyer*focal * (vertex.y-camera.y)/(vertex.z-camera.z)) + (window.height/2);
    float z = vertex.z-camera.z;

    CanvasPoint t(u, v, z);
    return t;
}

void drawObject(vector<ModelTriangle> triangles, DrawingWindow &window) {
	window.clearPixels();
    vector<vector<float>> depthBuffer(window.height, vector<float>(window.width, 0)); 

    for(int i=0; i<triangles.size();i++){
        ModelTriangle triangle = triangles[i];
        CanvasTriangle ct;
        for(int i=0; i<triangle.vertices.size(); i++){
            vec3 vertex = triangle.vertices[i];
            ct.vertices[i] = getCanvasIntersectionPoint(vertex, window);
        }
        filledTriangle(ct, triangle.colour, depthBuffer, window);
        // drawLine(ct.vertices[0], CanvasPoint(ct.vertices[0].x+1,ct.vertices[0].y+1),Colour(255,255,255), window);
        // drawLine(ct.vertices[1], CanvasPoint(ct.vertices[1].x+1,ct.vertices[1].y+1),Colour(255,255,255), window);
        // drawLine(ct.vertices[2], CanvasPoint(ct.vertices[2].x+1,ct.vertices[2].y+1),Colour(255,255,255), window);

        // strokedTriangle(ct, triangle.colour, depthBuffer, window);
    }
}

unordered_map<string, Colour> mlt_parser(string filename) {
	unordered_map<string, Colour> colourPallete;
	string colour_name;

	ifstream File(filename);
	string line;

	while(getline(File, line)) {
		if(line == "") continue;

		vector<string> parts = split(line, ' ');
		if(parts[0] == "newmtl") {
			colour_name = parts[1];
		} else if(parts[0] == "Kd") {
			Colour colour(int(stof(parts[1])*255),int(stof(parts[2])*255),int(stof(parts[3])*255));
			colourPallete.insert({colour_name, colour});
		}
	}
	File.close();
	return colourPallete;
}

vector<ModelTriangle> obj_parser(string filename, float scale, unordered_map<string, Colour> colourPallete){
    vector<ModelTriangle>triangles;
    vector<vec3> vertices;
    string colour;

    ifstream File(filename);
    string line;

    while(getline(File, line)){ //gets each line from file and store it in str line
        if(line== "") continue;

        vector<string> parts = split(line, ' ');
        if(parts[0] == "v"){
            vec3 vertex(
                stof(parts[1])*scale, 
                stof(parts[2])*scale, 
                stof(parts[3])*scale
            );
            vertices.push_back(vertex);
        }
        else if (parts[0] == "f"){
            ModelTriangle t( //connecting the faces to the correct vertices + colour
                vertices[stoi(parts[1])-1],
                vertices[stoi(parts[2])-1], 
                vertices[stoi(parts[3])-1], 
                colourPallete[colour]);
			triangles.push_back(t);
        }
        else if(parts[0] == "usemtl"){
            colour = parts[1];
        }
    }
    File.close();
    return triangles;
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

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

    vector<ModelTriangle> triangles = obj_parser("models/cornell-box.obj", 0.17, mlt_parser("models/cornell-box.mtl"));

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);

        drawObject(triangles, window);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
