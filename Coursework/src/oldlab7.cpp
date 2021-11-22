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
#include <RayTriangleIntersection.h>
#include <functional>

using namespace std;
using namespace glm;

#define WIDTH 640
#define HEIGHT 400
#define pi 3.14159265359

vec3 camera(0.0,0.0,4.0);
mat3 cameraOrientation(
	vec3(1.0,0.0,0.0),
	vec3(0.0,1.0,0.0),
	vec3(0.0,0.0,1.0)
);
vec3 light(0.0, 1.0, 0.0);
int renderMode = 2; //initially 2 = rasterised mode
bool orbiting = false;


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

vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
    float interval = (to - from)/(numberOfValues-1);
    vector<float> result(numberOfValues);
    for(size_t i=0 ; i< result.size(); i++){
        result[i] = from+(interval*i);
    }
    return result;
}

vector<CanvasPoint> interpolateRoundFloats(CanvasPoint from, CanvasPoint to, int numberOfValues){
    vector<CanvasPoint> points;
    vector<float> xs = interpolateSingleFloats(from.x, to.x, numberOfValues);
    vector<float> ys = interpolateSingleFloats(from.y, to.y, numberOfValues);
    vector<float> zs = interpolateSingleFloats(from.depth, to.depth, numberOfValues);
    for(int i=0; i<numberOfValues; i++){
        points.push_back(CanvasPoint(round(xs[i]), round(ys[i]), (zs[i])));
    }
    return points;
}

vector<CanvasPoint> interpolateRoundFloatsT(TexturePoint from, TexturePoint to, int numberOfValues){
    vector<CanvasPoint> points;
    vector<float> xs = interpolateSingleFloats(from.x, to.x, numberOfValues);
    vector<float> ys = interpolateSingleFloats(from.y, to.y, numberOfValues);
    for(int i=0; i<numberOfValues; i++){
        points.push_back(CanvasPoint(round(xs[i]), round(ys[i]))); 
    }
    return points;
}

void drawLine(CanvasPoint from, CanvasPoint to, Colour colour, vector<vector<float>> &depthBuffer, DrawingWindow &window){ //depth buffer ampersand (same for others)
    float numberOfSteps= fmax(fmax(abs(to.x - from.x), abs(to.y - from.y)), 1);
    vector<CanvasPoint> points = interpolateRoundFloats(from, to, numberOfSteps+1);

    for (float i=0.0; i<=numberOfSteps; i++){
        if(points[i].x >= 0 && points[i].x < window.width && points[i].y >= 0 && points[i].y < window.height) {
            float pointDepth = 1/-points[i].depth;
            if(pointDepth > depthBuffer[round(points[i].y)][round(points[i].x)]){
                depthBuffer[round(points[i].y)][round(points[i].x)] = pointDepth;
                uint32_t c = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
                window.setPixelColour(points[i].x, points[i].y, c);
            }
        }
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
		swap(triangle[0], triangle[1]);
	}
	if (triangle[1].y > triangle[2].y) {
		swap(triangle[1], triangle[2]);
		if (triangle[0].y > triangle[1].y) {
		    swap(triangle[0], triangle[1]);
		}	
	}
}

void filledTriangle(CanvasTriangle t, Colour colour,vector<vector<float>> &depthBuffer, DrawingWindow &window){
    sortTriangleVertices(t);

    vector<CanvasPoint> xStart = interpolateRoundFloats(t[0], t[1], t[1].y - t[0].y+1);
    if(t[2].y - t[1].y+1 >1){
        xStart.pop_back();
        std:: vector<CanvasPoint> start2 = interpolateRoundFloats(t[1], t[2], t[2].y - t[1].y+1); //t[1] not t[2]
        xStart.insert(xStart.end(), start2.begin(), start2.end());
    }
    vector<CanvasPoint> xEnd = interpolateRoundFloats(t[0], t[2], t[2].y - t[0].y+1);

    for(int i=0; i<=t[2].y - t[0].y; i++){
        drawLine(xStart[i],xEnd[i], colour, depthBuffer, window); 
    }
	strokedTriangle(t, colour, depthBuffer, window);
}


void textureTriangle(CanvasTriangle t, TextureMap texture, vector<vector<float>> &depthBuffer, DrawingWindow &window){
    sortTriangleVertices(t);

    vector<CanvasPoint> xStart = interpolateRoundFloats(t[0], t[1], t[1].y - t[0].y+1);
    xStart.pop_back();
    vector<CanvasPoint> xStart2 = interpolateRoundFloats(t[1], t[2], t[2].y - t[1].y+1);
    xStart.insert(xStart.end(), xStart2.begin(), xStart2.end());
    //begin() function is used to return an iterator pointing to the first element of the vector
    vector<CanvasPoint> xEnd = interpolateRoundFloats(t[0], t[2], t[2].y - t[0].y+1);

    vector<CanvasPoint> textureStart = interpolateRoundFloatsT(t[0].texturePoint, t[1].texturePoint, t[1].y - t[0].y+1);
    xStart.pop_back();
    vector<CanvasPoint> textureStart2 = interpolateRoundFloatsT(t[1].texturePoint, t[2].texturePoint, t[2].y - t[1].y+1);
    textureStart.insert(textureStart.end(), textureStart2.begin(), textureStart2.end());
    vector<CanvasPoint> textureEnd = interpolateRoundFloatsT(t[0].texturePoint, t[2].texturePoint, t[2].y - t[0].y+1);

    for(int i=0; i<=t[2].y - t[0].y; i++){
        float numberOfSteps = abs(xStart[i].x -xEnd[i].x);
        vector<CanvasPoint> texpoint = interpolateRoundFloats(textureStart[i], textureEnd[i],numberOfSteps+1);
        vector<CanvasPoint> points = interpolateRoundFloats(xStart[i], xEnd[i],numberOfSteps+1);

        for(float i=0.0; i<=numberOfSteps; i++){
            uint32_t col = texture.pixels[(round(texpoint[i].y)*texture.width) + round(texpoint[i].x)];
            window.setPixelColour(points[i].x, points[i].y, col);
        }  
    }
	strokedTriangle(t, Colour(255,255,255), depthBuffer, window);
}

void random_triangle(DrawingWindow &window, int event){
    //task 2: make random traingle
    CanvasPoint t_1(rand()%(window.width-1), rand()%(window.height-1));
    CanvasPoint t_2(rand()%(window.width-1), rand()%(window.height-1));
    CanvasPoint t_3(rand()%(window.width-1), rand()%(window.height-1));
    CanvasTriangle t(t_1,t_2,t_3);
    Colour colour = Colour(rand()%256, rand()%256, rand()%256);
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 vertexPosition, DrawingWindow &window){
    int planemultiplyer = 150;
    float focal = 2.0;
	vec3 vertex = (vertexPosition - camera)*cameraOrientation;

    float u = -round(planemultiplyer*focal * (vertex.x)/(vertex.z)) + (window.width/2);
    float v = round(planemultiplyer*focal * (vertex.y)/(vertex.z)) + (window.height/2);
    float z = vertex.z;

    CanvasPoint t(u, v, z);
    return t;
}

// DRAW WIREFRAME SCENE
void draw_wireframe(vector<ModelTriangle> triangles, float focal, int planeMultiplyer, DrawingWindow &window) {
	window.clearPixels();
    vector<vector<float>> depthBuffer(window.height, vector<float>(window.width, 0)); 

    for(int i=0; i<triangles.size();i++){
        ModelTriangle triangle = triangles[i];
        CanvasTriangle ct;
        for(int i=0; i<triangle.vertices.size(); i++){
            vec3 vertex = triangle.vertices[i];
            ct.vertices[i] = getCanvasIntersectionPoint(vertex, window);
        }
        strokedTriangle(ct, triangle.colour, depthBuffer, window);
    }
}

// DRAW RASTERISED SCENE
void draw_rasterise(vector<ModelTriangle> triangles, float focal, int planeMultiplyer, DrawingWindow &window) {
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
    }
}

bool inShadow(RayTriangleIntersection intersect, vector<ModelTriangle> faces) {
    vec3 shadowRay = light - intersect.intersectionPoint;

    for(int i = 0; i<faces.size(); i++){
        ModelTriangle triangle = faces[i];
        
        vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        vec3 SPVector = intersect.intersectionPoint - triangle.vertices[0];
        mat3 DEMatrix(-normalize(shadowRay), e0, e1);
        vec3 possibleSolution = inverse(DEMatrix) * SPVector;

        float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;

        if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
			if(t<length(shadowRay) && t> 0.01f && i != intersect.triangleIndex){
                return true;
            }
	    }
    } 
    return false;
}

void getPixelBrightness (){

}


RayTriangleIntersection getClosestIntersection(vec3 rayDirection, vector<ModelTriangle> faces) {
    RayTriangleIntersection intersection;
    intersection.distanceFromCamera = numeric_limits<float>::infinity();

    for(int i = 0; i<faces.size(); i++){
        ModelTriangle triangle = faces[i];
        
        vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        vec3 SPVector = camera - triangle.vertices[0];
        mat3 DEMatrix(-rayDirection, e0, e1);
        vec3 possibleSolution = inverse(DEMatrix) * SPVector;

        float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;

        if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
			if(intersection.distanceFromCamera > t && t > 0) {
				intersection.distanceFromCamera = t;
				intersection.intersectedTriangle = triangle;
				intersection.triangleIndex = i;

				vec3 intersect = triangle.vertices[0]+u*e0+v*e1;
				intersection.intersectionPoint = intersect;
			}
	    }
    } 
   return intersection;
}

// DRAW RAYTRACED SCENE
void draw_raytrace(vector<ModelTriangle> triangles, int focal, int PlaneMultiplyer, DrawingWindow &window){
	window.clearPixels();
    for (int x=0; x < window.width; x++){
        for(int y=0; y < window.height; y++){
            
            vec3 direction(-(float(window.width / 2) - x) / PlaneMultiplyer, (float(window.height / 2) - y) / PlaneMultiplyer, -2);
			glm::vec3 camDir = glm::normalize(cameraOrientation * direction);
			auto intersect = getClosestIntersection( camDir, triangles);

			if(!isinf(intersect.distanceFromCamera)){
                Colour colour = intersect.intersectedTriangle.colour;
				
                if (inShadow(intersect, triangles)==true){
                    uint32_t c = (255 << 24) + (int(colour.red/3) << 16) + (int(colour.green/3) << 8) + int(colour.blue/3);
                     window.setPixelColour(x, y, c);
                }else{
                    uint32_t c = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
                     window.setPixelColour(x, y, c);
                }
            }
        }
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
 vector<ModelTriangle> obj_parser(const std::string &filename, float scale, unordered_map<string, Colour> colourPallete) {
	cout << "[LOADING OBJ: " << filename << "]" << std::endl;
	vector<glm::vec3> vertices;
	vector<TexturePoint> textureVertices;
	vector<glm::vec3> normals;
	vector<ModelTriangle> faces;
    string colour;


	ifstream File(filename, ifstream::in);
	string nextLine;

	while (getline(File, nextLine)) {
		auto vector = split(nextLine, ' ');
		if (vector[0] == "usemtl") {
            colour = vector[1];
        }
		else if (vector[0] == "v") {
			vertices.push_back(vec3(
				stof(vector[1]) * scale,
				stof(vector[2]) * scale,
				stof(vector[3]) * scale
			));
		}
		else if (vector[0] == "vt") {
			textureVertices.push_back(TexturePoint(
				stof(vector[1]),
				stof(vector[2])
			));
		}
		else if (vector[0] == "vn") {
			normals.push_back(vec3(
				stof(vector[1]),
			    stof(vector[2]),
				stof(vector[3])
			));
		}
		else if (vector[0] == "f") {
			auto triangle = ModelTriangle();
			for (int i=0; i < 3; i++) {
				auto v = split(vector[i + 1], '/');
				triangle.vertices[i] = vertices[stoi(v[0]) - 1];
				if (v.size() > 1 && v[1] != "") triangle.texturePoints[i] = textureVertices[std::stoi(v[1]) - 1];
				if (v.size() > 2 && v[2] != "") triangle.vertexNormals[i] = normals[stoi(v[2]) - 1];
                
			}
			triangle.normal = normalize(cross(vec3(triangle.vertices[1] - triangle.vertices[0]), vec3(triangle.vertices[2] - triangle.vertices[0])));

            triangle.colour = colourPallete[colour];

			faces.push_back(triangle);
		}
	}
	cout << "[FINISHED LOADING OBJ]" << endl;
    File.close();
	return faces;
}

void look_at() {
	vec3 forward = normalize(camera - vec3(0.0,0.0,0.0));
	vec3 right = normalize(cross(vec3(0.0,1.0,0.0), forward));
	vec3 up = normalize(cross(forward, right));

	cameraOrientation[0] = right;
	cameraOrientation[1] = up;
	cameraOrientation[2] = forward;
}

mat3 rotation_x(float rads) {
	return {
        1.0f, 0.0f, 0.0f,
        0.0f, cos(rads), sin(rads), 
        0.0f, -sin(rads),  cos(rads), 
    };
}
mat3 rotation_y(float rads) {
	return {
        cos(rads), 0.0f, -sin(rads),
        0.0f, 1.0f, 0.0f,
        sin(rads), 0.0f, cos(rads), 
    };
}
mat3 rotation_z(float rads) {
	return {
        cos(rads), sin(rads), 0.0f, 
        -sin(rads), cos(rads), 0.0f, 
        0.0f, 0.0f, 1.0f, 
    };
}

void orbit(bool orb) {
	if(orb) {
        int rad = 180;
        if (renderMode == 1) rad = 50;
        if (renderMode ==3) rad = 400;
		camera = camera * rotation_y(-pi/rad);
		look_at();
	}
}

function<void(vector<ModelTriangle>,float, int , DrawingWindow &)> RenderMode = draw_rasterise;

void handleEvent(SDL_Event event, DrawingWindow &window) {
    float rotVal = pi/60;

	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_s) camera.y -= 0.1;
		else if (event.key.keysym.sym == SDLK_w) camera.y += 0.1;
		else if (event.key.keysym.sym == SDLK_a) camera.x -= 0.1;
		else if (event.key.keysym.sym == SDLK_d) camera.x += 0.1;
		else if (event.key.keysym.sym == SDLK_z) camera.z += 0.1;
        else if (event.key.keysym.sym == SDLK_x) camera.z -= 0.1;
        
        else if (event.key.keysym.sym == SDLK_RIGHT){
            camera = camera * rotation_y( -rotVal); 
            look_at();
        }
		else if (event.key.keysym.sym == SDLK_LEFT){
            camera = camera * rotation_y(rotVal);
            look_at();
        }
		else if (event.key.keysym.sym == SDLK_DOWN){
            camera = camera * rotation_x( -rotVal);
            look_at();
        }
		else if (event.key.keysym.sym == SDLK_UP) {
            camera = camera * rotation_x(rotVal);
            look_at();
        }
        else if (event.key.keysym.sym == SDLK_o) orbiting = (orbiting) ? false : true;
        else if (event.key.keysym.sym == SDLK_l) look_at();
        else if (event.key.keysym.sym == SDLK_1) { RenderMode = draw_raytrace; renderMode = 1; cout << "[Render Mode]: raytrace" << endl; }
		else if (event.key.keysym.sym == SDLK_2) { RenderMode = draw_rasterise; renderMode = 2;cout << "[Render Mode]: rasterise" << endl; }
		else if (event.key.keysym.sym == SDLK_3) { RenderMode = draw_wireframe; renderMode = 3; cout << "[Render Mode]: wireframe" << endl; }

	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

    vector<ModelTriangle> triangles = obj_parser("models/cornell-box.obj", 0.5, mlt_parser("models/cornell-box.mtl"));
    TextureMap texture("models/texture.ppm");

    int planemultiplyer = 150;
    float focal = 2.0;

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);

        orbit(orbiting);

        RenderMode(triangles,focal, planemultiplyer, window);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
