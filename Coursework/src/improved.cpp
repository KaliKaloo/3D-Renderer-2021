#include <CanvasTriangle.h>
#include <ModelTriangle.h>
#include <CanvasPoint.h>
#include <TextureMap.h>
#include <Colour.h>
#include <DrawingWindow.h>


#include <RayTriangleIntersection.h>
#include <Utils.h>
#include <utility>
#include <fstream>
#include <vector>
#include <functional>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <unordered_map>

using namespace std;
using namespace glm;

#define WIDTH 600
#define HEIGHT 600

#define pi 3.14159265359

vec3   cam(0.0, 0.0, 4.0); // cornell cam
vec3 light(0.0, 1.0, 2.0); // cornell light 

// vec3   cam(0.0, 1.3, 3.5); // sphere cam
// vec3 light(1.0, 2.0, 2.8); // sphere light

mat3 cam_orientation(
	vec3(1.0,0.0,0.0),
	vec3(0.0,1.0,0.0),
	vec3(0.0,0.0,1.0)
);
int renderMode = 2; //initially 2 = rasterised mode


bool orbiting = false;
bool proximity = true, angle_of = true, shadows = true, specular = true;
void draw(DrawingWindow &window) {
	window.clearPixels();
}

void update(DrawingWindow &window) {
	// Function for performing animation (shifting artifacts or moving the camera)
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

void drawLine(CanvasPoint from, CanvasPoint to, Colour colour, vector<vector<float>> &depthBuffer, DrawingWindow &window){ 
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
	vec3 vertex = (vertexPosition - cam)*cam_orientation;

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
        strokedTriangle(ct, Colour(255,255,255), depthBuffer, window);
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


bool is_shadow(RayTriangleIntersection intersect, vector<ModelTriangle> triangles) {

	vec3 shadow_ray = light - intersect.intersectionPoint;

	for(int i = 0; i < triangles.size(); i++) {
		ModelTriangle tri = triangles[i];

		vec3 e0 = tri.vertices[1] - tri.vertices[0];
		vec3 e1 = tri.vertices[2] - tri.vertices[0];
		vec3 sp_vector = intersect.intersectionPoint - tri.vertices[0];
		mat3 de_matrix(-normalize(shadow_ray), e0, e1);
		vec3 possible_s = inverse(de_matrix) * sp_vector;
		float t = possible_s.x, u = possible_s.y, v = possible_s.z;

		if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
			if(t < glm::length(shadow_ray) && t > 0.05 && i != intersect.triangleIndex) {
				return true;
			}
		}
	}
	return false;
}

float diffused(RayTriangleIntersection rt_int, int scale) {

	vec3 normal = normalize(rt_int.intersectedTriangle.normal);
	vec3 light_ray = light - rt_int.intersectionPoint;
	vec3 view_ray = normalize(cam - rt_int.intersectionPoint);
	vec3 reflection_ray = normalize(normalize(light_ray) - (normal * 2.0f * dot(normalize(light_ray), normal)));

	float scale_p = (proximity) ? 50/(4*pi*(pow(length(light_ray),2))) : 0;
	float scale_a = dot(normal, normalize(light_ray));
	float scale_s = pow(dot(reflection_ray, view_ray),scale);

	if(scale_a > 0 && angle_of) scale_p *= scale_a;
	if(scale_s > 0 && specular) scale_p += scale_s;
	return (scale_p < 1) ? scale_p : 1;
}

float gourad(RayTriangleIntersection rt_int, int scale) {
	ModelTriangle triangle = rt_int.intersectedTriangle;
    glm::vec3 lightRay = light - rt_int.intersectionPoint;
	glm::vec3 cameraRay = (cam * cam_orientation) - rt_int.intersectionPoint;
	
    vector<float> brightnesses;
	float length = glm::length(lightRay);

	for(int i = 0; i < triangle.normals.size(); i++) {
		float angleOfIncidence = (angle_of) ? dot(triangle.normals[i], normalize(lightRay)): 1;

		vec3 reflection = normalize(lightRay) - ((2.0f*triangle.normals[i])*dot(normalize(lightRay), triangle.normals[i]));

		float brightness = (proximity) ? 50/(4 * pi * length*length) : 0;

		if (angleOfIncidence > 0) {
			brightness *= angleOfIncidence;
		} else {
			brightness *= 0;
		}
        
		float spec = (specular) ? pow(dot(normalize(reflection), normalize(cameraRay)), 128) : 0;

		if (spec >= 0) {
			brightness += spec * 0.2;
		}

		brightnesses.push_back(brightness);
	}

	float finalBrightness = (1 - rt_int.u - rt_int.v) * brightnesses[0] + rt_int.u * brightnesses[1] + rt_int.v * brightnesses[2];

	if (finalBrightness > 1) {
		finalBrightness = 1;
	} 
	if (finalBrightness < 0.2) {
		finalBrightness = 0.2;
	}
	return finalBrightness;
}

float phong(RayTriangleIntersection rt_int, int scale) {
	ModelTriangle triangle = rt_int.intersectedTriangle;
    glm::vec3 lightRay = light - rt_int.intersectionPoint;
	glm::vec3 cameraRay = (cam * cam_orientation) - rt_int.intersectionPoint;

    float length = glm::length(lightRay);
	
	glm::vec3 interpolatedNormal = normalize((1 - rt_int.u - rt_int.v) * triangle.normals[0] + rt_int.u * triangle.normals[1] + rt_int.v * triangle.normals[2]);

	vec3 reflection_ray = normalize(lightRay) - (2.0f*interpolatedNormal*dot(glm::normalize(lightRay), interpolatedNormal));


	float brightness = (proximity) ? 50/(4 * pi * length*length) : 0;
    float angleOfIncidence = (angle_of) ? dot(glm::normalize(lightRay), interpolatedNormal) : 1;
	float spec = (specular) ? pow(dot(normalize(reflection_ray), normalize(cameraRay)), scale) : 0;

	if (angleOfIncidence > 0) {
		brightness *= angleOfIncidence;
	} else {
		brightness *= 0;
	}

	if (spec >= 0) {
		brightness += spec*0.5;
	}

	if (brightness > 1.0f) {
		brightness = 1;
	} 
	if (brightness < 0.2f) {
		brightness = 0.2;
	}
	return brightness;
}

RayTriangleIntersection get_closest_intersection(vec3 direction, vector<ModelTriangle> triangles) {
	RayTriangleIntersection rti;
	rti.distanceFromCamera = numeric_limits<float>::infinity();

	for(int i = 0; i < triangles.size(); i++) {
		ModelTriangle tri = triangles[i];

		vec3 e0 = tri.vertices[1] - tri.vertices[0];
		vec3 e1 = tri.vertices[2] - tri.vertices[0];
		vec3 sp_vector = cam - tri.vertices[0];
		mat3 de_matrix(-direction, e0, e1);
		vec3 possible_s = inverse(de_matrix) * sp_vector;
		float t = possible_s.x, u = possible_s.y, v = possible_s.z;

		if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
			if(rti.distanceFromCamera > t && t > 0) {
				rti.distanceFromCamera = t;
				rti.intersectedTriangle = tri;
				// rti.intersectedTriangle.normal = cross(e1,e0);
				rti.triangleIndex = i;
				rti.u = u;
				rti.v = v;

				vec3 intersect = tri.vertices[0]+u*e0+v*e1;
				rti.intersectionPoint = intersect;
			}
		}
	}
	return rti;
}

function<float(RayTriangleIntersection intersect, int scale)> brightness = phong;

void draw_raytrace(vector<ModelTriangle> triangles, float focal, int planeMultiplyer, DrawingWindow &window) {

	for(int x = 0; x < window.width; x++) {
		for(int y = 0; y < window.height; y++) {
            vec3 direction(-(float(window.width / 2) - x) / planeMultiplyer, (float(window.height / 2) - y) / planeMultiplyer, -2);
			glm::vec3 camDir = glm::normalize(cam_orientation * direction);
			RayTriangleIntersection rt_int = get_closest_intersection( camDir, triangles);

			float scale = brightness(rt_int, 64);
			scale = (scale > 0.15) ? scale : 0.15;

			if(!isinf(rt_int.distanceFromCamera)){
				Colour colour = rt_int.intersectedTriangle.colour;
				uint32_t c = (255 << 24) + (int(colour.red*scale) << 16) + (int(colour.green*scale) << 8) + int(colour.blue*scale);

				if(is_shadow(rt_int, triangles) && shadows) {
					float scale_s = 0.1;
					uint32_t s = (255 << 24) + (int(colour.red*scale_s) << 16) + (int(colour.green*scale_s) << 8) + int(colour.blue*scale_s);
					window.setPixelColour(x,y,s); 
				} else window.setPixelColour(x,y,c);
			}
		}
	}
}

void vertexNormals(std::vector<ModelTriangle> &triangles) {
    
    for(int i = 0; i < triangles.size(); i++) {
        ModelTriangle t = triangles[i];
        std::vector<glm::vec3> normals;
        for(int v = 0; v < t.vertices.size(); v++) {
            glm::vec3 vertex = t.normal;
            int count = 1;
            for(int j = 0; j < triangles.size(); j++) {
                ModelTriangle t_ = triangles[j];
                for(int u = 0; u < t_.vertices.size(); u++) {
                    if(i != j && t.vertices[v].x == t_.vertices[u].x && t.vertices[v].y == t_.vertices[u].y && t.vertices[v].z == t_.vertices[u].z) {
                        if (acos(dot(t.normal, t_.normal)/(length(t.normal)*length(t_.normal))) < pi/4) {
							vertex = vertex + t_.normal;
							count = count + 1;
						}
                    }
                }
            }
            vertex = vertex / float(count);
            triangles[i].normals[v] = normalize(vertex);
        }
    }
}

std::vector<ModelTriangle> parse_obj(std::string filename, float scale, std::unordered_map<std::string, Colour> colours) {
	std::vector<ModelTriangle> output;
	std::vector<glm::vec3> vertices;
	std::vector<TexturePoint> textureVertices;
	std::string colour;
	std::vector<glm::vec3> normalVecs;

	std::ifstream File(filename);
	std::string line;

	if (filename.compare("logo.obj") == 0) colour = "texture";

	// if (filename.compare("textured-cornell-box.obj") == 0) colour = "texture";

	// std::cout << colour << std::endl;
	
	while(std::getline(File, line)) {
		if(line == "") continue;

		std::vector<std::string> tokens = split(line, ' ');

		if (tokens[0] == "v") {
			glm::vec3 temp = glm::vec3(stof(tokens[1])*scale, stof(tokens[2])*scale, stof(tokens[3])*scale); 
			vertices.push_back(temp);

		} else if (tokens[0] == "f") {
			// when no texture map vertices, the second vector item is equal to ""
			std::vector<std::string> a = split(tokens[1],'/');
			std::vector<std::string> b = split(tokens[2],'/');
			std::vector<std::string> c = split(tokens[3],'/');
			
			ModelTriangle triangle(vertices[stoi(a[0])-1], vertices[stoi(b[0])-1], vertices[stoi(c[0])-1], colours[colour]);
			triangle.normal = glm::normalize(glm::cross(glm::vec3(triangle.vertices[1] - triangle.vertices[0]), glm::vec3(triangle.vertices[2] - triangle.vertices[0])));
			if (!normalVecs.empty()) {
				triangle.normals[0] = normalVecs[stoi(a[2])-1];
				triangle.normals[1] = normalVecs[stoi(b[2])-1];
				triangle.normals[2] = normalVecs[stoi(c[2])-1];
			} 
			if (colour.compare("Mirror") == 0)  {
				triangle.mirror = true;
			}
			if (colour.compare("Glass") == 0) {
				triangle.glass = true;
			} else {
				triangle.glass = false;
			}
			
			if(!textureVertices.empty() && a[1] != "") {
				// std::cout << colours[colour] << std::endl;
				// ModelTriangle triangle(vertices[stoi(a[0])-1], vertices[stoi(b[0])-1], vertices[stoi(c[0])-1], colours[colour]);
				triangle.texturePoints[0] = textureVertices[stoi(a[1])-1];
				triangle.texturePoints[1] = textureVertices[stoi(b[1])-1];
				triangle.texturePoints[2] = textureVertices[stoi(c[1])-1];
			}
			
			output.push_back(triangle);
			
		} else if (tokens[0] == "usemtl") {
			colour = tokens[1];
			if (colour[colour.size() - 1] == '\r') {
				colour.erase(colour.size() - 1);
			}
		} else if (tokens[0] == "vt") {
			TexturePoint temp = TexturePoint(stof(tokens[1]), stof(tokens[2]));
			textureVertices.push_back(temp);
		} 
		else if (tokens[0] == "vn") {
			glm::vec3 temp = glm::vec3(stof(tokens[1]), stof(tokens[2]), stof(tokens[3])); 
			normalVecs.push_back(temp);
		}
	}

	if (normalVecs.empty()) {
		vertexNormals(output);
	}

	File.close();
	return output;
}

unordered_map<string, Colour> parse_mtl(string filename) {
	unordered_map<string, Colour> colours;
	string colour_name;

	ifstream File(filename);
	string line;

	while(getline(File, line)) {
		if(line == "") continue;

		vector<string> tokens = split(line, ' ');
		if(tokens[0] == "newmtl") {
			colour_name = tokens[1];
		} else if(tokens[0] == "Kd") {
			Colour colour(int(stof(tokens[1])*255),int(stof(tokens[2])*255),int(stof(tokens[3])*255));
			colours.insert({colour_name, colour});
		} else if(tokens[0] == "map_Kd") {
			Colour colour = colours[colour_name];
			colour.name = tokens[1];
			colours[colour_name] = colour;
		}
	}
	File.close();
	return colours;
}

void reset_camera() {
	cam = vec3(0.0,0.0,4.0);
	light = vec3(0.0, 1.0, 2.0);
	cam_orientation = mat3(vec3(1.0,0.0,0.0),vec3(0.0,1.0,0.0),vec3(0.0,0.0,1.0));
}

void look_at() {
	vec3 forward = normalize(cam - vec3(0.0,0.0,0.0));
	vec3 right = normalize(cross(vec3(0.0,1.0,0.0), forward));
	vec3 up = normalize(cross(forward, right));

	cam_orientation[0] = right;
	cam_orientation[1] = up;
	cam_orientation[2] = forward;
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
		cam = cam * rotation_y(-pi/rad);
		look_at();
	}
}

function<void(vector<ModelTriangle>, float, int, DrawingWindow &)> drawing = draw_raytrace;

void handleEvent(SDL_Event event, DrawingWindow &window) {
    float rotVal = pi/60;

	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_s) cam.y -= 0.1;
		else if (event.key.keysym.sym == SDLK_w) cam.y += 0.1;
		else if (event.key.keysym.sym == SDLK_a) cam.x -= 0.1;
		else if (event.key.keysym.sym == SDLK_d) cam.x += 0.1;
		else if (event.key.keysym.sym == SDLK_x) cam.z -= 0.1; 
		else if (event.key.keysym.sym == SDLK_z) cam.z += 0.1;
		else if (event.key.keysym.sym == SDLK_RIGHT){
            cam = cam * rotation_y( -rotVal); 
            look_at();
        }
		else if (event.key.keysym.sym == SDLK_LEFT){
            cam = cam * rotation_y(rotVal);
            look_at();
        }
		else if (event.key.keysym.sym == SDLK_DOWN){
            cam = cam * rotation_x( -rotVal);
            look_at();
        }
		else if (event.key.keysym.sym == SDLK_UP) {
            cam = cam * rotation_x(rotVal);
            look_at();
        }
		else if (event.key.keysym.sym == SDLK_o) orbiting = (orbiting) ? false : true;
		else if (event.key.keysym.sym == SDLK_l) look_at();
		else if (event.key.keysym.sym == SDLK_r) reset_camera();
		else if (event.key.keysym.sym == SDLK_1) { drawing = draw_raytrace; cout << "[drawing]: raytrace" << endl; }
		else if (event.key.keysym.sym == SDLK_2) { drawing = draw_rasterise; cout << "[drawing]: rasterise" << endl; }
		else if (event.key.keysym.sym == SDLK_3) { drawing = draw_wireframe; cout << "[drawing]: wireframe" << endl; }
		else if (event.key.keysym.sym == SDLK_4) { brightness = diffused; cout << "[lighting]: diffused" << endl; }
		else if (event.key.keysym.sym == SDLK_5) { brightness = gourad; cout << "[lighting]: gourad" << endl; }
		else if (event.key.keysym.sym == SDLK_6) { brightness = phong; cout << "[lighting]: phong" << endl; }
		else if (event.key.keysym.sym == SDLK_KP_8) light.z -= 0.1;
		else if (event.key.keysym.sym == SDLK_KP_2) light.z += 0.1;
		else if (event.key.keysym.sym == SDLK_KP_6) light.x += 0.1;
		else if (event.key.keysym.sym == SDLK_KP_4) light.x -= 0.1;
		else if (event.key.keysym.sym == SDLK_KP_MINUS) light.y -= 0.1;
		else if (event.key.keysym.sym == SDLK_KP_PLUS) light.y += 0.1;

		else if (event.key.keysym.sym == SDLK_LEFTBRACKET)  { proximity = (proximity) ? false : true; cout << "[proximity]: " << proximity << endl; }
		else if (event.key.keysym.sym == SDLK_RIGHTBRACKET) { angle_of  = (angle_of)  ? false : true; cout << "[angle_of]: " << angle_of << endl; }
		else if (event.key.keysym.sym == SDLK_HASH)         { shadows   = (shadows)   ? false : true; cout << "[shadows]: " << shadows << endl; }
		else if (event.key.keysym.sym == SDLK_QUOTE)        { specular  = (specular)  ? false : true; cout << "[specular]: " << specular << endl; }
	}
}

int main(int argc, char *argv[]) {

	vector<ModelTriangle> t = parse_obj("models/cornell-box.obj", 0.5, parse_mtl("models/cornell-box.mtl"));
	// vector<ModelTriangle> t_2 = parse_obj("models/sphere.obj", 0.5, parse_mtl("models/cornell-box.mtl"));
	// t.insert(t.end(), t_2.begin(), t_2.end());

    float focal = 2.0;
    int planemultiplyer = 150;

	DrawingWindow window_grey = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window_grey.pollForInputEvents(event)) handleEvent(event, window_grey);
		orbit(orbiting);

		draw(window_grey);
		drawing(t, focal, planemultiplyer, window_grey); //cornell-box or cornell-box+sphere
		// drawing(t_2, focal, planemultiplyer, window_grey); //sphere only

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window_grey.renderFrame();
	}
}