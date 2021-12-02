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

// Function Headers
RayTriangleIntersection get_closest_refraction(vec3 point, vec3 dir_reflection, vector<ModelTriangle> triangles, int index, int recursion);

#define WIDTH 640
#define HEIGHT 480

#define pi 3.14159265359

vec3   cam(0.0, 0.0, 4.0); // cornell cam
vec3 light(0.0, 1.0, 1.0); // cornell light 

// vec3   cam(0.0, 0.0, 2.5); // sphere cam
// vec3 light(1.0, 2.0, 2.8); // sphere light

// vec3   cam(0.5, 0.5, 2.0); //logo cam
// vec3 light(1.0, 1.0, 2.0);  //logo light

mat3 cam_orientation(
	vec3(1.0,0.0,0.0),
	vec3(0.0,1.0,0.0),
	vec3(0.0,0.0,1.0)
);

int renderMode = 2; //initially 2 = rasterised mode
bool orbiting = false;
bool proximity = true, angleOfInc = true, softShadows = false, hardShadows = true, specular = true;

void draw(DrawingWindow &window) {
	window.clearPixels();
}

void update(DrawingWindow &window) {
}

float getNumberOfSteps(CanvasPoint from, CanvasPoint to) {
	return fmax(fmax(abs(to.x - from.x), abs(to.y - from.y)), 1);
}

vector<TexturePoint> interpolatePoints(TexturePoint from, TexturePoint to, int steps) {
	vector<TexturePoint> TexturePoints;
	float xs = (to.x - from.x) / (steps - 1);
	float ys = (to.y - from.y) / (steps - 1);
	for (int i=0; i<steps; i++) {
		TexturePoints.push_back(TexturePoint(from.x + (i * xs),  from.y + (i * ys)));
	}
	return TexturePoints;
}

vector<CanvasPoint> interpolatePoints(CanvasPoint from, CanvasPoint to, int numberOfValues){
    vector<CanvasPoint> points;
    float xs = (to.x - from.x)/(numberOfValues-1);
    float ys = (to.y - from.y)/(numberOfValues-1);
    float zs = (to.depth - from.depth) / (numberOfValues-1);

    for(int i=0; i<numberOfValues; i++){
        CanvasPoint cp(round(from.x+(i*xs)), round(from.y+(i*ys)), from.depth + (i*zs));
        points.push_back(cp); 
    }
    return points;
}

void drawLine(CanvasPoint from, CanvasPoint to, Colour colour, DrawingWindow &window){ 

    float numberOfSteps= getNumberOfSteps(from, to);
    vector<CanvasPoint> points = interpolatePoints(from, to, numberOfSteps+2);

    for (float i=0.0; i<numberOfSteps; i++){
        if(points[i].x >= 0 && points[i].x < window.width && points[i].y >= 0 && points[i].y < window.height) {
            uint32_t c = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
            window.setPixelColour(points[i].x, points[i].y, c);
        }
    }
}

void strokedTriangle(CanvasTriangle t, Colour colour, DrawingWindow &window){
    CanvasPoint p_1 = t[0];
    CanvasPoint p_2 = t[1];
    CanvasPoint p_3 = t[2];

    drawLine(p_1,p_2,colour, window);
    drawLine(p_1,p_3,colour, window);
    drawLine(p_2,p_3,colour,window);
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
    
void half_filledTriangle(CanvasPoint v1, CanvasPoint v2, CanvasPoint split, Colour colour, vector<vector<float>> &depthBuffer, DrawingWindow &window){
	v1.depth = -1/v1.depth;
	v2.depth = -1/v2.depth;
	split.depth = -1/split.depth;

    vector<CanvasPoint> left = interpolatePoints(v1, v2, abs(v2.y-v1.y)+2);

	vector<CanvasPoint> right = interpolatePoints(v1, split, abs(v2.y-v1.y)+2);

	for (int i = 0; i < left.size(); i++) {
		int steps = abs(left[i].x - right[i].x);
				
		vector<CanvasPoint> points = interpolatePoints(left[i], right[i], steps+2);
		
		for (int c = 0; c < points.size(); c++) {
			int newX = round(points[c].x);
			int newY = round(points[c].y);
			if (newX >= 0 && newX < window.width && newY >= 0 && newY < window.height){
				if (points[c].depth > depthBuffer[newX][newY]) {

					depthBuffer[newX][newY] = points[c].depth;
					uint32_t set = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
					window.setPixelColour(newX, newY, set);
				}
			}
		}
	}
}

void filledTriangle(CanvasTriangle triangle, Colour colour,vector<vector<float>> &depthBuffer, DrawingWindow &window){

    sortTriangleVertices(triangle);

    CanvasPoint top = triangle.vertices[0];
	CanvasPoint mid = triangle.vertices[1];
	CanvasPoint bot = triangle.vertices[2];

	CanvasPoint split;
	split.y = mid.y;

	split.x = round(top.x + ((mid.y - top.y)/(bot.y-top.y)) * (bot.x-top.x));

	split.depth = top.depth + ((mid.y - top.y)/(bot.y-top.y)) * (bot.depth-top.depth);

    half_filledTriangle(top, mid, split, colour, depthBuffer, window);
    half_filledTriangle(bot, mid, split, colour, depthBuffer, window);

}

void half_texturedTriangle(CanvasPoint v1, CanvasPoint v2, CanvasPoint midpoint, const TextureMap texture, vector<vector<float>> &depthBuffer, DrawingWindow &window){
	v1.depth = -1/v1.depth;
	v2.depth = -1/v2.depth;
	midpoint.depth = -1/midpoint.depth;

    vector<CanvasPoint> left = interpolatePoints(v1, v2, abs(v2.y-v1.y)+2);
	vector<TexturePoint> leftTexture = interpolatePoints(v1.texturePoint, v2.texturePoint, abs(v2.y-v1.y)+2);

	vector<CanvasPoint> right = interpolatePoints(v1, midpoint, abs(v2.y-v1.y)+2);
	vector<TexturePoint> rightTexture = interpolatePoints(v1.texturePoint, midpoint.texturePoint, abs(v2.y-v1.y)+2);

	for (int i = 0; i < left.size(); i++) {
		int steps = abs(left[i].x - right[i].x);
				
		vector<CanvasPoint> points = interpolatePoints(left[i], right[i], steps+2);

		vector<TexturePoint> texturePoints = interpolatePoints(leftTexture[i], rightTexture[i], steps+2);

		for (int c = 0; c < points.size(); c++) {
			int newX = round(points[c].x);
			int newY = round(points[c].y);
			if (newX >= 0 && newX < window.width && newY >= 0 && newY < window.height) {
				if (points[c].depth > depthBuffer[newX][newY]) {
					depthBuffer[newX][newY] = points[c].depth;
					int x_coord = round(texturePoints.at(c).x);
					int y_coord = round(texturePoints.at(c).y);
					uint32_t col = texture.pixels[y_coord*texture.width + x_coord];

					window.setPixelColour(newX, newY, col);	
				}			
			}
		}
	}
}

void texturedTriangle( CanvasTriangle triangle, const TextureMap texture, vector<vector<float>> &depthBuffer, DrawingWindow &window) {
    sortTriangleVertices(triangle);

    CanvasPoint top = triangle.vertices[0];
	CanvasPoint mid = triangle.vertices[1];
	CanvasPoint bot = triangle.vertices[2];

	CanvasPoint midpoint;
	midpoint.y = mid.y;
	midpoint.x = round(top.x + ((mid.y - top.y)/(bot.y-top.y)) * (bot.x-top.x));

	midpoint.depth = top.depth + ((mid.y - top.y)/(bot.y-top.y)) * (bot.depth-top.depth);

	float scale = (mid.y - top.y)/(bot.y-top.y);

	midpoint.texturePoint.x = top.texturePoint.x + scale * (bot.texturePoint.x - top.texturePoint.x);
	midpoint.texturePoint.y = top.texturePoint.y + scale * (bot.texturePoint.y - top.texturePoint.y);

    half_texturedTriangle(top, mid, midpoint, texture, depthBuffer, window);
    half_texturedTriangle(bot, mid, midpoint, texture, depthBuffer,window);
}



CanvasPoint getCanvasIntersectionPoint(glm::vec3 vertexPosition, float focal, int planeMultiplyer, DrawingWindow &window){
	vec3 vertex = (vertexPosition - cam)*cam_orientation;

    int u = -round(planeMultiplyer*focal * (vertex.x)/(vertex.z)) + (window.width/2);
    int v = round(planeMultiplyer*focal * (vertex.y)/(vertex.z)) + (window.height/2);

    CanvasPoint t(u, v, vertex.z);
    return t;
}

// DRAW WIREFRAME SCENE
void draw_wireframe(vector<ModelTriangle> triangles, float focal, int planeMultiplyer, unordered_map<string, TextureMap> textures,vector<glm::vec3> lightDirections, DrawingWindow &window) {
	window.clearPixels();

    for(int i=0; i<triangles.size();i++){
        ModelTriangle triangle = triangles[i];
        CanvasTriangle ct;
        for(int j=0; j<triangle.vertices.size(); j++){
            vec3 vertex = triangle.vertices[j];
            ct.vertices[j] = getCanvasIntersectionPoint(vertex, focal, planeMultiplyer, window);
        }
        strokedTriangle(ct, Colour(255,255,255), window);
    }
}

// DRAW RASTERISED SCENE
void draw_rasterise(vector<ModelTriangle> triangles, float focal, int planeMultiplyer, unordered_map<string, TextureMap> textures, vector<glm::vec3> lightDirections, DrawingWindow &window) {
	window.clearPixels();
    vector<vector<float>> depthBuffer(window.width, vector<float> (window.height,-numeric_limits<float>::infinity()));

    for(int i=0; i<triangles.size();i++){
        ModelTriangle triangle = triangles[i];
        CanvasTriangle ct;
        bool isTexture = false;
        TextureMap texture;

        if(triangle.colour.name != ""){
            texture = TextureMap(triangle.colour.name);
            isTexture = true;
        }
        
        for(int j=0; j<triangle.vertices.size(); j++){
            vec3 vertex = triangle.vertices[j];
            ct.vertices[j] = getCanvasIntersectionPoint(vertex,focal, planeMultiplyer, window);
            
            if(isTexture){
                ct.vertices[j].texturePoint = triangle.texturePoints[j];
                ct.vertices[j].texturePoint.x *= texture.width;
                ct.vertices[j].texturePoint.y *= texture.height;
            }
        }
        if(isTexture){
            texturedTriangle(ct, texture, depthBuffer, window);

        }else{
            filledTriangle(ct, triangle.colour, depthBuffer, window);
        }
    }
}

bool is_shadow(RayTriangleIntersection intersect, vector<ModelTriangle> triangles, vec3 l) {

	vec3 shadow_ray = l - intersect.intersectionPoint;

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

float diffused(RayTriangleIntersection intersect, vec3 l, int scale) {
	vec3 normal = normalize(intersect.intersectedTriangle.normal);
	vec3 light_ray = l - intersect.intersectionPoint;
	vec3 view_ray = normalize(cam - intersect.intersectionPoint);
	vec3 reflection_ray = normalize(normalize(light_ray) - (normal * 2.0f * dot(normalize(light_ray), normal)));

	float scale_p = (proximity) ? 50/(4*pi*(pow(length(light_ray),2))) : 0;
	float scale_a = dot(normal, normalize(light_ray));
	float scale_s = pow(dot(reflection_ray, view_ray),scale);

	if(scale_a > 0 && angleOfInc) {scale_p *= scale_a;
    }
    else{
        scale_p *= 0;
    }
	if(scale_s > 0 && specular) scale_p += scale_s;
	return (scale_p < 1) ? scale_p : 1;
}

float gouraurd(RayTriangleIntersection intersect, vec3 l,int scale) {
	ModelTriangle triangle = intersect.intersectedTriangle;
    glm::vec3 lightRay = l - intersect.intersectionPoint;
	glm::vec3 cameraRay = (cam - intersect.intersectionPoint) * cam_orientation ;
	
    vector<float> brightnesses;
	float length = glm::length(lightRay);

	for(int i = 0; i < triangle.normals.size(); i++) {
		float angleOfIncidence = (angleOfInc) ? dot(triangle.normals[i], normalize(lightRay)): 1;

		vec3 reflection = normalize(lightRay) - ((2.0f*triangle.normals[i])*dot(normalize(lightRay), triangle.normals[i]));

		float brightness = (proximity) ? 50/(4 * pi * length*length) : 0;

		if (angleOfIncidence > 0) {
			brightness *= angleOfIncidence;
		} else {
			brightness *= 0;
		}
        
		float spec = (specular) ? pow(dot(normalize(reflection), normalize(cameraRay)), 128) : 0;

		if (spec >= 0) {
			brightness += spec * 0.5;
		}

		brightnesses.push_back(brightness);
	}

	float finalBrightness = (1 - intersect.u - intersect.v) * brightnesses[0] + intersect.u * brightnesses[1] + intersect.v * brightnesses[2];

	if (finalBrightness > 1) {
		finalBrightness = 1;
	} 
	if (finalBrightness < 0.2) {
		finalBrightness = 0.2;
	}
	return finalBrightness;
}

float phong(RayTriangleIntersection intersect, vec3 l, int scale) {
	ModelTriangle triangle = intersect.intersectedTriangle;
    glm::vec3 lightRay = l - intersect.intersectionPoint;
	glm::vec3 cameraRay = (cam - intersect.intersectionPoint) * cam_orientation ;

    float length = glm::length(lightRay);
	
	glm::vec3 interpolatedNormal = normalize((1 - intersect.u - intersect.v) * triangle.normals[0] + intersect.u * triangle.normals[1] + intersect.v * triangle.normals[2]);

	vec3 reflection_ray = normalize(lightRay) - (2.0f*interpolatedNormal*dot(glm::normalize(lightRay), interpolatedNormal));

	float brightness = (proximity) ? 50/(4 * pi * length*length) : 1;
    float angleOfIncidence = (angleOfInc) ? dot(glm::normalize(lightRay), interpolatedNormal) : 1;
	float spec = (specular) ? pow(dot(normalize(reflection_ray), normalize(cameraRay)), scale) : 0;

	if (angleOfIncidence > 0) {
		brightness *= angleOfIncidence;
	} else {
		brightness *= 0;
	}
	if (spec >= 0) brightness += spec*0.2;
	if (brightness > 1.0f) brightness = 1;
	if (brightness < 0.2f) brightness = 0.2;
	return brightness;
}

vec3 refract(vec3 direction, vec3 normal, float refraction_index) {
	float a = refraction_index;
	float b = dot(direction, normal);

	if (b >= 0.0f) {
		a = 1.0f/a;
	}

	return normalize(direction * a - normal * (-b + a*b));
}

RayTriangleIntersection get_closest_reflection(vec3 point, vec3 dir_reflection, vector<ModelTriangle> triangles, int index, int recursion){
    RayTriangleIntersection intersect;
    intersect.distanceFromCamera = numeric_limits<float>::infinity();

    for(int i = 0; i < triangles.size(); i++) {
		ModelTriangle tri = triangles[i];

		vec3 e0 = tri.vertices[1] - tri.vertices[0];
		vec3 e1 = tri.vertices[2] - tri.vertices[0];
		vec3 sp_vector = point - tri.vertices[0];
		mat3 de_matrix(-dir_reflection, e0, e1);
		vec3 possible_s = inverse(de_matrix) * sp_vector;
		float t = possible_s.x, u = possible_s.y, v = possible_s.z;

		if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
            if(intersect.distanceFromCamera > t && t > 0 && i != index) {
				intersect.distanceFromCamera = t;
                vec3 point = tri.vertices[0]+u*e0+v*e1;
                intersect.intersectionPoint = point;

                intersect.intersectedTriangle = tri;
                intersect.triangleIndex = i;
                intersect.u = u;
                intersect.v = v;
            }
        }
    }

    if(recursion<5){
        if(intersect.intersectedTriangle.mirror){
            vec3 normal = normalize(intersect.intersectedTriangle.normal);
            vec3 reflection = normalize(dir_reflection - (normal * 2.0f * dot(dir_reflection, normal)));
            
            RayTriangleIntersection reflect_intersect = get_closest_reflection(intersect.intersectionPoint, reflection, triangles, intersect.triangleIndex, recursion++);
            intersect = reflect_intersect;
        }
        if(intersect.intersectedTriangle.glass){
            vec3 normal = normalize(intersect.intersectedTriangle.normal);
            vec3 refracted = refract(dir_reflection, normal, 1.5);

            RayTriangleIntersection refract_intersect = get_closest_refraction(intersect.intersectionPoint , refracted, triangles, intersect.triangleIndex, recursion++);
            intersect = refract_intersect;
        }
    }
    
    return intersect;
}

RayTriangleIntersection get_closest_refraction(vec3 point, vec3 dir_reflection, vector<ModelTriangle> triangles, int index, int recursion){
    RayTriangleIntersection intersect;
    intersect.distanceFromCamera = numeric_limits<float>::infinity();

    for(int i = 0; i < triangles.size(); i++) {
		ModelTriangle tri = triangles[i];

		vec3 e0 = tri.vertices[1] - tri.vertices[0];
		vec3 e1 = tri.vertices[2] - tri.vertices[0];
		vec3 sp_vector = point - tri.vertices[0];
		mat3 de_matrix(-dir_reflection, e0, e1);
		vec3 possible_s = inverse(de_matrix) * sp_vector;
		float t = possible_s.x, u = possible_s.y, v = possible_s.z;

		if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
            if(intersect.distanceFromCamera > t && t > 0 && i != index) {
				intersect.distanceFromCamera = t;
                vec3 point = tri.vertices[0]+u*e0+v*e1;
                intersect.intersectionPoint = point;

                intersect.intersectedTriangle = tri;
                intersect.triangleIndex = i;
                intersect.u = u;
                intersect.v = v;
            }
        }
    }
    if(recursion<5){
        if(intersect.intersectedTriangle.mirror){
            vec3 normal = normalize(intersect.intersectedTriangle.normal);
            vec3 reflection = normalize(dir_reflection - (normal * 2.0f * dot(dir_reflection, normal)));
            
            RayTriangleIntersection reflect_intersect = get_closest_reflection(intersect.intersectionPoint, reflection, triangles, intersect.triangleIndex, recursion++);
            intersect = reflect_intersect;
        }
        if(intersect.intersectedTriangle.glass){
            vec3 normal = normalize(intersect.intersectedTriangle.normal);
            vec3 refracted = refract(dir_reflection, normal, 1.5);

            RayTriangleIntersection refract_intersect = get_closest_refraction(intersect.intersectionPoint , refracted, triangles, intersect.triangleIndex, recursion++);
            intersect = refract_intersect;
        }
    }
    return intersect;
}

RayTriangleIntersection get_closest_intersection(vec3 direction, vector<ModelTriangle> triangles) {
	RayTriangleIntersection intersect;
	intersect.distanceFromCamera = numeric_limits<float>::infinity();
    float distance = numeric_limits<float>::infinity();

	for(int i = 0; i < triangles.size(); i++) {
		ModelTriangle tri = triangles[i];

		vec3 e0 = tri.vertices[1] - tri.vertices[0];
		vec3 e1 = tri.vertices[2] - tri.vertices[0];
		vec3 sp_vector = cam - tri.vertices[0];
		mat3 de_matrix(-direction, e0, e1);
		vec3 possible_s = inverse(de_matrix) * sp_vector;
		float t = possible_s.x, u = possible_s.y, v = possible_s.z;

		if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {

			if(distance > t && t > 0) {
                vec3 point = tri.vertices[0]+u*e0+v*e1;
                
                if(tri.mirror){
                    vec3 normal = normalize(tri.normal);
                    vec3 reflection = normalize(direction - (normal * 2.0f * dot(direction, normal)));
                    
                    RayTriangleIntersection reflect_intersect = get_closest_reflection(point, reflection, triangles, i, 1);
                    intersect = reflect_intersect;
                    if(isinf(intersect.distanceFromCamera)) intersect.isInfinity = true;

                }else if(tri.glass){
                    vec3 normal = normalize(tri.normal);
                    vec3 refracted = refract(direction, normal, 1.5);

                    RayTriangleIntersection refract_intersect = get_closest_refraction(point, refracted, triangles, i, 1);
                    intersect = refract_intersect;
                }
                else{
                    intersect.intersectionPoint = point;

                    intersect.distanceFromCamera = t;
                    intersect.intersectedTriangle = tri;
                    intersect.triangleIndex = i;
                    intersect.u = u;
                    intersect.v = v;
                }
                distance = t;
			}
		}
	}
	return intersect;
}

function<float(RayTriangleIntersection intersect, vec3 l,  int scale)> pixelBrightness = phong;

void draw_raytrace(vector<ModelTriangle> triangles, float focal, int planeMultiplyer, unordered_map<string, TextureMap> textures, vector<glm::vec3> lightDirections, DrawingWindow &window) {

    for(int x = 0; x < window.width; x++) {
		for(int y = 0; y < window.height; y++) {
            vec3 direction(-(float(window.width / 2) - x) / planeMultiplyer, (float(window.height / 2) - y) / planeMultiplyer, -focal);
			glm::vec3 camDir = normalize(cam_orientation * (direction));
			RayTriangleIntersection intersect = get_closest_intersection( camDir, triangles);
            Colour colour = intersect.intersectedTriangle.colour;
            float brightness = 0.0;

			if(!isinf(intersect.distanceFromCamera)){
                ModelTriangle t = intersect.intersectedTriangle;

                if(softShadows){
                    for(int l=0; l<lightDirections.size(); l++){
                        float temp = pixelBrightness(intersect,lightDirections[l], 256);
                        bool shadow = is_shadow(intersect, triangles, lightDirections[l]);
                        if(shadow) temp = 0.18;
                        brightness += temp;
                        }
                    brightness = brightness/lightDirections.size();
                }else if(hardShadows){
                    float temp = pixelBrightness(intersect,lightDirections[0], 256);
                    bool shadow = is_shadow(intersect, triangles, lightDirections[0]);
                    if(shadow) temp = 0.18;
                    brightness += temp;
                }
                else{
                    float temp = pixelBrightness(intersect, lightDirections[0], 256);
                    brightness += temp;
                }

                if(intersect.intersectedTriangle.colour.name != ""){
                    //get texture map
                    TextureMap texture = textures[t.colour.name];

                    float t_x = ((1 - intersect.u - intersect.v) * t.texturePoints[0].x + intersect.u * t.texturePoints[1].x + intersect.v * t.texturePoints[2].x);

					float t_y = ((1 - intersect.u - intersect.v) * t.texturePoints[0].y + intersect.u * t.texturePoints[1].y + intersect.v * t.texturePoints[2].y);

                    int tt_x = t_x * texture.width;
                    int tt_y = t_y * texture.height;
                    
                    uint32_t  t_value = texture.pixels[tt_y*texture.width + tt_x];

                    int r = (t_value >> 16) & 0xff;
                    int g = (t_value >>  8) & 0xff;
                    int b =        t_value  & 0xff;
                    colour = Colour(r,g,b);
                }
                colour.red *=  brightness;
                colour.green *= brightness;
                colour.blue *=  brightness;
                
                if(intersect.isInfinity){
                    colour.red = 0;
                    colour.green = 0;
                    colour.blue =0;
                }
			}
            uint32_t col = (255 << 24) + (int(colour.red) <<16) + (int(colour.green) << 8) + int(colour.blue);
            window.setPixelColour(x,y,col);
		}
	}
}

void vertexNormals(vector<ModelTriangle> &triangles) {
    
    for(int i = 0; i < triangles.size(); i++) {
        ModelTriangle t = triangles[i];
        vector<glm::vec3> normals;
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

vector<ModelTriangle> load_obj(string filename, float scale, unordered_map<string, Colour> colours) {
	vector<ModelTriangle> output;
	vector<glm::vec3> vertices;
	vector<TexturePoint> textureVertices;
	string colour;
	vector<glm::vec3> normalVecs;

	ifstream File(filename);
	string line;

	while(getline(File, line)) {
		if(line == "") continue;

		vector<string> tokens = split(line, ' ');

		if (tokens[0] == "v") {
			glm::vec3 temp = glm::vec3(stof(tokens[1])*scale, stof(tokens[2])*scale, stof(tokens[3])*scale); 
			vertices.push_back(temp);
        } else if (tokens[0] == "vt") {
			TexturePoint temp = TexturePoint(stof(tokens[1]), stof(tokens[2]));
			textureVertices.push_back(temp);
		} 
		else if (tokens[0] == "vn") {
			glm::vec3 temp = glm::vec3(stof(tokens[1]), stof(tokens[2]), stof(tokens[3])); 
			normalVecs.push_back(temp);
		}
        else if (tokens[0] == "vt") {
			TexturePoint temp = TexturePoint(stof(tokens[1]), stof(tokens[2]));
			textureVertices.push_back(temp);
		} 
        else if (tokens[0] == "f") {
			vector<string> a = split(tokens[1],'/');
			vector<string> b = split(tokens[2],'/');
			vector<string> c = split(tokens[3],'/');
			
			ModelTriangle triangle(
                vertices[stoi(a[0])-1], 
                vertices[stoi(b[0])-1], 
                vertices[stoi(c[0])-1], 
                colours[colour]);
			
            triangle.normal = normalize(cross(vec3(triangle.vertices[1] - triangle.vertices[0]), vec3(triangle.vertices[2] - triangle.vertices[0])));
			
            if (!normalVecs.empty()) {
				triangle.normals[0] = normalVecs[stoi(a[2])-1];
				triangle.normals[1] = normalVecs[stoi(b[2])-1];
				triangle.normals[2] = normalVecs[stoi(c[2])-1];
			}
            if (colour.compare("Glass") == 0)  {
				triangle.glass = true;
			} else{triangle.glass = false;}
			if (colour.compare("Mirror") == 0)  {
				triangle.mirror = true;
			}else{triangle.mirror = false;}
			
			if(!textureVertices.empty() && a[1] != "") {
		
				triangle.texturePoints[0] = textureVertices[stoi(a[1])-1];
				triangle.texturePoints[1] = textureVertices[stoi(b[1])-1];
				triangle.texturePoints[2] = textureVertices[stoi(c[1])-1];
			}
			output.push_back(triangle);

        }else if (tokens[0] == "usemtl") {
			colour = tokens[1];
		}
	}

	if (normalVecs.empty()) {
		vertexNormals(output);
	}

	File.close();
	return output;
}

unordered_map<string, Colour> load_mtl(string filename, unordered_map<string, TextureMap> &textures) {
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
			textures.insert({tokens[1], TextureMap(tokens[1])});
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
        int rad = 60;
        if (renderMode == 1) rad = 400;
        if (renderMode ==3) rad = 50;
		cam = cam * rotation_y(-pi/rad);
		look_at();
	}
}

function<void(vector<ModelTriangle>, float, int, unordered_map<string, TextureMap>,vector<glm::vec3> lightDirections, DrawingWindow &)> drawing = draw_raytrace;

void handleEvent(SDL_Event event, DrawingWindow &window) {
    float rotVal = (2*pi)/60;

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
        else if (event.key.keysym.sym == SDLK_l)  {
            cam_orientation =  rotation_y(-rotVal)*cam_orientation;
        }
		else if (event.key.keysym.sym == SDLK_j) {
            cam_orientation = rotation_y( rotVal)*cam_orientation ;
        }
		else if (event.key.keysym.sym == SDLK_k) {
            cam_orientation = rotation_x(-rotVal)*cam_orientation ;
        }
		else if (event.key.keysym.sym == SDLK_i)  {
            cam_orientation = rotation_x( rotVal)* cam_orientation ;
        }

		else if (event.key.keysym.sym == SDLK_o) orbiting = (orbiting) ? false : true;
		else if (event.key.keysym.sym == SDLK_l) look_at();
		else if (event.key.keysym.sym == SDLK_r) reset_camera();
		else if (event.key.keysym.sym == SDLK_1) { drawing = draw_wireframe; cout << "[Rendering Mode]: wireframe" << endl;}
		else if (event.key.keysym.sym == SDLK_2) { drawing = draw_rasterise; cout << "[Rendering Mode]: rasterise" << endl; }
		else if (event.key.keysym.sym == SDLK_3) { drawing = 
        draw_raytrace; cout << "[Rendering Mode]: raytrace" << endl; }
		else if (event.key.keysym.sym == SDLK_4) { pixelBrightness = diffused; cout << "[Lighting]: diffused" << endl; }
		else if (event.key.keysym.sym == SDLK_5) { pixelBrightness = gouraurd; cout << "[Lighting]: gouraurd" << endl; }
		else if (event.key.keysym.sym == SDLK_6) { pixelBrightness = phong; cout << "[Lighting]: phong" << endl; }
		else if (event.key.keysym.sym == SDLK_KP_8) light.z -= 0.1;
		else if (event.key.keysym.sym == SDLK_KP_2) light.z += 0.1;
		else if (event.key.keysym.sym == SDLK_KP_6) light.x += 0.1;
		else if (event.key.keysym.sym == SDLK_KP_4) light.x -= 0.1;
		else if (event.key.keysym.sym == SDLK_KP_MINUS) light.y -= 0.1;
		else if (event.key.keysym.sym == SDLK_KP_PLUS) light.y += 0.1;

		else if (event.key.keysym.sym == SDLK_v)  { proximity = (proximity) ? false : true; cout << "[Proximity]: " << proximity << endl; }
		else if (event.key.keysym.sym == SDLK_b) { angleOfInc  = (angleOfInc)  ? false : true; cout << "[Angle of Inc]: " << angleOfInc << endl; }
		else if (event.key.keysym.sym == SDLK_n)        { specular  = (specular)  ? false : true; cout << "[Specular]: " << specular << endl; }
        else if (event.key.keysym.sym == SDLK_m)         { softShadows   = (softShadows)   ? false : true; cout << "[Shadows]: " << softShadows << endl; }
	}
}

int main(int argc, char *argv[]) {
    float focal = 2.8;
    int planemultiplyer = 150;
    vector<glm::vec3> lightDirections;
    lightDirections.push_back(light);
    for(int n=0; n < 7; n++){
        double i = (n+1)*0.02;
        lightDirections.push_back(light + vec3(-i,-i,i));
        lightDirections.push_back(light + vec3(-i,-i,-i));
        lightDirections.push_back(light + vec3(-i,i,i));
        lightDirections.push_back(light + vec3(-i,i,-i));
        lightDirections.push_back(light + vec3(i,-i,i));
        lightDirections.push_back(light + vec3(-i,-i,i));
        lightDirections.push_back(light + vec3(i,i,-i));
        lightDirections.push_back(light + vec3(i,i,i));
    }

	unordered_map<string, TextureMap> textures;
	vector<ModelTriangle> t = load_obj("models/cornell-box-mirror.obj", 0.5, load_mtl("models/cornell-box-mirror.mtl", textures));

    // vector<ModelTriangle> t = load_obj("models/cornell-box-bunny.obj", 0.5, load_mtl("models/cornell-box.mtl",textures));
	// vector<ModelTriangle> t_sphere = load_obj("models/newestsphere.obj", 0.5, load_mtl("models/cornell-box.mtl",textures));
	// vector<ModelTriangle> t = load_obj("models/textured-cornell-box.obj", 0.5, load_mtl("models/textured-cornell-box.mtl", textures));
	// vector<ModelTriangle> t_logo = load_obj("models/logo.obj", 0.002, load_mtl("models/materials.mtl", textures));
	// vector<ModelTriangle> t_sphere = load_obj("models/high-res-sphere.obj", 0.4, load_mtl("models/cornell-box.mtl",textures));

	// t.insert(t.end(),t_logo.begin(), t_logo.end());
	// t.insert(t.end(),t_sphere.begin(), t_sphere.end());
    
    cout <<""<<endl;
    cout << "[Rendering Mode]: raytrace"  << endl;
    cout << "[Lighting]: phong" << endl;
    cout << "[Proximity]: " << proximity << endl;
    cout << "[Angle of Inc]: " << angleOfInc << endl;
    cout << "[Specular]: " << specular << endl;
    cout << "[Soft Shadows]: " << softShadows << endl;
    
    cout <<""<<endl;
    cout <<"Settings:"<<endl;
    cout <<"1 - Wireframe, 2 - Rasterise, 3 - Raytrace"<<endl;
    cout <<"4 - Diffused, 5 - Gourard, 6 - Phong"<<endl;
    cout <<"v - Proximity, b - Angle Of Inc, n - specular"<<endl;
    cout <<"m - Soft shadow"<<endl;
    cout <<""<<endl;

	DrawingWindow window_grey = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window_grey.pollForInputEvents(event)) handleEvent(event, window_grey);
		orbit(orbiting);

		draw(window_grey);
        drawing(t, focal, planemultiplyer, textures, lightDirections, window_grey);
		// drawing(t_logo, focal, planemultiplyer, textures, lightDirections, window_grey); //logo only
		// drawing(t_sphere, focal, planemultiplyer, textures, lightDirections, window_grey); //sphere only

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window_grey.renderFrame();
	}
}