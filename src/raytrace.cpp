#include <iostream>
#include <tira/parser.h>
#include <tira/image.h>
#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <vector>
#include <string>
#include <chrono>
using namespace std;



struct Ray {
    glm::vec3 origin;
    glm::vec3 direction;

    Ray(const glm::vec3& o, const glm::vec3& d) : origin(o), direction(glm::normalize(d)) {}

    glm::vec3 at(float t) const {
        return origin + t * direction;
    }
};

struct Hit {
    float t;
    glm::vec3 color;
    glm::vec3 normal;
    glm::vec3 center;

    Hit() : t(0), color(0.0f) {}
};

struct Light {
    glm::vec3 position;
    glm::vec3 color;
};

struct Sphere {
    glm::vec3 center;
    float radius;
    glm::vec3 color;

    Sphere(const glm::vec3& c, float r, const glm::vec3& col) : center(c), radius(r), color(col) {}

    bool intersect(const Ray& ray, Hit& hit) const {
        glm::vec3 oc = ray.origin - center;
        float a = glm::dot(ray.direction, ray.direction);
        float b = 2.0 * glm::dot(oc, ray.direction);
        float c = glm::dot(oc, oc) - radius * radius;
        float discriminant = b * b - 4 * a * c;

        if (discriminant > 0) {
            float t = (-b - sqrt(discriminant)) / (2.0 * a);
            if (t < 0) t = (-b + sqrt(discriminant)) / (2.0 * a);

            if (t >= 0) {
                hit.t = t;
                hit.color = color;
                glm::vec3 point = ray.at(t);
                hit.normal = glm::normalize(point - center);
                hit.center = center;
                return true;
            }
        }
        return false;
    }
};

struct Triangle {
    glm::vec3 v0, v1, v2;
    glm::vec3 normal;
    glm::vec3 color;
	Sphere bounding_sphere;

    Triangle(const glm::vec3& _v0, const glm::vec3& _v1, const glm::vec3& _v2, const glm::vec3& col)
        : v0(_v0), v1(_v1), v2(_v2), color(col), bounding_sphere(glm::vec3(0.0f), 0.0f, col) {
        normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
		compute_bounding_sphere();

    }

    void compute_bounding_sphere() {
        glm::vec3 center = (v0 + v1 + v2) / 3.0f;
        float radius = glm::length(glm::max(glm::max(glm::distance(center, v0), glm::distance(center, v1)), glm::distance(center, v2)));
        bounding_sphere = Sphere(center, radius, color);
    }


    bool intersect(const Ray& ray, Hit& hit) const {
		// Step 1: Check if the ray intersects the bounding sphere
        if (!bounding_sphere.intersect(ray, hit)) {
			return false;
		}

        // Step 2: Check if the ray is pointing away from the plane
		float denom = glm::dot(ray.direction, normal);
        if (denom > 0) {
            return false;
        }

		// Step 3: Check if the ray intersects the triangle
        const float EPSILON = 1e-8;
        glm::vec3 edge1 = v1 - v0;
        glm::vec3 edge2 = v2 - v0;
        glm::vec3 h = glm::cross(ray.direction, edge2);
        float a = glm::dot(edge1, h);
        if (fabs(a) < EPSILON) {
            return false;
        }
        float f = 1.0f / a;
        glm::vec3 s = ray.origin - v0;
        float u = f * glm::dot(s, h);
        if (u < 0.0 || u > 1.0) {
            return false;
        }
        glm::vec3 q = glm::cross(s, edge1);
        float v = f * glm::dot(ray.direction, q);
        if (v < 0.0 || u + v > 1.0) {
            return false;
        }
        float t = f * glm::dot(edge2, q);
        if (t > EPSILON) {
            hit.t = t;
            hit.normal = normal;
            hit.color = color;
            return true;
        }
        else {
            return false;
        }
    }
};


class Camera {
public:
    glm::vec3 position;
    glm::vec3 look;
    glm::vec3 up;
    float fov;

    Camera(const glm::vec3& pos, const glm::vec3& look_dir, const glm::vec3& up_dir, float field_of_view)
        : position(pos), look(glm::normalize(look_dir - pos)), up(glm::normalize(up_dir)), fov(field_of_view) {}

    Ray get_ray(float u, float v) const {
        float aspect_ratio = 1.0f;
        float theta = glm::radians(fov);
        float h = glm::tan(theta / 2.0f);
        float viewport_height = 2.0f * h;
        float viewport_width = aspect_ratio * viewport_height;

        glm::vec3 w = glm::normalize(position - look);
        glm::vec3 u_vec = glm::normalize(glm::cross(up, w));
        glm::vec3 v_vec = glm::cross(w, u_vec);

        glm::vec3 lower_left = position - (viewport_width / 2.0f) * u_vec - (viewport_height / 2.0f) * v_vec - w;
        glm::vec3 horizontal = viewport_width * u_vec;
        glm::vec3 vertical = viewport_height * v_vec;

        glm::vec3 direction = glm::normalize(lower_left + u * horizontal + v * vertical - position);
        return Ray(position, direction);
    }
};

class Plane {
public:
    glm::vec3 point;
    glm::vec3 normal;
    glm::vec3 color;

    Plane(const glm::vec3& p, const glm::vec3& n, const glm::vec3& col)
        : point(p), normal(glm::normalize(n)), color(col) {}

    bool intersect(const Ray& ray, Hit& hit) const {
        float denom = glm::dot(ray.direction, normal);
        if (fabs(denom) > 1e-6) {
            float t = glm::dot(point - ray.origin, normal) / denom;
            if (t >= 0) {
                hit.t = t;
                hit.color = color;
                hit.normal = normal;
                return true;
            }
        }
        return false;
    }
};


bool load_scene(const std::string& scene_file, std::vector<Sphere>& spheres, std::vector<Plane>& planes, Camera& camera, std::vector<Light>& lights, glm::vec3& background) {
    tira::parser scene_parser(scene_file);

    // Load camera settings
    camera.position = glm::vec3(scene_parser.get<float>("camera_position", 0, 0),
        scene_parser.get<float>("camera_position", 0, 1),
        scene_parser.get<float>("camera_position", 0, 2));

    camera.look = glm::vec3(scene_parser.get<float>("camera_look", 0, 0),
        scene_parser.get<float>("camera_look", 0, 1),
        scene_parser.get<float>("camera_look", 0, 2));

    camera.up = glm::vec3(scene_parser.get<float>("camera_up", 0, 0),
        scene_parser.get<float>("camera_up", 0, 1),
        scene_parser.get<float>("camera_up", 0, 2));

    camera.fov = scene_parser.get<float>("camera_fov", 0);

    // Load spheres
    size_t num_spheres = scene_parser.count("sphere");
    for (size_t i = 0; i < num_spheres; ++i) {
        float radius = scene_parser.get<float>("sphere", i, 0);
        float x = scene_parser.get<float>("sphere", i, 1);
        float y = scene_parser.get<float>("sphere", i, 2);
        float z = scene_parser.get<float>("sphere", i, 3);
        float r = scene_parser.get<float>("sphere", i, 4);
        float g = scene_parser.get<float>("sphere", i, 5);
        float b = scene_parser.get<float>("sphere", i, 6);

        spheres.push_back(Sphere(glm::vec3(x, y, z), radius, glm::vec3(r, g, b)));
    }

    // Load lights
    size_t num_lights = scene_parser.count("light");
    for (size_t i = 0; i < num_lights; ++i) {
        float x = scene_parser.get<float>("light", i, 0);
        float y = scene_parser.get<float>("light", i, 1);
        float z = scene_parser.get<float>("light", i, 2);
        float r = scene_parser.get<float>("light", i, 3);
        float g = scene_parser.get<float>("light", i, 4);
        float b = scene_parser.get<float>("light", i, 5);
        lights.push_back(Light{ glm::vec3(x, y, z), glm::vec3(r, g, b) });
    }

    // Load background color
    background = glm::vec3(scene_parser.get<float>("background", 0, 0),
        scene_parser.get<float>("background", 0, 1),
        scene_parser.get<float>("background", 0, 2));


    // Load planes
    size_t num_planes = scene_parser.count("plane");
    for (size_t i = 0; i < num_planes; ++i) {
        float px = scene_parser.get<float>("plane", i, 0);
        float py = scene_parser.get<float>("plane", i, 1);
        float pz = scene_parser.get<float>("plane", i, 2);
        glm::vec3 point(px, py, pz);

        float nx = scene_parser.get<float>("plane", i, 3);
        float ny = scene_parser.get<float>("plane", i, 4);
        float nz = scene_parser.get<float>("plane", i, 5);
        glm::vec3 normal(nx, ny, nz);

        float r = scene_parser.get<float>("plane", i, 6);
        float g = scene_parser.get<float>("plane", i, 7);
        float b = scene_parser.get<float>("plane", i, 8);
        glm::vec3 color(r, g, b);

        planes.push_back(Plane(point, normal, color));
    

    }



    return true;
}

bool load_obj(const std::string& file_path, std::vector<Triangle>& triangles, const glm::vec3& color) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Failed to open OBJ file: " << file_path << std::endl;
        return false;
    }

    std::vector<glm::vec3> vertices;
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;
        if (prefix == "v") {
            glm::vec3 vertex;
            iss >> vertex.x >> vertex.y >> vertex.z;
            vertices.push_back(vertex);
        }
        else if (prefix == "f") {
            int idx0, idx1, idx2;
            iss >> idx0 >> idx1 >> idx2;
            triangles.push_back(Triangle(vertices[idx0 - 1], vertices[idx1 - 1], vertices[idx2 - 1], color));
        }
    }
    return true;
}


bool check_shadow(const glm::vec3& point, const glm::vec3& normal, const Light& light, const std::vector<Sphere>& spheres,const std::vector<Triangle>& triangles) {
    const float epsilon = 0.001f;
    glm::vec3 light_dir = glm::normalize(light.position - point);
    glm::vec3 shadow_origin = point + epsilon * normal;
    Ray shadow_ray(shadow_origin, light_dir);
    float light_distance = glm::length(light.position - point);

    // Check for intersection with spheres
    for (const Sphere& sphere : spheres) {
        Hit hit_record;
        if (sphere.intersect(shadow_ray, hit_record)) {
            if (hit_record.t < light_distance) {
                return true;
            }
        }
    }


    // Check for intersection with triangles
    for (const Triangle& triangle : triangles) {
        Hit hit_record;
        if (triangle.intersect(shadow_ray, hit_record)) {
            if (hit_record.t < light_distance) {
                return true;
            }
        }
    }

    return false;
}

glm::vec3 compute_lighting(const glm::vec3& point, const glm::vec3& normal, const std::vector<Light>& lights, const std::vector<Sphere>& spheres, const std::vector<Triangle>& triangles) {
    glm::vec3 light_color(0.0f);
    for (const auto& light : lights) {
        if (check_shadow(point, normal, light, spheres, triangles)) {
            continue;
        }
        glm::vec3 light_dir = glm::normalize(light.position - point);
        float diffuse_intensity = std::max(glm::dot(normal, light_dir), 0.0f);
        light_color += diffuse_intensity * light.color;
    }
    return glm::clamp(light_color, 0.0f, 1.0f);
}




string get_file_extension(const std::string& filename) {
    size_t dot_pos = filename.find_last_of(".");
    if (dot_pos == std::string::npos) return ""; // No extension
    return filename.substr(dot_pos + 1);
}

string output_name(string filename) {
    string name = "";
    for (int i = 0; i < filename.length(); i++) {
        if (filename.at(i) != '.') {
            name += filename.at(i);
        }
        else {
            break;
        }
    }
    return name;
}

tira::image<unsigned char> render_image(int height, int width, const Camera& camera, const glm::vec3& background, const std::vector<Sphere>& spheres, const std::vector<Plane>& planes, const std::vector<Triangle>& triangles,const std::vector<Light>& lights) {
    tira::image<unsigned char> image(width, height, 3);

    for (int j = 0; j < height; ++j) {
        std::cout << "Rendering row " << j + 1 << " out of " << height << std::endl;

        for (int i = 0; i < width; ++i) {
            float u = float(width - 1 - i) / float(width - 1);
            float v = float(height - 1 - j) / float(height - 1);
            
            Ray ray = camera.get_ray(u, v);
            glm::vec3 pixel_color = background;

            Hit closest_hit_record;
            bool hit_boolean = false;
            float closest_t = std::numeric_limits<float>::max();

            for (const Sphere& sphere : spheres) {
                Hit hit_record;
                if (sphere.intersect(ray, hit_record) && hit_record.t < closest_t) {
                    closest_t = hit_record.t;
                    closest_hit_record = hit_record;
                    hit_boolean = true;
                   
                }
            }

            for (const Plane& plane : planes) {
                Hit hit_record;
                if (plane.intersect(ray, hit_record) && hit_record.t < closest_t) {
                    closest_t = hit_record.t;
                    closest_hit_record = hit_record;
                    hit_boolean = true;
            
                }
            }


            for (const Triangle& triangle : triangles) {
                Hit hit_record;
                if (triangle.intersect(ray, hit_record) && hit_record.t < closest_t) {
                    closest_t = hit_record.t;
                    closest_hit_record = hit_record;
                    hit_boolean = true;
                    
                }
            }
            

            if (hit_boolean) {
                glm::vec3 point = ray.at(closest_hit_record.t);
                glm::vec3 lighting = compute_lighting(point, closest_hit_record.normal, lights, spheres, triangles);
                pixel_color = glm::clamp(closest_hit_record.color * lighting, 0.0f, 1.0f);
            }
            else {
                pixel_color = background;
            }

            image(i, j, 0) = static_cast<unsigned char>(255 * pixel_color.r);
            image(i, j, 1) = static_cast<unsigned char>(255 * pixel_color.g);
            image(i, j, 2) = static_cast<unsigned char>(255 * pixel_color.b);


        }
    }
    return image;
}



int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 2) {
        cout << "Error: there was no filename. Use the following command: ./program <filename>" << endl;
        return -1;
    }

    vector<Sphere> spheres; vector<Plane> planes; vector<Triangle> triangles; vector<Light> lights;
    Camera camera(glm::vec3(0, 0, 0), glm::vec3(0, 0, 0), glm::vec3(0, 0, 0), 0.0f);
    glm::vec3 background(0.0f, 0.0f, 0.0f);
    string filename = argv[1]; string file_type = "scene";

    if (!load_scene(filename, spheres, planes, camera, lights, background)) {
        cout << "Error loading scene!\t" << filename << endl;
        return -1;
    }
    
    if (filename == "mesh.scene") {
        if (!load_obj("subdiv.obj", triangles, glm::vec3(1.0, 1.0, 1.0))) {
            cout << "Error loading OBJ file!\t" << filename << endl;
            return -1;
		}
	}

    tira::parser p(filename);

    int height = p.get<int>("resolution", 0); int width = p.get<int>("resolution", 1);
    
	// for testing, 100x100 resolution
    //int height = 100; int width = 100;

    tira::image<unsigned char> image(height, width, 3);

    auto parse_start = chrono::high_resolution_clock::now();


    auto parse_end = chrono::high_resolution_clock::now();
    chrono::duration<double> parse_elapsed = parse_end - parse_start;

    auto image_start = chrono::high_resolution_clock::now();

    image = render_image(height, width, camera, background, spheres, planes, triangles, lights);

    auto image_end = chrono::high_resolution_clock::now();
    chrono::duration<double> image_elapsed = image_end - image_start;

    cout << "Execution time of parsing scene: " << parse_elapsed << " seconds" << endl;

    cout << "Execution time of rendering image: " << image_elapsed << " seconds" << endl;
    
    string output = output_name(filename) + ".bmp";

    image.save(output);

    return 0;
}