#define _CRT_SECURE_NO_WARNINGS
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "simulation.h"
#include "particle_type.h"
#include "config.h"

using namespace std;
namespace fs = std::filesystem;

namespace {
constexpr double EPS = 1e-9;

string trim(const string& s) {
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == string::npos) return "";
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

string lower_copy(string s) {
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

string remove_comment(const string& line) {
    size_t pos = line.find('#');
    return trim(pos == string::npos ? line : line.substr(0, pos));
}

bool nearly_same(const Point& a, const Point& b) {
    return PointMag(a - b) < EPS;
}

double cross(const Point& a, const Point& b) {
    return a.x * b.y - a.y * b.x;
}

double signed_polygon_area(const vector<Point>& vertices) {
    double sum = 0.0;
    for (size_t i = 0; i < vertices.size(); ++i) {
        size_t j = (i + 1) % vertices.size();
        sum += cross(vertices[i], vertices[j]);
    }
    return 0.5 * sum;
}

Point polygon_centroid(const vector<Point>& vertices, double signed_area) {
    double cx = 0.0;
    double cy = 0.0;
    for (size_t i = 0; i < vertices.size(); ++i) {
        size_t j = (i + 1) % vertices.size();
        double cr = cross(vertices[i], vertices[j]);
        cx += (vertices[i].x + vertices[j].x) * cr;
        cy += (vertices[i].y + vertices[j].y) * cr;
    }
    double factor = 1.0 / (6.0 * signed_area);
    return Point(cx * factor, cy * factor);
}

bool is_convex_polygon(const vector<Point>& vertices) {
    if (vertices.size() < 3) return false;
    int sign = 0;
    const size_t n = vertices.size();
    for (size_t i = 0; i < n; ++i) {
        Point e1 = vertices[(i + 1) % n] - vertices[i];
        Point e2 = vertices[(i + 2) % n] - vertices[(i + 1) % n];
        double cr = cross(e1, e2);
        if (fabs(cr) <= EPS) continue;
        int s = cr > 0 ? 1 : -1;
        if (sign == 0) sign = s;
        else if (sign != s) return false;
    }
    return sign != 0;
}

Point parse_point_line(string line, const fs::path& path, int line_no) {
    replace(line.begin(), line.end(), ',', ' ');
    istringstream iss(line);
    double x, y;
    if (!(iss >> x >> y)) {
        cerr << "Invalid point in " << path << " at line " << line_no << ": " << line << endl;
        exit(1);
    }
    return Point(x, y);
}

ParticleType load_polygon_shape(const fs::path& shapes_dir, int shape_id) {
    fs::path shape_file = shapes_dir / (to_string(shape_id) + ".txt");
    ifstream file(shape_file);
    if (!file.is_open()) {
        cerr << "Failed to open shape file: " << shape_file << endl;
        exit(1);
    }

    vector<Point> vertices;
    string line;
    int line_no = 0;
    while (getline(file, line)) {
        line_no++;
        line = remove_comment(line);
        if (line.empty()) continue;
        vertices.push_back(parse_point_line(line, shape_file, line_no));
    }

    if (vertices.size() >= 2 && nearly_same(vertices.front(), vertices.back())) {
        vertices.pop_back();
    }
    if (vertices.size() < 3) {
        cerr << "Shape " << shape_file << " must contain at least 3 distinct vertices." << endl;
        exit(1);
    }
    if (!is_convex_polygon(vertices)) {
        cerr << "Shape " << shape_file << " is not a strictly convex polygon." << endl;
        exit(1);
    }

    double signed_area = signed_polygon_area(vertices);
    if (fabs(signed_area) <= EPS) {
        cerr << "Shape " << shape_file << " has near-zero area." << endl;
        exit(1);
    }

    // The input is expected to be clockwise.  If it is counterclockwise,
    // normalize it to clockwise so outward normals are consistent.
    if (signed_area > 0.0) {
        reverse(vertices.begin(), vertices.end());
        signed_area = signed_polygon_area(vertices);
    }

    Point centroid = polygon_centroid(vertices, signed_area);
    vector<Point> local_vertices;
    local_vertices.reserve(vertices.size());
    for (const auto& v : vertices) local_vertices.emplace_back(v.x - centroid.x, v.y - centroid.y);

    double outer_r = 0.0;
    for (const auto& v : local_vertices) outer_r = max(outer_r, PointMag(v));

    double inner_r = numeric_limits<double>::max();
    Point origin(0.0, 0.0);
    for (size_t i = 0; i < local_vertices.size(); ++i) {
        size_t j = (i + 1) % local_vertices.size();
        inner_r = min(inner_r, distanceFromPointToSegment(origin, local_vertices[i], local_vertices[j]));
    }

    ParticleType t;
    t.shape = ParticleShape::POLYGON;
    t.shape_id = shape_id;
    t.radius = outer_r;
    t.vertices = local_vertices;
    t.area = fabs(signed_area);
    t.inner_r = inner_r;
    t.character_L = 2.2 * outer_r;
    t.label = "S" + to_string(shape_id);
    return t;
}

ParticleType make_circle_type(double radius) {
    if (radius <= 0.0) {
        cerr << "Circle radius must be positive." << endl;
        exit(1);
    }
    ParticleType t;
    t.shape = ParticleShape::CIRCLE;
    t.shape_id = 0;
    t.radius = radius;
    t.area = PI * radius * radius;
    t.inner_r = radius;
    t.character_L = 2.2 * radius;
    ostringstream label;
    label << "C" << radius;
    t.label = label.str();
    return t;
}

// Supported config lines:
//   circle  <radius>   [probability]
//   polygon <shape_id> [probability]
// Also accepts c/circ and p/poly.  If all probabilities are omitted, equal
// probabilities are assigned.  If any probability is provided, all must be provided.
vector<ParticleType> read_particle_config(const fs::path& config_path, const fs::path& shapes_dir) {
    vector<ParticleType> types;
    map<int, ParticleType> polygon_cache;
    vector<bool> has_prob;

    ifstream file(config_path);
    if (!file.is_open()) {
        cerr << "Failed to open config file: " << config_path << endl;
        exit(1);
    }

    string line;
    int line_no = 0;
    while (getline(file, line)) {
        line_no++;
        line = remove_comment(line);
        if (line.empty()) continue;

        istringstream iss(line);
        string shape_str;
        iss >> shape_str;
        shape_str = lower_copy(shape_str);

        ParticleType t;
        if (shape_str == "circle" || shape_str == "circ" || shape_str == "c") {
            double radius;
            if (!(iss >> radius)) {
                cerr << "Missing circle radius in " << config_path << " at line " << line_no << endl;
                exit(1);
            }
            t = make_circle_type(radius);
        } else if (shape_str == "polygon" || shape_str == "poly" || shape_str == "p") {
            int shape_id;
            if (!(iss >> shape_id)) {
                cerr << "Missing polygon shape id in " << config_path << " at line " << line_no << endl;
                exit(1);
            }
            auto it = polygon_cache.find(shape_id);
            if (it == polygon_cache.end()) {
                it = polygon_cache.emplace(shape_id, load_polygon_shape(shapes_dir, shape_id)).first;
            }
            t = it->second;
        } else {
            cerr << "Invalid particle type in " << config_path << " at line " << line_no << ": " << shape_str << endl;
            exit(1);
        }

        double prob;
        if (iss >> prob) {
            if (prob < 0.0) {
                cerr << "Probability must be non-negative in " << config_path << " at line " << line_no << endl;
                exit(1);
            }
            t.prob = prob;
            has_prob.push_back(true);
        } else {
            t.prob = 0.0;
            has_prob.push_back(false);
        }
        types.push_back(t);
    }

    if (types.empty()) {
        cerr << "Config file has no particle types: " << config_path << endl;
        exit(1);
    }

    bool any_prob = any_of(has_prob.begin(), has_prob.end(), [](bool v) { return v; });
    bool all_prob = all_of(has_prob.begin(), has_prob.end(), [](bool v) { return v; });
    if (!any_prob) {
        double p = 1.0 / static_cast<double>(types.size());
        for (auto& t : types) t.prob = p;
    } else if (!all_prob) {
        cerr << "In " << config_path << ", either provide probabilities on every line or omit all of them." << endl;
        exit(1);
    } else {
        double total_prob = 0.0;
        for (const auto& t : types) total_prob += t.prob;
        if (fabs(total_prob - 1.0) > 1e-6) {
            cerr << "Error in " << config_path << ": total probability is " << total_prob << " (should be 1.0)." << endl;
            exit(1);
        }
    }

    return types;
}

string safe_component(string s) {
    for (char& ch : s) {
        if (!(isalnum(static_cast<unsigned char>(ch)) || ch == '-' || ch == '_')) ch = '_';
    }
    return s;
}

string generate_folder_name(const fs::path& config_path, const vector<ParticleType>& types) {
    ostringstream name;
    name << safe_component(config_path.stem().string());
    for (const auto& t : types) {
        name << "_" << t.label << "_p" << fixed << setprecision(3) << t.prob;
    }
    return safe_component(name.str());
}

fs::path resolve_config_path(const string& config_arg) {
    fs::path direct(config_arg);
    const bool has_separator = config_arg.find('/') != string::npos || config_arg.find('\\') != string::npos;

    // Explicit paths are honored directly. If only "1.txt" is provided and it
    // is not in the current directory, fall back to conf/1.txt.
    if (has_separator || direct.has_extension()) {
        if (fs::exists(direct)) return direct;
        if (!direct.has_parent_path()) {
            fs::path in_conf = fs::path("conf") / direct;
            if (fs::exists(in_conf)) return in_conf;
        }
        return direct;
    }

    // Required usage: .\main 1 means read conf/1.txt.
    return fs::path("conf") / (config_arg + ".txt");
}
}

int main(int argc, char* argv[]) {
    srand(static_cast<unsigned int>(time(NULL)));

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <conf_id|conf_path> [shapes_dir]" << endl;
        cerr << "Example: " << argv[0] << " 1   # reads conf/1.txt" << endl;
        return 1;
    }

    fs::path config_path = resolve_config_path(argv[1]);
    fs::path shapes_dir = (argc >= 3) ? fs::path(argv[2]) : fs::path("shapes");

    if (!fs::is_regular_file(config_path)) {
        cerr << "Config file does not exist: " << config_path << endl;
        return 1;
    }
    if (!fs::is_directory(shapes_dir)) {
        cerr << "Shapes directory does not exist: " << shapes_dir << endl;
        return 1;
    }

    const long long area = static_cast<long long>(AREA);
    const long long rsa_step = static_cast<long long>(RSA_STEP_COEFFICIENT * AREA);
    const int cycle = CYCLE_NUM;

    if (area <= 0 || rsa_step <= 0 || cycle <= 0) {
        cerr << "Invalid values in config.h. AREA, RSA_STEP_COEFFICIENT, and CYCLE_NUM must be positive." << endl;
        return 1;
    }

    fs::path root = fs::current_path();
    ofstream outfile(root / "results.txt", ios_base::app);
    if (!outfile.is_open()) {
        cerr << "Unable to open results.txt" << endl;
        return 1;
    }
    outfile << fixed << setprecision(8);

    cout << "Running config: " << config_path << endl;
    cout << "AREA=" << AREA
         << ", RSA_STEP_COEFFICIENT=" << RSA_STEP_COEFFICIENT
         << ", FINAL_ADDNUM_MIN=" << FINAL_ADDNUM_MIN
         << ", CYCLE_NUM=" << CYCLE_NUM << endl;

    vector<ParticleType> particle_types = read_particle_config(config_path, shapes_dir);

    string folder_name = generate_folder_name(config_path, particle_types);
    fs::create_directories(root / folder_name);
    fs::current_path(root / folder_name);

    RsaSummary summary = rsa_cycle(particle_types, area, rsa_step, cycle);

    fs::current_path(root);
    outfile << config_path.filename().string() << "\t";
    for (const auto& t : particle_types) outfile << t.label << "\t" << t.prob << "\t";
    outfile << "avg_counts\t";
    for (double val : summary.avg_counts) outfile << val << "\t";
    outfile << "avg_packing_fraction\t" << summary.avg_packing_fraction << endl;

    outfile.close();
    return 0;
}
