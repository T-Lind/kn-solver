#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <random>

using namespace std;

struct VectorOperations {
    static vector<vector<double>> normalize(const vector<vector<double>>& vectors) {
        vector<vector<double>> normalized(vectors.size(), vector<double>(vectors[0].size()));
        for (size_t i = 0; i < vectors.size(); ++i) {
            double norm = 0.0;
            for (double val : vectors[i]) {
                norm += val * val;
            }
            norm = sqrt(norm);
            for (size_t j = 0; j < vectors[i].size(); ++j) {
                normalized[i][j] = vectors[i][j] / norm;
            }
        }
        return normalized;
    }

    static vector<vector<double>> calculate_distances(const vector<vector<double>>& positions) {
        size_t n = positions.size();
        vector<vector<double>> dist_matrix(n, vector<double>(n, 0.0));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                double dist = 0.0;
                for (size_t k = 0; k < positions[i].size(); ++k) {
                    double diff = positions[i][k] - positions[j][k];
                    dist += diff * diff;
                }
                dist = sqrt(dist);
                dist_matrix[i][j] = dist;
                dist_matrix[j][i] = dist;
            }
        }
        return dist_matrix;
    }
};

class SphereDistributor {
public:
    SphereDistributor(int num_points, int dimensions, int iterations = 10000, double learning_rate = 0.01)
        : num_points(num_points), dimensions(dimensions), iterations(iterations), learning_rate(learning_rate) {}

    vector<vector<double>> distribute_points() {
        random_device rd;
        mt19937 gen(rd());
        normal_distribution<> d(0, 1);

        vector<vector<double>> points(num_points, vector<double>(dimensions));
        for (int i = 0; i < num_points; ++i) {
            for (int j = 0; j < dimensions; ++j) {
                points[i][j] = d(gen);
            }
        }
        points = VectorOperations::normalize(points);

        int update_interval = iterations / 20;

        for (int iter = 0; iter < iterations; ++iter) {
            if (iter % update_interval == 0) {
                cout << "Progress: " << (iter * 100 / iterations) << "%\n";
            }

            vector<vector<double>> forces = calculate_repulsive_forces(points);
            for (int i = 1; i < num_points; ++i) {
                for (int j = 0; j < dimensions; ++j) {
                    points[i][j] += learning_rate * forces[i][j];
                }
            }
            points = VectorOperations::normalize(points);
        }

        return points;
    }

private:
    int num_points;
    int dimensions;
    int iterations;
    double learning_rate;

    vector<vector<double>> calculate_repulsive_forces(const vector<vector<double>>& positions) {
        size_t n = positions.size();
        size_t d = positions[0].size();
        vector<vector<double>> forces(n, vector<double>(d, 0.0));
        vector<vector<double>> dist_matrix = VectorOperations::calculate_distances(positions);

        for (size_t i = 0; i < n; ++i) {
            dist_matrix[i][i] = numeric_limits<double>::infinity();
        }

        for (size_t i = 1; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                vector<double> diff(d);
                for (size_t k = 0; k < d; ++k) {
                    diff[k] = positions[i][k] - positions[j][k];
                }
                double dist = dist_matrix[i][j];
                double force_magnitude = 1.0 / (dist * dist);
                for (size_t k = 0; k < d; ++k) {
                    double force_vector = force_magnitude * diff[k] / dist;
                    forces[i][k] += force_vector;
                    forces[j][k] -= force_vector;
                }
            }
        }

        return forces;
    }
};

struct SphereOperations {
    static vector<vector<double>> generate_outer_points(const vector<vector<double>>& inner_points, double scale_factor) {
        vector<vector<double>> outer_points(inner_points.size(), vector<double>(inner_points[0].size()));
        for (size_t i = 0; i < inner_points.size(); ++i) {
            for (size_t j = 0; j < inner_points[i].size(); ++j) {
                outer_points[i][j] = inner_points[i][j] * scale_factor;
            }
        }
        return outer_points;
    }

    static int count_overlapping_spheres(const vector<vector<double>>& points, double radius) {
        int overlap_count = 0;
        for (size_t i = 0; i < points.size(); ++i) {
            for (size_t j = i + 1; j < points.size(); ++j) {
                double dist = 0.0;
                for (size_t k = 0; k < points[i].size(); ++k) {
                    double diff = points[i][k] - points[j][k];
                    dist += diff * diff;
                }
                dist = sqrt(dist);
                if (dist < 2 * radius) {
                    overlap_count++;
                }
            }
        }
        return overlap_count;
    }

    static void print_distance_matrix(const vector<vector<double>>& points) {
        vector<vector<double>> dist_matrix = VectorOperations::calculate_distances(points);

        cout << "Distance Matrix:\n";
        for (const auto& row : dist_matrix) {
            for (double val : row) {
                cout << val << " ";
            }
            cout << endl;
        }
    }

    static void print_points(const vector<vector<double>>& points) {
        for (const auto& point : points) {
            for (double val : point) {
                cout << val << " ";
            }
            cout << endl;
        }
    }
};

int main(int argc, char* argv[]) {
    int num_points = 11;
    int dimensions = 3;
    int n_iterations = 1000000;
    double radius = 1.0;
    double scale_factor = 2.0;

    SphereDistributor distributor(num_points, dimensions, n_iterations);
    vector<vector<double>> points = distributor.distribute_points();

    vector<vector<double>> outer_points = SphereOperations::generate_outer_points(points, scale_factor);

    int overlap_count = SphereOperations::count_overlapping_spheres(outer_points, radius);
    cout << "Number of overlapping spheres: " << overlap_count << endl;

    cout << "Inner Points:" << endl;
    SphereOperations::print_points(points);

//    cout << "Inner Points Distance Matrix:" << endl;
//    SphereOperations::print_distance_matrix(points);
//    cout << "Outer Points Distance Matrix:" << endl;
//    SphereOperations::print_distance_matrix(outer_points);

    return 0;
}