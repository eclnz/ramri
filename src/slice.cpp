#include <Rcpp.h>
#include <queue>
#include <cmath>
#include <chrono>
#include <unordered_set>
using namespace Rcpp;
// using namespace std::chrono;  // Comment out chrono namespace

// Define constants
// Maximum theoretical number of neighbors (all 26 surrounding cells)
const int MAX_THEORETICAL_NEIGHBORS = 26;
// Define a very good neighbor count that would allow early termination
const int EARLY_TERMINATION_THRESHOLD = 18;
// Define the radius for orthogonal plane checking
const int PLANE_CHECK_RADIUS = 10;
// Add hard limit for total surface/trajectory checks to avoid infinite loops
const int MAX_SURFACE_CHECKS = 1000000;

// Helper struct for 3D coordinates with neighbor count
struct Point3D {
    int x, y, z, h;
    int neighbor_count;
    Point3D(int x_, int y_, int z_, int h_, int n_) : 
        x(x_), y(y_), z(z_), h(h_), neighbor_count(n_) {}
    
    // Comparison operator for priority queue (max heap)
    bool operator<(const Point3D& other) const {
        return neighbor_count < other.neighbor_count;
    }
};

// Hash function for points to use in unordered_set
struct PointHash {
    std::size_t operator()(const std::tuple<int, int, int, int>& p) const {
        auto [x, y, z, h] = p;
        // Use a more efficient hash combining function
        // This is a simplified version of Boost's hash_combine
        std::size_t seed = x;
        seed ^= y + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= z + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

// Helper function to check if a point is within bounds
inline bool is_valid_point(int x, int y, int z, int D0, int D1, int D2) {
    return x >= 0 && x < D0 && y >= 0 && y < D1 && z >= 0 && z < D2;
}

// Count neighbors in 3D space including diagonals
inline int count_neighbors_3d(const NumericVector& array, const std::vector<int>& strides,
                            int x, int y, int z, int h, int D0, int D1, int D2,
                            const IntegerVector& dimensions, int array_size, double threshold) {
    int count = 0;
    int ndim = dimensions.length();
    int nonzero_checks = 0;
    int total_checks = 0;
    
    // Check all 26 neighbors (6 face, 12 edge, 8 corner)
    for(int dx = -1; dx <= 1; dx++) {
        for(int dy = -1; dy <= 1; dy++) {
            for(int dz = -1; dz <= 1; dz++) {
                if(dx == 0 && dy == 0 && dz == 0) continue;  // Skip self
                
                int nx = x + dx;
                int ny = y + dy;
                int nz = z + dz;
                
                if(!is_valid_point(nx, ny, nz, D0, D1, D2)) continue;
                
                // For this neighbor position, check if any higher dimension has a non-zero value
                bool has_nonzero = false;
                
                // First check the specified h coordinate since it's most likely non-zero
                int idx_base = nx * strides[0] + ny * strides[1] + nz * strides[2];
                int idx = idx_base;
                
                // Add offset for higher dimensions from the current h
                if(ndim > 3) {
                    int remaining_h = h;
                    for(int d = ndim-1; d >= 3; d--) {
                        int dim_size = dimensions[d];
                        int dim_idx = remaining_h % dim_size;
                        idx += dim_idx * strides[d];
                        remaining_h /= dim_size;
                    }
                    
                    total_checks++;
                    // Check if this index has a non-zero value
                    if(idx >= 0 && idx < array_size && std::abs(array[idx]) > threshold) {
                        has_nonzero = true;
                        nonzero_checks++;
                    } else if(h > 0) {
                        total_checks++;
                        // If not, try h=0 which often has non-zero values
                        idx = idx_base;
                        if(idx >= 0 && idx < array_size && std::abs(array[idx]) > threshold) {
                            has_nonzero = true;
                            nonzero_checks++;
                        }
                    }
                    
                    // If still not found and there are many higher dimensions, we need a more efficient approach
                    if(!has_nonzero && ndim > 4) {
                        // Instead of checking all combinations, sample intelligently
                        // First check the axes-aligned positions in higher dimensions
                        for(int d = 3; d < ndim && !has_nonzero; d++) {
                            for(int val = 0; val < dimensions[d] && !has_nonzero; val++) {
                                total_checks++;
                                // Create an index where only this dimension has value 'val'
                                idx = idx_base;
                                for(int d2 = 3; d2 < ndim; d2++) {
                                    idx += (d2 == d ? val : 0) * strides[d2];
                                }
                                
                                if(idx >= 0 && idx < array_size && std::abs(array[idx]) > threshold) {
                                    has_nonzero = true;
                                    nonzero_checks++;
                                }
                            }
                        }
                    }
                } else {
                    total_checks++;
                    // For 3D arrays, just check the single index
                    if(idx >= 0 && idx < array_size && std::abs(array[idx]) > threshold) {
                        has_nonzero = true;
                        nonzero_checks++;
                    }
                }
                
                if(has_nonzero) count++;
            }
        }
    }
    
    if(count > 0) {
        Rcout << "  Neighbor check at (" << x << "," << y << "," << z << "," << h << "): " 
              << count << " neighbors, checked " << total_checks 
              << " positions, found " << nonzero_checks << " non-zero values\n";
    }
    
    return count;
}

// Function to check orthogonal planes extending from a center point
// Returns the percentage of non-zero voxels in the three orthogonal planes
inline double check_orthogonal_planes(const NumericVector& array, const std::vector<int>& strides,
                                     int x, int y, int z, int h, int D0, int D1, int D2,
                                     const IntegerVector& dimensions, int array_size, int radius, double threshold) {
    int total_points = 0;
    int nonzero_points = 0;
    int ndim = dimensions.length();
    
    // Function to calculate array index
    auto calc_index = [&](int px, int py, int pz) -> int {
        if(!is_valid_point(px, py, pz, D0, D1, D2)) return -1;
        
        int idx = px * strides[0] + py * strides[1] + pz * strides[2];
        int remaining_h = h;
        for(int d = ndim-1; d >= 3; d--) {
            int dim_size = dimensions[d];
            int dim_idx = remaining_h % dim_size;
            idx += dim_idx * strides[d];
            remaining_h /= dim_size;
        }
        
        return (idx >= 0 && idx < array_size) ? idx : -1;
    };
    
    // Check YZ plane (fixed X)
    for(int dy = -radius; dy <= radius; dy++) {
        for(int dz = -radius; dz <= radius; dz++) {
            if(dy == 0 && dz == 0) continue; // Skip center point
            
            int idx = calc_index(x, y + dy, z + dz);
            if(idx >= 0) {
                total_points++;
                if(std::abs(array[idx]) > threshold) {
                    nonzero_points++;
                }
            }
        }
    }
    
    // Check XZ plane (fixed Y)
    for(int dx = -radius; dx <= radius; dx++) {
        for(int dz = -radius; dz <= radius; dz++) {
            if(dx == 0 && dz == 0) continue; // Skip center point
            
            int idx = calc_index(x + dx, y, z + dz);
            if(idx >= 0) {
                total_points++;
                if(std::abs(array[idx]) > threshold) {
                    nonzero_points++;
                }
            }
        }
    }
    
    // Check XY plane (fixed Z)
    for(int dx = -radius; dx <= radius; dx++) {
        for(int dy = -radius; dy <= radius; dy++) {
            if(dx == 0 && dy == 0) continue; // Skip center point
            
            int idx = calc_index(x + dx, y + dy, z);
            if(idx >= 0) {
                total_points++;
                if(std::abs(array[idx]) > threshold) {
                    nonzero_points++;
                }
            }
        }
    }
    
    // Return percentage of non-zero points
    return total_points > 0 ? (double)nonzero_points / total_points : 0.0;
}

// Function to optimize center point based on orthogonal plane checks
Point3D optimize_center_point(const NumericVector& array, const std::vector<int>& strides,
                             Point3D initial_point, int D0, int D1, int D2,
                             const IntegerVector& dimensions, int array_size, int higher_dims_size, double threshold) {
    Point3D best_point = initial_point;
    double best_score = check_orthogonal_planes(array, strides, initial_point.x, initial_point.y, 
                                              initial_point.z, initial_point.h, D0, D1, D2, 
                                              dimensions, array_size, PLANE_CHECK_RADIUS, threshold);
    
    Rcout << "Initial center point (" << initial_point.x << "," << initial_point.y << "," 
          << initial_point.z << ") with orthogonal plane score: " << best_score * 100 << "%\n";
    
    // If score is already perfect (1.0), return immediately
    if(best_score >= 0.99) {
        Rcout << "Center point is optimal with 100% non-zero orthogonal planes\n";
        return best_point;
    }
    
    // Use a more aggressive search approach to find a perfect center
    bool improved = true;
    int iteration = 0;
    const int MAX_ITERATIONS = 200; // Allow more iterations for finding the perfect center
    
    // Try increasingly larger search radii to find a better center
    const int MAX_SEARCH_RADIUS = 10;
    
    while(improved && iteration < MAX_ITERATIONS && best_score < 0.99) {
        improved = false;
        iteration++;
        
        // Start with a small search radius and gradually increase it
        int search_radius = std::min(iteration / 5 + 1, MAX_SEARCH_RADIUS);
        
        // Use a spiral search pattern to explore points more efficiently
        // Create an array of candidate points sorted by distance to current best point
        std::vector<std::tuple<int, int, int, double>> candidates;
        
        for(int dx = -search_radius; dx <= search_radius; dx++) {
            for(int dy = -search_radius; dy <= search_radius; dy++) {
                for(int dz = -search_radius; dz <= search_radius; dz++) {
                    // Skip current point
                    if(dx == 0 && dy == 0 && dz == 0) continue;
                    
                    // Skip points that are too far (use Manhattan distance for efficiency)
                    if(std::abs(dx) + std::abs(dy) + std::abs(dz) > search_radius * 2) continue;
                    
                    int nx = best_point.x + dx;
                    int ny = best_point.y + dy;
                    int nz = best_point.z + dz;
                    
                    if(!is_valid_point(nx, ny, nz, D0, D1, D2)) continue;
                    
                    // Calculate Euclidean distance for sorting
                    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                    candidates.push_back(std::make_tuple(nx, ny, nz, distance));
                }
            }
        }
        
        // Sort candidates by distance (closest first)
        std::sort(candidates.begin(), candidates.end(), 
                 [](const auto& a, const auto& b) { 
                     return std::get<3>(a) < std::get<3>(b); 
                 });
        
        // Try each candidate
        for(const auto& candidate : candidates) {
            int nx = std::get<0>(candidate);
            int ny = std::get<1>(candidate);
            int nz = std::get<2>(candidate);
            
            // Calculate flat index
            int idx = nx * strides[0] + ny * strides[1] + nz * strides[2];
            int remaining_h = best_point.h;
            for(int d = dimensions.length()-1; d >= 3; d--) {
                int dim_size = dimensions[d];
                int dim_idx = remaining_h % dim_size;
                idx += dim_idx * strides[d];
                remaining_h /= dim_size;
            }
            
            // Check if point is non-zero
            if(idx >= 0 && idx < array_size && std::abs(array[idx]) > threshold) {
                double score = check_orthogonal_planes(array, strides, nx, ny, nz, 
                                                     best_point.h, D0, D1, D2, 
                                                     dimensions, array_size, PLANE_CHECK_RADIUS, threshold);
                
                if(score > best_score) {
                    best_score = score;
                    best_point = Point3D(nx, ny, nz, best_point.h, 
                                        // Recalculate neighbor count for new point
                                        count_neighbors_3d(array, strides, nx, ny, nz, 
                                                         best_point.h, D0, D1, D2, 
                                                         dimensions, array_size, threshold));
                    improved = true;
                    
                    Rcout << "Improved center to (" << best_point.x << "," << best_point.y << "," 
                          << best_point.z << ") with orthogonal plane score: " << best_score * 100 << "%\n";
                    
                    // Early termination if we found a perfect center
                    if(best_score >= 0.99) {
                        Rcout << "Found optimal center with 100% non-zero orthogonal planes\n";
                        return best_point;
                    }
                    
                    // Once we find an improvement, break and start a new iteration
                    // with this point as the new center
                    break;
                }
            }
        }
        
        // If we've been searching for a while and still haven't found 100%, try other h dimensions
        if(!improved && iteration > 50 && best_score < 0.95) {
            Rcout << "Trying different h dimensions to find better center...\n";
            
            // Store the best h we find
            int best_h = best_point.h;
            double best_h_score = best_score;
            
            // Try all possible h values
            for(int h = 0; h < higher_dims_size; h++) {
                if(h == best_point.h) continue; // Skip current h
                
                // Calculate index for current point with new h
                int idx = best_point.x * strides[0] + best_point.y * strides[1] + best_point.z * strides[2];
                int remaining_h = h;
                for (int d = dimensions.length()-1; d >= 3; d--) {
                    int dim_size = dimensions[d];
                    int dim_idx = remaining_h % dim_size;
                    idx += dim_idx * strides[d];
                    remaining_h /= dim_size;
                }
                // Check if point is non-zero with this h
                if(idx >= 0 && idx < array_size && std::abs(array[idx]) > threshold) {
                    double score = check_orthogonal_planes(array, strides, best_point.x, best_point.y, best_point.z, 
                                                         h, D0, D1, D2, dimensions, array_size, PLANE_CHECK_RADIUS, threshold);
                    if(score > best_h_score) {
                        best_h_score = score;
                        best_h = h;
                        if(score >= 0.99) {
                            break; // Found a perfect h dimension
                        }
                    }
                }
            }
            
            // If we found a better h dimension, update our best point
            if(best_h != best_point.h) {
                best_score = best_h_score;
                best_point = Point3D(best_point.x, best_point.y, best_point.z, best_h,
                                    count_neighbors_3d(array, strides, best_point.x, best_point.y, best_point.z,
                                                     best_h, D0, D1, D2, dimensions, array_size, threshold));
                improved = true;
                
                Rcout << "Better h dimension found: h=" << best_h << " with score: " << best_score * 100 << "%\n";
                
                if(best_score >= 0.99) {
                    Rcout << "Found optimal center with 100% non-zero orthogonal planes\n";
                    return best_point;
                }
            }
        }
    }
    
    if(iteration >= MAX_ITERATIONS) {
        Rcout << "Reached maximum iterations (" << MAX_ITERATIONS 
              << ") during center optimization\n";
    } else if(!improved) {
        Rcout << "No further improvements possible, best score: " << best_score * 100 << "%\n";
    } else {
        Rcout << "Center optimization complete after " << iteration 
              << " iterations with final score: " << best_score * 100 << "%\n";
    }
    
    return best_point;
}

// [[Rcpp::export]]
IntegerVector extract_nonzero_slices_cpp(NumericVector array, IntegerVector dimensions, double threshold = -1.0) {
    // Uncommented timing code
    auto start_time = std::chrono::high_resolution_clock::now();
    auto last_time = start_time;
    auto log_time = [&](const std::string& message) {
        auto current_time = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - last_time).count();
        auto total = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
        Rcout << "[Time: " << elapsed << "ms, Total: " << total << "ms] " << message << "\n";
        last_time = current_time;
    };
    
    // Input validation
    int ndim = dimensions.length();
    if (ndim < 3) {
        stop("Array must have at least 3 dimensions");
    }
    
    // Get dimensions
    int D0 = dimensions[0];
    int D1 = dimensions[1];
    int D2 = dimensions[2];
    
    int higher_dims_size = 1;
    for (int i = 3; i < ndim; i++) {
        higher_dims_size *= dimensions[i];
    }
    
    // Calculate strides (R uses column-major order)
    std::vector<int> strides(ndim);
    strides[0] = 1;  // First dimension has stride 1 in column-major
    for (int i = 1; i < ndim; i++) {
        strides[i] = strides[i-1] * dimensions[i-1];
    }
    
    int array_size = array.length();
    
    // Calculate threshold by sampling corners if not provided
    if (threshold < 0) {
        double corner_sum = 0.0;
        int corner_count = 0;
        
        // Helper function to calculate flat index for a corner
        auto calc_corner_index = [&](int x, int y, int z, int h) -> int {
            if (!is_valid_point(x, y, z, D0, D1, D2)) return -1;
            
            int idx = x * strides[0] + y * strides[1] + z * strides[2];
            int remaining_h = h;
            for (int d = ndim-1; d >= 3; d--) {
                int dim_size = dimensions[d];
                int dim_idx = remaining_h % dim_size;
                idx += dim_idx * strides[d];
                remaining_h /= dim_size;
            }
            
            return (idx >= 0 && idx < array_size) ? idx : -1;
        };
        
        // Sample all 8 corners across all higher dimensions
        std::vector<int> corner_xs = {0, D0-1};
        std::vector<int> corner_ys = {0, D1-1};
        std::vector<int> corner_zs = {0, D2-1};
        
        for (int x_idx = 0; x_idx < corner_xs.size(); x_idx++) {
            for (int y_idx = 0; y_idx < corner_ys.size(); y_idx++) {
                for (int z_idx = 0; z_idx < corner_zs.size(); z_idx++) {
                    int x = corner_xs[x_idx];
                    int y = corner_ys[y_idx];
                    int z = corner_zs[z_idx];
                    
                    // Sample each corner in all higher dimensions
                    for (int h = 0; h < higher_dims_size; h++) {
                        int idx = calc_corner_index(x, y, z, h);
                        if (idx >= 0) {
                            corner_sum += std::abs(array[idx]);
                            corner_count++;
                        }
                    }
                }
            }
        }
        
        // Calculate threshold as double the average corner value
        if (corner_count > 0) {
            double corner_avg = corner_sum / corner_count;
            threshold = corner_avg * 2.0;
            Rcout << "Automatically calculated threshold from " << corner_count << " corners: " 
                  << threshold << " (2x average of " << corner_avg << ")\n";
        } else {
            // Fallback to default if sampling fails
            threshold = 1.0;
            Rcout << "Corner sampling failed, using default threshold: " << threshold << "\n";
        }
    } else {
        Rcout << "Using provided threshold: " << threshold << "\n";
    }
    
    Rcout << "Array size: " << array_size << "\n";
    Rcout << "Expected size: " << D0 * D1 * D2 * higher_dims_size << "\n";
    Rcout << "Dimensions: ";
    for(int i = 0; i < ndim; i++) Rcout << dimensions[i] << " ";
    Rcout << "\n";
    Rcout << "Strides: ";
    for(int i = 0; i < ndim; i++) Rcout << strides[i] << " ";
    Rcout << "\n";
    
    log_time("Initialization complete");
    
    // Start from center
    int center_x = D0 / 2;
    int center_y = D1 / 2;
    int center_z = D2 / 2;
    
    // Priority queue to store points with their neighbor counts (max heap)
    std::priority_queue<Point3D> pq;
    
    // Replace 4D visited array with a hash-based solution
    // Using unordered_set with a tuple of (x,y,z,h) as key
    // Start with a reasonably small size and let it grow as needed
    const size_t INITIAL_BUCKETS = 128;
    std::unordered_set<std::tuple<int, int, int, int>, PointHash> visited(INITIAL_BUCKETS);
    
    log_time("Data structures initialized");
    
    int center_checks = 0;
    // Start from center
    for(int h = 0; h < higher_dims_size; h++) {
        // Calculate flat index for center point using strides
        int idx = center_x * strides[0] + center_y * strides[1] + center_z * strides[2];
        // Add offset for higher dimensions by unraveling h
        int remaining_h = h;
        for (int d = ndim-1; d >= 3; d--) {
            int dim_size = dimensions[d];
            int dim_idx = remaining_h % dim_size;
            idx += dim_idx * strides[d];
            remaining_h /= dim_size;
        }
        center_checks++;
        // Strict bounds checking
        if(idx >= 0 && idx < array_size) {
            if(std::abs(array[idx]) > threshold) {
                Rcout << "Checking center point at (" << center_x << "," << center_y << "," << center_z << "," << h << ")\n";
                int neighbors = count_neighbors_3d(array, strides, center_x, center_y, center_z, h, D0, D1, D2, dimensions, array_size, threshold);
                pq.push(Point3D(center_x, center_y, center_z, h, neighbors));
                visited.insert(std::make_tuple(center_x, center_y, center_z, h));
                // Early termination - if we find a good center point, stop searching other higher dimensions
                if(neighbors >= EARLY_TERMINATION_THRESHOLD) {
                    Rcout << "Found optimal center point with " << neighbors << " neighbors. Early termination.\n";
                }
                break; // Stop searching h for this center point after first nonzero
            }
        }
    }
    
    log_time("Center point checks complete. Checked " + std::to_string(center_checks) + " center points");
    
    int surface_checks = 0;
    // If center is zero, expand outward in layers until we find a non-zero point
    if(pq.empty()) {
        Rcout << "Center points empty, expanding search...\n";
        int max_radius = std::min({D0/2, D1/2, D2/2});
        bool found_good_point = false;
        
        // After finding a non-zero point, focus search on areas with non-zero data
        for(int radius = 1; radius <= max_radius && !found_good_point; radius++) {
            // Only check a subset of points at larger radii to avoid explosion
            int step = 1;
            if(radius > 50) {
                step = radius / 5;
            } else if(radius > 10) {
                step = 2;
            }
            
            // Track whether we found any non-zero points at this radius
            bool found_nonzero_at_radius = false;
            
            for(int dx = -radius; dx <= radius && !found_good_point; dx += step) {
                for(int dy = -radius; dy <= radius && !found_good_point; dy += step) {
                    for(int dz = -radius; dz <= radius && !found_good_point; dz += step) {
                        // Only check points on the surface of the current radius
                        if(std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) != radius) continue;
                        
                        int x = center_x + dx;
                        int y = center_y + dy;
                        int z = center_z + dz;
                        
                        if(!is_valid_point(x, y, z, D0, D1, D2)) continue;
                        
                        // Look for non-zero points in higher dimensions
                        for(int h = 0; h < higher_dims_size; h++) {
                            // Skip already visited points
                            auto point_tuple = std::make_tuple(x, y, z, h);
                            if(visited.find(point_tuple) != visited.end()) continue;
                            // Calculate index
                            int idx = x * strides[0] + y * strides[1] + z * strides[2];
                            int remaining_h = h;
                            for (int d = ndim-1; d >= 3; d--) {
                                idx += (remaining_h % dimensions[d]) * strides[d];
                                remaining_h /= dimensions[d];
                            }
                            surface_checks++;
                            if(idx < 0 || idx >= array_size) continue;
                            if(std::abs(array[idx]) > threshold) {
                                found_nonzero_at_radius = true;
                                Rcout << "Checking surface point at (" << x << "," << y << "," << z << "," << h << ")\n";
                                int neighbors = count_neighbors_3d(array, strides, x, y, z, h, D0, D1, D2, dimensions, array_size, threshold);
                                pq.push(Point3D(x, y, z, h, neighbors));
                                visited.insert(point_tuple);
                                // Early termination for very good points
                                if(neighbors >= EARLY_TERMINATION_THRESHOLD) {
                                    found_good_point = true;
                                    Rcout << "Found optimal point at surface with " << neighbors << " neighbors. Early termination.\n";
                                }
                                break; // Stop searching h for this (x, y, z) after first nonzero
                            }
                        }
                    }
                }
            }
            
            // After finding the first non-zero point, search along the line to center and beyond
            if(found_nonzero_at_radius && pq.size() > 0) {
                Rcout << "Found non-zero points at radius " << radius << ", searching along center trajectory...\n";
                
                // Make a copy of current queue to avoid modifying it while iterating
                std::vector<Point3D> current_points;
                while(!pq.empty()) {
                    current_points.push_back(pq.top());
                    pq.pop();
                }
                
                // Maximum neighbors found so far in trajectory search
                int best_trajectory_neighbors = 0;
                
                // For each point found at this radius, search along trajectory to/through center
                for(const auto& point : current_points) {
                    // Add back to priority queue for later BFS
                    pq.push(point);
                    
                    // Calculate vector from this point to center
                    float dx = center_x - point.x;
                    float dy = center_y - point.y;
                    float dz = center_z - point.z;
                    
                    // Normalize to create a unit vector
                    float length = std::sqrt(dx*dx + dy*dy + dz*dz);
                    if(length < 0.001f) continue; // Skip if very close to center
                    
                    dx /= length;
                    dy /= length;
                    dz /= length;
                    
                    // Search along trajectory in both directions (toward and away from center)
                    // Start from original point and take steps both toward and past center
                    bool found_nonzero_on_trajectory = false;
                    for(float t = -2.0f * length; t <= 2.0f * length && !found_nonzero_on_trajectory; t += 1.0f) {
                        int x = std::round(point.x + dx * t);
                        int y = std::round(point.y + dy * t);
                        int z = std::round(point.z + dz * t); 
                        
                        if(!is_valid_point(x, y, z, D0, D1, D2)) continue;
                        
                        // Check if already visited
                        auto point_tuple = std::make_tuple(x, y, z, point.h);
                        if(visited.find(point_tuple) != visited.end()) continue;
                        
                        // Check if point has non-zero value
                        int idx = x * strides[0] + y * strides[1] + z * strides[2];
                        int remaining_h = point.h;
                        for (int d = ndim-1; d >= 3; d--) {
                            idx += (remaining_h % dimensions[d]) * strides[d];
                            remaining_h /= dimensions[d];
                        }
                        
                        if(idx < 0 || idx >= array_size) continue;
                        
                        if(std::abs(array[idx]) > threshold) {
                            Rcout << "Checking trajectory point at (" << x << "," << y << "," << z << "," << point.h << ")\n";
                            int neighbors = count_neighbors_3d(array, strides, x, y, z, point.h, D0, D1, D2, dimensions, array_size, threshold);
                            pq.push(Point3D(x, y, z, point.h, neighbors));
                            visited.insert(point_tuple);
                            surface_checks++;
                            
                            if(neighbors > best_trajectory_neighbors) {
                                best_trajectory_neighbors = neighbors;
                            }
                            
                            // If we found a very good point, stop searching this trajectory
                            if(neighbors >= EARLY_TERMINATION_THRESHOLD) {
                                Rcout << "Found optimal point along trajectory with " << neighbors << " neighbors. Early termination.\n";
                                found_good_point = true;
                                found_nonzero_on_trajectory = true;
                            }
                            
                            // Even if not optimal, break to move to next trajectory - no need to check every point
                            found_nonzero_on_trajectory = true;
                        }
                        
                        // If we've exceeded max checks, abort the search
                        if(surface_checks > MAX_SURFACE_CHECKS) {
                            Rcout << "Exceeded maximum surface checks (" << MAX_SURFACE_CHECKS << "). Terminating search.\n";
                            found_good_point = true;
                            break;
                        }
                    }
                    
                    // If we found a good point on any trajectory (high neighbor count), stop searching other trajectories
                    if(best_trajectory_neighbors >= EARLY_TERMINATION_THRESHOLD / 2) {
                        Rcout << "Found good trajectory point with " << best_trajectory_neighbors << " neighbors. Moving to BFS phase.\n";
                        break;
                    }
                }
            }
        }
    }
    
    log_time("Surface point checks complete. Checked " + std::to_string(surface_checks) + " surface points");
    
    // Search using priority queue to always explore best points first
    Point3D best_point(-1, -1, -1, -1, -1);
    int max_neighbors = -1;
    int neighbor_checks = 0;
    
    Rcout << "Queue size after initial search: " << pq.size() << "\n";
    Rcout << "Visited set size: " << visited.size() << "\n";
    
    while(!pq.empty()) {
        Point3D current = pq.top();
        pq.pop();
        
        if(current.neighbor_count > max_neighbors) {
            max_neighbors = current.neighbor_count;
            best_point = current;
            Rcout << "Found better point at (" << current.x << "," << current.y << "," << current.z 
                  << ") in volume " << current.h << " with " << current.neighbor_count << " neighbors\n";
            
            // Early termination if we've found the theoretical maximum or a very good point (18+ neighbors)
            if(max_neighbors == MAX_THEORETICAL_NEIGHBORS || max_neighbors >= EARLY_TERMINATION_THRESHOLD) {
                Rcout << "Found optimal point with " << max_neighbors << " neighbors. Early termination.\n";
                break; // Exit the search loop immediately
            }
        }
        
        // Add unvisited neighbors to queue
        for(int dx = -1; dx <= 1; dx++) {
            for(int dy = -1; dy <= 1; dy++) {
                for(int dz = -1; dz <= 1; dz++) {
                    if(dx == 0 && dy == 0 && dz == 0) continue;  // Skip self
                    
                    int nx = current.x + dx;
                    int ny = current.y + dy;
                    int nz = current.z + dz;
                    
                    if(!is_valid_point(nx, ny, nz, D0, D1, D2)) continue;
                    
                    // Check if already visited using the hash set
                    auto point_tuple = std::make_tuple(nx, ny, nz, current.h);
                    if(visited.find(point_tuple) != visited.end()) continue;
                    
                    // Calculate flat index using strides
                    int idx = nx * strides[0] + ny * strides[1] + nz * strides[2];
                    
                    // Add offset for higher dimensions by unraveling h
                    int remaining_h = current.h;
                    for (int d = ndim-1; d >= 3; d--) {
                        int dim_size = dimensions[d];
                        int dim_idx = remaining_h % dim_size;
                        idx += dim_idx * strides[d];
                        remaining_h /= dim_size;
                    }
                    
                    // Strict bounds checking
                    if(idx < 0 || idx >= array_size) continue;
                    
                    if(std::abs(array[idx]) > threshold) {
                        neighbor_checks++;
                        int neighbors = count_neighbors_3d(array, strides, nx, ny, nz, current.h, D0, D1, D2, dimensions, array_size, threshold);
                        pq.push(Point3D(nx, ny, nz, current.h, neighbors));
                        visited.insert(point_tuple);
                    }
                }
            }
        }
    }
    
    log_time("BFS search complete. Checked " + std::to_string(neighbor_checks) + " neighbor points");
    
    Rcout << "Best point found: (" << best_point.x << "," << best_point.y << "," << best_point.z 
          << ") in volume " << best_point.h << " with " << max_neighbors << " neighbors\n";
          
    // ADDED: Optimize center point based on orthogonal plane checks
    if(best_point.x != -1) {
        log_time("Starting orthogonal plane validation and center optimization");
        
        // Check the current best point's orthogonal planes
        double initial_score = check_orthogonal_planes(array, strides, best_point.x, best_point.y, 
                                                     best_point.z, best_point.h, D0, D1, D2, 
                                                     dimensions, array_size, PLANE_CHECK_RADIUS, threshold);
        
        Rcout << "Initial orthogonal plane check: " << (initial_score * 100) << "% non-zero voxels\n";
        
        if(initial_score < 0.99) {  // Only optimize if needed (less than 99% non-zero)
            // Optimize the center point
            best_point = optimize_center_point(array, strides, best_point, D0, D1, D2, 
                                              dimensions, array_size, higher_dims_size, threshold);
            
            log_time("Center optimization complete");
        } else {
            Rcout << "Center point has optimal orthogonal planes, no optimization needed\n";
        }
    }

    // Return just the x, y, z coordinates as a vector of length 3
    IntegerVector result;
    if(best_point.x != -1) {
        // Create a vector with the coordinates (converting to 1-based indexing for R)
        result = IntegerVector::create(best_point.x + 1, best_point.y + 1, best_point.z + 1);
    } else {
        // If no point found, return zeros
        result = IntegerVector::create(0, 0, 0);
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    Rcout << "Total execution time: " << total_time << "ms\n";
    
    return result;
}