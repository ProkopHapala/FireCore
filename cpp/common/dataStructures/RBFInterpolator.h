#ifndef RBFInterpolator_h
#define RBFInterpolator_h

// sketched using OpenAI-o1 https://chatgpt.com/share/672cd6e1-5e44-8012-a103-e80bdbb77c8d

#include "Buckets2D.h"
#include "Vec2.h"
#include <cmath>
#include <algorithm> // For std::min, std::max

/// @brief Class for performing Radial Basis Function interpolation in 2D
class RBFInterpolator {
public:
    // All members are public
    Buckets2D buckets2D;         // Spatial partitioning map
    Buckets& map;                // Reference to Buckets2D::map
    double Rcut;
    double epsilon;
    const Vec2d* reference_points;
    const double* values;
    int num_reference;

    int* nearby_indices_buffer;  // Preallocated buffer for nearby indices
    int buffer_size;

    // Constructor
    RBFInterpolator(const Vec2d* reference_points_, int num_reference_,
                    const double* values_, double Rcut_, double epsilon_)
        : buckets2D(), map(buckets2D.map),
          Rcut(Rcut_), epsilon(epsilon_),
          reference_points(reference_points_), values(values_),
          num_reference(num_reference_),
          nearby_indices_buffer(nullptr), buffer_size(0)
    {
        // Initialize Buckets2D with reference points
        buckets2D.initialize(reference_points_, num_reference_, Rcut, epsilon);

        // Allocate buffer for nearby_indices
        buffer_size = map.maxInBucket * 9;
        nearby_indices_buffer = new int[buffer_size];
    }

    // Destructor
    ~RBFInterpolator() {
        if (nearby_indices_buffer) {
            delete[] nearby_indices_buffer;
            nearby_indices_buffer = nullptr;
        }
    }

    /// @brief Performs interpolation on an array of sample points
    /// @param sample_points Array of sample points
    /// @param num_sample Number of sample points
    /// @param interpolated_values Output array for interpolated values
    void interpolate(const Vec2d* sample_points, int num_sample_,
                    double* interpolated_values) const
    {
        for (int m = 0; m < num_sample_; ++m) {
            interpolated_values[m] = interpolatePoint(sample_points[m]);
        }
    }

    /// @brief Interpolates the value at a single sample point
    double interpolatePoint(const Vec2d& point) const {
        // Determine grid cell
        int cell_idx = buckets2D.getCellIndex(point);

        // Retrieve neighboring cells
        int neighbor_cells[9];
        buckets2D.getNeighborCells(cell_idx, neighbor_cells);

        // Collect nearby reference points using preallocated buffer
        int count = 0;

        for (int k = 0; k < 9; ++k) {
            int neighbor_idx = neighbor_cells[k];
            if (neighbor_idx == -1) continue;

            int num_in_cell = map.cellNs[neighbor_idx];
            int start_idx = map.cellI0s[neighbor_idx];
            for (int m = 0; m < num_in_cell; ++m) {
                nearby_indices_buffer[count++] = map.cell2obj[start_idx + m];
            }
        }

        // Compute RBF weights
        double weight_sum = 0.0;
        double interpolated_value = 0.0;

        for (int m = 0; m < count; ++m) {
            int ref_idx = nearby_indices_buffer[m];
            const Vec2d& ref_pt = reference_points[ref_idx];
            double dist = (ref_pt - point).norm();

            if (dist <= Rcut) {
                double w = rbf(dist);
                interpolated_value += w * values[ref_idx];
                weight_sum += w;
            }
        }

        if (weight_sum > 0.0) {
            interpolated_value /= weight_sum;
            return interpolated_value;
        } else {
            return 0.0; // Default value if no neighbors within cutoff
        }
    }

    /// @brief RBF kernel (Gaussian)
    inline double rbf(double r) const {
        return std::exp(-(epsilon * r) * (epsilon * r));
    }
};

#endif // RBFInterpolator_h
