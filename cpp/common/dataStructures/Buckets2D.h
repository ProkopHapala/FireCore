#ifndef Buckets2D_h
#define Buckets2D_h

// sketched using OpenAI-o1 https://chatgpt.com/share/672cd6e1-5e44-8012-a103-e80bdbb77c8d

#include "Buckets.h"    // Include your existing Buckets class
#include "Vec2.h"      // 2D vector class
#include <cmath>
#include <cstring>      // For memset
#include <algorithm>    // For std::min, std::max

/// @brief Class for managing 2D spatial partitioning using Buckets
class Buckets2D {
public:
    // Parameters
    double Rcut;
    double epsilon;

    // Reference points
    const Vec2d* reference_points;
    int num_reference;

    // Grid parameters
    double x_min, y_min;
    double x_max, y_max;
    int n_x, n_y;
    int ncell;

    // Buckets instance, renamed to 'map'
    Buckets map;

    // Array for object to cell mapping
    int* obj2cell_array;

    // Constructor
    Buckets2D()
        : Rcut(5.0), epsilon(1.0),
          x_min(0.0), y_min(0.0),
          x_max(0.0), y_max(0.0),
          n_x(0), n_y(0), ncell(0),
          obj2cell_array(nullptr)
    {}

    // Destructor
    ~Buckets2D() {
        if (obj2cell_array) {
            delete[] obj2cell_array;
            obj2cell_array = nullptr;
        }
        // Clean up Buckets memory
        delete[] map.cellNs;
        delete[] map.cellI0s;
        delete[] map.cell2obj;
        delete[] map.obj2cell;
    }

    /// @brief Initializes the Buckets2D with reference points and parameters
    void initialize(const Vec2d* reference_points_, int num_reference_,  double Rcut_, double epsilon_) {
        reference_points = reference_points_;
        num_reference = num_reference_;
        Rcut = Rcut_;
        epsilon = epsilon_;

        computeBoundingBox();
        computeGridDimensions();
        setupMap();
        assignPointsToMap();
    }

    /// @brief Computes the bounding box of the reference points
    void computeBoundingBox() {
        if (num_reference <= 0) return;

        x_min = reference_points[0].x;
        y_min = reference_points[0].y;
        x_max = reference_points[0].x;
        y_max = reference_points[0].y;

        for (int i = 1; i < num_reference; ++i) {
            if (reference_points[i].x < x_min) x_min = reference_points[i].x;
            if (reference_points[i].y < y_min) y_min = reference_points[i].y;
            if (reference_points[i].x > x_max) x_max = reference_points[i].x;
            if (reference_points[i].y > y_max) y_max = reference_points[i].y;
        }

        // Expand the bounding box slightly to include edge points
        x_min -= Rcut;
        y_min -= Rcut;
        x_max += Rcut;
        y_max += Rcut;
    }

    /// @brief Computes the number of grid cells in x and y directions
    void computeGridDimensions() {
        n_x = static_cast<int>(std::ceil((x_max - x_min) / Rcut));
        n_y = static_cast<int>(std::ceil((y_max - y_min) / Rcut));
        ncell = n_x * n_y;
    }

    /// @brief Sets up the Buckets instance for 2D grid
    void setupMap() {
        map.ncell = ncell;
        map.nobj = num_reference;

        // Allocate memory for Buckets arrays
        map.cellNs = new int[ncell];
        map.cellI0s = new int[ncell];
        map.cell2obj = new int[num_reference];
        map.obj2cell = new int[num_reference];

        // Initialize cellNs and cellI0s to zero
        std::memset(map.cellNs, 0, sizeof(int) * ncell);
        std::memset(map.cellI0s, 0, sizeof(int) * ncell);
        std::memset(map.cell2obj, -1, sizeof(int) * num_reference);
        std::memset(map.obj2cell, -1, sizeof(int) * num_reference);
    }

    /// @brief Assigns reference points to grid cells
    void assignPointsToMap() {
        // Allocate temporary obj2cell_array
        obj2cell_array = new int[num_reference];

        // Determine cell for each reference point
        for (int i = 0; i < num_reference; ++i) {
            int cell_idx = getCellIndex(reference_points[i]);
            obj2cell_array[i] = cell_idx;
        }

        // Bind the obj2cell mapping to Buckets
        map.bindObjs(num_reference, obj2cell_array);

        // Update Buckets with the object-to-cell mapping
        map.updateCells(num_reference, obj2cell_array);
    }

    /// @brief Retrieves cell index for a given point
    inline int getCellIndex(const Vec2d& point) const {
        int i = static_cast<int>(std::floor((point.x - x_min) / Rcut));
        int j = static_cast<int>(std::floor((point.y - y_min) / Rcut));

        // Clamp indices to be within grid bounds
        i = std::max(0, std::min(i, n_x - 1));
        j = std::max(0, std::min(j, n_y - 1));

        return i + j * n_x;
    }

    /// @brief Retrieves neighboring cell indices (including the current cell)
    void getNeighborCells(int cell_idx, int* neighbor_cells) const {
        int i = cell_idx % n_x;
        int j = cell_idx / n_x;

        int count = 0;
        for (int dj = -1; dj <= 1; ++dj) {
            for (int di = -1; di <= 1; ++di) {
                int ni = i + di;
                int nj = j + dj;

                // Check boundaries
                if (ni >= 0 && ni < n_x && nj >= 0 && nj < n_y) {
                    neighbor_cells[count++] = ni + nj * n_x;
                }
            }
        }

        // Fill remaining with -1 if less than 9 neighbors
        while (count < 9) {
            neighbor_cells[count++] = -1;
        }
    }
};

#endif // Buckets2D_h
