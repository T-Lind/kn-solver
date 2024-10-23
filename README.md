# Kissing Number Solver

This project is a C++ implementation of a Kissing Number solver. A kissing number is the number of non-overlapping unit
spheres that can be arranged to touch another unit sphere in a given space. In two dimensions, the kissing number is 6,
while in three dimensions, it is 12. The kissing number is unknown for higher dimensions, like 5D and 6D. This program
aims to find the kissing number for higher dimensions.

The program distributes points on the surface of a sphere using electrostatic repulsion and checks for overlapping
spheres. In this way, it can check whether or not a certain number of unit spheres can be arranged to touch another unit
sphere in a given space.

## Features

- Distributes points on a sphere using repulsive forces.
- Generates outer points scaled by a factor.
- Counts overlapping spheres.
- Prints distance matrices and points.

## Requirements

- C++20
- CMake 3.22 or higher

## Build and Run

1. **Clone the repository:**
   ```sh
   git clone <repository-url>
   cd kn_solver
   ```

2. **Build the project:**
   ```sh
   mkdir build
   cd build
   cmake ..
   make
   ```

3. **Run the executable:**
   ```sh
   ./kn_solver
   ```

## Usage

The main parameters can be adjusted in the `main.cpp` file:

- `num_points`: Number of points to distribute.
- `dimensions`: Dimensionality of the space.
- `n_iterations`: Number of iterations for the distribution algorithm.
- `radius`: Radius of the spheres.
- `scale_factor`: Scaling factor for generating outer points.