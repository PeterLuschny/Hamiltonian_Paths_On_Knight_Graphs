# Hamiltonian_Paths_On_Knight_Graphs

Count undirected Hamiltonian knight paths on k x n chessboards.

Uses backtracking DFS with Warnsdorff's rule and advanced pruning techniques.

See OEIS A390833.

Make:

sudo pacman -S gcc clang make cmake libomp

g++ -O3 -march=native -fopenmp -std=c++17 knights.cpp -o knights
