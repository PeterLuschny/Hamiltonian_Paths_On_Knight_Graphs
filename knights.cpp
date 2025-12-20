// knights.cpp
// Count undirected Hamiltonian knight paths on k x n chessboards.
// Uses backtracking DFS with Warnsdorff's rule and advanced pruning techniques.
// See OEIS A390833.
//
// Peter Luschny, compiled from various AI-sources, December 2025

/*
**Key Optimization Strategies**

1. **Bit-Level Pruning (Degree & Connectivity):**
-- **Degree Heuristic:** Before attempting moves, we analyze the subgraph of unvisited nodes (`rem`). If *any* node is isolated (degree 0) or if *more than two* nodes are dead ends (degree 1), the branch is immediately pruned. This is mathematically necessary for a Hamiltonian path (which can have at most two endpoints) and cuts down the search space drastically.
-- **Degree Cutoff:** The condition `if (ends > 2) return false;` is the "killer heuristic". In Hamiltonian path problems, most branches fail because the path "splinters" the remaining graph into disconnected segments or segments with too many endpoints. Detecting this with `popcount` before running the expensive BFS saves massive time.
-- **Bit-Parallel Flood Fill:** Instead of an array-based queue for connectivity checks, implement a bitwise flood-fill. This leverages bitwise OR operations to propagate connectivity, which is significantly faster for small dense grids (up to 64 nodes) than managing queue pointers. Uses register-based bit-flood (`wavefront |= neighbors[v]`). This stays entirely inside the CPU's general-purpose registers, avoiding L1 cache latency.

2. **Symmetry Reduction (Group Theory):**
-- **Symmetry (Square Boards):** The board has 8 symmetries (rotations and reflections). The code identifies the "canonical" start nodes (the smallest lexicographical index in their orbit) and only runs the expensive DFS for these. The result is then multiplied by the orbit size (1, 2, 4, or 8). This yields an **~8x speedup** for square boards.
-- **Symmetry (Rectangular Boards):** It correctly identifies the 4 symmetries (Identity, 180° rotation, Horizontal/Vertical flips), yielding a **~4x speedup**.
-- **Parity Pruning:** On odd boards the start and end of a Hamiltonian path *must* be on the majority color. This instantly prunes half of the possible start nodes.

3. **Memory & Structure:**
-- **Stack Allocation:** Critical arrays (neighbors, candidates) are allocated on the stack or in fixed-size structs to prevent heap fragmentation and pointer chasing.
-- **Lookahead Sorting (Warnsdorff’s Rule):** Moves are dynamically sorted so that the most constrained neighbors (lowest degree in `rem`) are visited first. This forces "hard" decisions early, leading to faster cutoffs.

4. **Parallelization:**
-- **OpenMP:** The outer loop over canonical start positions is parallelized using OpenMP. This is an embarrassingly parallel problem since each start position's DFS is independent. The reduction clause efficiently aggregates results.
-- **Dynamic Scheduling:** Using dynamic scheduling in OpenMP helps balance the workload, as some start positions may lead to significantly deeper searches than others.
*/
 
#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <chrono>

using namespace std;
using u64 = unsigned long long;

// Compiler intrinsics for speed (GCC/Clang)
#define POPCOUNT(x) __builtin_popcountll(x)
#define CTZ(x) __builtin_ctzll(x)

struct Solver {
    int K, N, V;
    u64 ALL_MASK;
    u64 neighbors[64];

    Solver(int k, int n) : K(k), N(n), V(k*n) {
        ALL_MASK = (V == 64) ? ~0ULL : ((1ULL << V) - 1ULL);
        
        // Precompute neighbor masks
        int moves[8][2] = {
            {1,2}, {1,-2}, {-1,2}, {-1,-2},
            {2,1}, {2,-1}, {-2,1}, {-2,-1}
        };

        for (int r = 0; r < K; r++) {
            for (int c = 0; c < N; c++) {
                int u = r * N + c;
                u64 m = 0;
                for (auto &mv : moves) {
                    int rr = r + mv[0];
                    int cc = c + mv[1];
                    if (rr >= 0 && rr < K && cc >= 0 && cc < N) {
                        m |= (1ULL << (rr * N + cc));
                    }
                }
                neighbors[u] = m;
            }
        }
    }

    // --- Heuristic Pruning ---
    // Checks if the unvisited subgraph 'rem' can possibly support a path extension.
    inline bool is_valid_state(u64 rem) const {
        if (POPCOUNT(rem) <= 1) return true;

        // 1. Degree Check within the unvisited subgraph
        int ends = 0;
        u64 temp_rem = rem;
        
        while (temp_rem) {
            int v = CTZ(temp_rem);
            temp_rem &= temp_rem - 1; 

            // Intersection of neighbors and remaining unvisited nodes
            u64 n_rem = neighbors[v] & rem;
            int d = POPCOUNT(n_rem);

            // If a node has no unvisited neighbors, it is isolated.
            // Since we established popcount > 1, this is a dead end.
            if (d == 0) return false; 
            
            // Nodes with degree 1 in the unvisited set are potential endpoints.
            // A simple path can have at most 2 endpoints.
            // Note: The entry point from the 'visited' set effectively "saves" one endpoint.
            if (d == 1) {
                ends++;
                if (ends > 2) return false; 
            }
        }

        // 2. Connectivity Check (Bitwise Flood Fill)
        int start_node = CTZ(rem);
        u64 connected = (1ULL << start_node);
        u64 wavefront = connected;

        while (wavefront) {
            u64 new_wavefront = 0;
            while (wavefront) {
                int v = CTZ(wavefront);
                wavefront &= wavefront - 1;
                new_wavefront |= (neighbors[v] & rem);
            }
            new_wavefront &= ~connected;
            connected |= new_wavefront;
            wavefront = new_wavefront;
        }

        return (connected & rem) == rem;
    }

    // Backtracking DFS
    u64 dfs(int current, u64 visited) {
        if (visited == ALL_MASK) return 1;

        u64 moves_mask = neighbors[current] & ~visited;
        if (!moves_mask) return 0;

        // Lookahead pruning
        u64 rem = ALL_MASK & ~visited;
        if (!is_valid_state(rem)) return 0;

        // Warnsdorff's Rule Sorting
        int candidates[8];
        int degrees[8];
        int count = 0;

        while (moves_mask) {
            int v = CTZ(moves_mask);
            moves_mask &= moves_mask - 1;

            // Sort by degree in the *remaining* graph minus the node v itself
            int deg = POPCOUNT(neighbors[v] & rem) - 1; 

            int i = count;
            while (i > 0 && degrees[i - 1] > deg) {
                degrees[i] = degrees[i - 1];
                candidates[i] = candidates[i - 1];
                i--;
            }
            candidates[i] = v;
            degrees[i] = deg;
            count++;
        }

        u64 total_paths = 0;
        for (int i = 0; i < count; i++) {
            total_paths += dfs(candidates[i], visited | (1ULL << candidates[i]));
        }

        return total_paths;
    }
};

// --- Symmetry & Canonical Start Logic ---

struct BoardSym {
    int r, c;
    bool operator<(const BoardSym& o) const {
        if (r != o.r) return r < o.r;
        return c < o.c;
    }
    bool operator==(const BoardSym& o) const {
        return r == o.r && c == o.c;
    }
};

void get_canonical(int r, int c, int K, int N, int& out_r, int& out_c, int& orbit) {
    vector<BoardSym> syms;
    syms.reserve(8);
    
    // D2 Symmetries (Rectangle)
    syms.push_back({r, c});
    syms.push_back({K - 1 - r, c});
    syms.push_back({r, N - 1 - c});
    syms.push_back({K - 1 - r, N - 1 - c});

    if (K == N) {
        // D4 Symmetries (Square)
        syms.push_back({c, r});
        syms.push_back({N - 1 - c, K - 1 - r});
        syms.push_back({c, K - 1 - r});
        syms.push_back({N - 1 - c, r});
    }

    BoardSym min_s = syms[0];
    for (auto& s : syms) {
        if (s < min_s) min_s = s;
    }
    
    out_r = min_s.r;
    out_c = min_s.c;

    sort(syms.begin(), syms.end());
    orbit = unique(syms.begin(), syms.end()) - syms.begin();
}

u64 knight_hamiltonian_paths(int k, int n) {
    if (k > n) std::swap(k, n);
    if (k * n > 64) return 0;

    bool odd_board = (k * n) % 2 != 0;
    Solver solver(k, n);
    u64 total_directed = 0;

    vector<pair<int, int>> tasks;

    for (int r = 0; r < k; r++) {
        for (int c = 0; c < n; c++) {
            int can_r, can_c, orbit;
            get_canonical(r, c, k, n, can_r, can_c, orbit);
            
            if (r == can_r && c == can_c) {
                // Parity Pruning: On odd boards, path must start on majority color.
                // Assuming (0,0) is white/majority.
                if (odd_board) {
                    bool is_white = (r + c) % 2 == 0;
                    if (!is_white) continue; 
                }
                tasks.push_back({r * n + c, orbit});
            }
        }
    }

    #pragma omp parallel for schedule(dynamic) reduction(+:total_directed)
    for (size_t i = 0; i < tasks.size(); i++) {
        int start_node = tasks[i].first;
        int weight = tasks[i].second;

        Solver local_solver(k, n); 
        u64 paths = local_solver.dfs(start_node, 1ULL << start_node);
        total_directed += (paths * weight);
    }

    return total_directed / 2;
}

void A(int k, int n) {
    auto start = chrono::steady_clock::now();
    u64 result = knight_hamiltonian_paths(k, n);
    auto end = chrono::steady_clock::now();
    double seconds = chrono::duration<double>(end - start).count();

    cout << "A(" << k << "," << n << ") = " << result 
         << "   \tTime: " << seconds << "s" << endl;
}

void benchmark() {
        for (int k = 3; k < 12; k++) {
            for (int n = 3; n < 12; n++) {
                if (n * k <= 40) { // Limit to manageable sizes
                    A(k, n);
                }
            }
        }
}

int main() {
    omp_set_num_threads(omp_get_max_threads());

    auto start = chrono::steady_clock::now();
    benchmark();
    auto end = chrono::steady_clock::now();
    double total_seconds = chrono::duration<double>(end - start).count();
    cout << "Total Time: " << total_seconds << "s" << endl;

    return 0;
}

/*
A(3,3) = 0      Time: 0.0040712s
A(3,4) = 8      Time: 3.25e-05s
A(3,5) = 0      Time: 1.14e-05s
A(3,6) = 0      Time: 2.31e-05s
A(3,7) = 52     Time: 8.47e-05s
A(3,8) = 396    Time: 0.000403s
A(3,9) = 560    Time: 0.0009808s
A(3,10) = 3048  Time: 0.0041568s
A(3,11) = 10672 Time: 0.0189267s
A(4,3) = 8      Time: 0.000235s
A(4,4) = 0      Time: 3.51e-05s
A(4,5) = 82     Time: 0.0027262s
A(4,6) = 744    Time: 0.0119047s
A(4,7) = 6378   Time: 0.0185291s
A(4,8) = 31088  Time: 0.140034s
A(4,9) = 189688 Time: 1.00172s
A(4,10) = 1213112 Time: 7.46396s
A(5,3) = 0      Time: 0.0033393s
A(5,4) = 82     Time: 0.0020133s
A(5,5) = 864    Time: 0.005744s
A(5,6) = 18784  Time: 0.0384103s
A(5,7) = 622868 Time: 1.06086s
A(5,8) = 18061054 Time: 27.3598s
A(6,3) = 0      Time: 0.0001351s
A(6,4) = 744    Time: 0.0015218s
A(6,5) = 18784  Time: 0.0281101s
A(6,6) = 3318960 Time: 2.82377s
A(7,3) = 52     Time: 9.22e-05s
A(7,4) = 6378   Time: 0.0135029s
A(7,5) = 622868 Time: 1.01951s
A(8,3) = 396    Time: 0.0003808s
A(8,4) = 31088  Time: 0.100319s
A(8,5) = 18061054 Time: 27.9002s
A(9,3) = 560    Time: 0.0008884s
A(9,4) = 189688 Time: 0.831437s
A(10,3) = 3048  Time: 0.004737s
A(10,4) = 1213112 Time: 7.21754s
A(11,3) = 10672 Time: 0.0163133s
Total Time: 77.0983s
*/

/*
sudo pacman -Syu
sudo pacman -S gcc clang make cmake libomp
g++ -O3 -march=native -fopenmp -std=c++17 knights.cpp -o knights
./knights
*/
