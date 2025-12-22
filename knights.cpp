// Count undirected Hamiltonian knight paths on k x n chessboards.
// Uses backtracking DFS with Warnsdorff's rule and advanced pruning techniques.
//
// Compiled from various AI-sources (vibe coding), Peter Luschny, December 2025
// See OEIS A390833, but note that this code was **not** used there.

/*
**Key Optimization Strategies**

1. **Bit-Level Pruning (Degree & Connectivity):**

-- **Degree Heuristic:** Before attempting moves, we analyze the subgraph of 
unvisited nodes (`rem`). If *any* node is isolated (degree 0) or if *more than 
two* nodes are dead ends (degree 1), the branch is immediately pruned. This is
mathematically necessary for a Hamiltonian path (which can have at most two 
endpoints) and cuts down the search space drastically.

-- **Degree Cutoff:** The condition `if (ends > 2) return false;` is the "killer 
heuristic". In Hamiltonian path problems, most branches fail because the path 
"splinters" the remaining graph into disconnected segments or segments with too 
many endpoints. Detecting this with `popcount` before running the expensive BFS 
saves massive time.

-- **Bit-Parallel Flood Fill:** Instead of an array-based queue for connectivity 
checks, implement a bitwise flood-fill. This leverages bitwise OR operations to 
propagate connectivity, which is significantly faster for small dense grids (up 
to 64 nodes) than managing queue pointers. Uses register-based bit-flood 
(`wavefront |= neighbors[v]`). This stays entirely inside the CPU's general-purpose
registers, avoiding L1 cache latency.

2. **Symmetry Reduction (Group Theory):**

-- **Symmetry (Square Boards):** The board has 8 symmetries (rotations and 
reflections). The code identifies the "canonical" start nodes (the smallest 
lexicographical index in their orbit) and only runs the expensive DFS for these. 
The result is then multiplied by the orbit size (1, 2, 4, or 8). This yields an 
**~8x speedup** for square boards.

-- **Symmetry (Rectangular Boards):** It correctly identifies the 4 symmetries 
(Identity, 180° rotation, Horizontal/Vertical flips), yielding a **~4x speedup**.

-- **Parity Pruning:** On odd boards the start and end of a Hamiltonian path *must* 
be on the majority color. This instantly prunes half of the possible start nodes.

3. **Memory & Structure:**

-- **Stack Allocation:** Critical arrays (neighbors, candidates) are allocated on
the stack or in fixed-size structs to prevent heap fragmentation and pointer chasing.

-- **Lookahead Sorting (Warnsdorff’s Rule):** Moves are dynamically sorted so 
that the most constrained neighbors (lowest degree in `rem`) are visited first. 
This forces "hard" decisions early, leading to faster cutoffs.

4. **Parallelization:**

-- **OpenMP:** The outer loop over canonical start positions is parallelized using 
OpenMP. This is an embarrassingly parallel problem since each start position's 
DFS is independent. The reduction clause efficiently aggregates results.

-- **Dynamic Scheduling:** Using dynamic scheduling in OpenMP helps balance the 
workload, as some start positions may lead to significantly deeper searches than others.
*/
 
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>
#include <omp.h>
using u64 = std::uint64_t;

// Compiler intrinsics for speed (GCC/Clang)
static inline int popcount64(u64 x) { return __builtin_popcountll(x); }
static inline int ctz64(u64 x) { return __builtin_ctzll(x); }
static inline bool has_at_most_one_bit(u64 x) {
    return (x & (x - 1)) == 0;
}

struct Solver {
    int K, N, V;
    u64 ALL_MASK;
    // Align to 64 bytes to match cache line size, preventing cache line splits.
    alignas(64) std::array<u64, 64> neighbors{};

    Solver(int k, int n) : K(k), N(n), V(k*n) {
        ALL_MASK = (V == 64) ? ~0ULL : ((1ULL << V) - 1ULL);
        
        static constexpr int moves[8][2] = {
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

    // Heuristic Pruning: Degree checks & Connectivity
    inline bool is_valid_state(u64 rem) const {
        if (rem == 0 || has_at_most_one_bit(rem)) return true;

        int ends = 0;
        u64 temp_rem = rem;

        while (temp_rem) {
            int v = ctz64(temp_rem);
            temp_rem &= temp_rem - 1; 
            
            // Check neighbors within unvisited set
            u64 n_rem = neighbors[v] & rem;
            
            // Optimization: Zero neighbors = Isolated = Fail
            if (n_rem == 0) return false;

            // Optimization: One neighbor = Endpoint
            if (has_at_most_one_bit(n_rem)) {
                ends++;
                if (ends > 2) return false; 
            }
        }

        // Connectivity Check (Bit-Flood)
        int start_node = ctz64(rem);
        u64 connected = (1ULL << start_node);
        u64 wavefront = connected;

        while (wavefront) {
            u64 new_wavefront = 0;
            while (wavefront) {
                int v = ctz64(wavefront);
                wavefront &= wavefront - 1;
                new_wavefront |= (neighbors[v] & rem);
            }
            new_wavefront &= ~connected;
            connected |= new_wavefront;
            wavefront = new_wavefront;
        }

        return (connected & rem) == rem;
    }

    // DFS with Forced Move Pruning
    u64 dfs(int current, u64 visited) const {
        if (visited == ALL_MASK) return 1;

        u64 rem = ALL_MASK & ~visited;
        u64 moves_mask = neighbors[current] & rem;
        
        // 1. Immediate Dead-End Check
        if (!moves_mask) return 0;
        
        // 2. Global Validity Check
        if (!is_valid_state(rem)) return 0;

        int candidates[8];
        int degrees[8];
        int count = 0;
        int forced_move = -1;

        // 3. Move Generation & Forced Move Detection
        while (moves_mask) {
            int v = ctz64(moves_mask);
            moves_mask &= moves_mask - 1;

            // Calculate degree in the remaining graph (lookahead)
            int deg = popcount64(neighbors[v] & rem);

            // Critical optimization: Forced Move
            // If 'v' has 0 neighbors in 'rem', it is a leaf attached only to 'current'.
            // We *must* visit it now, otherwise it becomes isolated.
            if (deg == 0) {
                if (forced_move != -1) return 0; // Two forced moves -> Impossible (fork)
                forced_move = v;
            }

            // Insertion Sort (Warnsdorff's Rule)
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

        // 4. Execution
        // If a forced move exists, ignore all other candidates.
        if (forced_move != -1) {
            return dfs(forced_move, visited | (1ULL << forced_move));
        }

        // Otherwise, recurse on sorted candidates
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
    // Keep this allocation-free: at most 8 symmetries.
    BoardSym syms[8];
    int count = 0;

    // D2 symmetries (rectangle)
    syms[count++] = {r, c};                  // 1. Identity
    syms[count++] = {K - 1 - r, c};          // 2. Reflect Vertical (across mid-row)
    syms[count++] = {r, N - 1 - c};          // 3. Reflect Horizontal (across mid-col)
    syms[count++] = {K - 1 - r, N - 1 - c};  // 4. Rotate 180

    if (K == N) {  // D4 symmetries (square)
        syms[count++] = {c, r};                  // 5. Transpose (Main Diag)
        syms[count++] = {N - 1 - c, K - 1 - r};  // 6. Anti-Transpose (Anti Diag)
        syms[count++] = {c, K - 1 - r};          // 7. Rot 90 CW
        syms[count++] = {N - 1 - c, r};          // 8. Rot 270 CW
    }

    // Find canonical (min)
    BoardSym min_s = syms[0];
    for (int i = 1; i < count; i++) {
        if (syms[i] < min_s) min_s = syms[i];
    }

    out_r = min_s.r;
    out_c = min_s.c;

    // Count distinct symmetries (Orbit size)
    // orbit = number of unique symmetries
    int uniq = 0;
    for (int i = 0; i < count; i++) {
        bool seen = false;
        for (int j = 0; j < uniq; j++) {
            if (syms[i] == syms[j]) {
                seen = true;
                break;
            }
        }
        if (!seen) syms[uniq++] = syms[i];
    }
    orbit = uniq;
}

u64 knight_hamiltonian_paths(int k, int n) {
    if (k > n) std::swap(k, n); // Ensure k <= n
    if (k * n > 64) {
        std::cerr << "Error: Board too large for 64-bit mask implementation.\n";
        return 0;
    }

    // Pruning: Parity Check
    // On an odd-sized board, start and end must be the SAME color (Majority).
    // On an even-sized board, start and end must be DIFFERENT colors.
    // If the board is odd, we can only start on the majority color.
    // However, undirected paths allow A->B (start min, end maj).
    // But Hamiltonian path on odd graph MUST start/end on Majority.
    // (Majority count = Minority + 1. Path alternates. M-m-M-m...-M).
    // So if we pick a start node on Minority color on an Odd board, count is 0.
    bool odd_board = (k * n) % 2 != 0;
    const Solver solver(k, n);
    u64 total_directed = 0;

    // We collect tasks: (start_node, orbit_size)
    // to parallelize over canonical start nodes.
    struct Task {
        int start_node;
        int weight;
    };
    std::vector<Task> tasks;
    tasks.reserve(static_cast<size_t>(k) * static_cast<size_t>(n));

    for (int r = 0; r < k; r++) {
        for (int c = 0; c < n; c++) {
            int can_r, can_c, orbit;
            get_canonical(r, c, k, n, can_r, can_c, orbit);

            // Only process if current is the canonical representative
            if (r == can_r && c == can_c) {
                // Parity Pruning: On odd boards, path must start on majority color.
                // Assuming (0,0) is white/majority.
                if (odd_board) {
                    // We just need to know if it's the minority color.
                    // White squares = ceil(V/2), Black = floor(V/2).
                    // If V is odd, White is Majority.
                    // If (r+c)%2 != 0, it's Black (Minority). Impossible to start path here.
                    bool is_white = (r + c) % 2 == 0;
                    if (!is_white) continue; 
                }
                tasks.push_back(Task{r * n + c, orbit});
            }
        }
    }

    #pragma omp parallel for schedule(dynamic) reduction(+:total_directed)
    for (size_t i = 0; i < tasks.size(); i++) {
        int start_node = tasks[i].start_node;
        int weight = tasks[i].weight;

        u64 paths = solver.dfs(start_node, 1ULL << start_node);
        total_directed += (paths * static_cast<u64>(weight));
    }

    return total_directed / 2;
}

// --- Driver & Benchmarking ---

void A390833(int k, int n) {
    auto start = std::chrono::steady_clock::now();
    u64 result = knight_hamiltonian_paths(k, n);
    auto end = std::chrono::steady_clock::now();
    double seconds = std::chrono::duration<double>(end - start).count();

    std::cout << "A(" << k << "," << n << ") = " << result 
              << "  \tTime: " << seconds << "s" << std::endl;
}

void benchmark() {
    for (int k = 3; k < 12; k++) {
        for (int n = 3; n < 12; n++) {
            if (n * k <= 40) { // Limit to manageable sizes
                A390833(k, n);
            }
        }
    }
}

int main() {
    omp_set_num_threads(omp_get_max_threads());

    auto start = std::chrono::steady_clock::now();
    benchmark();
    auto end = std::chrono::steady_clock::now();
    double total_seconds = std::chrono::duration<double>(end - start).count();
    std::cout << "Total Time: " << total_seconds << "s" << std::endl;

    return 0;
}

/*  Makefile commands to compile and run:
sudo pacman -Syu
sudo pacman -S gcc clang make cmake libomp
g++ -O3 -march=native -fopenmp -std=c++20 knights.cpp -o knights
./knights
*/
