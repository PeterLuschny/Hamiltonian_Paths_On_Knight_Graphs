// Count undirected Hamiltonian knight paths on k x n chessboards.
// Uses backtracking DFS with Warnsdorff's rule and advanced pruning techniques.
// Compiled from various AI-sources (vibe coding), Peter Luschny, December 2025
// OEIS A390833, but note that this code was **not** used there.

/*
   Optimizations:
   1. Forced Move Pruning (Leaf Detection)
   2. Small-Stride Orientation (Canonicalization)
   3. Bit-Parallel Connectivity Checks
   4. D2/D4 Symmetry Reduction
   5. Parallelization with OpenMP

1. **Forced Move Pruning (Leaf Detection):**
   During DFS, if a move leads to a position where a neighbor has no further
   moves (a leaf), that move must be taken immediately. If multiple such forced
   moves exist, the path is invalid. This drastically reduces the search space
   by eliminating paths that would inevitably fail.

2. **Small-Stride Orientation (Canonicalization):**
   By ensuring the board is oriented such that the smaller dimension is treated
   as the width (N), we minimize the bit-distance between knight moves in memory.
   This improves cache locality and the effectiveness of heuristic move ordering.

3. **Bit-Parallel Connectivity Checks:**
   Using bitwise operations to perform connectivity checks allows for rapid
   determination of whether the remaining unvisited nodes form a connected
   subgraph. This is crucial for pruning invalid states early in the search.

4. **D2/D4 Symmetry Reduction:**
   By identifying canonical starting positions under the board's symmetry group,
   we reduce redundant computations. Each unique starting position is weighted
   by the size of its orbit, allowing us to count all equivalent paths without
   explicitly exploring them.

5. **Parallelization:**
   The outer loop over canonical starting positions is parallelized, allowing
   multiple DFS searches to occur simultaneously. This takes advantage of multi-
   core processors to significantly speed up the overall computation.
*/

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>
#include <omp.h>

using u64 = std::uint64_t;

// Compiler intrinsics for fast bit manipulation
static inline int popcount64(u64 x) { return __builtin_popcountll(x); }
static inline int ctz64(u64 x) { return __builtin_ctzll(x); }
static inline bool has_at_most_one_bit(u64 x) { return (x & (x - 1)) == 0; }

struct Solver
{
    int K, N;
    u64 ALL_MASK;
    // Align to 64 bytes to match cache lines
    alignas(64) std::array<u64, 64> neighbors{};

    Solver(int k, int n) : K(k), N(n)
    {
        int V = k * n;
        ALL_MASK = (V == 64) ? ~0ULL : ((1ULL << V) - 1ULL);

        static constexpr int moves[8][2] = {
            {1, 2}, {1, -2}, {-1, 2}, {-1, -2}, {2, 1}, {2, -1}, {-2, 1}, {-2, -1}};

        for (int r = 0; r < K; r++)
        {
            for (int c = 0; c < N; c++)
            {
                int u = r * N + c;
                u64 m = 0;
                for (auto &mv : moves)
                {
                    int rr = r + mv[0];
                    int cc = c + mv[1];
                    if (rr >= 0 && rr < K && cc >= 0 && cc < N)
                    {
                        m |= (1ULL << (rr * N + cc));
                    }
                }
                neighbors[u] = m;
            }
        }
    }

    // --- Heuristic Pruning ---
    inline bool is_valid_state(u64 rem) const
    {
        if (rem == 0 || has_at_most_one_bit(rem))
            return true;

        int ends = 0;
        u64 temp_rem = rem;

        // 1. Degree Check
        while (temp_rem)
        {
            int v = ctz64(temp_rem);
            temp_rem &= temp_rem - 1;

            u64 n_rem = neighbors[v] & rem;

            // Isolated node -> Impossible
            if (n_rem == 0)
                return false;

            // Degree 1 node -> Endpoint
            if (has_at_most_one_bit(n_rem))
            {
                ends++;
                if (ends > 2)
                    return false;
            }
        }

        // 2. Connectivity Check (Bit-Flood)
        int start_node = ctz64(rem);
        u64 connected = (1ULL << start_node);
        u64 wavefront = connected;

        while (wavefront)
        {
            u64 new_wavefront = 0;
            while (wavefront)
            {
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

    // --- DFS with Forced Move Optimization ---
    u64 dfs(int current, u64 visited) const
    {
        if (visited == ALL_MASK)
            return 1;

        u64 rem = ALL_MASK & ~visited;
        u64 moves_mask = neighbors[current] & rem;

        if (!moves_mask)
            return 0;
        if (!is_valid_state(rem))
            return 0;

        int candidates[8];
        int degrees[8];
        int count = 0;
        int forced_move = -1;

        // Move Generation & Sorting
        while (moves_mask)
        {
            int v = ctz64(moves_mask);
            moves_mask &= moves_mask - 1;

            // Lookahead degree
            int deg = popcount64(neighbors[v] & rem);

            // Leaf Detection
            // If a neighbor has degree 0 in 'rem', it is a leaf.
            // We must visit it immediately or it becomes unreachable.
            if (deg == 0)
            {
                if (forced_move != -1)
                    return 0; // Fork -> Fail
                forced_move = v;
            }

            // Insertion Sort (Warnsdorff's Rule)
            int i = count;
            while (i > 0 && degrees[i - 1] > deg)
            {
                degrees[i] = degrees[i - 1];
                candidates[i] = candidates[i - 1];
                i--;
            }
            candidates[i] = v;
            degrees[i] = deg;
            count++;
        }

        // Execution
        if (forced_move != -1)
        {
            return dfs(forced_move, visited | (1ULL << forced_move));
        }

        u64 total_paths = 0;
        for (int i = 0; i < count; i++)
        {
            total_paths += dfs(candidates[i], visited | (1ULL << candidates[i]));
        }

        return total_paths;
    }
};

// --- Symmetry Handling ---

struct BoardSym
{
    int r, c;
    bool operator<(const BoardSym &o) const
    {
        if (r != o.r)
            return r < o.r;
        return c < o.c;
    }
    bool operator==(const BoardSym &o) const { return r == o.r && c == o.c; }
};

// Map (r,c) to canonical representative and get orbit size
void get_canonical(int r, int c, int K, int N, int &out_r, int &out_c, int &orbit)
{
    BoardSym syms[8];
    int count = 0;

    // D2 Symmetries
    // 1. Identity
    syms[count++] = {r, c};
    // 2. Reflect Vertical (across mid-row)
    syms[count++] = {K - 1 - r, c};
    // 3. Reflect Horizontal (across mid-col)
    syms[count++] = {r, N - 1 - c};
    // 4. Rotate 180
    syms[count++] = {K - 1 - r, N - 1 - c};

    // D4 Symmetries (Square only)
    if (K == N)
    {
        // 5. Transpose (Main Diag)
        syms[count++] = {c, r};
        // 6. Anti-Transpose (Anti Diag)
        syms[count++] = {N - 1 - c, K - 1 - r};
        // 7. Rot 90 CW
        syms[count++] = {c, K - 1 - r};
        // 8. Rot 270 CW
        syms[count++] = {N - 1 - c, r};
    }

    // Find canonical (min)
    BoardSym min_s = syms[0];
    for (int i = 1; i < count; i++)
    {
        if (syms[i] < min_s)
            min_s = syms[i];
    }
    out_r = min_s.r;
    out_c = min_s.c;

    // Orbit Size
    int uniq = 0;
    for (int i = 0; i < count; i++)
    {
        bool seen = false;
        for (int j = 0; j < uniq; j++)
        {
            if (syms[i] == syms[j])
            {
                seen = true;
                break;
            }
        }
        if (!seen)
            syms[uniq++] = syms[i];
    }
    orbit = uniq;
}

u64 knight_hamiltonian_paths(int k, int n)
{
    // Enforce "tall" orientation (small stride)
    // Minimizing N minimizes the bit-distance between neighbors.
    // This standardizes the heuristic search tree stability.
    if (n > k)
        std::swap(k, n);

    if (k * n > 64)
    {
        std::cerr << "Error: Board > 64 squares.\n";
        return 0;
    }

    bool odd_board = (k * n) % 2 != 0;
    const Solver solver(k, n);
    u64 total_directed = 0;

    // To parallelize over canonical start nodes.
    struct Task
    {
        int start_node;
        int weight;
    };
    std::vector<Task> tasks;
    tasks.reserve(k * n);

    for (int r = 0; r < k; r++)
    {
        for (int c = 0; c < n; c++)
        {
            int can_r, can_c, orbit;
            get_canonical(r, c, k, n, can_r, can_c, orbit);

            if (r == can_r && c == can_c)
            {
                // Parity Pruning: Odd boards must start on Majority Color (White)
                if (odd_board)
                {
                    bool is_white = (r + c) % 2 == 0;
                    if (!is_white)
                        continue;
                }
                tasks.push_back(Task{r * n + c, orbit});
            }
        }
    }

#pragma omp parallel for schedule(dynamic) reduction(+ : total_directed)
    for (size_t i = 0; i < tasks.size(); i++)
    {
        total_directed += solver.dfs(tasks[i].start_node, 1ULL << tasks[i].start_node) * tasks[i].weight;
    }

    return total_directed / 2;
}

// --- Driver & Benchmarking ---

void A390833(int k, int n)
{
    auto start = std::chrono::steady_clock::now();
    u64 result = knight_hamiltonian_paths(k, n);
    auto end = std::chrono::steady_clock::now();
    double seconds = std::chrono::duration<double>(end - start).count();

    std::cout << "A(" << k << "," << n << ") = " << result
              << "  \tTime: " << seconds << "s" << std::endl;
}

void benchmark()
{
    for (int k = 3; k <= 12; k++)
    {
        for (int n = 3; n <= 12; n++)
        {
            if (n * k <= 40)
            { // Limit to manageable sizes
                A390833(k, n);
            }
        }
    }
}

int main()
{
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
