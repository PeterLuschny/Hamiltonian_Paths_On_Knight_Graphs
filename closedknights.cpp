// Count undirected Hamiltonian knight paths on k x n chessboards.
// Uses backtracking DFS with Warnsdorff's rule and advanced pruning techniques.
// Compiled from various AI-sources (vibe coding), Peter Luschny, December 2025
// OEIS A392000, but note that this code was **not** used there.

/**
 * Count undirected closed Hamiltonian Knight's tours (cycles) on k x n boards.
 * 1. A cycle must visit (0,0). (0,0) has exactly two neighbors: v1 and v2.
 * 2. We find all Hamiltonian paths in (Board - {(0,0)}) that start at v1 and end at v2.
 * 3. This count corresponds exactly to the number of unique undirected cycles.
 */

#include <iostream>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <chrono>
#include <omp.h>

using u64 = std::uint64_t;

// Compiler intrinsics for fast bit manipulation
static inline int popcount64(u64 x) { return __builtin_popcountll(x); }
static inline int ctz64(u64 x) { return __builtin_ctzll(x); }
static inline bool has_at_most_one_bit(u64 x) { return (x & (x - 1)) == 0; }

struct Solver
{
    int K, N, V;
    u64 ALL_MASK;
    u64 neighbors[64];

    Solver(int k, int n) : K(k), N(n), V(k * n)
    {
        ALL_MASK = (V == 64) ? ~0ULL : ((1ULL << V) - 1ULL);
        int dr[] = {2, 2, 1, 1, -2, -2, -1, -1};
        int dc[] = {1, -1, 2, -2, 1, -1, 2, -2};

        for (int r = 0; r < K; r++)
        {
            for (int c = 0; c < N; c++)
            {
                u64 mask = 0;
                for (int i = 0; i < 8; i++)
                {
                    int nr = r + dr[i], nc = c + dc[i];
                    if (nr >= 0 && nr < K && nc >= 0 && nc < N)
                        mask |= (1ULL << (nr * N + nc));
                }
                neighbors[r * N + c] = mask;
            }
        }
    }

    /**
     * Advanced Backtracking DFS
     * @param curr The knight's current position.
     * @param visited Bitmask of visited squares.
     * @param target The neighbor of (0,0) we must end at.
     */
    u64 dfs(int curr, u64 visited, int target) const
    {
        // Base Case: All nodes except (0,0) visited.
        // We are at 'curr'. If 'curr' can jump to 'target', we've closed the cycle.
        if (popcount64(visited) == V - 1)
        {
            return (neighbors[curr] & (1ULL << target)) ? 1 : 0;
        }

        u64 rem = ALL_MASK & ~visited;

        // --- PRUNING ---
        int forced_move = -1;
        u64 t_rem = rem;
        while (t_rem)
        {
            int v = ctz64(t_rem);
            t_rem &= t_rem - 1;

            u64 v_neighs_in_rem = neighbors[v] & rem;
            int d = popcount64(v_neighs_in_rem);
            bool is_reachable_from_curr = (neighbors[v] & (1ULL << curr));

            // Effective degree = degree in unvisited subgraph + connection to current
            int eff_d = d + (is_reachable_from_curr ? 1 : 0);

            if (v == target)
            {
                // Target must have at least 1 way to be reached (exit to (0,0) is implicit)
                if (eff_d < 1)
                    return 0;
            }
            else
            {
                // Every other node must have degree 2 (one in, one out)
                if (eff_d < 2)
                    return 0;
                // If a neighbor of curr has only one exit left in rem, we MUST take it now
                if (eff_d == 2 && is_reachable_from_curr)
                {
                    if (forced_move != -1 && forced_move != v)
                        return 0; // Contradiction
                    forced_move = v;
                }
            }
        }

        u64 moves = neighbors[curr] & rem;
        if (forced_move != -1)
        {
            if (!(moves & (1ULL << forced_move)))
                return 0;
            moves = (1ULL << forced_move);
        }

        if (!moves)
            return 0;

        // --- WARNSDORFF'S RULE ---
        struct Cand
        {
            int v, d;
        };
        Cand cands[8];
        int count = 0;
        while (moves)
        {
            int v = ctz64(moves);
            moves &= moves - 1;
            cands[count++] = {v, (int)popcount64(neighbors[v] & rem)};
        }

        // Sort moves by accessibility (fewest onward options first)
        for (int i = 0; i < count; i++)
        {
            for (int j = i + 1; j < count; j++)
            {
                if (cands[j].d < cands[i].d)
                    std::swap(cands[i], cands[j]);
            }
        }

        u64 total = 0;
        for (int i = 0; i < count; i++)
        {
            total += dfs(cands[i].v, visited | (1ULL << cands[i].v), target);
        }
        return total;
    }

    u64 solve()
    {
        // Parity: Cycles impossible on odd-sized boards
        if (V % 2 != 0)
            return 0;

        int start_node = 0; // Fix (0,0)
        u64 start_neighs = neighbors[start_node];
        std::vector<int> adj;
        u64 temp = start_neighs;
        while (temp)
        {
            adj.push_back(ctz64(temp));
            temp &= temp - 1;
        }

        if (adj.size() < 2)
            return 0;

        // To parallelize, we take the first step from one of the neighbors
        int v1 = adj[0];
        int target = adj[1];
        u64 visited_init = (1ULL << start_node) | (1ULL << v1);

        u64 moves_from_v1 = neighbors[v1] & (ALL_MASK & ~visited_init);
        std::vector<int> first_steps;
        while (moves_from_v1)
        {
            first_steps.push_back(ctz64(moves_from_v1));
            moves_from_v1 &= moves_from_v1 - 1;
        }

        u64 total_cycles = 0;
#pragma omp parallel for reduction(+ : total_cycles) schedule(dynamic)
        for (size_t i = 0; i < first_steps.size(); i++)
        {
            int next_v = first_steps[i];
            total_cycles += dfs(next_v, visited_init | (1ULL << next_v), target);
        }

        return total_cycles;
    }
};

// --- Driver & Benchmarking ---

void A392000(int k, int n)
{
    auto start = std::chrono::steady_clock::now();
    Solver s(k, n);
    u64 result = s.solve();
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
            // Limit to manageable sizes
            if (n * k <= 40)
            {
                A392000(k, n);
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
g++ -O3 -march=native -fopenmp -std=c++20 closedknights.cpp -o closedknights
./closedknights
*/
