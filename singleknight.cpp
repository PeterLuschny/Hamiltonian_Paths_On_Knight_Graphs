// Find and print the FIRST Closed Hamiltonian Knight's Tour.
// Uses backtracking DFS with Warnsdorff's rule and advanced pruning techniques.
// Compiled from various AI-sources (vibe coding), Peter Luschny, December 2025

#include <iostream>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <iomanip>

using u64 = std::uint64_t;

// Compiler intrinsics for fast bit manipulation
static inline int popcount64(u64 x) { return __builtin_popcountll(x); }
static inline int ctz64(u64 x) { return __builtin_ctzll(x); }

struct Solver
{
    int K, N, V;
    u64 ALL_MASK;
    std::vector<u64> neighbors;
    std::vector<int> path; // Stores the solution sequence

    Solver(int k, int n) : K(k), N(n), V(k * n)
    {
        // Safe resizing based on actual volume
        neighbors.resize(V, 0);

        // shift by 64 is undefined behavior, handle gracefully
        ALL_MASK = (V >= 64) ? ~0ULL : ((1ULL << V) - 1ULL);

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
                    {
                        int neighbor_idx = nr * N + nc;
                        // Only add neighbor if it fits in 64-bit mask
                        if (neighbor_idx < 64)
                        {
                            mask |= (1ULL << neighbor_idx);
                        }
                    }
                }
                neighbors[r * N + c] = mask;
            }
        }
    }

    // Returns true if a solution is found
    bool dfs(int curr, u64 visited, int target)
    {
        path.push_back(curr);

        // Base Case: All nodes visited except (0,0)
        // We are currently at 'curr'. If we can jump to 'target', the cycle is closed.
        if (popcount64(visited) == V - 1)
        {
            if (neighbors[curr] & (1ULL << target))
            {
                // Complete the path to the target neighbor
                path.push_back(target);
                return true;
            }
            path.pop_back();
            return false;
        }

        u64 rem = ALL_MASK & ~visited;
        if (rem == 0)
        {
            path.pop_back();
            return false;
        }

        // --- PRUNING ---
        int forced_move = -1;
        u64 t_rem = rem;
        while (t_rem)
        {
            int v = ctz64(t_rem);
            t_rem &= t_rem - 1;

            int d = popcount64(neighbors[v] & rem);
            bool is_reachable = (neighbors[v] & (1ULL << curr));
            int eff_d = d + (is_reachable ? 1 : 0);

            if (v == target)
            {
                if (eff_d < 1)
                {
                    path.pop_back();
                    return false;
                }
            }
            else
            {
                if (eff_d < 2)
                {
                    path.pop_back();
                    return false;
                }
                if (eff_d == 2 && is_reachable)
                {
                    if (forced_move != -1 && forced_move != v)
                    {
                        path.pop_back();
                        return false;
                    }
                    forced_move = v;
                }
            }
        }

        u64 moves = neighbors[curr] & rem;
        if (forced_move != -1)
        {
            if (!(moves & (1ULL << forced_move)))
            {
                path.pop_back();
                return false;
            }
            moves = (1ULL << forced_move);
        }

        if (!moves)
        {
            path.pop_back();
            return false;
        }

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

        // Sort: fewest onward moves first
        for (int i = 0; i < count; i++)
        {
            for (int j = i + 1; j < count; j++)
            {
                if (cands[j].d < cands[i].d)
                    std::swap(cands[i], cands[j]);
            }
        }

        // --- RECURSION ---
        for (int i = 0; i < count; i++)
        {
            if (dfs(cands[i].v, visited | (1ULL << cands[i].v), target))
                return true; // Found! Bubble up success
        }

        path.pop_back(); // Backtrack
        return false;
    }

    void solve_first()
    {
        // Setup start node (0, 0) and its two neighbors
        int start_node = 0;
        u64 start_neighs = neighbors[start_node];
        std::vector<int> adj;
        while (start_neighs)
        {
            adj.push_back(ctz64(start_neighs));
            start_neighs &= start_neighs - 1;
        }

        if (adj.size() < 2)
        {
            std::cout << "No solution (Corner blocked).\n";
            return;
        }

        // Path must start at (0,0) -> adj[0] -> ... -> adj[1] -> (0,0)
        // We start the DFS at adj[0] and aim for adj[1].
        int first_step = adj[0];
        int target_step = adj[1];

        path.clear();
        path.push_back(start_node);

        u64 initial_visited = (1ULL << start_node) | (1ULL << first_step);

        if (dfs(first_step, initial_visited, target_step))
        {
            print_solution_algebraic();
        }
        else
        {
            std::cout << "No solution found.\n";
        }
    }

    void print_solution_algebraic()
    {
        std::cout << "Closed Knight's Tour found on " << K << "x" << N << " board:\n\n";
        std::cout << "Algebraic notation:\n";
        for (size_t i = 0; i < path.size(); i++)
        {
            int u = path[i];
            int r = u / N;
            int c = u % N;
            std::cout << (char)('a' + c) << (r + 1);
            if (i < path.size() - 1)
                std::cout << " -> ";
        }
        std::cout << " -> a1 (Loop closed)\n\n";

        std::cout << "Board Visualization:\n";
        std::vector<int> grid(V);
        for (size_t i = 0; i < path.size(); ++i)
            grid[path[i]] = i + 1;

        std::cout << "   ";
        for (int c = 0; c < N; c++)
            std::cout << " " << (char)('a' + c) << " ";
        std::cout << "\n";

        for (int r = 0; r < K; r++)
        {
            std::cout << std::setw(2) << (r + 1) << " ";
            for (int c = 0; c < N; c++)
            {
                std::cout << std::setw(2) << grid[r * N + c] << " ";
            }
            std::cout << "\n";
        }
    }
};

int findtour(int k, int n)
{
    // Check constraints BEFORE creating Solver ---
    // Parity Check
    if ((k * n) % 2 != 0)
    {
        std::cout << "No solution (Odd board size: " << k * n << " squares).\n";
        return 0;
    }

    // Size Limit Check (Implementation limitation)
    if (k * n > 64)
    {
        std::cout << "Error: Board size (" << k * n << ") exceeds the 64-square limit of this bitmask solver.\n";
        return 0;
    }

    Solver s(n, k);
    s.solve_first();
    return 1;
}

int main()
{
    // Example: 9x6 board
    int k = 9;
    int n = 6;
    return findtour(k, n);
}
