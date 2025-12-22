# Hamiltonian Paths on Knight Graphs

Count undirected Hamiltonian knight paths on k x n chessboards. See OEIS A390833.

Uses backtracking DFS with Warnsdorff's rule and advanced pruning techniques. 

Make: g++ -O3 -march=native -fopenmp -std=c++20 knights.cpp -o knights

# The challenge
 
Let's turn this into a little challenge! 


![Challenge](pictures/FightingKnights.pgn)


The benchmark is:
```
    // Limit to manageable sizes
    for (k = 3; k < 12; k++) 
       for (n = 3; n < 12; n++) 
          if (n * k <= 40) A390833(k, n);
```
To keep the code simple and clear, we also require < 200 lines of code (not counting comment lines).

Can you improve the following runtimes? (Times measured on a standard laptop).

| Board | Paths | Time (s) |
|---|---:|---:|
| A(3,3) | 0 | 0.0037753 |
| A(3,4) | 8 | 2.36e-05 |
| A(3,5) | 0 | 1.1e-05 |
| A(3,6) | 0 | 2.19e-05 |
| A(3,7) | 52 | 9.2e-05 |
| A(3,8) | 396 | 0.0003373 |
| A(3,9) | 560 | 0.000828 |
| A(3,10) | 3048 | 0.0049929 |
| A(3,11) | 10672 | 0.0148218 |
| A(4,3) | 8 | 0.0001479 |
| A(4,4) | 0 | 3.66e-05 |
| A(4,5) | 82 | 0.0001478 |
| A(4,6) | 744 | 0.0012874 |
| A(4,7) | 6378 | 0.0128137 |
| A(4,8) | 31088 | 0.0944202 |
| A(4,9) | 189688 | 0.854279 |
| A(4,10) | 1213112 | 6.90553 |
| A(5,3) | 0 | 0.0035265 |
| A(5,4) | 82 | 0.0001744 |
| A(5,5) | 864 | 0.0013642 |
| A(5,6) | 18784 | 0.030783 |
| A(5,7) | 622868 | 0.948344 |
| A(5,8) | 18061054 | 25.603 |
| A(6,3) | 0 | 2.71e-05 |
| A(6,4) | 744 | 0.0012968 |
| A(6,5) | 18784 | 0.0268025 |
| A(6,6) | 3318960 | 2.56409 |
| A(7,3) | 52 | 8.75e-05 |
| A(7,4) | 6378 | 0.0126201 |
| A(7,5) | 622868 | 1.01799 |
| A(8,3) | 396 | 0.0047561 |
| A(8,4) | 31088 | 0.104557 |
| A(8,5) | 18061054 | 25.8258 |
| A(9,3) | 560 | 0.0008302 |
| A(9,4) | 189688 | 0.829167 |
| A(10,3) | 3048 | 0.0080581 |
| A(10,4) | 1213112 | 6.91208 |
| A(11,3) | 10672 | 0.0187905 |
|   | Total Time | 71.8105 |
