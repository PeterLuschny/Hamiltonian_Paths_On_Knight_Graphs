# Hamiltonian Paths on Knight Graphs

Count undirected Hamiltonian knight paths on k x n chessboards. See OEIS A390833.

Uses backtracking DFS with Warnsdorff's rule and advanced pruning techniques. 

Make: g++ -O3 -march=native -fopenmp -std=c++20 knights.cpp -o knights

# The challenge
 
Let's turn this into a little challenge! The benchmark is:
```
    // Limit to manageable sizes
    for (k = 3; k < 12; k++) 
       for (n = 3; n < 12; n++) 
          if (n * k <= 40) A390833(k, n);
```
To keep the code simple and clear, we also require < 200 lines of code (not counting comment lines).

Who can improve the following runtimes? (Measured on a standard 2016 laptop).


| Board | Paths | Time (s) |
|---|---:|---:|
| A(3,3) | 0 | 0.0040712 |
| A(3,4) | 8 | 3.25e-05 |
| A(3,5) | 0 | 1.14e-05 |
| A(3,6) | 0 | 2.31e-05 |
| A(3,7) | 52 | 8.47e-05 |
| A(3,8) | 396 | 0.000403 |
| A(3,9) | 560 | 0.0009808 |
| A(3,10) | 3048 | 0.0041568 |
| A(3,11) | 10672 | 0.0189267 |
| A(4,3) | 8 | 0.000235 |
| A(4,4) | 0 | 3.51e-05 |
| A(4,5) | 82 | 0.0027262 |
| A(4,6) | 744 | 0.0119047 |
| A(4,7) | 6378 | 0.0185291 |
| A(4,8) | 31088 | 0.140034 |
| A(4,9) | 189688 | 1.00172 |
| A(4,10) | 1213112 | 7.46396 |
| A(5,3) | 0 | 0.0033393 |
| A(5,4) | 82 | 0.0020133 |
| A(5,5) | 864 | 0.005744 |
| A(5,6) | 18784 | 0.0384103 |
| A(5,7) | 622868 | 1.06086 |
| A(5,8) | 18061054 | 27.3598 |
| A(6,3) | 0 | 0.0001351 |
| A(6,4) | 744 | 0.0015218 |
| A(6,5) | 18784 | 0.0281101 |
| A(6,6) | 3318960 | 2.82377 |
| A(7,3) | 52 | 9.22e-05 |
| A(7,4) | 6378 | 0.0135029 |
| A(7,5) | 622868 | 1.01951 |
| A(8,3) | 396 | 0.0003808 |
| A(8,4) | 31088 | 0.100319 |
| A(8,5) | 18061054 | 27.9002 |
| A(9,3) | 560 | 0.0008884 |
| A(9,4) | 189688 | 0.831437 |
| A(10,3) | 3048 | 0.004737 |
| A(10,4) | 1213112 | 7.21754 |
| A(11,3) | 10672 | 0.0163133 |
|   | Total Time | 77.0983 |

