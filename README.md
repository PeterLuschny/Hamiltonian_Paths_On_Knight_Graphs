# Hamiltonian Paths on Knight Graphs

Count undirected Hamiltonian knight paths on k x n chessboards. See OEIS A390833.

Uses backtracking DFS with Warnsdorff's rule and advanced pruning techniques.

Make: g++ -O3 -march=native -fopenmp -std=c++20 knights.cpp -o knights

# The challenge
 
Let's turn this into a little challenge! 


![Challenge](https://github.com/PeterLuschny/Hamiltonian_Paths_On_Knight_Graphs/blob/main/FightingKnights.png)


The benchmark is:
```
    // Limit to manageable sizes
    for (k = 3; k <= 12; k++) 
       for (n = 3; n <= 12; n++) 
          if (n * k <= 40) A390833(k, n);
```
To keep the code simple and clear, we also require < 300 lines of code (not counting comment lines).

Can you improve the following runtimes? (Times measured on a standard laptop).

| Board | Paths  | Time (s)  |
|---    |---:    |---:       |
| A(3,3)| 0      | 0.0051301 |
| A(3,4)| 8      | 3.2e-05   |
| A(3,5)| 0      | 1.14e-05  |
| A(3,6)| 0      | 2.07e-05  |
| A(3,7)| 52     | 8.82e-05  |
| A(3,8)| 396    | 0.000374  |
| A(3,9)| 560    | 0.0007934 |
| A(3,10)| 3048  | 0.0041759 |
| A(3,11)| 10672 | 0.0197895 |
| A(3,12)| 57248 | 0.0739956 |
| A(4,3)| 8      | 8.87e-05  |
| A(4,4)| 0      | 0.0012474 |
| A(4,5)| 82     | 0.0003997 |
| A(4,6)| 744    | 0.0015467 |
| A(4,7)| 6378   | 0.0143789 |
| A(4,8)| 31088  | 0.095957  |
| A(4,9)| 189688 | 0.820108  |
| A(4,10)| 1213112 | 6.57979 |
| A(5,3)| 0      | 0.0002406 |
| A(5,4)| 82     | 0.0001471 |
| A(5,5)| 864    | 0.0013945 |
| A(5,6)| 18784  | 0.028394  |
| A(5,7)| 622868 | 0.979126  |
| A(5,8)| 18061054 | 26.1308 |
| A(6,3)| 0      | 2.65e-05  |
| A(6,4)| 744    | 0.0013012 |
| A(6,5)| 18784  | 0.0264616 |
| A(6,6)| 3318960 | 2.63667  |
| A(7,3)| 52     | 8.28e-05  |
| A(7,4)| 6378   | 0.0136763 |
| A(7,5)| 622868 | 0.983143  |
| A(8,3)| 396    | 0.0004183 |
| A(8,4)| 31088  | 0.093957  |
| A(8,5)| 18061054 | 25.4271 |
| A(9,3)| 560    | 0.0009313 |
| A(9,4)| 189688 | 0.780557  |
| A(10,3)| 3048  | 0.0037522 |
| A(10,4)| 1213112 | 6.56309 |
| A(11,3)| 10672 | 0.0140439 |
| A(12,3)| 57248 | 0.065328  |
|        | Total | 71.3712   |

If you have found a more efficient solution, share it with us and send us a pull request.
