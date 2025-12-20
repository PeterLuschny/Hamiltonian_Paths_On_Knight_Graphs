# Hamiltonian_Paths_On_Knight_Graphs

Count undirected Hamiltonian knight paths on k x n chessboards.
Uses backtracking DFS with Warnsdorff's rule and advanced pruning techniques.
See OEIS A390833.

Make:
sudo pacman -S gcc clang make cmake libomp
g++ -O3 -march=native -fopenmp -std=c++17 knights.cpp -o knights

Benchmark:
A(3,3) = 0        Time: 0.0040712s
A(3,4) = 8        Time: 3.25e-05s
A(3,5) = 0        Time: 1.14e-05s
A(3,6) = 0        Time: 2.31e-05s
A(3,7) = 52       Time: 8.47e-05s
A(3,8) = 396      Time: 0.000403s
A(3,9) = 560      Time: 0.0009808s
A(3,10) = 3048    Time: 0.0041568s
A(3,11) = 10672   Time: 0.0189267s
A(4,3) = 8        Time: 0.000235s
A(4,4) = 0        Time: 3.51e-05s
A(4,5) = 82       Time: 0.0027262s
A(4,6) = 744      Time: 0.0119047s
A(4,7) = 6378     Time: 0.0185291s
A(4,8) = 31088    Time: 0.140034s
A(4,9) = 189688   Time: 1.00172s
A(4,10) = 1213112 Time: 7.46396s
A(5,3) = 0        Time: 0.0033393s
A(5,4) = 82       Time: 0.0020133s
A(5,5) = 864      Time: 0.005744s
A(5,6) = 18784    Time: 0.0384103s
A(5,7) = 622868   Time: 1.06086s
A(5,8) = 18061054 Time: 27.3598s
A(6,3) = 0        Time: 0.0001351s
A(6,4) = 744      Time: 0.0015218s
A(6,5) = 18784    Time: 0.0281101s
A(6,6) = 3318960  Time: 2.82377s
A(7,3) = 52       Time: 9.22e-05s
A(7,4) = 6378     Time: 0.0135029s
A(7,5) = 622868   Time: 1.01951s
A(8,3) = 396      Time: 0.0003808s
A(8,4) = 31088    Time: 0.100319s
A(8,5) = 18061054 Time: 27.9002s
A(9,3) = 560      Time: 0.0008884s
A(9,4) = 189688   Time: 0.831437s
A(10,3) = 3048    Time: 0.004737s
A(10,4) = 1213112 Time: 7.21754s
A(11,3) = 10672   Time: 0.0163133s
Total Time: 77.0983s
