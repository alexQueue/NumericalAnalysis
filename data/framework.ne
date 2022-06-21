# Node-Edges-File .ne
NODES
0.1 0.1
0.0 1.0
1.0 0.0
1.0 1.0
2.0 0.0
2.0 1.0
EDGES
1 2
2 3
1 3
2 4
3 4
5 6
4 6
4 5
TYPE
1 FIXED
2 FORCE [100 10] [0] # Q first M second
3 FREE
4 FREE
5 MOVABLE [1 0]
6 MOVABLE [0 1]
PARAMETERS
E 2.0
I 2.0x + 3
A 3.0x + x^2
mu 1.0