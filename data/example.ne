NODES
0.0 0.0
1.0 0.0
1.0 1.0
0.0 1.0
EDGES
1 2
2 3 params=[x,1,2,3] gp=5
3 4 gp=12
1 3 params=[1.5,1,1,x-1]
1 4
TYPE
1 FIXED
2 FREE
3 FORCE [-1.0 -1.0] [0]
4 MOVABLE [1 0]
PARAMETERS
E 10x^2-sin(x)
I 1
A 1
mu 1
