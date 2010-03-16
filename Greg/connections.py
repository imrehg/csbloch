from scipy import *

conn = array([[0, 0, 1],
              [0, 0, 1],
              [-1, -1, 0]])

def getconnxjj(i):

    print "x(%d, %d):"%(i, i)
    for l in range(1, 4):
        if ((conn[j-1, l-1] <> 0) & (l <> j)):
            print "  M(%d,%d) y(%d, %d)" %(l, j, j, l)

def getconnx(j, k):
    print "x(%d, %d):"%(j, k)
    for l in range(1, 4):
        if ((conn[j-1, l-1] <> 0) & (l not in (j, k))):
            print "- M(%d,%d) y(%d, %d)" %(j, l, l, k)
        if ((conn[l-1, k-1] <> 0) & (l not in (j, k))):
            print "  M(%d,%d) y(%d, %d)" %(l, k, j, l)

def getconny(j, k):
    print "y(%d, %d):"%(j, k)
    for l in range(1, 4):
        if (conn[j-1, l-1] <> 0):
            print "  M(%d,%d) x(%d, %d)" %(j, l, l, k)
        if (conn[l-1, k-1] <> 0):
            print "- M(%d,%d) x(%d, %d)" %(l, k, j, l)


print "#Note: x(i, j) = x(j, i)"
print "#Note y(i, j) = -y(j, i); i <> j"


for j in range(1,4):
    getconnxjj(j)

for j in range(1,3):
    for k in range(j+1,4):
        getconnx(j, k)

for j in range(1,3):
    for k in range(j+1,4):
        getconny(j, k)
