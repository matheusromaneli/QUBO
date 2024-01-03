def funcP(A,B,C):
    return (5*A*B) + (5*A*C) + (5*B*C) - (15*A*B*C)

def funcN(A,B,C):
    return 2*A + 2*B + 2*C - A*B - A*C - C*B - A*B*C -2

for a in range(2):
    for b in range(2):
        for c in range(2):
            print(a,b,c,"=",funcN(a,b,c), funcP(a,b,c))