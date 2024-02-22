def funcP(A,B,C):
    return (5*A*B) + (5*A*C) + (5*B*C) - (15*A*B*C)

def funcN(A,B,C):
    return 2*A + 2*B + 2*C - A*B - A*C - C*B - A*B*C -2

def r1(a,b,c,x,y,z):
    return 2*a*b + 2*a*c + 2*b*c - 2*a*x - 2*b*x - 2*a*y - 2*c*y - 2*b*z - 2*c*z + 3*x + 3*y + 3*z - a*z - b*y - c*x

for a in range(2):
    for b in range(2):
        for c in range(2):
            print()
            for x in range(2):
                for y in range(2):
                    for z in range(2):
                        print(a,b,c,x,y,z,"=", r1(a,b,c,x,y,z))