import matplotlib.pyplot as plt

def f(x):
    x1 = x[0]
    x2 = x[1]
    return 100*(x2**2-2*x2*x1**2+x1**4)+(1-2*x2+x2**2)

def gradient_f(x):
    x1 = x[0]
    x2 = x[1]
    return [400*x1**3-400*x1*x2, 202*x2-200*x1**2-2]

def H_f(x):
    x1 = x[0]
    x2 = x[1]
    return [ [1200*x1*x1-400*x2, -400*x1], [-400*x1, 202]]

# 2x2 matrix inverse
def inverse(f): 
    a, b, c, d = f[0][0], f[0][1], f[1][0], f[1][1]
    fac = 1 / (a*d - b*c)
    return [ [fac*d, -fac*b], [-fac*c, fac*a] ]

def definite_pt(pt):
    return pt[0]*pt[0]+pt[1]*pt[1]

# For solve H * gradient
def multiple_2x2_2x1(H, f): 
    a = H[0][0]*f[0] + H[0][1]*f[1]
    b = H[1][0]*f[0] + H[1][1]*f[1]
    return [a, b]

# Backtracking Line Search Algorithm
def bls(x, a, t, c, pt):
    gradient = gradient_f(x)
    while f([ x[0] + a*pt[0] , x[1] + a*pt[1] ]) > f(x) + c*a*(gradient[0] * pt[0] + gradient[1] * pt[1]):
        a = t*a
    return a

# Vectors
x = [1.8, 1.8]
ft = f(x)
gradient = gradient_f(x)

# Original vectors
ft_o = ft
gradient_o = gradient

# Factors
alpha = 0.99
tau = 0.7
c = 0.7
T = 0.000000001

# For matplotlib
plotx = []
ploty = []
plotx1 = []
plotx2 = []

i = 0 # iterator

while definite_pt(gradient)/definite_pt(gradient_o) > T: # origin : ||gradient f(x)|| / ||gradient f(x0)|| < T for contant T
    plotx.append(i)
    ploty.append(f(x))
    plotx1.append(x[0])
    plotx2.append(x[1])

    gradient = gradient_f(x)
    H = H_f(x)
    H_inverse = inverse(H)
    print ( "iter: " + str(i) + "    x1: " + str(round(x[0],4)) + "    x2: " + str(round(x[1],4)))
    
    x[0] = x[0] - alpha*multiple_2x2_2x1(H_inverse, gradient)[0]
    x[1] = x[1] - alpha*multiple_2x2_2x1(H_inverse, gradient)[1]
    alpha = bls(x, alpha, tau, c, [-gradient_f(x)[0], -gradient_f(x)[1]])
    
    i += 1

plt.plot(plotx, ploty)
plt.savefig('nm.png')
plt.clf()
plt.scatter(plotx1, plotx2)
plt.savefig('nm_scatter.png')