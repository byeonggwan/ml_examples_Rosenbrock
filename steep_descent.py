import matplotlib.pyplot as plt
import time

def f(x):
    x1 = x[0]
    x2 = x[1]
    return 100*(x2**2-2*x2*x1**2+x1**4)+(1-2*x2+x2**2)

def pt_gen(x, positive_sign):
    x1 = x[0]
    x2 = x[1]
    if positive_sign == 1:
        return [400*x1**3-400*x1*x2, 202*x2-200*x1**2-2]
    return [-400*x1**3+400*x1*x2, -202*x2+200*x1**2+2]


def definite_pt(pt):
    return pt[0]*pt[0]+pt[1]*pt[1]

# Backtracking Line Search Algorithm
def bls(x, a, t, c, pt):
    gradient = pt_gen(x, 1)
    while f([ x[0] + a*pt[0] , x[1] + a*pt[1] ]) > f(x) + c*a*(gradient[0] * pt[0] + gradient[1] * pt[1]):
        a = t*a
    return a

# Vectors
x = [1.8, 1.8]
pt = pt_gen(x, -1)
pt_o = pt

# Factors. You may change this part.
c = 0.5
alpha = 0.5
tau = 0.5
T = 0.0000000001

# For matplotlib
plotx = []
ploty = []
plotx1 = []
plotx2 = []

i = 0 # iterator

print("Timer start.")
start = time.time()
while definite_pt(pt)/definite_pt(pt_o) > T: # origin : ||gradient f(x)|| / ||gradient f(x0)|| < T for contant T
    plotx.append(i)
    ploty.append(f(x))
    plotx1.append(x[0])
    plotx2.append(x[1])

    pt = pt_gen(x, -1)
    #print ( "iter: " + str(i) + "    x1: " + str(round(x[0],4)) + "    x2: " + str(round(x[1],4)) + "    f(x): " + str(round(f(x),4)))

    alpha = bls(x, alpha, tau, c, pt_gen(x, -1))
    x[0] = x[0] + alpha*pt[0] # x_1 = x_0 + alpha_0 * p_0
    x[1] = x[1] + alpha*pt[1]
    

    i += 1

print("Timer end. Each iteration time is : ", (time.time() - start)/i)

plt.plot(plotx[:50], ploty[:50])
plt.savefig('gd.png')
plt.clf()
plt.scatter(plotx1, plotx2)
plt.savefig('gd_scatter.png')