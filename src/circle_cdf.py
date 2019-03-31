from sympy import symbols, limit, diff, pi, asin, sqrt

w = symbols("w")

n = 11

print("g[{}]={}".format(n, limit(diff((pi*w/(asin(w) + w*sqrt(1 - w**2)))**n, w, n - 1), w, 0)))

