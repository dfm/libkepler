import numpy as np
from numpy.polynomial import Polynomial
from scipy.optimize import minimize_scalar, brentq


def polynomial_approximation(
    scalar_func, min_x, max_x, poly_order, maxiter=100, tol=5e-13, init_points=None
):
    if init_points is None:
        maxerr, points, coeffs = init_chebyshev(scalar_func, poly_order, (min_x, max_x))
        if maxerr < tol:
            return maxerr, coeffs
    else:
        points = init_points
    for i in range(maxiter):
        maxerr, points, coeffs = iter_remez(scalar_func, points, (min_x, max_x))
        if maxerr < tol:
            break
    return maxerr, points, coeffs


def init_chebyshev(func, order, domain):
    mn, mx = sorted(domain)
    unknowns = order + 2
    terms = order + 1
    zeros = np.zeros(unknowns + 1)
    for n in range(terms):
        cheb = 0.5 * (1 + np.cos((2 * terms - 1 - 2 * n) * np.pi / (2 * terms)))
        cheb *= mx - mn
        cheb += mn
        zeros[n + 1] = cheb
    zeros[0] = mn
    zeros[-1] = mx

    A = np.vander(zeros[1:-1], order + 1, True)
    b = func(zeros[1:-1])
    soln = np.linalg.solve(A, b)
    poly = Polynomial(soln)

    def resid(x):
        return -np.abs(func(x) - poly(x))

    maxerr = 0.0
    points = np.zeros(unknowns)
    for n in range(unknowns):
        result = minimize_scalar(
            resid, bounds=(zeros[n], zeros[n + 1]), method="bounded"
        )
        points[n] = result.x
        if -result.fun > maxerr:
            maxerr = -result.fun

    return maxerr, points, soln


def iter_remez(func, points, domain):
    mn, mx = sorted(domain)
    b = func(points)
    A = np.vander(points, len(points), True)
    A[:, -1] = 1 - 2 * (np.arange(len(points)) % 2)

    soln = np.linalg.solve(A, b)
    poly = Polynomial(soln[:-1])

    def resid(x):
        return func(x) - poly(x)

    zeros = np.zeros(len(points) + 1)
    zeros[0] = mn
    zeros[-1] = mx
    for n in range(len(points) - 1):
        result = brentq(resid, points[n], points[n + 1])
        zeros[n + 1] = result

    def resid(x):
        return -np.abs(func(x) - poly(x))

    maxerr = 0.0
    points = np.zeros_like(points)
    for n in range(len(points)):
        result = minimize_scalar(
            resid, bounds=(zeros[n], zeros[n + 1]), method="bounded"
        )
        points[n] = result.x
        if -result.fun > maxerr:
            maxerr = -result.fun

    return maxerr, points, soln[:-1]
