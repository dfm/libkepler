from functools import partial
import kepler
import numpy as np
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt

from remez import polynomial_approximation


def approximate(name, target_func, transform, domain, max_order=10, tol=5e-15):
    for order in range(3, max_order + 1):
        err, _, coeffs = polynomial_approximation(
            target_func, *transform(domain), order, tol=tol
        )
        if err < tol:
            break

    print(name, order, err)
    poly = Polynomial(coeffs)
    approx = lambda x: poly(transform(x))

    fig, axes = plt.subplots(2, 1, sharex=True)
    x = np.linspace(*domain, 1000)
    ax = axes[0]
    ax.plot(x, target_func(transform(x)), label=name)
    ax.plot(x, approx(x), "--", label="approx")
    ax.legend()
    ax.set_ylabel("y")
    ax.set_title(f"{name}; order = {order}; max error = {err:.2e}")

    ax = axes[1]
    ax.plot(x, approx(x) - target_func(transform(x)))
    ax.set_ylabel("error")
    ax.set_xlabel("x")

    return fig, err


#
# Approximate sin(x) in the range [0, pi/4)
#
# In this case it's useful to factor out a linear term and then reparameterize
# in terms of x^2
#
def target_func(x_sq):
    x = np.sqrt(x_sq) / np.pi
    return np.sinc(x)


def transform(x):
    return np.square(x)


domain = [0.0, 0.25 * np.pi]
fig, _ = approximate("sin(x) / x", target_func, transform, domain, tol=1e-15)
fig.savefig("sin.png", bbox_inches="tight")
plt.close(fig)

#
# Approximate cos(x) in the range [0, pi/4)
#
# Reparameterize in terms of x^2
#
def target_func(x_sq):
    x = np.sqrt(x_sq)
    return np.cos(x)


def transform(x):
    return np.square(x)


domain = [0.0, 0.25 * np.pi]
fig, _ = approximate("cos(x)", target_func, transform, domain, tol=5e-15)
fig.savefig("cos.png", bbox_inches="tight")
plt.close(fig)


#
# Approximate Kepler's equation for a set of eccentricities
#
# Reparameterize in terms of x^2
#
def mikkola_starter(M, ecc, correct=True):
    factor = 1.0 / (4.0 * ecc + 0.5)
    alpha = (1.0 - ecc) * factor
    alpha3 = alpha * alpha * alpha

    beta = 0.5 * M * factor
    z = np.cbrt(beta + np.sqrt(beta * beta + alpha3))
    s = z - alpha / z
    if correct:
        s -= 0.078 * s**5.0 / (1.0 + ecc)
    return M + ecc * s * (3.0 - 4.0 * s * s)


def target_func(ecc, x, starter=False, correct=False):
    if starter:
        return kepler.solve(np.pi * x, ecc) - mikkola_starter(
            np.pi * x, ecc, correct=correct
        )
    else:
        return kepler.solve(np.pi * x, ecc)


transform = lambda x: x

eccs = np.linspace(0.0, 1.0, 25)[1:-1]

for starter, correct in [(False, False), (True, False)]:  # , (True, True)
    suffix = ""
    if starter:
        suffix += "_starter"
        if correct:
            suffix += "_correct"
    print(suffix)
    errs = np.zeros_like(eccs)
    for n, ecc in enumerate(eccs):
        domain = [0.0, 1.0]
        fig, err = approximate(
            f"kepler; e = {ecc:.2f}",
            partial(target_func, ecc, starter=starter, correct=correct),
            transform,
            domain,
            max_order=15,
            tol=5e-8,
        )
        fig.savefig(f"kepler{suffix}_{ecc:.2f}.png", bbox_inches="tight")
        plt.close(fig)
        errs[n] = err

    fig = plt.figure()
    plt.semilogy(eccs, errs)
    fig.savefig(f"kepler{suffix}_summary.png", bbox_inches="tight")
    plt.close(fig)
