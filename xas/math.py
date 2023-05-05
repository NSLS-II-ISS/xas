import numpy as np


def gauss(x, *p):
    A, mu, sigma = p
    return A * np.exp(-((x - mu) ** 2) / (2.0 * sigma**2))
