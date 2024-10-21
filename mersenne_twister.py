import random


def mersenne_twister(n):
    L = []
    for _ in range(n):
        L.append(random.randint(0, 1))
    return L
