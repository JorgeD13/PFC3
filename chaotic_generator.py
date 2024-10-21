import math
import time


# 3.9
def mapa_logistico(r=3.9, x0=0.5, num_iteraciones=100):
    resultados = []
    x = x0

    # semillas
    seeds = [3.66, 3.73, 3.81]

    for i in range(num_iteraciones):
        x = r * x * (1 - x)

        if x < 0.5:
            resultados.append(0)
        else:
            resultados.append(1)

    return resultados


def variante_mapa_logistico(x0=0.5, num_iteraciones=100):
    resultados = []
    x = x0

    # semillas
    seeds = [3.879, 3.901, 3.914, 3.899, 3.921]

    for i in range(num_iteraciones):
        r = seeds[i % 5]

        x = r * x * (1 - x)

        if x < 0.5:
            resultados.append(0)
        else:
            resultados.append(1)

    return resultados