import numpy as np
from scipy.stats import gaussian_kde
from sklearn.neighbors import NearestNeighbors
from collections import Counter


def mle_entropy_binary(data):
    # Calcular la frecuencia de 0s y 1s
    p_ones = np.mean(data)  # Probabilidad de 1
    p_zeros = 1 - p_ones  # Probabilidad de 0

    # Evitar log(0) usando un condicional para las probabilidades
    entropy = 0
    if p_zeros > 0:
        entropy -= p_zeros * np.log2(p_zeros)
    if p_ones > 0:
        entropy -= p_ones * np.log2(p_ones)

    return entropy


def miller_madow_entropy_binary(data):
    # Calcular la frecuencia de 0s y 1s
    p_ones = np.mean(data)  # Probabilidad de 1
    p_zeros = 1 - p_ones  # Probabilidad de 0

    # Evitar log(0) usando un condicional para las probabilidades
    entropy_mle = 0
    if p_zeros > 0:
        entropy_mle -= p_zeros * np.log2(p_zeros)
    if p_ones > 0:
        entropy_mle -= p_ones * np.log2(p_ones)

    # Aplicar corrección de Miller-Madow
    N = len(data)  # Número total de elementos
    K = 2  # Número de categorías (0 y 1)
    correction = (K - 1) / (2 * N)

    entropy_mm = entropy_mle + correction
    return entropy_mm


def chao_shen_entropy_binary(data):
    # Calcular la frecuencia de 0s y 1s
    p_ones = np.mean(data)  # Probabilidad de 1
    p_zeros = 1 - p_ones  # Probabilidad de 0

    # Calcular la entropía MLE
    entropy_mle = 0
    if p_zeros > 0:
        entropy_mle -= p_zeros * np.log2(p_zeros)
    if p_ones > 0:
        entropy_mle -= p_ones * np.log2(p_ones)

    # Aplicar corrección de Chao-Shen
    N = len(data)  # Número total de elementos
    K = 2  # Número de categorías (0 y 1)

    correction = (K - 1) / N * np.log(N / (K - 1))

    entropy_cs = entropy_mle + correction
    return entropy_cs


def bayesian_entropy_binary(data):
    # Calcular la frecuencia de 0s y 1s
    p_ones = np.mean(data)  # Probabilidad de 1
    p_zeros = 1 - p_ones  # Probabilidad de 0

    # Calcular la entropía MLE
    entropy_mle = 0
    if p_zeros > 0:
        entropy_mle -= p_zeros * np.log2(p_zeros)
    if p_ones > 0:
        entropy_mle -= p_ones * np.log2(p_ones)

    # Aplicar corrección Bayesiana
    N = len(data)  # Número total de elementos
    K = 2  # Número de categorías (0 y 1)

    correction = (1 / N) * (np.log2(N + 1) - np.log2(K + 1))

    entropy_bayes = entropy_mle + correction
    return entropy_bayes


def james_stein_entropy_binary(data):
    # Calcular la frecuencia de 0s y 1s
    p_ones = np.mean(data)  # Probabilidad de 1
    p_zeros = 1 - p_ones  # Probabilidad de 0

    # Calcular la entropía MLE
    entropy_mle = 0
    if p_zeros > 0:
        entropy_mle -= p_zeros * np.log2(p_zeros)
    if p_ones > 0:
        entropy_mle -= p_ones * np.log2(p_ones)

    # Aplicar corrección James-Stein
    N = len(data)  # Número total de elementos
    K = 2  # Número de categorías (0 y 1)

    correction = (K - 1) / N * ((N - K) / (N + 1))

    entropy_js = entropy_mle + correction
    return entropy_js


def jackknife_entropy_binary(data):
    N = len(data)

    # Calcular la entropía MLE para el conjunto completo
    p_ones = np.mean(data)
    p_zeros = 1 - p_ones

    # Calcular la entropía MLE
    entropy_mle = 0
    if p_zeros > 0:
        entropy_mle -= p_zeros * np.log2(p_zeros)
    if p_ones > 0:
        entropy_mle -= p_ones * np.log2(p_ones)

    # Inicializar la suma para la entropía jackknife
    entropy_jk = 0

    # Calcular la entropía MLE excluyendo cada elemento
    for i in range(N):
        # Crear un nuevo conjunto de datos excluyendo el i-ésimo elemento
        data_excluded = np.delete(data, i)

        p_ones_excluded = np.mean(data_excluded)
        p_zeros_excluded = 1 - p_ones_excluded

        # Calcular la entropía MLE excluyendo el i-ésimo elemento
        entropy_mle_excluded = 0
        if p_zeros_excluded > 0:
            entropy_mle_excluded -= p_zeros_excluded * np.log2(p_zeros_excluded)
        if p_ones_excluded > 0:
            entropy_mle_excluded -= p_ones_excluded * np.log2(p_ones_excluded)

        # Sumar la entropía excluida
        entropy_jk += entropy_mle_excluded

    # Calcular el estimador Jackknife
    entropy_jk = N * entropy_mle - (N - 1) * (entropy_jk / N)

    return entropy_jk


def best_upper_bound_entropy_binary(data):
    # Calcular la frecuencia de 0s y 1s
    p_ones = np.mean(data)  # Probabilidad de 1
    p_zeros = 1 - p_ones  # Probabilidad de 0

    # Evitar problemas con log(0) y calcular el valor de H BUB
    H_BUB = np.log2(2)  # Para secuencias binarias, K = 2
    print("hbub: ", H_BUB)
    print("pzeros: ", p_zeros)
    print("pones: ", p_ones)

    # Calcular la entropía usando la probabilidad de 0s y 1s
    if p_zeros > 0:
        H_BUB -= p_zeros * np.log2(p_zeros)
    if p_ones > 0:
        H_BUB -= p_ones * np.log2(p_ones)

    return H_BUB


def kernel_density_entropy_binary(data):
    # Asegurarse de que los datos son numpy array
    data = np.asarray(data)

    # Estimar la densidad usando KDE
    kde = gaussian_kde(data, bw_method='scott')  # Puedes ajustar bw_method según sea necesario
    x = np.linspace(0, 1, 1000)  # Crear puntos en el rango [0, 1]
    pdf = kde(x)

    # Calcular la entropía usando la estimación de densidad
    # Evitar log(0) añadiendo un pequeño valor
    pdf = np.maximum(pdf, 1e-10)  # Para evitar log(0)

    entropy = -np.sum(pdf * np.log(pdf)) * (x[1] - x[0])  # Integración numérica aproximada

    return entropy


def knn_entropy_binary(data, k=1):
    # Asegurarse de que los datos son numpy array y reshaped para K-NN
    data = np.asarray(data).reshape(-1, 1)  # Convertir a formato 2D para K-NN

    # Calcular las distancias usando K-NN
    nbrs = NearestNeighbors(n_neighbors=k + 1).fit(data)  # +1 para incluir el punto mismo
    distances, _ = nbrs.kneighbors(data)

    # Tomar las distancias a los k vecinos más cercanos (sin incluir el punto mismo)
    distances = distances[:, 1:]  # Excluir la primera columna que es la distancia al mismo punto
    avg_distances = np.mean(distances, axis=1)  # Calcular la distancia promedio

    # Calcular la entropía
    entropy = np.log(k) + np.mean(np.log(avg_distances + 1e-10))  # Añadir pequeño valor para evitar log(0)

    return entropy


def zhang_entropy_binary(data):
    # Asegurarse de que los datos son numpy array
    data = np.asarray(data)

    # Contar las frecuencias de los elementos
    freq = Counter(data)
    n = len(data)  # Tamaño total de la muestra
    m = len(freq)  # Número de eventos únicos

    # Calcular las probabilidades
    probabilities = np.array([count / n for count in freq.values()])

    # Calcular la entropía
    entropy = -np.sum(probabilities * np.log(probabilities + 1e-10))  # Añadir pequeño valor para evitar log(0)
    correction_term = (m - 1) / (2 * n) * (np.sum(1 / probabilities) - 1)

    return entropy + correction_term


def nsb_entropy_binary(data):
    # Asegurarse de que los datos son numpy array
    data = np.asarray(data)

    # Contar las frecuencias de los elementos
    freq = Counter(data)
    n = len(data)  # Tamaño total de la muestra
    m = len(freq)  # Número de eventos únicos

    # Calcular las probabilidades
    probabilities = np.array([count / n for count in freq.values()])
    frequencies = np.array([count for count in freq.values()])

    # Calcular la entropía
    entropy = -np.sum(probabilities * np.log(probabilities + 1e-10))  # Añadir pequeño valor para evitar log(0)

    # Término de corrección NSB
    correction_term = (1 / n) * np.sum(np.log(n / frequencies))

    return entropy + correction_term

