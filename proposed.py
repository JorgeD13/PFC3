import hashlib
import time

from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
from cryptography.hazmat.backends import default_backend
import math
import os

def variante_mapa_logistico(x0=0.5, num_iteraciones=100):
    x = x0

    # semillas
    seeds = [3.66, 3.73, 3.81, 3.94]

    # resultado
    res = 0

    for i in range(num_iteraciones):
        r = seeds[i % 4]

        x = r * x * (1 - x)

        if x >= 0.5:
            res += 1

    # retorna una "16-byte word"
    return res.to_bytes(16, 'big')

def encrypt(key, counter):
    # Usar counter como plaintext convirtiendolo a bytes
    c = counter.to_bytes(16, 'big')

    # Cifrado AES
    iv = os.urandom(16)
    cipher = Cipher(algorithms.AES(key), modes.CBC(iv), backend=default_backend())

    encryptor = cipher.encryptor()
    datos_cifrados = encryptor.update(c) + encryptor.finalize()

    return datos_cifrados


class ProposedAccumulator:
    def __init__(self):
        # Entropy pool:
        self.P = [b''] * 32
        self.reseed_counter = 0
        self.generator = ProposedGenerator()
        self.last_reseed = 0
        self.MINPOOLSIZE = 64
        for i in range(32):
            output = variante_mapa_logistico(num_iteraciones=i*10)
            self.P[i] = output

    def random_data(self, n):
        s = b''
        for i in range(32):
            if time.time()-self.last_reseed > 100:
                s += hashlib.sha256(self.P[i]).digest()
                self.P[i] = variante_mapa_logistico(num_iteraciones=i*100)
        self.last_reseed = self.generator.reseed(s)
        secuencia = self.generator.generate_random_data(n)
        secuencia = list(secuencia)
        secuencia = [int(x >= 128) for x in secuencia]
        return secuencia


class ProposedGenerator:
    def __init__(self):
        # G:
        self.key = int(time.time()).to_bytes(16, 'big')
        # self.key = b'0123456789012345'
        self.counter = 0

    def reseed(self, s):
        t = time.time()

        entrophy = os.urandom(32)

        self.key = hashlib.sha256(self.key + s + entrophy).digest()
        self.counter += 1
        return t

    def generate_blocks(self, k):
        r = b''
        for i in range(k):
            r += encrypt(self.key, self.counter)
            self.counter += 1
        return r

    def generate_random_data(self, n):
        assert 0 <= n <= 2 ** 20
        r = self.generate_blocks(math.ceil(n / 16))[:n]
        self.key = self.generate_blocks(2)
        return r




