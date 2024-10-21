import hashlib
import time

from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
from cryptography.hazmat.backends import default_backend
import math
import os

def encrypt(key, counter):
    # Usar counter como plaintext convirtiendolo a bytes
    c = counter.to_bytes(16, 'big')

    # Cifrado AES
    iv = os.urandom(16)
    cipher = Cipher(algorithms.AES(key), modes.CBC(iv), backend=default_backend())

    encryptor = cipher.encryptor()
    datos_cifrados = encryptor.update(c) + encryptor.finalize()

    return datos_cifrados


class Accumulator:
    def __init__(self):
        # Entropy pool:
        self.P = [b''] * 32
        self.reseed_counter = 0
        self.generator = Generator()
        self.last_reseed = 0
        self.MINPOOLSIZE = 64
        for i in range(32):
            entropy = os.urandom(32)
            self.P[i] = entropy

    def random_data(self, n):
        # print(time.time() - self.last_reseed)
        if len(self.P[0]) >= self.MINPOOLSIZE and time.time()-self.last_reseed > 100:
            self.reseed_counter += 1
            s = b''
            for i in range(32):
                if 2 ** i % self.reseed_counter == 0:
                    s += hashlib.sha256(self.P[i]).digest()
                    # self.P[i] = b''
                    new_entropy = os.urandom(32)
                    self.P[i] = new_entropy
            self.last_reseed = self.generator.reseed(s)

        secuencia = self.generator.generate_random_data(n)
        secuencia = list(secuencia)
        secuencia = [int(x >= 128) for x in secuencia]
        return secuencia


class Generator:
    def __init__(self):
        # G:
        self.key = int(time.time()).to_bytes(16, 'big')
        # self.key = b'0123456789012345'
        self.counter = 0

    def reseed(self, s):
        t = time.time()
        self.key = hashlib.sha256(self.key + s)
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




