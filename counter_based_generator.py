import time

from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
from cryptography.hazmat.backends import default_backend

class PRNG:
    def __init__(self):
        self.key = int(time.time())
        # self.key = b'0123456789012345'
        self.c = 0

    def generate_random(self, cantidad):
        # Rellena los bytes del contador para que tenga 16 bytes
        self.key = int(self.key * 1000000)
        if self.key > 9999999999999999:
            self.key = int(time.time() * 1000000)
        contador = self.c.to_bytes(16, 'big')

        cipher = Cipher(algorithms.AES(self.key.to_bytes(16, 'big')), modes.CTR(contador), backend=default_backend())

        encryptor = cipher.encryptor()
        secuencia = encryptor.update(b'\x00' * cantidad) + encryptor.finalize()
        self.c += 1

        secuencia = list(secuencia)
        return secuencia
