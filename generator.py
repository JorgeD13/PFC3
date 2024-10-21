import hashlib
import os
import time


class Fortuna:
    def __init__(self):
        self.pool = bytearray()  # Pool de acumulación de datos
        self.counter = 0  # Contador de bytes añadidos al pool

    def _get_entropy(self):
        # Generar entropía utilizando fuentes del sistema operativo
        entropy = os.urandom(64)
        self._add_entropy(entropy)

    def _add_entropy(self, data):
        # Añadir datos al pool de acumulación de datos
        self.pool += data
        self.counter += len(data)

    def _reseed(self):
        # Generar nueva semilla utilizando datos del pool
        seed = hashlib.sha256(self.pool).digest()
        self._add_entropy(seed)

    def random_data(self, length):
        # Generar una secuencia de datos aleatorios de la longitud especificada
        if self.counter > 64:
            self._reseed()
        elif time.time() % 3 == 0:
            self._get_entropy()

        data = bytearray()
        while len(data) < length:
            if self.counter > 64:
                self._reseed()

            self.counter += 1
            output = hashlib.sha256(self.pool).digest()
            self.pool = hashlib.sha256(self.pool + output).digest()
            data += output

        return list(data[:length])
