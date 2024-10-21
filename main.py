# import tests_run
#
# import os
# import time
# from chaotic_generator import variante_mapa_logistico, mapa_logistico
# from proposed import ProposedAccumulator

import requests

# URL de la API del ANU QRNG
url = "https://qrng.anu.edu.au/API/jsonI.php"

# Parámetros de la solicitud
params = {
    "length": 10,  # Número de elementos a generar
    "type": "uint8"  # Tipo de dato (0-255)
}

try:
    # Realizar la solicitud GET a la API
    response = requests.get(url, params=params)
    print(response)

    # Verificar si la solicitud fue exitosa
    if response.status_code == 200:
        data = response.json()["data"]
        print("Números aleatorios cuánticos:", data)
    else:
        print("Error en la solicitud:", response.status_code)
except Exception as e:
    print("Error durante la conexión:", e)
