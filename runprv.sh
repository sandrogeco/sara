#!/bin/bash

# Percorso ambiente virtuale
VENV_PATH="venv1"

# Attiva l'ambiente virtuale
source venv1/bin/activate

# Lancia lo script Python con tutti i parametri ricevuti
python3 envelopeprv.py "$@"

# Disattiva l'ambiente virtuale (facoltativo)
deactivate
