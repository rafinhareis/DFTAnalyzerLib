import os
os.system('python pyqe.py')
import matplotlib.pyplot as plt
import numpy as np

from pyqe import *


file = File()
plot = Plot()


names = {}
names['files_folder'] = 'Examples/GeS/ml/normal/pdos'
file.Set_files_atributes(names)
file.Load()

file.Dos_files('projwfc.x')
plot.dos(file.dos, E_max = 3, E_min=  -3)


file.Pdos_files()
plot.pdos_atoms(file.pdos_per_atoms,subplotsize=(8,6),dos_max= 10, E_max = 3, E_min=  -3)



names = {}
names['files_folder'] = 'Examples/GeS/ml/trocado/pdos'
file.Set_files_atributes(names)
file.Load()

file.Dos_files('projwfc.x')
plot.dos(file.dos, E_max = 3, E_min=  -3)


file.Pdos_files()
plot.pdos_atoms(file.pdos_per_atoms,subplotsize=(8,6),dos_max= 10, E_max = 3, E_min=  -3)