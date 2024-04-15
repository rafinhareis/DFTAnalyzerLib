import os
os.system('python pyqe.py')
import matplotlib.pyplot as plt
import numpy as np

from pyqe import *


file = File()
plot = Plot()

names = {}
names['files_folder'] = 'Examples/Titu'
names['bands_in_name'] = 'scf-bands.in'
names['bands_out_name'] = 'scf-bands.out'
names['bandsx_name'] = 'bands.out'

file.Set_files_atributes(names)
file.Load()
file.Bands_files()
plot.bands(file.bands, file.k_points_letter, file.k_path, E_min = -.5, E_max = .5)
