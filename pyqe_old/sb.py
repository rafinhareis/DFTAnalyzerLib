import os
os.system('python pyqe.py')
import matplotlib.pyplot as plt
import numpy as np

from pyqe import *


file = File()
plot = Plot()

#ML
names = {}
names['files_folder'] = 'Examples/sb2/ml/bands_k'
file.Set_files_atributes(names)
file.Load()
file.Bands_pdos()
file.orbitals_numbers()
plot.bands_project(file.bands_proj,states = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],E_min= -10, E_max= 10, dE = 2, cmap = 'bone_r')


file.Pdos_kresolved()
plot.kresolved(file.pdos_kresolved_plot,states =[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],E_min= -6, E_max= 6, dE = 2,cmap = 'bwr')

#ml2

names = {}
names['files_folder'] = 'Examples/sb2/ml2/bands_k'
file.Set_files_atributes(names)
file.Load()
file.Bands_pdos()
file.orbitals_numbers()
plot.bands_project(file.bands_proj,states = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],E_min= -10, E_max= 10, dE = 2, cmap = 'bone_r')


file.Pdos_kresolved()
plot.kresolved(file.pdos_kresolved_plot,states =[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],E_min= -6, E_max= 6, dE = 2,cmap = 'bwr')


#BL
names = {}
names['files_folder'] = 'Examples/sb2/bl/bands_k'
file.Set_files_atributes(names)
file.Load()
file.Bands_pdos()
file.orbitals_numbers()
plot.bands_project(file.bands_proj,states = range(1,33),E_min= -10, E_max= 10, dE = 2, cmap = 'bone_r')


file.Pdos_kresolved()
plot.kresolved(file.pdos_kresolved_plot,states =[3,4,11,12,19,20,27,28],E_min= -4, E_max= 4, dE = 2,cmap = 'bwr',thresh = True,tmax = 0.5)
