import os
os.system('python pyqe.py')
import matplotlib.pyplot as plt
import numpy as np

from pyqe import *


file = File()
plot = Plot()


#zn
names = {}
names['files_folder'] = 'Examples/Zn/band_proj'
file.Set_files_atributes(names)
file.Load()
file.Bands_pdos()
file.orbitals_numbers()
plot.bands_project(file.bands_proj,states = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],E_min= -10, E_max= 10, dE = 2, cmap = 'bone_r')

plot.bands_project(file.bands_proj,states = [2,3,4,11,12,13],E_min= -20, E_max= 20, dE = 5, cmap = 'bone_r')
#4s
plot.bands_project(file.bands_proj,states = [1,10],E_min= -20, E_max= 20, dE = 4, cmap = 'bone_r')
#4p
plot.bands_project(file.bands_proj,states = [2,3,4,11,12,13],E_min= -10, E_max= 10, dE = 2, cmap = 'bone_r')
#3d
plot.bands_project(file.bands_proj,states = [5,6,7,8,9,14,15,16,17,18],E_min= -10, E_max= 10, dE = 2, cmap = 'bone_r')


plot.bands_project(file.bands_proj,states = [1,10],E_min= -20, E_max= 20, dE = 4, cmap = 'viridis', linewidth=8)

p4 = [2,3,4,11,12,13]
p3_s4 =  [1,10]+[5,6,7,8,9,14,15,16,17,18]

plot.bands_project(file.bands_proj, corretalion = True, states1 = p4, 
                   states2 = p3_s4,E_min= -10, E_max= 10, dE = 2, cbar_label=['4d','3s4p'],cbar_tick= [1,-1]
                   , linewidth=8)


#cu
names = {}
names['files_folder'] = 'Examples/Cu/band_proj'
file.Set_files_atributes(names)
file.Load()
file.Bands_pdos()
file.orbitals_numbers()
plot.bands_project(file.bands_proj,states = [1,2,3,4,5,6,7,8,9],E_min= -10, E_max= 10, dE = 2, cmap = 'bone_r')

#4s
plot.bands_project(file.bands_proj,states = [1],E_min= -5, E_max= 5, dE = 1, cmap = 'bone_r')

#4p
plot.bands_project(file.bands_proj,states = [2,3,4],E_min= -5, E_max= 5, dE = 1, cmap = 'bone_r')

#3d
plot.bands_project(file.bands_proj,states = [5,6,7,8,9],E_min= -5, E_max= 5, dE = 1, cmap = 'bone_r')



p4 = [2,3,4]
p3_s4 =  [1]+[2,3,4]

plot.bands_project(file.bands_proj, corretalion = True, states1 = p4, 
                   states2 = p3_s4,E_min= -10, E_max= 10, dE = 2, cbar_label=['4d','3s4p'],cbar_tick= [1,-1]
                   , linewidth=8)
