import pyprocar


###Standard bands###
pyprocar.bandsplot('PROCAR', outcar='OUTCAR', elimit=[-20,20], mode='plain', color='black', code='vasp', kpointsfile='KPOINTS', savefig='bands.png')

###Projection on atoms/orbitals###

###B and N###
pyprocar.bandsplot('PROCAR', outcar='OUTCAR', elimit=[-20,20], kpointsfile='KPOINTS', cmap='bone_r' ,mode='parametric', atoms=[0,1], savefig='b_pdos.png', vmin=0, vmax=1)

pyprocar.bandsplot('PROCAR', outcar='OUTCAR', elimit=[-20,20], kpointsfile='KPOINTS', cmap='bone_r' ,mode='parametric', atoms=[2,3], savefig='n_pdos.png', vmin=0, vmax=1)

pyprocar.bandsplot('PROCAR', outcar='OUTCAR', elimit=[-20,20], kpointsfile='KPOINTS', cmap='bone_r' ,mode='parametric', atoms=[0], orbitals=[0], savefig='b1_s_pdos.png', vmin=0, vmax=0.5)

pyprocar.bandsplot('PROCAR', outcar='OUTCAR', elimit=[-20,20], kpointsfile='KPOINTS', cmap='bone_r' ,mode='parametric', atoms=[0], orbitals=[1,2,3], savefig='b1_p_pdos.png', vmin=0, vmax=0.5)

pyprocar.bandsplot('PROCAR', outcar='OUTCAR', elimit=[-20,20], kpointsfile='KPOINTS', cmap='bone_r' ,mode='parametric', atoms=[2], orbitals=[0], savefig='n1_s_pdos.png', vmin=0, vmax=0.5)

pyprocar.bandsplot('PROCAR', outcar='OUTCAR', elimit=[-20,20], kpointsfile='KPOINTS', cmap='bone_r' ,mode='parametric', atoms=[2], orbitals=[1,2,3], savefig='n1_p_pdos.png', vmin=0, vmax=0.5)
