from pandas import DataFrame, read_csv
from numpy import array,arange, sqrt, concatenate
from scipy import interpolate
from scipy.signal import savgol_filter as sv
import os
import subprocess
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import warnings
from .base import *
from .siesta_base import *
warnings.filterwarnings('ignore')


def Load(names = {'A':0,'B':0}, code = 'qe', prefix = ''):
    if code == 'qe':
        def Set_files_atributes(dict2):
            load_dict =  {
            'files_folder' : './',
            'scf_in_name' : 'scf.in',
            'scf_out_name' : 'scf.out',
            'nscf_in_name' : 'nscf.in',
            'nscf_out_name': 'nscf.out',
            # must have ! in front of k letter, like 0.3333333 0.3333333 0.00000000 !K , for gamma point use G
            'bands_in_name' : 'bands.in',
            'bands_out_name' : 'bands.out',
            'bandsx_name' : 'bands.x.out',
            'bands_dat_name' : 'bands.dat.gnu',
            'dos_dat_name' : 'dos',
            'pdos_dat_name' : 'pdos.dat.pdos_tot',
            'pdos_prefix' : '',
            'fermi_energy' : 0,
            'projwfc_out_name' : 'projwfc.out',
            'atomic_proj' : 'atomic_proj.xml',
            'ppin_name' : 'pp.in',
            'charge2D_dat_name' : 'charge2D.out'
            }
            keys = dict2.keys()
            for key in keys:
                if key == 'files_folder':
                    load_dict['files_folder'] = dict2[key]
                elif key == 'scf_in_name':
                    load_dict['scf_in_name'] = dict2[key]
                elif key == 'scf_out_name':
                    load_dict['scf_out_name'] = dict2[key]
                elif key == 'nscf_in_name':
                    load_dict['nscf_in_name'] = dict2[key]
                elif key == 'nscf_out_name':
                    load_dict['nscf_out_name'] = dict2[key]
                elif key == 'bands_in_name':
                    load_dict['bands_in_name'] = dict2[key]
                elif key == 'bands_out_name':
                    load_dict['bands_out_name'] = dict2[key]
                elif key == 'bandsx_name':
                    load_dict['bandsx_name'] = dict2[key]
                elif key == 'bands_dat_name':
                    load_dict['bands_dat_name'] = dict2[key]
                elif key == 'dos_dat_name':
                    load_dict['dos_dat_name'] = dict2[key]
                elif key == 'pdos_dat_name':
                    load_dict['pdos_dat_name'] = dict2[key]
                elif key == 'pdos_prefix':
                    load_dict['pdos_prefix'] = dict2[key]
                elif key == 'projwfc_out_name':
                    load_dict['projwfc_out_name'] = dict2[key]
                elif key == 'ppin_name':
                    load_dict['ppin_name'] = dict2[key]
                elif key == 'charge2D_dat_name':
                    load_dict['charge2D_dat_name'] = dict2[key]

            return load_dict
        
        load_dict = Set_files_atributes(names)

        scfin_file_path = os.path.join(load_dict['files_folder'], load_dict['scf_in_name'])
        scfin_file =  open_file(scfin_file_path)

        scfout_file_path = os.path.join(load_dict['files_folder'],load_dict['scf_out_name'])
        scfout_file = open_file(scfout_file_path)

        version_QE = version_quantum(scfout_file)

        try:
            load_dict['fermi_energy'] = E_fermi(scfout_file,load_dict['scf_out_name'],version_QE)
        except UnboundLocalError:
            try:
                load_dict['fermi_energy'] = E_fermi(scfout_file,load_dict['nscf_out_name'],version_QE)
            except UnboundLocalError:
                    load_dict['fermi_energy'] = E_fermi(scfout_file,load_dict['bands_out_name'] ,version_QE)

        cell_volumn = volumn(scfout_file)
        alatice = alat(scfout_file)

        return [scfin_file,scfout_file,load_dict,[cell_volumn,alatice]]

    elif code == 'siesta':
        def Set_files_atributes2(dict2, prefix =prefix):
            load_dict =  {
            'files_folder' : './',
            'fdf_name' : prefix+str('.fdf'),
            'out_name' : prefix+str('.out'),
            'bands_name' : prefix+str('.bands'),
            'dos_name': prefix+str('.DOS'),
            }
            keys = dict2.keys()
            for key in keys:
                if key == 'files_folder':
                    load_dict['files_folder'] = dict2[key]
                elif key == 'fdf_name':
                    load_dict['fdf_name'] = dict2[key]
                elif key == 'out_name':
                    load_dict['out_name'] = dict2[key]
                elif key == 'bands_name':
                    load_dict['bands_name'] = dict2[key]
                elif key == 'dos_name':
                    load_dict['dos_name'] = dict2[key]
            return load_dict
        
        load_dict = Set_files_atributes2(names,prefix=prefix)
        fdf_file_path = os.path.join(load_dict['files_folder'], load_dict['fdf_name'])
        fdf_file =  open_file(fdf_file_path)

        out_file_path = os.path.join(load_dict['files_folder'],load_dict['out_name'])
        out_file = open_file(out_file_path)



        return [load_dict,fdf_file,out_file]



def Bands_files(load_dict,fermi_energy_corr=True, bands_project = False, from_bandsout = False, code = 'qe'):
        
    if code == 'qe':
        bandout_file_path= os.path.join(load_dict['files_folder'],load_dict['bands_out_name'])
        bandout_file = open_file(bandout_file_path)

        bandin_file_path = os.path.join( load_dict['files_folder'],load_dict['bands_in_name'])

        bandin_file = open_file(bandin_file_path)

        nspin = n_spin(bandin_file)
        k_points_letter, k_points_bands, nk = k_point(bandin_file)

        verbosity_value = verbosity(bandin_file)
                    
        version_QE = version_quantum(bandout_file)
        fermi_energy = load_dict['fermi_energy']
                
        if version_QE == '7.1':
                    fermi_energy = E_fermi(bandout_file,'bands.out',version_QE)
        elif version_QE == '7.2' and verbosity_value == 'low':
                    fermi_energy = E_fermi(bandout_file,'bands.out',version_QE)
        elif version_QE == '6.7MaX':
                    fermi_energy = load_dict['fermi_energy']
        if bands_project == False:
            if from_bandsout == False:
                if nspin == 1:
                    band_file_path = os.path.join(load_dict['files_folder'],load_dict['bands_dat_name'] )
                    bands_files = bandas_df(band_file_path)
                    bands, E_valence, E_conduction, E_valence_point,E_conduction_point = orbitais( bands_files, fermi_energy_corr, fermi_energy)
                    gap = E_conduction - E_valence

                #implementar depis
        #        elif self.nspin == 2:
        #            self.bandup_file_path = os.path.join(self.file_dir, self.bandup_file_name)
        #            bands_up = bandas_df(self.bandup_file_path)
        #            self.bands_up, self.E_valence_up, self.E_conduction_up, self.E_valence_point_up,self.E_conduction_point_up = orbitais(
        #                bands_up, fermi_energy_corr, self.fermi_energy)#
        #
        #            self.banddown_file_path = os.path.join(self.file_dir, self.banddown_file_name)
        #            bands_down = bandas_df(self.banddown_file_path)
        #            self.bands_down, self.E_valence_down, self.E_conduction_down, self.E_valence_point_down,self.E_conduction_point_down = orbitais(
        #                bands_down, fermi_energy_corr, self.fermi_energy)#
        #
        #            columns = self.bands_up.columns
        #            self.bands = DataFrame(self.bands_up[columns[0]])
        #            for i in range(1,len(columns)):
        #                self.bands[columns[i]] = self.bands_up[columns[i]] + self.bands_down[columns[i]]#
        #
        #            if self.E_conduction_down <= self.E_conduction_up:
        #                self.E_conduction = self.E_conduction_down
        #                self.E_conduction_point = self.E_conduction_down
        #           else:
        #                self.E_conduction = self.E_conduction_up
        #                self.E_conduction_point = self.E_conduction_up
        #            
        #            if self.E_valence_down >= self.E_valence_up:
        #                self.E_valence= self.E_valence_down
        #                self.E_valence_point = self.E_valence_down
        #            else:
        #                self.E_valence = self.E_valence_up
        #                self.E_valence_point = self.E_valence_up
        #
        #
        #            self.gap = self.E_conduction - self.E_valence
                    
                else:
                    print('Worng value for nspin atribute.')

                bandxout_file_path = os.path.join(load_dict['files_folder'], load_dict['bandsx_name'])
                bandxout_file = open_file(bandxout_file_path)
                k_path, k_points_bandsx = k_points_path(bandxout_file, len(k_points_letter))

                return [bands,k_points_letter,k_path, E_valence, E_conduction, E_valence_point,E_conduction_point]

            if from_bandsout == True:
                if verbosity_value== 'high':
                    band_df= band_structure_from_bandsout(bandout_file, fermi_energy_corr, fermi_energy , version_QE)

                    points_break = [0]
                    ct = 0
                    for n in nk:
                        ct +=n
                        points_break.append(ct)
                    k_axis = []
                    k = band_df[band_df.columns[0]]
                    for i in range(len(band_df[band_df.columns[0]])):
                        if (i in points_break) == True:
                            k_axis.append(k[i] )
                    k_points_PAT = array(k_axis)

                    return [band_df,k_points_letter, k_points_PAT]
        elif bands_project ==True:
            
            if verbosity_value== 'high':
                atomic_xml_file_path = os.path.join(load_dict['files_folder'],load_dict['atomic_proj'])
                atomic_xml_file = open_file(atomic_xml_file_path)
                proj,n_band,nspin,n_wfc = atomic_proj(atomic_xml_file)
                

                projwfc_out_file_path = os.path.join(load_dict['files_folder'],  load_dict['projwfc_out_name'])      
                projwfc_out_file = open_file(projwfc_out_file_path)
                atom_number , atom_orb_number = States(projwfc_out_file)



                band_df= band_structure_from_bandsout(bandout_file, fermi_energy_corr, fermi_energy , version_QE)

                points_break = [0]
                ct = 0
                for n in nk:
                    ct +=n
                    points_break.append(ct)
                k_axis = []
                k = band_df[band_df.columns[0]]
                for i in range(len(band_df[band_df.columns[0]])):
                    if (i in points_break) == True:
                        k_axis.append(k[i] )
                k_points_PAT = array(k_axis)
                print(atom_number)
                print(atom_orb_number)
                return  [band_df,proj,k_points_letter,k_points_PAT,atom_number,atom_orb_number]
    elif code == 'siesta':
        bands_file_path = os.path.join(load_dict['files_folder'],load_dict['bands_name'])
        bands_file = open_file(bands_file_path)
        band_df,k_points_letter, k_points_PAT = bands_from_siesta(bands_file)
        return [band_df,k_points_letter, k_points_PAT]
    
def Bands_plot(data, subplotsize=(10, 8), subplot=False, ax=None, vline=True, vline_linewidth=1, E_min=-5, E_max=5, dE=1, font_axis=20, font_label=20, ylabel='$E -E_{F} (eV)$', xlabel='', title='Band Structure',
              linewidth=3, legend=False, loc="upper right", color='black', color_vline='gray', vline_linestyle='dashed', label_band='eletronic strucuture', fermi_level_line=False, fermi_level=0, fermi_level_color='red', 
              fermi_level_linewidth=1, fermi_level_linestyle='dashed', occupied_and_unoccupied_states_dot = False, occupied_and_unoccupied_states_points = ((0,0),(0,0)), scatter_type = 'o',occupied_and_unoccupied_states_color_dot = 'blue',
              occupied_and_unoccupied_states_line=False, occupied_and_unoccupied_states=(-1, 1), occupied_and_unoccupied_states_color=('blue', 'blue'), occupied_and_unoccupied_states_linewidth=1,occupied_and_unoccupied_states_size_dot=80
              , dpi = 150, fig = None, shift_E = 0,fill = False,fill_x0 = [0],fill_x1= [1],colorfill = ['blue'], alpha = [1],exclude_bands = [], highligth_bands = [],highligth_bands_color = []
              ,combine = False,combine_bands = [1,2], combine_point = 0, combine_delta = [-0.1,0.1],combine_color = ['red','blue'],combine_exclude = -1, grad_color = False, color1= "#8A5AC2", color2="#3575D5",grad_n = 1
              , smooth = False, smooth_par = [5,1],dy = 0, nbdn_dy = [0], linestyle = 'solid', which_kpoints = []  ) :
        # flag subplot
        if subplot == False:
            del (ax)
            fig, ax = plt.subplots(figsize=subplotsize)#, dpi = dpi)
        # -----------------------------------------------------------------------------------------------------
        bands, k_points_letter2, k_path2 = data
        if len(which_kpoints)>0:
            k_path = [];k_points_letter = []
            ct = 0
            for i in range(len(which_kpoints)):
                for j in range(ct,len(k_points_letter2)):
                    if which_kpoints[i] == 'G' or which_kpoints[i] == '$\\Gamma $' or which_kpoints[i] == 'Gamma'  :
                        if '$\\Gamma $' == k_points_letter2[j]:
                            k_points_letter.append(k_points_letter2[j])
                            k_path.append(k_path2[j])
                            ct =j
                            break
                    else:
                        if which_kpoints[i] == k_points_letter2[j]:
                            k_points_letter.append(k_points_letter2[j])
                            k_path.append(k_path2[j])
                            ct =j
                            break
        else:
            k_points_letter = k_points_letter2.copy()
            k_path = k_path2.copy()


        # plot data
        columns = bands.columns
        k = bands[columns[0]]
        for i in range(1, len(columns)):
            if (i not in exclude_bands):
                if i not in combine_bands:
                    if i in nbdn_dy:
                        energy = bands[columns[i]]+shift_E+dy
                    else:
                        energy = bands[columns[i]]+shift_E
                    if smooth:
                        energy = sv(energy,smooth_par[0],smooth_par[1])

                    if legend == True:
                        if i == 1:
                            ax.plot(k, energy, color=color,
                                    linewidth=linewidth, label=label_band,linestyle = linestyle)
                        else:
                            ax.plot(k, energy, color=color, linewidth=linewidth,linestyle = linestyle)
                    else:
                        ax.plot(k, energy, color=color, linewidth=linewidth,linestyle = linestyle)
        # -----------------------------------------------------------------------------------------------------
        #combine bands
        if combine:
            conca = []
            conca2 = []
            if combine_bands[0] in nbdn_dy:
                energy1 = bands[columns[combine_bands[0]]]+shift_E+dy
            else:
                energy1 = bands[columns[combine_bands[0]]]+shift_E
            if combine_bands[1] in nbdn_dy:   
                energy2 = bands[columns[combine_bands[1]]]+shift_E+dy
            else:
                energy2 = bands[columns[combine_bands[1]]]+shift_E
            for i in range(len(k)):
                if k[i] < combine_point-combine_delta[0] or k[i]> combine_point+combine_delta[1]:
                    conca.append(energy1[i])
                    conca2.append(energy2[i])
                else:
                    conca.append(energy2[i])
                    conca2.append(energy1[i])
            
            conca = array(conca)
            conca2 = array(conca2)
            if smooth:
                conca = sv(conca,smooth_par[0],smooth_par[1])
                conca2 = sv(conca2,smooth_par[0],smooth_par[1])
            if combine_exclude ==-1:
                ax.plot(k, conca, color=combine_color[0], linewidth=linewidth,linestyle = linestyle)
                ax.plot(k, conca2, color=combine_color[1], linewidth=linewidth,linestyle = linestyle)
            elif combine_exclude == 0:
                ax.plot(k, conca2, color=combine_color[1], linewidth=linewidth,linestyle = linestyle)
            elif combine_exclude == 1:
                ax.plot(k, conca, color=combine_color[0], linewidth=linewidth,linestyle = linestyle)

        # -----------------------------------------------------------------------------------------------------
        #highlidted bands
        for i in range(len(highligth_bands)):
            if highligth_bands[i] in nbdn_dy:
                energy = bands[columns[highligth_bands[i]]]+shift_E+dy
            else:
                energy = bands[columns[highligth_bands[i]]]+shift_E
            if smooth:
                energy = sv(energy,smooth_par[0],smooth_par[1])
            ax.plot(k,energy,color = highligth_bands_color[i],linewidth=linewidth,linestyle = linestyle)


        # -----------------------------------------------------------------------------------------------------
        #fill_between

        if fill:
            for i in range(len(fill_x1)):
                if i in nbdn_dy:
                    energy1 = bands[columns[fill_x0[i]]]+shift_E+dy
                else:
                    energy1 = bands[columns[fill_x0[i]]]+shift_E
                if i in nbdn_dy:     
                    energy2 = bands[columns[fill_x1[i]]]+shift_E+dy
                else:
                    energy2 = bands[columns[fill_x1[i]]]+shift_E
                if smooth:
                    energy1 = sv(energy1,smooth_par[0],smooth_par[1])
                    energy2 = sv(energy2,smooth_par[0],smooth_par[1])
                if grad_color:
                    color_grad = get_color_gradient(color1,color2,grad_n)
                    grad_n_temp = int(grad_n/(len(fill_x1)))
                    delta = (energy2 - energy1)/grad_n_temp
                    for j in range(grad_n_temp):
                        ax.fill_between(k,energy1+j*delta,energy1+(j+1)*delta, color = color_grad[j+i*grad_n_temp],alpha = alpha[i])    
                else:
                    ax.fill_between(k,energy1,energy2, color = colorfill[i],alpha = alpha[i])    


        # -----------------------------------------------------------------------------------------------------
        # flag vline
        if vline == True:
            for i in range(1, len(k_path)-1):
                ax.vlines(k_path[i], E_min, E_max, linestyles=vline_linestyle,
                          colors=color_vline, linewidth=vline_linewidth)
        # -----------------------------------------------------------------------------------------------------
        # flag fermy_level_line
        if fermi_level_line == True:
            ax.hlines(fermi_level, k_path[0], k_path[len(
                k_path)-1], linestyles=fermi_level_linestyle, color=fermi_level_color, linewidth=fermi_level_linewidth)
        # -----------------------------------------------------------------------------------------------------
        # flag occupied_and_unoccupied_states_line
        if occupied_and_unoccupied_states_line == True and occupied_and_unoccupied_states_dot == False:
            for i in range(len(occupied_and_unoccupied_states)):
                item = occupied_and_unoccupied_states[i]
                ax.hlines(item, k_path[0], k_path[len(k_path)-1], linestyles=fermi_level_linestyle,
                          linewidth=occupied_and_unoccupied_states_linewidth, color=occupied_and_unoccupied_states_color[i])
        elif occupied_and_unoccupied_states_line == False and occupied_and_unoccupied_states_dot == True:
            occupied_and_unoccupied_states_points = []
            columns = bands.columns
            k = bands[columns[0]]
            for i in range(1, len(columns)):
                if (i not in exclude_bands):
                    if i not in combine_bands:
                        break
                        if i in nbdn_dy:
                            energy = bands[columns[i]]+shift_E+dy
                        else:
                            energy = bands[columns[i]]+shift_E
                        if smooth:
                            energy = sv(energy,smooth_par[0],smooth_par[1])
                        



            x = []; y = []
            for item in occupied_and_unoccupied_states_points:
                x.append(item[0]); y.append(item[1])
            ax.scatter(x,y, marker = scatter_type,s= occupied_and_unoccupied_states_size_dot, color = occupied_and_unoccupied_states_color_dot)
        # -----------------------------------------------------------------------------------------------------

        # flag legend
        if legend == True:
            ax.legend(loc=loc, fontsize=font_label)
        # -----------------------------------------------------------------------------------------------------
        # globar parameters
        ax.set_xlim(k_path[0], k_path[len(k_path)-1])
        ax.set_xticks(k_path, k_points_letter, fontsize=font_axis)
        ax.set_xlabel(xlabel)

        ax.set_ylim(E_min, E_max)
        ytick = array(list(map( lambda x: round(x,2),arange(E_min, E_max+dE, dE))))
        ax.set_yticks(ytick,
                      ytick, fontsize=font_axis)
        ax.set_ylabel(ylabel)
        ax.xaxis.label.set_size(font_label)
        ax.yaxis.label.set_size(font_label)

        ax.set_title(title, fontsize=font_label)

        return ax

def Bands_project(bands_proj, subplotsize=(10, 8), subplot = False, ax = None,dpi = 150,
    legend = False, font_label = 20, font_axis = 20, xlabel = '', ylabel = 'E - E_f (eV)', E_min = -5, E_max = 5, dE = 1,loc = 'upper right',title = 'Band projected',
    vline = True, color_vline='gray', vline_linestyle='dashed', vline_linewidth = 1.5, backgorund = 'white', cmap = 'viridis', norm = (0,1)
    , index = False, states = [1,2], corretalion = False, states1 = [1,2,3], states2 = [4,5,6],norm_corr = (-1,1), cbar_label = ['0','1'], cbar_tick = [0,1]
    ,linewidth = 3, facecolor = 'white',fig = None):


        if subplot == False:
            del (ax)
            fig, ax = plt.subplots(figsize=subplotsize)#, dpi = dpi)
        bands = bands_proj[0]; k_letter = bands_proj[2]; k_path = bands_proj[3]
        columns = bands.columns
        x_interp = bands[columns[0]]
        x = arange(x_interp.min(),x_interp.max(),0.001)
        proj = bands_proj[1]
        for i in range(1,len(columns)):
            y_int = bands[columns[i]]
            if corretalion == False:
                marker = False
                for state in states:
                    if marker == False:
                        z_int = array(proj[state][i])
                        marker = True
                    else:
                        z_int= z_int  + array(proj[state][i])
                        
            elif corretalion == True:
                marker = False
                for state1 in states1:
                    if marker == False:
                        z1 = array(proj[state1][i])
                        marker = True
                    else:
                        z1= z1  + array(proj[state1][i])
                marker = False
                for state2 in states2:
                    if marker == False:
                        z2 = array(proj[state1][i])
                        marker = True
                    else:
                        z2= z2  + array(proj[state1][i])
                z_int = z2 - z1
            else:
                print('Wrong value for correlation, only True and False are accepted. Default is False')
            
                
            f = interpolate.interp1d(x_interp, y_int)
            g = interpolate.interp1d(x_interp, z_int)
            y = f(x)
            z = g(x)
            points = array([x, y]).T.reshape(-1, 1, 2)
            segments = concatenate([points[:-1], points[1:]], axis=1)
            if corretalion == False:
                norma = plt.Normalize(norm[0],norm[1])
            elif corretalion == True:
                norma = plt.Normalize(norm_corr[0],norm_corr[1])
            #norm = plt.Normalize(0, z.max())
            #norm = plt.Normalize(mini,maxi)
            lc = LineCollection(segments, cmap=cmap, norm=norma, linewidth = linewidth) #cmap='viridis'
            lc.set_array(z)
            lc.set_linewidth( linewidth)
            line = ax.add_collection(lc)

        ax.set_facecolor(backgorund)
        cbar = fig.colorbar(line, ticks = cbar_tick)
        cbar.ax.set_yticklabels(cbar_label,fontsize = 20)


        # flag vline
        if vline == True:
            for i in range(1, len(k_path)-1):
                ax.vlines(k_path[i], E_min, E_max, linestyles=vline_linestyle,
                          colors=color_vline, linewidth=vline_linewidth)
        # -----------------------------------------------------------------------------------------------------

                # flag legend
        if legend == True:
            ax.legend(loc=loc, fontsize=font_label)
        # -----------------------------------------------------------------------------------------------------
        # globar parameters
        ax.set_xlim(k_path[0], k_path[len(k_path)-1])
        ax.set_xticks(k_path, k_letter, fontsize=font_axis)
        ax.set_xlabel(xlabel)

        ax.set_ylim(E_min, E_max)
        ax.set_yticks(arange(E_min, E_max+dE, dE),
                      arange(E_min, E_max+dE, dE), fontsize=font_axis)
        ax.set_ylabel(ylabel)
        ax.xaxis.label.set_size(font_label)
        ax.yaxis.label.set_size(font_label)
        ax.set_facecolor(facecolor)
        ax.set_title(title, fontsize=font_label)

        return ax 

def Load_bands_files(names = {'A':0,'B':0}, type_bands = 'default', code = 'qe',fermi_energy_corr=True, from_bandsout = False,prefix = ''):
    if code == 'qe':
        if type_bands == 'default':
            scfin_file,scfout_file,load_dict , parameters = Load(names = names)
            bands_project = False
            if from_bandsout==False:
                bands,k_points_letter,k_path, E_valence, E_conduction, E_valence_point,E_conduction_point = Bands_files(
                    load_dict= load_dict,fermi_energy_corr=fermi_energy_corr, bands_project = bands_project, from_bandsout = from_bandsout)
                data = [bands,k_points_letter,k_path]
                return data #E_valence, E_conduction, E_valence_point,E_conduction_point]
            elif from_bandsout == True:
                bands,k_points_letter,k_path = Bands_files(
                    load_dict= load_dict,fermi_energy_corr=fermi_energy_corr, bands_project = bands_project, from_bandsout = from_bandsout)
                data = [bands,k_points_letter,k_path]
                return data #E_valence, E_conduction, E_valence_point,E_conduction_point]
        elif type_bands == 'projected':
            scfin_file,scfout_file,load_dict , parameters = Load(names = names)
            bands_project = True
            data = Bands_files(
            load_dict= load_dict,fermi_energy_corr=fermi_energy_corr, bands_project = bands_project, from_bandsout = from_bandsout)
            return data #E_valence, E_conduction, E_valence_point,E_conduction_point]
    elif code == 'siesta':
        if type_bands == 'default':
            load_dict,fdf,out = Load(names = names,code = code,prefix=prefix)
            k_letters = get_kpath(fdf)
            bands =Bands_files(load_dict,code = code)
            return bands
