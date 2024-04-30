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
warnings.filterwarnings('ignore')

def get_kpath(fdf_file):
    n_linha=0
    n_linha_BandLines=[]
    for line in fdf_file:
        #busca os pontos da BZ uzados
        for item in line:
            if item=="BandLines":
                n_linha_BandLines.append(n_linha)
        n_linha+=1
    #grava a linha na lista linhas
    #salva os nomes K_points usados no calculo
    K_points=[]
    for i in range(n_linha_BandLines[0]+1,n_linha_BandLines[1]):
        k_point=fdf_file[i][len(fdf_file[i])-1]
        if k_point == "\Gamma" or k_point=='G':
            k_point="$\\Gamma$"
        K_points.append(k_point)
    return K_points
def bands_from_siesta(BANDS):
    ef = float(BANDS[0][0])
    kmin = float(BANDS[1][0]);kmax = float(BANDS[1][1])
    emin = float(BANDS[2][0]);emax = float(BANDS[2][1])
    nbnd = int(BANDS[3][0]);nkpoints = int(BANDS[3][2])
    kmin_index = -1
    k_letter = []; k_path = []
    for i in range(3,len(BANDS)):
        if len(BANDS[i])>0:
            try:
                if float(BANDS[i][0])==kmin:
                    if kmin_index==-1:
                        kmin_index = i
                    else:
                        kmin_index_second = i
            except TypeError:
                print('erro')
                print(BANDS[i])
        if len(BANDS[i])==2:
            try:
                float(BANDS[i][1])
            except ValueError:           
                k_path.append(float(BANDS[i][0]))
                if BANDS[i][1][1:len(BANDS[i][1])-1] == "\Gamma" or BANDS[i][1][1:len(BANDS[i][1])-1] =='G' or BANDS[i][1][1:len(BANDS[i][1])-1] =='Gamma':
                    k_letter.append("$\Gamma $")
                else:
                    k_letter.append(BANDS[i][1][1:len(BANDS[i][1])-1] )
    k_points = [float(BANDS[kmin_index][0])]
    cont = -1
    data = {'k':[]}
    for i in range(nbnd):
        data['band_'+str(i+1)] = []
    e = []
    list_e = []
    for i in range(kmin_index,kmin_index_second):
        if cont==nbnd:
            #print(i)
            #print(BANDS[i][0])
            if len(k_points)==nkpoints:
                pass
            else:
                k_points.append(float(BANDS[i][0]))
            cont = -1+len(BANDS[i])
            if len(list_e)==nkpoints:
                break
            else:
                list_e.append(e)
            e = []
            for j in range(1,len(BANDS[i])):
                e.append(BANDS[i][j])      

        else:
            if i == kmin_index:
                for j in range(1,len(BANDS[i])):
                    e.append(BANDS[i][j])
            else:
                for j in range(len(BANDS[i])):
                    e.append(BANDS[i][j])
            cont+=len(BANDS[i])
    data['k'] = k_points
    keys = list(data.keys())
    for j in range(len(list_e[0])):
        e_bnd = []
        for i in range(len(list_e)):
            #print(i,j)
            e_bnd.append(float(list_e[i][j]))
        data[keys[j+1]] = array(e_bnd)-ef
    
    return  [DataFrame(data),k_letter,k_path]