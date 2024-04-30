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

####global functions
#--------------------------------------------------------------------------------------------------------------

#general functions

#find function, use to find a word in a file. input reg: string, input linhas: list of list from the orifinal file
def find(reg, linhas):  
    # cria a lista que vai ser retornada como resultado da funcao com a palavra e seu valor
    busca_resultado = []
    # for para varrer todas as linhas, o valor de i indica qual á linha em linhas[i]
    for i in range(len(linhas)):
        # for para varrer cada linha a procura da palavra buscada
        for item in linhas[i]:
            # verifica se o registro esta na linha
            if item == reg:
                # este for e para varrer novamente a linha se achar o registro para gravar todos os valores da linha em uma lista
                for item in linhas[i]:
                    if item != reg:
                        busca_resultado.append(item)
    return busca_resultado

#find n line function, use to find a line number of a word in a file. input reg: string; input linhas: list of list from the orifinal file
def find_n_line(reg, linhas):
    # cria a lista que vai ser retornada como resultado da funcao com a palavra e seu valor
    busca_resultado = []
    # for para varrer todas as linhas, o valor de i indica qual á linha em linhas[i]
    j = 0
    for i in range(0, len(linhas)):
        # for para varrer cada linha a procura da palavra buscada
        for item in linhas[i]:
            # verifica se o registro esta na linha
            if item == reg:
                # este for e para varrer novamente a linha se achar o registro para gravar todos os valores da linha em uma lista
                busca_resultado.append(j)
        j += 1
    return busca_resultado


# function that transform file in to list. input name: path and the name of the file.
def open_file(name):
    file = open(name, 'r', encoding='utf-8',
                errors='ignore')  # open the file.
    # creat a list that will be full filed with the file with a line as a list element.
    file_list = []
    for line in file:  # open a loop that cover the file.
        line = line.strip('\n')  # drop out all '\n' contained in every line.
        # change the spaces for a element of a list, ex: 'the energy is' --> ['the','energy','is'].
        line = line.split()
        file_list.append(line)  # add the line in the list file_list.
    file.close()  # close de file.
    return file_list

def number_identify(string):
        number = ''
        number_list = list(map(lambda x:str(x),range(10)))
        for caracter in string:
            if (caracter in number_list) == True:
                number+=caracter
            elif caracter == '.':
                number+=caracter
            elif caracter == '+':
                old_number = number
                number = ''
        return float(number)

#--------------------------------------------------------------------------------------------------------------

#functions to get information from scf.out, bands.out and nscf.out

#function to get the quantum expresso version from output. input: file, list of list from original output file
def version_quantum(file):
    for line in file:
        for item in line:
            if item == 'Program':
                version = line[2][2:]
    return version

#function to get the fermy energy from outpu. input: file, list of list from original output file; input type: string for whitch ouput, scf.out,nscf.out or bands.out; input version: string with the version from quantum expresso 
def E_fermi(file,type,version):
    i =0
    for line in file:
        for item in line:
            if version == '7.1':
                if type == 'bands.out' or type == 'nscf.out':
                    if item == "Fermi":
                        count = i
                    elif item == "occupied," or item == "occupied":
                        count = i
                elif type == 'scf.out':
                    if item == "Fermi":
                        count = i
                    elif item == "occupied," or item == "occupied":
                        count = i
            elif version == '7.2':
                if type == 'bands.out' or type == 'nscf.out' or 'scf.out':
                    if item == "Fermi":
                        count = i
                    elif item == "occupied," or item == "occupied":
                        count = i
            elif version == '7.0':
                if type == 'bands.out' or type == 'nscf.out':
                    if item == "Fermi":
                        count = i
                    elif item == "occupied," or item == "occupied":
                        count = i
                elif type == 'scf.out':
                    if item == "Fermi":
                        count = i
                    elif item == "occupied," or item == "occupied":
                        count = i                      
                        
            elif version == '6.7MaX':
                if type == 'scf.out':
                    if item == "Fermi":
                        count = i
                    elif item == "occupied," or  item ==  "occupied":
                        count = i    
                
        i+=1
    for item in file[count]:
        try:
            E_f = round(float(item), 4)
            break
        except ValueError:
                pass
    return E_f

#function to get the volume of cell. input: file, list of list from original output file
def volumn(file):
    for line in file:
        for i in range(len(line)):
            if line[i] == 'volume':
                vol =  float(line[i+2])
                return vol


#function to get the a latice. input: file, list of list from original output file
def alat(file):
    for line in file:
        for i in range(len(line)):
            if line[0]== "lattice" and line[1] == 'parameter':
                return float(line[4])

#--------------------------------------------------------------------------------------------------------------

#functions for pdos files

#function to correct the head for pdos files. input arquivo , file with the path from the original pdos file
def change_coluns_pdos(arquivo):
    # abre o arquivo scf.out
    arq = open(arquivo, 'r')
    # cria a lista que contenhaa cada linha do arquvio cif como um elemento da lista
    new_arq = []
    ct = 0
    for linha in arq:
        # retira o \n de quebra de linnha da linha
        line = linha.strip('\n')
        # retira os espaços entre os elementos da linha, assim criando uma lista
        line = line.split()
        if len(line)>0:
            if ct == 0 and line[0] == '#':
                line.remove('#')
                new_line = []
                for item in line:
                    if item != '(eV)':
                        new_line.append(item)
            else:
                new_line = line
                new_arq.append(new_line)
        ct += 1

    # fecha o arquivo
    arq.close()

    arq = open(arquivo,'w')
    for item in new_arq:
        line = ''
        for word in item:
            line+= word + '   '
        line += '\n'
        arq.writelines(line)
    arq.close()
    return 0

#function to get the atoms from scf.in file; input scf_in_file: list os list from original scf.in file
def atoms(scf_in_file):
    i = 0;atom = []
    for line in scf_in_file:
        for item in line:
            if item == 'ATOMIC_SPECIES' or item == 'atomic_species':
                count1 = i+1
            elif item == 'ATOMIC_POSITIONS' or item == 'atomic_positions':
                count2 = i
        i+=1
    for j in range(count1,count2):
        atom.append(scf_in_file[j][0])
    return atom


#function to get all of pdos files and covert in to a dicitionary delimted for orbiatis( s,p,etc) and for atoms. input: path , path from the pdos files are; input prefix, the prefix if exist in the pdos files; input atoms, atoms get from the functions atoms
def pdos_files_get(path,prefix,atoms):
    pdos_files = []
    for dir,subdir,files in os.walk(path):
        for file in files:
            if file[:len(prefix)] == prefix and file != prefix+'.pdos.dat.pdos_tot' and file != prefix + '.pdos-proj.dat.projwfc_up' and file != prefix+'.pdos_tot':
                pdos_files.append(os.path.join(path,file))

    pdos_atom_dict = {};pdos_orb_dict = {}
    for a in range(len(atoms)):
        list = [];orb_list = []
        for name in pdos_files:
            for i in range(1,len(name)):
                if len(atoms[a]) == 1:
                    if name[i] == '(' and name[i+1] == atoms[a]:
                        list.append(name)
                    elif name[i] == '(' and name[i+1] == 's':
                        if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])
                    elif name[i] == '(' and name[i+1] == 'p':
                        if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])
                    elif name[i] == '(' and name[i+1] == 'd':
                       if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])
                    elif name[i] == '(' and name[i+1] == 'f':
                        if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])                       
                else:
                    if name[i-1] == '(' and name[i] == atoms[a][0] and name[i+1] == atoms[a][1]:
                        list.append(name)
                    elif name[i] == '(' and name[i+1] == 's':
                        if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])
                    elif name[i] == '(' and name[i+1] == 'p':
                        if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])
                    elif name[i] == '(' and name[i+1] == 'd':
                       if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])
                    elif name[i] == '(' and name[i+1] == 'f':
                        if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])     

        pdos_atom_dict[atoms[a]] = list
        pdos_orb_dict[atoms[a]] = orb_list

    pdos_atom_orb_dict = {}
    for atom in atoms:
        pdos_atom_orb_dict[atom] = {}
        for orb in pdos_orb_dict[atom]:
            pdos_atom_orb_dict[atom][orb] = []
            list = []
            for file in pdos_atom_dict[atom]:
                for i in range(len(file)):
                    if file[i] == '(' and file[i+1] == orb:
                        list.append(file)
            pdos_atom_orb_dict[atom][orb] = list

    return [pdos_atom_dict,pdos_atom_orb_dict]


#--------------------------------------------------------------------------------------------------------------

#functions for bands files
#function to transform the bands.dat.gnu file in to a data frame. input arquivo, is the orifinal aquivo bands.dat.gnu
def bandas_df(arquivo):
    bandas = read_csv(arquivo, delim_whitespace=True,
                         names=['k', 'E'], dtype=float)
    # subtrai todos os valores pelo nivel de fermi
    return bandas

#function to transform the original dataframe from function bandas_df for a new dataframe separeted for orbitals. input bandas, dataframme obtained from function bandas_df; input ef, bool varable, if u want to correct the fermi level or not; input Ef: float, the value of fermyenergy
def orbitais(bandas, ef, Ef):
    # nesse bloco sera feito a separaçao do data frame em varias colunas, onde cada coluna é referente a uma banda
    # neste passo vamos obter do dataFrame Bandas original apenas os valores de momento que nao são repetidos
    # cria a lista para armazenaros valroes
    if ef == True:
        bandas['E'] -= Ef
    momentos = []
    # for para varrer todos os indices do dataframe
    for i in bandas.index:
        # pega o valor do momento para o indice i, no caso a linha i do data frame
        mom = bandas.loc[i][0]
        # cria uma condicao para adicionar na lista apenas itens nao repetidos, ou seja se mom nao estiver na  lista momentos, entao ele sera adicionado a mesma, caso contrario nada e feito
        condicao = mom in momentos
        if condicao == False:
            momentos.append(mom)
    # variavel para fazer a contagem total do numero de bandas
    n_bnd = 0
    # para fazer a contagem vamso pegar todos os valores de k na coluna "k" do dataframe Bandas, sempre que k=0, quer dize que começamos uma nova banda, portanto a
    # condicao de troca de bandas sera quando k=0, sempre que isso acontecer, adicionamos +1 no valor de n_bnd(numero de bandas)
    cond = bandas["k"] == 0
    # começa o for varrendo todos os valores de k contidos na coluna bandas["k"]
    for item in cond:
        if item == True:
            n_bnd += 1
        else:
            pass
    # cria o novo dataFrame com a coluna contendo apenas os valores de k que nao sao repetidos
    df = DataFrame(momentos, columns=["K"], dtype=float)
    # começa com o menor valor posivel para a energia da banda de valencia
    E_valencia = min(bandas['E'])
    # começa com o maior valor posivel para a energia da banda de conduçao
    E_conducao = max(bandas['E'])
    gap = True
    # loop para varrer todas as bandas do data frame
    for j in range(n_bnd):
        # cria a lista que vai armazear os valores de energia para a banda j
        energia = []
        # loop para varrer todos os valores de energia contidos no dataframe original que coresponem a banda j
        for i in range(len(momentos)):
            e = bandas.loc[i+j*len(momentos)][1]
            energia.append(e)
        # cria uma nova coluna no novo dataframe adicionando todos os valores de energia da banda j
        df["orb "+str(j)] = energia
        #verifica se é metal
        e_min = min(energia); e_max = max(energia)
        if ef == False:
            if e_min <Ef and e_max>Ef:
                gap = False
        else:
            if e_min <0 and e_max>0:
                gap = False

        # atualiza o valor da energia de valencia
        if ef == False:
            if max(energia) >= E_valencia and max(energia) <= Ef:
                E_valencia = max(energia)
                for k in range(len(energia)):
                    if energia[k] == E_valencia:
                        E_valencia_point = [momentos[k],energia[k]]
            # atualiza o valor da energia de conduçao
            elif min(energia) <= E_conducao and min(energia) >= Ef:
                E_conducao = min(energia)
                for k in range(len(energia)):
                    if energia[k] == E_conducao:
                        E_conducao_point = [momentos[k],energia[k]]
        else:
            if max(energia) >= E_valencia and max(energia) <= 0:
                E_valencia = max(energia)
                for k in range(len(energia)):
                    if energia[k] == E_valencia:
                        E_valencia_point = [momentos[k],energia[k]]                
            # atualiza o valor da energia de conduçao
            elif min(energia) <= E_conducao and min(energia) >= 0:
                E_conducao = min(energia)
                for k in range(len(energia)):
                    if energia[k] == E_conducao:
                        E_conducao_point = [momentos[k],energia[k]]
    if gap == True:
        return [df, round(E_valencia,3), round(E_conducao,3),E_valencia_point,E_conducao_point]
    else:
        if ef == False:
            return [df, Ef, Ef,0,0]
        else:
            return [df,0, 0,0,0]

#function to get the kpoints letter from the bands.in. input bandin_file, list of list from the original bands.in file
def k_point(bandin_file):
    k_points_band = []
    k_points_letter = []
    nk = []
    i = 0
    for line in bandin_file:
        if len(line) >= 1:
            if line[0] == 'K_POINTS':
                n_line = i
        i += 1
    n_kpoints = int(bandin_file[n_line+1][0])
    for i in range(n_line+2, n_line+2+n_kpoints):
        k_points_band.append(array(bandin_file[i][:3], dtype=float))
        nk.append(int(bandin_file[i][3]))
        if bandin_file[i][4][1:] == 'G':
            k_points_letter.append('$\Gamma$')
        else:
            k_points_letter.append(bandin_file[i][4][1:])
    return [k_points_letter, k_points_band,nk]

#function to que n_spin from bands.in.
def n_spin(bandsin_file):
    marker = False
    for line in bandsin_file:
        if len(line) >=1:
            if line[0] == 'nspin=':
                marker = True
                nspin = int(line[1])
            elif line[0] == 'nspin':
                marker = True
                nspin = int(line[2])
            elif line[0][:-1] == 'nspin=':
                marker = True
                nspin = int(line[0][-1:])
                
    if marker == False:
        return 1
    else:
        return nspin

def verbosity(bandsin_file):
    marker = False
    ct = 0
    for i in range(len(bandsin_file)):
        for item in bandsin_file[i]:
            if ('verbosity' in item) == True:
                ct = i
    for item in  bandsin_file[ct]:
        if ('high' in item) == True:
            marker = True
            return 'high'
    if marker == False:
        return 'low'

#funtion to get the momentum axis values from bands.x.out file. input file, list of list from the original bands.x.out file. input n_k, number of high simmetric points in brilluion zone, geted from bands.in file.
def k_points_path(file, n_k):
    i = 0
    k_points_bandsx = []
    k_path = []
    for line in file:
        for item in line:
            if item == 'wavefunctions':
                n_line = i+1
        i += 1
    for j in range(n_line, n_line + n_k):
        
        
        try:
            n1 = file[j][2]
            n2 = file[j][3]
            n3 = file[j][4]
            k_points_bandsx.append(array(file[j][2:5], dtype=float))
            k_path.append(float(file[j][7]))
        except ValueError:
            try:
                n1 = float( file[j][2] )
                for k in range(len(file[j][3])):
                    if file[j][3][k] == '-':
                        count = k 
                        break
                n2 = float(file[j][3][:count])
                n3 = float(file[j][3][count:])
            except ValueError:
                for k in range(1,len(file[j][2])):
                    if file[j][2][k] == '-':
                        count = k 
                        break
                n1 = float(file[j][2][:count])
                try:
                    n2 = float(file[j][2][count:])
                    n3  = float(file[j][3])
                    
                except ValueError:
                    for k in range(count,len(file[j][2])):
                        if file[j][2][k] == '-':
                            count2 = k 
                            break
                        n2 = float(file[j][2][count:count2])
                        n3  = float(file[j][2][count2:])
            k_points_bandsx.append([n1,n2,n3] )
            l = 0
            for item in file[j]:
                if item == 'coordinate':
                    k_path.append(float(file[j][l+1]))
                    break
                l+=1
    return [k_path, k_points_bandsx]


#--------------------------------------------------------------------------------------------------------------
#functions for projected bands

#function to get projeceted contribuitions from every state on atomitc_proj.xml. Input: File, list of list from the original file
def atomic_proj(file):
    def number_identify(string):
        number = ''
        number_list = list(map(lambda x:str(x),range(10)))
        for caracter in string:
            if (caracter in number_list) == True:
                number+=caracter
            elif caracter == '.':
                number+=caracter
            elif caracter == '+':
                old_number = number
                number = ''
        return float(number)

    for line in file:
        if line[0] == '<HEADER':
            n_band = int(number_identify(line[1]))
            n_spin = int(number_identify(line[3]))
            n_wfc = int(number_identify(line[4]))


    n = 0
    for i in range(len(file)):
            if file[i][0] == '<PROJS>':
                ct = i
            elif file[i][0] == '</PROJS>':
                    n = i-ct
            elif n != 0:
                break
    wfc_proj = {}
    for i in range(1,n_wfc+1):
        wfc_proj[i] = {}
        for j in range(1,n_band+1):
            wfc_proj[i][j] = []
    for wfc in range(1,n_wfc+1):
            for i in range(len(file)):
                if file[i][0] == '<ATOMIC_WFC':
                    #print(file[i])
                    index = int(number_identify(file[i][1]))
                    if index == wfc:
                        ct = 1
                        for j in range(i+1,i+1+n_band):
                            psi_real = float(file[j][0]);psi_imag = float(file[j][1])
                            psi_2 = pow(psi_real,2) +pow(psi_imag,2)
                            wfc_proj[wfc][ct].append(psi_2)
                            ct+=1          
    return [wfc_proj,n_band,n_spin,n_wfc]

#function to get states from bn.prowfc.out. Input:File, list of list from the original file
def States(file):
    states_line = []; atom_number = {};atom_number2 = {}
    for line in file:
        if len(line) >=1:
            if line[0] == 'state':
                states_line.append(line)

    for line in states_line:
        atom_number2[line[5][1:]] = []
    for line in states_line:
        if (int(line[4])in atom_number2[line[5][1:]]) == False:
            atom_number2[line[5][1:]].append(int(line[4]))
        atom_number[int(line[4])] = line[5][1:]
    atom_orb = {}
    for key in atom_number.keys():
        marc_s = False; marc_p = False; marc_d = False
        for line in states_line:
            if int(line[4]) == key:
                if line[8] == '1':
                    if marc_s == False:
                        atom_orb[key] = {'s':[line[2][:-1]]}
                        marc_s = True
                    else:
                        atom_orb[key]['s'].append(line[2][:-1])
                elif line[8] == '2':
                    if marc_p == False:
                        atom_orb[key].update({'p':[line[2][:-1]]})
                        marc_p = True
                    else:
                        atom_orb[key]['p'].append(line[2][:-1])
                elif line[8] == '3':
                    if marc_d == False:
                        atom_orb[key].update({'d':[line[2][:-1]]})
                        marc_d = True
                    else:
                        atom_orb[key]['d'].append(line[2][:-1])
        
    return [atom_number2,atom_orb]

#function to transform all information in to dictonary. Input: atom_number, dict with atom and the correspondent wfc number;Input: atom_number_orb, dict with atom and the correspondent state, Input: proj, geted from atomic_proj function, Input: nbands, number of bands
def projcted_atoms(atom_number,atom_orb_number,proj,n_band):
    atom_orb = {}
    for atom in atom_number.keys():
        atom_orb[atom] = {}
        for number in atom_number[atom]:
            for orb in atom_orb_number[number].keys():
                s=0;p=0;d=0;marker_s = False;marker_p = False;marker_d = False
                for item in atom_orb_number[number][orb]:
                    if orb == 's':                      
                        atom_orb[atom][orb] = proj[int(item)]
                    elif orb == 'p':
                        if marker_p == False:
                            atom_orb[atom][orb] = {}
                            marker_p = True
                        if p == 0:
                            atom_orb[atom][orb].update({'pz':proj[int(item)]})
                        elif p == 1:
                            atom_orb[atom][orb].update({'px':proj[int(item)]})
                        elif p == 2:
                            atom_orb[atom][orb].update({'py':proj[int(item)]})
                        p+=1
                    elif orb == 'd':
                        if marker_d == False:
                            atom_orb[atom][orb] = {}
                            marker_d = True
                        if d == 0:
                            atom_orb[atom][orb].update({'dz2':proj[int(item)]})
                        elif d == 1:
                            atom_orb[atom][orb].update({'dzx':proj[int(item)]})
                        elif d == 2:
                            atom_orb[atom][orb].update({'dzy':proj[int(item)]})
                        elif d == 3:
                            atom_orb[atom][orb].update({'dx2-y2':proj[int(item)]})
                        elif d == 4:
                            atom_orb[atom][orb].update({'dxy':proj[int(item)]})
                        d+=1                        

    for atom in atom_orb.keys():    
        temp = {};temp2 = {}; temp3 = {}
        for i in range(1,n_band+1):
            for orb in atom_orb[atom]:
                marker = False
                if orb == 's':
                    temp2[i] = array(atom_orb[atom][orb][i])
                    temp3[i] = array(atom_orb[atom][orb][i])
                elif orb == 'p':
                    for item in atom_orb[atom][orb]:
                        if marker == False:
                            temp[i] = array(atom_orb[atom][orb][item][i])
                            marker = True
                        else:
                            temp[i] = temp[i] + array(atom_orb[atom][orb][item][i])
                        temp2[i] = temp2[i] + temp[i]
                        temp3[i] = temp3[i] + temp[i]
                elif orb == 'd':
                    for item in atom_orb[atom][orb]:
                        if marker == False:
                            temp[i] = array(atom_orb[atom][orb][item][i])
                            marker = True
                        else:
                            temp[i] = temp[i] + array(atom_orb[atom][orb][item][i])
                        temp2[i] = temp2[i] + temp[i]
        atom_orb[atom].update({'s+p':temp3})                        
        atom_orb[atom].update({'tot':temp2})            
        atom_orb[atom]['p'].update({'tot':temp})
    return atom_orb

#function to get band data from bands.out. Input: file, list of list from the original bands.out file
def band_structure_from_bandsout(file, fermy_energ_cor, fermy_energy,version):
    n_band = 0
    for line in file:
        if len(line)>=3:
            if line[2] == 'Kohn-Sham':
                n_band = int(line[4])
    bands = {}
    for e in range(1,n_band+1):
        bands[e] = []
    k_value = [];k_tol = 0
    for i in range(len(file)):
        if len(file[i])>=1:
            if file[i][0] == 'k':
                line = file[i]
                if len(k_value)>=1:
                    if ( '-' in line[1]) == True:
                        try:
                            kx = float(line[1][1:])
                            ky = float(line[2])
                            kz = float(line[3])
                        except ValueError:
                            try:
                                kx = float(line[1][1:])
                                for k in range(len(line[2])):
                                    if line[2][k] == '-':
                                        count = k 
                                        break
                                ky = float(line[2][:count])
                                kz = float(line[2][count:])
                            except ValueError:
                                for k in range(2,len(line[1])):
                                    if line[1][k] == '-':
                                        count = k 
                                        break
                                kx = float(line[1][1:count])
                                try:
                                    ky = float(line[1][count:])
                                    kz  = float(line[2])
                                except ValueError:
                                    for k in range(count+1,len(line[1])):
                                        if line[1][k] == '-':
                                            count2 = k 
                                            break
                                    ky = float(line[1][count:count2])
                                    kz  = float(line[1][count2:])
                    
                    else:
                    
                        try:
                            kx = float(line[2])
                            ky =  float(line[3])
                            kz = float(line[4])
                        except ValueError:
                            try:
                                kx = float(line[2])
                                for k in range(len(line[3])):
                                    if line[3][k] == '-':
                                        count = k 
                                        break
                                ky = float(line[3][:count])
                                kz = float(line[3][count:])
                            except ValueError:
                                for k in range(1,len(line[2])):
                                    if line[2][k] == '-':
                                        count = k 
                                        break
                                kx = float(line[2][:count])
                                try:
                                    ky = float(line[2][count:])
                                    kz  = float(line[3])
                                    
                                except ValueError:
                                    for k in range(count+1,len(line[2])):
                                        if line[2][k] == '-':
                                            count2 = k 
                                            break
                                    ky = float(line[2][count:count2])
                                    kz  = float(line[2][count2:])
                         
    
                    k_mod = sqrt((kx-kx_old)**2+(ky-ky_old)**2+(kz-kz_old)**2)
                    k_tol+=k_mod
                    k_value.append(round(k_tol,4))
                    del(kx_old);del(ky_old);del(kz_old)
                    kx_old = kx; ky_old = ky; kz_old = kz
                else:
                    try:
                        kx = float(line[2])
                        ky =  float(line[3])
                        kz = float(line[4])
                    except ValueError:
                        try:
                            kx = float(line[2])
                            for k in range(len(line[3])):
                                if line[3][k] == '-':
                                    count = k 
                                    break
                            ky = float(line[3][:count])
                            kz = float(line[3][count:])
                        except ValueError:
                            for k in range(1,len(line[2])):
                                if line[2][k] == '-':
                                    count = k 
                                    break
                            kx = float(line[2][:count])
                            try:
                                ky = float(line[2][count:])
                                kz  = float(line[3])
                                
                            except ValueError:
                                for k in range(count,len(line[2])):
                                    if line[2][k] == '-':
                                        count2 = k 
                                        break
                                    ky = float(line[2][count:count2])
                                    kz  = float(line[2][count2:])
                         
            

                    k_mod = sqrt(kx**2+ky**2+kz**2)
                    k_tol+=k_mod
                    k_value.append(round(k_tol,4))
                    kx_old = kx; ky_old = ky; kz_old = kz
                    #print(kx,ky,kz)
                    #print(round(k_tol,4))
                bands_line = []
                for k in range(i+1,i+len(file)):
                    if len(file[k]) >=1: 
                        if file[k][0] == 'k' or file[k][0] == 'Writing' or file[k][0] == 'highest':
                            break
                        else:
                            if version == '7.1' or version == '7.0' or version == '7.2':
                                if file[k][0] != ' ':
                                    for item in file[k]:
                                        bands_line.append(float(item))
                            elif version == '6.7MaX' or version == '6.7':
                                if file[k][0] != ' ' and file[k][0] != 'occupation' and file[k][0] != '1.0000' and file[k][0] != '0.0000':
                                    for item in file[k]:
                                        bands_line.append(float(item))
                for j in range(n_band):
                    bands[j+1].append(bands_line[j])
    df = DataFrame({'k':k_value})
    #print(df)
    for e in range(1,n_band+1):
        if fermy_energ_cor == True:
            df[str(e)]= array(bands[e]) - fermy_energy
        elif fermy_energ_cor == False:
            df[str(e)]= bands[e] 
        else:
            print('Worng value for fermy energy correction parameter, only True and False are accepted. Default: True')
    return df

#--------------------------------------------------------------------------------------------------------------
#functions for change density

def get_vectors(file):
    e11 = 'e1(1)';e12 = 'e1(2)';e13 = 'e1(3)'
    e1 = [e11,e12,e13]
    e21 = 'e2(1)';e22 = 'e2(2)';e23 = 'e2(3)'
    e2 = [e21,e22,e23]
    x1 = 'x0(1)';x2 = 'x0(2)';x3= 'x0(3)'
    x_list = [x1,x2,x3]
    e1_v = []; e2_v = []; x_0 = []
    for e in e1:
        for line in file:
            for i in range(len(line)):
                if line[i] == e:
                    if (',' in line[i+2]):
                        e1_v.append(float(line[i+2][:-1]))
                    else:                  
                        e1_v.append(float(line[i+2]))
                        
    for e in e2:
        for line in file:
            for i in range(len(line)):
                if line[i] == e:
                    if (',' in line[i+2]):
                        e2_v.append(float(line[i+2][:-1]))
                    else:                  
                        e2_v.append(float(line[i+2]))
                        
    for x in x_list:
        for line in file:
            for i in range(len(line)):
                if line[i] == x:
                    if (',' in line[i+2]):
                        x_0.append(float(line[i+2][:-1]))
                    else:                  
                        x_0.append(float(line[i+2]))
                        
                        
    for line in file:
       for i in range(len(line)):
           if line[i] == 'nx':
               nx = line[i+2]
           elif line[i] == 'ny':
               ny = line[i+2]
    return [e1_v,e2_v,x_0,int(nx),int(ny)]

def hex_to_RGB(hex_str):
    """ #FFFFFF -> [255,255,255]"""
    #Pass 16 to the integer function for change of base
    return [int(hex_str[i:i+2], 16) for i in range(1,6,2)]

def get_color_gradient(c1, c2, n):
    """
    Given two hex colors, returns a color gradient
    with n colors.
    """
    assert n > 1
    c1_rgb = array(hex_to_RGB(c1))/255
    c2_rgb = array(hex_to_RGB(c2))/255
    mix_pcts = [x/(n-1) for x in range(n)]
    rgb_colors = [((1-mix)*c1_rgb + (mix*c2_rgb)) for mix in mix_pcts]
    return ["#" + "".join([format(int(round(val*255)), "02x") for val in item]) for item in rgb_colors]
