import numpy as np
import read as rd

# Create padding file (binary read by FastS) from multiple Shift_partXX text files created by compute_shift.py

# Les fichiers Shift_partXX sont obtenus a l'issue de runs de compute_shift.py (sortie Shift) auquel on a enlevé les lignes d'entête et l'éventuelle dimension incomplète en I (en fin de fichier)

filename = 'padding_BDW.bin'

data = {}
data['Part0']   = rd.readdatashift_final('./Shift_part0')
data['Part1']   = rd.readdatashift_final('./Shift_part1')
data['Part2']   = rd.readdatashift_final('./Shift_part2')
data['Part3']   = rd.readdatashift_final('./Shift_part3')
data['Part4']   = rd.readdatashift_final('./Shift_part4')
data['Part5']   = rd.readdatashift_final('./Shift_part5')
data['Part6']   = rd.readdatashift_final('./Shift_part6')
data['Part7']   = rd.readdatashift_final('./Shift_part7')
data['Part8']   = rd.readdatashift_final('./Shift_part8')
data['Part9']   = rd.readdatashift_final('./Shift_part9')
data['Part10']   = rd.readdatashift_final('./Shift_part10')
data['Part11']   = rd.readdatashift_final('./Shift_part11')
data['Part12']   = rd.readdatashift_final('./Shift_part12')
data['Part13']   = rd.readdatashift_final('./Shift_part13')
data['Part14']   = rd.readdatashift_final('./Shift_part14')

Data = np.empty(1000*1000+2,dtype='int32')
Data2D = np.empty([1000,1000],dtype='int32')

for k in data.keys():
    print("Case:",k)
    i1 = int(data[k][1,0,0])
    i2 = int(data[k][1,-1,-1])
    j1 = int(data[k][2,0,0])
    j2 = int(data[k][2,-1,-1])
    print("i1, i2",i1,i2)
    print("j1, j2",j1,j2)
    Data2D[i1:i2+1,j1:j2+1] = data[k][5,:,:].astype('int32')

Data[0] = 1000
Data[1] = 1000

Data[2:] = np.reshape(Data2D,1000000)

Data.astype('int32').tofile(filename)
