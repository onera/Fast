from array import array
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cl
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import read as rd
import write as wt
import sys
import matplotlib.cm as cmx

majorLocator = MultipleLocator(10)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(1)

lab="HSW"

# Usage python plotData.py 0 0 => graph 2D
# Usage python plotData.py I 0 => graph 1D I constant TpsCPU(J)
# Usage python plotData.py 0 J => graph 1D J constant TpsCPU(I)

ref=0.085#0.03#0.15
maxval=2.0
Iplt=0
Jplt=0
shifted=0

Iplt=int(sys.argv[1])
Jplt=int(sys.argv[2])

#------------------------------- Read database ---------------------------
data = {}
data['Part0']   = rd.readdatashift_final('./Shift_part0')
#data['Part1']   = rd.readdatashift_final('./Shift_part1')
#data['Part2']   = rd.readdatashift_final('./Shift_part2')
#data['Part3']   = rd.readdatashift_final('./Shift_part3')
#data['Part4']   = rd.readdatashift_final('./Shift_part4')
#data['Part5']   = rd.readdatashift_final('./Shift_part5')
#data['Part6']   = rd.readdatashift_final('./Shift_part6')
#data['Part7']   = rd.readdatashift_final('./Shift_part7')
#data['Part8']   = rd.readdatashift_final('./Shift_part8')
#data['Part9']   = rd.readdatashift_final('./Shift_part9')
#data['Part10']   = rd.readdatashift_final('./Shift_part10')
#data['Part11']   = rd.readdatashift_final('./Shift_part11')
#data['Part12']   = rd.readdatashift_final('./Shift_part12')
#data['Part13']   = rd.readdatashift_final('./Shift_part13')
#data['Part14']   = rd.readdatashift_final('./Shift_part14')

minI = 1000
maxI = 0
minJ = 1000
maxJ = 0

for k in data.keys():
    print("Case:",k)
    print("Shape ",data[k].shape)
    print("imin,imax ",data[k][1,0,0],data[k][1,-1,-1])
    print("jmin,jmax ",data[k][2,0,0],data[k][2,-1,-1])
    if data[k][1,0,0]<=minI:  minI=data[k][1,0,0]
    if data[k][1,-1,-1]>=maxI:maxI=data[k][1,-1,-1]
    if data[k][2,0,0]<=minJ:  minJ=data[k][2,0,0]
    if data[k][2,-1,-1]>=maxJ:maxJ=data[k][2,-1,-1]

print("Min Max I",minI,maxI)
print("Min Max J",minJ,maxJ)

#------------------------------- Set figure properties ---------------------

fig, ax = plt.subplots(figsize=(25.0,16.0))
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=25)
params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

##-------------------------------- Set colors ------------------------------

#cmp = plt.get_cmap('Blues')
cmp = plt.get_cmap('gist_ncar')
cNorm  = cl.Normalize(vmin=0, vmax=100)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmp)

cmps = plt.get_cmap('Blues')
cNorms  = cl.Normalize(vmin=0, vmax=255)
scalarMaps = cmx.ScalarMappable(norm=cNorms, cmap=cmps)

color=np.empty([257,3],dtype=float)
for i in xrange(256):
    color[i,0:3]=scalarMaps.to_rgba(i)[0:3]

color[256,0:3]=[102.0/255.0,0.0,51.0/255.0]

blue_cstm=LinearSegmentedColormap.from_list('blue_cstm', color)

origin = 'lower'

X={}
Y={}
map2d={}

for k in data.keys():   
#   ni       = data[k].shape[2]
#   nj       = data[k].shape[3]
    ni       = data[k].shape[1]
    nj       = data[k].shape[2]
    X[k]     = np.empty([ni,nj],dtype=float)
    Y[k]     = np.empty([ni,nj],dtype=float)
    map2d[k] = np.empty([ni,nj],dtype=float)

lev=np.arange(1.0,maxval,0.01)

#================= Create data to plot ======================

for k in data.keys():
    for i in xrange(0,data[k].shape[1]):
        for j in xrange(0,data[k].shape[2]):
#            map2d[k][i,j]=data[k][1,0,i,j,0]/ref*10**6*data[k][2,0,i,j,0]
            map2d[k][i,j]=data[k][4,i,j]/ref*10**6
            if(map2d[k][i,j]>=lev[-1]):map2d[k][i,j]=lev[-1]
            X[k][i,j]    =data[k][1,i,j]
            Y[k][i,j]    =data[k][2,i,j]

# 2D Map
ctrue=True
if(Iplt+Jplt==0):
    if(shifted==0):
        for k in data.keys():
            Perf = plt.contourf(X[k], Y[k], map2d[k], 50,
                                #[-1, -0.1, 0, 0.1],
                                #alpha=0.5,
                                levels=lev,
                                cmap=blue_cstm,#plt.get_cmap('Blues'),
                                origin=origin)
            Perf2 = plt.contour(Perf, levels=[1.5,2.0,2.5,3.0,4.0,5.0],#Perf.levels[::10],
                                colors='r',
                                hold='on')
            if(ctrue):
                cbar = plt.colorbar(Perf)
                cbar.ax.set_ylabel('Effective time 'r'($\mu s$)')
                cbar.ax.set_ylabel("Performance ratio / "+str(ref)+"$\mu s$")
                ctrue=False
##
            axes = plt.gca()
            axes.set_ylim([6,260])
            axes.set_xlim([6,300])

            cbar.add_lines(Perf2)
            plt.xlabel('Size I',fontsize=25)
            plt.ylabel('Size J',fontsize=25)

            ax.xaxis.set_major_locator(majorLocator)
            ax.xaxis.set_major_formatter(majorFormatter)
            ax.xaxis.set_minor_locator(minorLocator)
        
            ax.yaxis.set_major_locator(majorLocator)
            ax.yaxis.set_major_formatter(majorFormatter)
            ax.yaxis.set_minor_locator(minorLocator)
        
            ax.tick_params('both', length=10, width=1, which='major',direction='out')
            ax.tick_params('both', length=5, width=1, which='minor',direction='out')
            plt.savefig('effectivesize_2D_noshift_'+k+'.png')
            plt.cla()
    else:
        print("NRY !!! ")

else:
    kcolors = 5
    if(Iplt!=0):        
        for k,value in sorted(data.items()):
            kcolors = kcolors + 20
            ipos = np.argwhere(X[k][:,0]==Iplt)
            if len(ipos) == 0: 
                print("I position not found for case "+k+", X is ",X[k][:,0])
                sys.exit()                
            ipos = np.squeeze(ipos)
            print(map2d[k][ipos,:])
            plt.plot(Y[k][ipos,:],map2d[k][ipos,:],                                                                                                             
                 marker='o',                                                                                                                    
                 markersize=10,                                                                                                                 
                 linewidth=2,                                                                                                                    
                 linestyle='-',     
                 label=k+", Constant Size I ="+str(Iplt),
                 color=scalarMap.to_rgba(kcolors)[0:3])  
        plt.xlabel('Size J',fontsize=25)
    else:
        for k,value in sorted(data.items()):
            kcolors = kcolors + 20
            jpos = np.argwhere(Y[k][0,:]==Jplt)
            if len(jpos) == 0: 
                print("J position not found for case "+k+", Y is ",Y[k][0,:])
                sys.exit()                
            jpos = np.squeeze(jpos)
            plt.plot(X[k][:,jpos],map2d[k][:,jpos],                                                                                                             
                     marker='o',                                                                                                                    
                     markersize=10,                                                                                                                 
                     linewidth=2,                                                                                                                    
                     linestyle='-',     
                     label=k+", Constant Size J ="+str(Jplt),
                     color=scalarMap.to_rgba(kcolors)[0:3])          
        plt.xlabel('Size I',fontsize=25)
    axes = plt.gca()
#    axes.set_ylim([1,maxval])
#    axes.set_xlim([6,505])       

    plt.ylabel("Performance ratio / "+str(ref)+"$\mu s$",fontsize=25)

    plt.legend(fontsize=25,frameon=False,loc=0)
    plt.savefig('effectivesize_1D_noshift_'+lab+'.png')
##---------------------------------- Save figure -------------------------------------
#

#
plt.show()
