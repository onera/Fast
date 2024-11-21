import numpy as np
import sys

def writedatashifttot(filename,jkstep):
    ishift=jkstep
    i=0
    with file(filename) as input:
        for count, line in enumerate(input):
            i=i+1
    data =np.empty(i-4,dtype=object)
    i=0
    with file(filename) as input:
        for count, line in enumerate(input):
            if i>4:
                data[i-5] = line.strip().split(',')    
            i       = i+1

    shifttab=np.asarray([0,64,128,256,384,512,640,768])
    ishft = 2

    f = open('tempread','a')
    for  n in xrange(1,data.shape[0]-2):
        out = ''
        for count,x in enumerate(data[n]):
            out = out+str(',')+str(x)
        f.write( out[1:]+'\n')
        indexShift = np.argwhere(shifttab==(int(data[n][9])))
        if(indexShift != 7):
            if( int(data[n+1][9]) != shifttab[indexShift+1]):
                for i in range(indexShift[0][0]+1,len(shifttab)):
                    f.write('1,1.0,1.0,1,1,1,1,1,1,'+str(shifttab[i])+',1.0'+'\n')
    f.close()

    return [];
