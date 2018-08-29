import numpy as np
import write as wt
import os

def readdata(filename):
  ishift=2
  i=0
  with file(filename) as input:
   for count, line in enumerate(input):
      i=i+1
  data =np.empty(i,dtype=object)
  i=0
  with file(filename) as input:
    for count, line in enumerate(input):
      data[i] = line.strip().split(',')      
      i       = i+1
  i  = 1 
  while ((data[i-1][6]==data[i][6]) and (i<data.shape[0]-ishift)):
    i = i + 1 
  ni=i

  j  = 1 
  while ((int(data[j-1][7])==int(data[j][7])-1) and (j<data.shape[0]-ishift)):
    j = j + 1 
  ni=i
  nj=j
  nij=ni*nj
  nsize=nij**2
  nproc=1
  dataout = np.empty([9,nproc,ni,nj] , dtype=float)

  for k in xrange(nproc):
    for i in xrange(ni):
      for j in xrange(nj):
        dataout[0][k][i][j] = float(data[k*nij**2+j*ni+i][1])
        dataout[1][k][i][j] = float(data[k*nij**2+j*ni+i][1])
        dataout[2][k][i][j] = float(data[k*nij**2+j*ni+i][1])/float(data[k*nij**2+j*ni+i][1])
        dataout[3][k][i][j] = float(data[k*nij**2+j*ni+i][2])
        dataout[4][k][i][j] = float(data[k*nij**2+j*ni+i][0])
        dataout[5][k][i][j] = float(data[k*nij**2+j*ni+i][6])
        dataout[6][k][i][j] = float(data[k*nij**2+j*ni+i][7])
        dataout[7][k][i][j] = float(data[k*nij**2+j*ni+i][8])        
        dataout[8][k][i][j] = float(data[k*nij**2+j*ni+i][-1])
  return dataout;

def readdatashift(filename):
  i=0
  with file(filename) as input:
   for count, line in enumerate(input):
      i=i+1
  data =np.empty(i,dtype=object)
  i=0
  with file(filename) as input:
    for count, line in enumerate(input):
      data[i] = line.strip().split(',')      
      i       = i+1
   
  i  = 1 
  for n in xrange(1,data.shape[0]-1):
    if(int(data[n-1][1])==int(data[n][1])-1): 
      i=i+1
  ni=i

  j  = 1 
  while ((int(data[j-1][2])==int(data[j][2])-1) and (j<data.shape[0]-1)):
    j = j + 1 
  ni=i
  nj=j
  nij=ni*nj
  nsize=nij
  nproc=1
  dataout = np.empty([6,nproc,ni,nj] , dtype=float)

  for k in xrange(nproc):
    for i in xrange(ni):
      for j in xrange(nj):
        dataout[0][k][i][j] = float(data[k*nij**2+i*nj+j][0])
        dataout[1][k][i][j] = float(data[k*nij**2+i*nj+j][1])
        dataout[2][k][i][j] = float(data[k*nij**2+i*nj+j][2])
        dataout[3][k][i][j] = float(data[k*nij**2+i*nj+j][3])
        dataout[4][k][i][j] = float(data[k*nij**2+i*nj+j][4])
        dataout[5][k][i][j] = float(data[k*nij**2+i*nj+j][5])
  return dataout;


def readdatashifttot(filename,jkstep):
  wt.writedatashifttot(filename,jkstep)
  ishift=jkstep
  i=0
  with file('tempread') as input:
   for count, line in enumerate(input):
      i=i+1
  data =np.empty(i-4,dtype=object)
  i=0
  with file('tempread') as input:
    for count, line in enumerate(input):
      if i>4:
        data[i-5] = line.strip().split(',')    
      i       = i+1

  ni  = 1 
  nj  = 1 

  for n in xrange(1,data.shape[0]-1):
    if(int(data[n-1][6])==int(data[n][6])-ishift): 
      ni=ni+1

  n=1
  while(n<data.shape[0] and (int(data[n-1][7])<=int(data[n][7]))):
    if(int(data[n-1][7])==int(data[n][7])-ishift):
      nj=nj+1 
    n=n+1

  n=1
  nproc=1  
  while(n<data.shape[0] and (int(data[n-1][0])<=int(data[n][0]))):
    if(int(data[n-1][0])!=int(data[n][0])):
      nproc=nproc+1 
    n=n+1

  n=1
  nshift=1  
  while(n<data.shape[0] and (int(data[n-1][9])<=int(data[n][9]))):
    if(int(data[n-1][9])!=int(data[n][9])):
      nshift=nshift+1 
    n=n+1

  nij=ni*nj
  nsize=nij
#
#  print "file length?",data.shape[0]
#
#  for n in xrange(1,data.shape[0]-1):
#    if( (int(data[n-1][6])!=int(data[n][6]))and (int(data[n-1][6])<int(data[n][6])) ): 
#      ni=ni+1
#
#  njfinal = 1
#  for n in xrange(1,data.shape[0]-1):
#    if( (int(data[n-1][7])!=int(data[n][7]))and (int(data[n-1][7])<int(data[n][7])) ): 
#      nj=nj+1
#    if(int(data[n-1][7])>int(data[n][7])):
#      if(nj>=njfinal):
#        njfinal=nj
#        nj=0
#      
#
#  nprocfinal = 1
#  nproc=1
#  for n in xrange(1,data.shape[0]-1):
#    if( (int(data[n-1][0])!=int(data[n][0]))and (int(data[n-1][0])<int(data[n][0])) ): 
#      nproc=nproc+1
#    if(int(data[n-1][0])>int(data[n][0])):
#      if(nproc>=nprocfinal):
#        nprocfinal=nproc
#        nproc=0
#
#  nshiftfinal = 0
#  nshift=1
#  for n in xrange(1,data.shape[0]-1):
#    if (int(data[n][9])==0):nshift=1
#    if( (int(data[n-1][9])<int(data[n][9])) ): 
#      nshift=nshift+1
#    if(nshift>=nshiftfinal):
#        nshiftfinal=nshift
#
#  if(nshift>=nshiftfinal):nshiftfinal=nshift
#
#  print "ni = ",ni
#  print "nj = ",njfinal
#  print "nproc = ",nprocfinal
#  print "nshift = ",nshiftfinal,nshift
#

  nij=ni*nj
  nsize=nij


  dataout = np.empty([11,nproc,ni,nj,nshift] , dtype=float)

  for k in xrange(nproc):
    for i in xrange(ni):
      for j in xrange(nj):
        for n in xrange(nshift):
          dataout[0][k][i][j][n] = float(data[k*nij**2+i*nj*nshift+j*nshift+n][0])
          dataout[1][k][i][j][n] = float(data[k*nij**2+i*nj*nshift+j*nshift+n][1])
          dataout[2][k][i][j][n] = float(data[k*nij**2+i*nj*nshift+j*nshift+n][2])
          dataout[3][k][i][j][n] = float(data[k*nij**2+i*nj*nshift+j*nshift+n][3])
          dataout[4][k][i][j][n] = float(data[k*nij**2+i*nj*nshift+j*nshift+n][4])
          dataout[5][k][i][j][n] = float(data[k*nij**2+i*nj*nshift+j*nshift+n][5])
          dataout[6][k][i][j][n] = float(data[k*nij**2+i*nj*nshift+j*nshift+n][6])
          dataout[7][k][i][j][n] = float(data[k*nij**2+i*nj*nshift+j*nshift+n][7])
          dataout[8][k][i][j][n] = float(data[k*nij**2+i*nj*nshift+j*nshift+n][8])
          dataout[9][k][i][j][n] = float(data[k*nij**2+i*nj*nshift+j*nshift+n][9])
 #         dataout[10][k][i][j][n] = float(data[k*nij**2+i*nj*nshift+j*nshift+n][10])
  os.remove('./tempread')
  return dataout;


def readdatashift_final(filename):
  i=0
  with file(filename) as input:
   for count, line in enumerate(input):
      i=i+1
  data =np.empty(i,dtype=object)
  i=0
  with file(filename) as input:
    for count, line in enumerate(input):
      data[i] = line.strip().split(',')    
      i       = i+1

  ni  = 1 
  nj  = 1 
  for n in xrange(1,data.shape[0]-1):
    if(int(data[n-1][1])==int(data[n][1])-1): 
      ni=ni+1

  n=1
  while(n<data.shape[0] and (int(data[n-1][2])<=int(data[n][2]))):
    if(int(data[n-1][2])==int(data[n][2])-1):
      nj=nj+1 
    n=n+1

  nij=ni*nj
  nsize=nij

  dataout = np.empty([6,ni,nj] , dtype=float)

  for i in xrange(ni):
    for j in xrange(nj):
      dataout[0][i][j] = float(data[i*nj+j][0])
      dataout[1][i][j] = float(data[i*nj+j][1])
      dataout[2][i][j] = float(data[i*nj+j][2])
      dataout[3][i][j] = float(data[i*nj+j][3])
      dataout[4][i][j] = float(data[i*nj+j][4])
      dataout[5][i][j] = float(data[i*nj+j][5])

  return dataout;
