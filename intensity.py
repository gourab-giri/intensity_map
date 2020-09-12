import pyPLUTO as pp
import numpy as np
from matplotlib import pyplot as plt
import math
import time

I = pp.Tools()
######################################################################################################
def viewangle(phi, theta, imx, imy, sx, sy, sz):
  
  pl=np.zeros([imx, imy, 5])
  nx=np.sin(theta)*np.cos(phi)
  ny=np.sin(theta)*np.sin(phi)
  nz=np.cos(theta)

  if nx >= 0:
    xi=0.0
    xf=sx
  else:
    xi=sx
    xf=0.0

  if ny >= 0:
    yi=0.0
    yf=sy
  else:
    yi=sy
    yf=0.0

  if nz >= 0:
    zi=0.0
    zf=sz
  else:
    zi=sz
    zf=0.0
     
  for xn in range((-imx/2), (imx/2)): #Step through x pixels of image frame
    for yn in range((-imy/2) ,(imy/2)): #Step through y pixels of image frame

      xp=sx/2+np.cos(theta)*np.cos(phi)*yn-np.sin(phi)*xn
      yp=sy/2+np.cos(theta)*np.sin(phi)*yn+np.cos(phi)*xn
      zp=sz/2-np.sin(theta)*yn

      t0_x=(xi-xp)/nx
      tm_x=(xf-xp)/nx   

      t0_y=(yi-yp)/ny
      tm_y=(yf-yp)/ny
      
      t0_z=(zi-zp)/nz
      tm_z=(zf-zp)/nz   

      t0=np.array([t0_x, t0_y, t0_z])
      tf=np.array([tm_x, tm_y, tm_z]) 
      ti=min(abs(t0[0]), abs(t0[1]), abs(t0[2]))
      tm=min(abs(tm_x), abs(tm_y), abs(tm_z))

      for i in range(0,3):
        if ti==abs(t0[i]):
          tii=t0[i]  
        if tm==abs(tf[i]):
          tmm=tf[i]
 
  #print(t0)
      pl[(xn+imx/2),(yn+imy/2),:]=[round(xp), round(yp), round(zp), round(tii), round(tmm)]
  #print(np.shape(pl))
  return pl

####################################################################################################
def cells(phi, theta, pl, sx, sy, sz, f, Iv, gs, length_unit):
  
  #xp=np.array(pl[:,:,0])
  #yp=np.array(pl[:,:,1])
  #zp=np.array(pl[:,:,2])
  #t0=np.array(pl[:,:,3])
  tmin=np.amin(pl[:,:,3])
  tmax=np.amax(pl[:,:,4])

  #print(tmin)
  dIv=np.zeros(np.shape(Iv))
  j=np.zeros([sr, sz/2])
  for t in range(int(tmin), int(tmax)): #Numerical integration of the emisivity along the normal
    
    xc=np.array(pl[:,:,0]+t*np.sin(theta)*np.cos(phi), dtype=int)-sx/2
    yc=np.array(pl[:,:,1]+t*np.sin(theta)*np.sin(phi), dtype=int)-sy/2
    zc=abs(np.array(pl[:,:,2]+t*np.cos(theta), dtype=int)-sz/2)
    rc=np.array(np.sqrt(xc**2 + yc**2), dtype=int)
    ind=(t-tmin)*f
    w=ind%1
    ind=int(ind//1)
    j[:,:]=w*(jc[ind+1,:, :]-jc[ind, :, :])+jc[ind, :, :]

    for m in range(0, gs[1]): #Step through y pixels of image frame
      for n in range(0, gs[0]): #Step through x pixels of image frame
        if 0 <= rc[n, m] < sr and 0 <= zc[n, m] < sz/2:
          #print(round((t-tmin)*unitl/f))
          dIv[n,m]=j[rc[n, m],zc[n, m]]*length_unit
        else:
          dIv[n,m]=0

        Iv[n,m]=Iv[n,m]+dIv[n,m]
        log.write("{:^6},{:^6}{:^10}{:^8},{:^8},{:^8}{:^20}\n".format(n, m, int(round((t-tmin)*f)), xc[n, m],yc[n, m],zc[n, m], Iv[n,m]))
    print("Execution time: %s" % (time.time() - start_time))
  return Iv
  



############################################################################################################
def gridsize(sx, sy, sz, phi, theta):

  ty1=int(sx/(np.cos(theta)*np.cos(phi)))
  ty2=int(sy/(np.cos(theta)*np.sin(phi)))
  ty3=int(sz/(-np.sin(theta)))

  max_imy=min(abs(ty1), abs(ty2), abs(ty3))
  
  tx1=sx/(-np.sin(phi))
  tx2=sy/(np.cos(phi))

  max_imx=min(abs(tx1), abs(tx2))

  print("Image size; ")
  print(max_imx, max_imy)
  return (int(max_imx),int( max_imy))


###########################################################################################################
print("################################################################### \n Integration of emission for PLUTO data \nVersion 5.1\n ##########################################################################")

nr=pp.nlast_info(datatype="vtk") #Reading simulation information
out=[] #Creating an array to store coordinates for integration
sx=0
sy=0
sz=0 

#Read parameters from input file
fp = open('paramsee.dat', 'r')
buffer=fp.read()
bufarr=buffer.splitlines()

buffer1,fr=bufarr[0].split('=')
buffer1,phi=bufarr[1].split('=')
buffer1,theta=bufarr[2].split('=')
buffer1,d_name=bufarr[3].split('=')
buffer1,lengthunit=bufarr[4].split('=')
buffer1,m=bufarr[5].split('=')
buffer1,n=bufarr[6].split('=')

fr =int(fr)	#User input for frame number
phi=float(phi)	#User input for the Azimuth angle
theta=float(theta)  #User input for the polar angle
lengthunit=float(lengthunit) #User input for the length of a cell
m=int(m)
n=int(n)

#Avoid division by zero
if theta==0:
  theta=1e-40

if phi==0:
  phi=1e-40

#Loading dataframe
start_time=time.time()
D = pp.pload(fr,datatype="vtk")

sqz_m = int(len(D.x1)/m)
sqz_n = int(len(D.x2)/n)
lengthunit = lengthunit*sqz_m
#Determining the dimensions of the loaded data
sx=2*len(D.x1)/sqz_m
sy=2*len(D.x1)/sqz_m
sz=2*len(D.x2)/sqz_n
sr=len(D.x1)/sqz_m

print("Grid size")
print("x1: %f" %(sx))
print("x2: %f" %(sy))
print("x3: %f" %(sz))

#photon travel distance between saved files
crstme=(D.x1[1]-D.x1[0])/round((nr['time'])/(nr['nlast']))
crstme = crstme*sqz_m

#Calculating the size of the image
gs=gridsize(sx,sy,sz,phi,theta)
Iv=np.zeros([gs[0],gs[1]])

#Deteriming the position of image pixel in 3D grid environment
pl=viewangle(phi, theta, gs[0], gs[1], sx, sy, sz) 

#Determine the number of files necessary to account for light travel time
tmin=np.amin(pl[:,:,3])
tmax=np.amax(pl[:,:,4])
totalfiles=int(math.ceil((tmax-tmin)*crstme)+1)
print("Number of files that will be loaded for time of arival calculation: %i" %(totalfiles))

jc=np.zeros([(totalfiles), sr, sz/2])
#ac=np.zeros([(totalfiles), sx, sy, sz])

j=eval("D." + d_name)
jc[totalfiles-1, :, :]=I.congrid(j,(m,n))
#ac[totalfiles, :, :, :]=D.av

del(D)

for f in range(fr-totalfiles+1, fr):
  if f!=fr:
    D = pp.pload(f,datatype="vtk")
    
    findex=f-fr+totalfiles-1

    #Calculating emission and absorption coefficients
    j = eval("D." + d_name)
    jc[findex, :, :] = I.congrid(j,(m,n))
#    ac[findex, :, :, :] = D.av
    
    del(D)


log=open("intemis.log", 'w')
log.write("Calculating intensity map for theta={} and phi={}".format(theta,phi) +"\n########################################\n{:^12} {:^10} {:^24} {:^20}".format("IMx,IMy", "Time step", "x,y,z", "Intensity\n")) 

print("Calculating intensity map for frame " + str(fr))

#Integration of emissivity
Map=cells(phi, theta, pl, sx, sy, sz, crstme, Iv, gs, lengthunit)

#Writing output    
fname="imap"+pp.get_nstepstr(fr)+"_"+d_name+"_{0}_{1}".format(round(180/np.pi*theta),round(180/np.pi*phi))
print("writing file " + fname)
np.save(fname +".npy",Map)

#plt.imshow(np.transpose(Iv), origin='image')
plt.imshow(np.log10(Iv).T, origin='image')
plt.colorbar()
plt.savefig(fname+".png", dpi=600)

print("Execution time: %s" % (time.time() - start_time))

log.close()
######################################################################################################


