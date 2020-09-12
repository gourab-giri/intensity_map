from pylab import *
import nice_plots as nplt
nplt.nplts(fs=24)
import math

#l = 10
I = np.load('imap0495_JSyf3_70_0.npy')
Q = np.load('imap0495_QSyf3_70_0.npy')
U = np.load('imap0495_USyf3_70_0.npy')
I = I[:1024,544:1089]
Q = Q[:1024,544:1089]
U = U[:1024,544:1089]
I = I
#x,y = np.shape(I)
#for i in range(0,10,1):
#    for j in range(0,10,1):
#        theta = 0.5*np.arctan(U[j,i]/Q[j,i])
#        x0 = i
#        y0 = j
#        x1 = x0 + l*np.cos(theta)
#        y1 = y0 + l*np.sin(theta)
#        x_vals = [x0,x1]
#        y_vals = [y0,y1]
#mindex1 = np.where(U < 1e-18)
#mindex2 = np.where(Q < 1e-3)
#U[mindex1[0],mindex1[1]] = 0.0
#Q[mindex2[0],mindex2[1]] = 1.0e-18
theta = 0.5*np.arctan(U/Q)
length = 7.0*np.sqrt(np.multiply(Q,Q)+np.multiply(U,U))/I
x,y = np.shape(I)
for i in range(0,549,20):
    for j in range(0,1024,30):
            l = length[j,i]
            #print(l)
            c = theta[j,i] + math.pi/2
            x0 = i*0.0377
            x1 = x0 + l*np.cos(c)*0.037
            if j<512:
               y0 = j*0.0377-19.37
               y1 = y0 + l*np.sin(c)*0.037
            else:
               y0 = j*0.037-19.37
               y1 = y0 + l*np.sin(c)*0.037
            x_vals = [x0,x1]
            y_vals = [y0,y1]
            plt.plot(y_vals, x_vals,'k',linewidth=1.5)
#        if Q[j,i] > -0.1:
#            x0 = 0.0
#            y0 = 0.0
#            x1 = 0.0
#            y1 = 0.0
#            x_vals = [x0,x1]
#            y_vals = [y0,y1]
#            plt.plot(x_vals, y_vals)
#        else:
#            c = theta[j,i]+90
#            x0 = i
#            y0 = j
#            x1 = x0 + l*np.cos(c)
#            y1 = y0 + l*np.sin(c)
#            x_vals = [x0,x1]
#            y_vals = [y0,y1]
#            plt.plot(x_vals, y_vals)


#mindex = np.where(I1 < 1e-18)
#I1[mindex[0],mindex[1]] = 1.0e-18
#mindex = np.where(I2 < 1e-18)
#I2[mindex[0],mindex[1]] = 1.0e-18
#frac_pol = 100*np.sqrt(np.multiply(Q,Q)+np.multiply(U,U))/I
#plt.plot(I.T)
#plt.colorbar()
#plt.imshow(np.log10(I).T,origin='image', extent = [-19.37,19.37,0.0,20.55])

plt.imshow(np.log10(I/(4.0*math.pi)).T,origin='image',vmax = -11,vmin = -18,extent = [-19.37,19.37,0.0,20.55])
plt.colorbar()

#plt.title('percentage fractional polarization-map')
#plt.title('B field map')
plt.xlabel('$\\Delta$RA in arcsec')
plt.ylabel('$\\Delta$Dec in arcsec')
plt.savefig('polarization_map')
plt.show()

