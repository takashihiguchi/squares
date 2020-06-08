#!/usr/bin/python

# Fri May 24 10:11:45 CDT 2019 Jeff added this line.

# Tue Feb 11 13:43:43 CST 2020 Jeff taking original patch.py and
# updating to solve the zero mode issue.  Will now update to use the
# patchlib submodule.

# Fri Feb 14 11:45:17 CST 2020 Jeff speeding up code by improving the
# sensorarray class to use numpy structures.

from scipy.constants import mu_0, pi
import numpy as np
from patchlib.patch import *
from Pis.Pislib import *

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
mpl.rcParams['legend.fontsize'] = 10

from scipy.optimize import curve_fit




class coilcube:
    def __init__(self,xdim,ydim,zdim,corners):
        self.xdim = xdim
        self.ydim = ydim
        self.zdim = zdim
        self.corners = corners
        x = corners[1]-corners[0]
        y = corners[2]-corners[0]
        z = corners[3]-corners[0]
        self.face = []
        thesecorners=(corners[0],corners[1],corners[2])
        self.face.append(face(xdim,ydim,thesecorners))
        thesecorners=thesecorners+z
        self.face.append(face(xdim,ydim,thesecorners))
        thesecorners=(corners[0],corners[1],corners[3])
        self.face.append(face(xdim,zdim,thesecorners))
        thesecorners=thesecorners+y
        self.face.append(face(xdim,zdim,thesecorners))
        thesecorners=(corners[0],corners[2],corners[3])
        self.face.append(face(ydim,zdim,thesecorners))
        thesecorners=thesecorners+x
        self.face.append(face(ydim,zdim,thesecorners))
        self.numcoils=(xdim*ydim+xdim*zdim+ydim*zdim)*2

    def coil(self,number):
        xdim=self.xdim
        ydim=self.ydim
        zdim=self.zdim
        if(number<xdim*ydim):
            return self.face[0].coil[number]
        elif(number<xdim*ydim*2):
            return self.face[1].coil[number-xdim*ydim]
        elif(number<xdim*ydim*2+xdim*zdim):
            return self.face[2].coil[number-xdim*ydim*2]
        elif(number<xdim*ydim*2+xdim*zdim*2):
            return self.face[3].coil[number-xdim*ydim*2-xdim*zdim]
        elif(number<xdim*ydim*2+xdim*zdim*2+ydim*zdim):
            return self.face[4].coil[number-xdim*ydim*2-xdim*zdim*2]
        else:
            return self.face[5].coil[number-xdim*ydim*2-xdim*zdim*2-ydim*zdim]

    def set_independent_current(self,number,current):
        # # wire the last two coils together
        # # only works if xdim=ydim=zdim=1
        # made it back to not wired together Jeff
        xdim=self.xdim
        ydim=self.ydim
        zdim=self.zdim
        if(number<xdim*ydim):
            self.face[0].coil[number].set_current(current)
        elif(number<xdim*ydim*2):
            self.face[1].coil[number-xdim*ydim].set_current(current)
        elif(number<xdim*ydim*2+xdim*zdim):
            self.face[2].coil[number-xdim*ydim*2].set_current(current)
        elif(number<xdim*ydim*2+xdim*zdim*2):
            self.face[3].coil[number-xdim*ydim*2-xdim*zdim].set_current(current)
        elif(number<xdim*ydim*2+xdim*zdim*2+ydim*zdim):
            self.face[4].coil[number-xdim*ydim*2-xdim*zdim*2].set_current(current)
        else:
            self.face[5].coil[number-xdim*ydim*2-xdim*zdim*2-ydim*zdim].set_current(current)

    def set_currents(self,vec_i):
        # set the currents to the vector given as the argument
        for i in range(self.numcoils):
            self.set_independent_current(i,vec_i[i])
            
    def zero_currents(self):
        # set all currents to zero
        for i in range(self.numcoils):
            self.set_independent_current(i,0.0)
            
    def set_common_current(self,curr):
        # set all currents to curr
        for i in range(self.numcoils):
            self.set_independent_current(i,curr)

    def draw_coil(self,number,ax):
        coil = self.coil(number)
        points = coil.points + (coil.points[0],)
        x = ([p[0] for p in points])
        y = ([p[1] for p in points])
        z = ([p[2] for p in points])
        ax.plot(x,y,z,label='coil')
        #ax.plot(x,y,z,label='coil')

    def draw_coils(self,ax):
        for number in range(self.numcoils):
            self.draw_coil(number,ax)

    def b(self,r):
        b_total = 0.0
        for number in range(self.numcoils):
            b_total = b_total + self.coil(number).b(r)
        return b_total

    def b_prime(self,x,y,z):
        b_total_x=0.*x
        b_total_y=0.*y
        b_total_z=0.*z
        for coilnum in range(self.numcoils):
            b_coil_x,b_coil_y,b_coil_z=self.coil(coilnum).b_prime(x,y,z)
            b_total_x=b_total_x+b_coil_x
            b_total_y=b_total_y+b_coil_y
            b_total_z=b_total_z+b_coil_z
        return b_total_x,b_total_y,b_total_z

class face:
    def __init__(self,xdim,ydim,corners):
        self.xdim = xdim
        self.ydim = ydim
        self.corners = corners
        x = corners[1]-corners[0]
        xstep = x/xdim
        y = corners[2]-corners[0]
        ystep = y/ydim
        coilnum = 0
        self.coil = []
        for i in range(xdim):
            for j in range(ydim):
                p0 = corners[0]+xstep*i+ystep*j
                p1 = p0+xstep
                p2 = p1+ystep
                p3 = p2-xstep
                points = (p0,p1,p2,p3)
                self.coil.append(coil(points,0.0))
                coilnum = coilnum + 1
        self.coilnum = coilnum

        
class sensor:
    def __init__(self,pos):
        self.pos = pos

class sensorarray:
    def __init__(self,xdim,ydim,zdim,corners):
        x = corners[1]-corners[0]
        y = corners[2]-corners[0]
        z = corners[3]-corners[0]
        #self.sensorgrid=np.mgrid[-a:a:xdim*1j,-a:a:ydim*1j,-a:a:zdim*1j]
        #print(self.sensorgrid)
        self.sensors = []
        if(xdim==1 and ydim==1 and zdim==1):
            pos = corners[0]+x/2+y/2
            self.sensors.append(sensor(pos))
            pos = corners[0]+x/2+y/2+z
            self.sensors.append(sensor(pos))
            pos = corners[0]+y/2+z/2
            self.sensors.append(sensor(pos))
            pos = corners[0]+y/2+z/2+x
            self.sensors.append(sensor(pos))
            pos = corners[0]+x/2+z/2
            self.sensors.append(sensor(pos))
            pos = corners[0]+x/2+z/2+y
            self.sensors.append(sensor(pos))
        else:
            for i in range(xdim):
                for j in range(ydim):
                    for k in range(zdim):
                        pos = corners[0]+x*i/(xdim-1)+y*j/(ydim-1)+z*k/(zdim-1)
                        self.sensors.append(sensor(pos))
        self.numsensors = len(self.sensors)
    def draw_sensor(self,number,ax):
        x = self.sensors[number].pos[0]
        y = self.sensors[number].pos[1]
        z = self.sensors[number].pos[2]
        c = 'r'
        m = 'o'
        ax.scatter(x,y,z,c=c,marker=m)
    def draw_sensors(self,ax):
        for number in range(self.numsensors):
            self.draw_sensor(number,ax)
    def vec_b(self, sp):
        # makes a vector of magnetic fields in the same ordering as
        # the_matrix class below
        the_vector=np.zeros((self.numsensors*3))
        for j in range(self.numsensors): #editted by TH June 8, 2020
            r = self.sensors[j].pos  #editted by TH June 8, 2020
            b=np.array([sp.fPix(r[0],r[1],r[2]),
                        sp.fPiy(r[0],r[1],r[2]),
                        sp.fPiz(r[0],r[1],r[2])])
            for k in range(3):
                the_vector[j*3+k]=b[k]
        return the_vector


class the_matrix:
    def __init__(self,mycube,myarray):
        self.m=np.zeros((mycube.numcoils,myarray.numsensors*3))
        #self.fill(mycube,myarray)
        self.fillspeed(mycube,myarray)
        self.condition = np.linalg.cond(self.m)

        # for some reason I chose to create the transpose of the usual
        # convention, when I first wrote the fill method
        self.capital_M=self.m.T # M=s*c=sensors*coils Matrix

        # Do the svd
        self.U,self.s,self.VT=np.linalg.svd(self.capital_M)

        # s is just a list of the diagonal elements, rather than a matrix
        # You can make the matrix this way:
        self.S=np.zeros(self.capital_M.shape)
        self.S[:self.capital_M.shape[1],:self.capital_M.shape[1]]=np.diag(self.s)
        # Or, I've seen people use "full_matrices=True" in the svd command

        # Start to calculate the inverse explicitly
        # list of reciprocals
        d=1./self.s
        self.D=np.zeros(self.capital_M.shape)
        # matrix of reciprocals
        self.D[:self.capital_M.shape[1],:self.capital_M.shape[1]]=np.diag(d)

        # inverse of capital_M
        self.Minv=self.VT.T.dot(self.D.T).dot(self.U.T)
        #self.Minv=np.linalg.pinv(self.capital_M)
        
        # now gets to fixin'
        # remove just the last mode
        n_elements=mycube.numcoils-1

        self.Dp=self.D[:,:n_elements]
        self.VTp=self.VT[:n_elements,:]
        self.Minvp=self.VTp.T.dot(self.Dp.T).dot(self.U.T)
        
    def fill(self,mycube,myarray):
        for i in range(mycube.numcoils):
            mycube.set_independent_current(i,1.0)
            for j in range(myarray.numsensors):
                r = myarray.sensors[j].pos
                b = mycube.b(r)
                for k in range(3):
                    self.m[i,j*3+k]=b[k]
            mycube.set_independent_current(i,0.0)

    def fillspeed(self,mycube,myarray):
        mycube.set_common_current(1.0)
        for i in range(mycube.numcoils):
            print("Doing coil %d"%i)
            for j in range(myarray.numsensors):
                r = myarray.sensors[j].pos
                bx,by,bz=mycube.coil(i).b_prime(r[0],r[1],r[2])
                b=[bx,by,bz]
                for k in range(3):
                    self.m[i,j*3+k]=b[k]
        mycube.zero_currents()
            
    def check_field_graphically(self,mycube,myarray):
        # test each coil by graphing field at each sensor
        for i in range(mycube.numcoils):
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            mycube.draw_coil(i,ax)
            mycube.coil(i).set_current(1.0)
            for j in range(myarray.numsensors):
                r = myarray.sensors[j].pos
                b=mycube.b(r)
                bhat=b*5.e4
                points = []
                points.append(r)
                points.append(r+bhat)
                xs = ([p[0] for p in points])
                ys = ([p[1] for p in points])
                zs = ([p[2] for p in points])
                ax.plot(xs,ys,zs)
            mycube.coil(i).set_current(0.0)
            ax.legend()
            plt.show()

    def show_matrices(self):
        fig1,ax1=plt.subplots()
        fig2,ax2=plt.subplots()
        fig3,ax3=plt.subplots()
        fig4,ax4=plt.subplots()
        fig5,ax5=plt.subplots()
        fig6,ax6=plt.subplots()
        fig7,ax7=plt.subplots()
        fig8,ax8=plt.subplots()
        fig9,ax9=plt.subplots()

        ax1.imshow(self.capital_M,cmap=cm.bwr)
        ax2.imshow(self.U,cmap=cm.bwr)
        ax3.imshow(self.S,cmap=cm.bwr)
        ax4.imshow(self.VT,cmap=cm.bwr)
        ax5.imshow(self.D,cmap=cm.bwr)
        ax6.imshow(self.Minv,cmap=cm.bwr)
        ax7.imshow(self.Dp,cmap=cm.bwr)
        ax8.imshow(self.VTp,cmap=cm.bwr)
        ax9.imshow(self.Minvp,cmap=cm.bwr)
        '''
        ax1.set_xticks(np.arange(self.m.shape[1]))
        ax1.set_yticks(np.arange(self.m.shape[0]))
        ax1.set_xlabel('Fluxgate positions')
        ax1.set_ylabel('Coils')
        ax1.set_title('Matrix M* (nT/A) ('+str(len(np.arange(self.m.shape[0])))+'coils * '+str(len(np.arange(self.m.shape[1])))+'sensors)')

        ax2.set_xticks(np.arange(self.capital_M1.shape[1]))
        ax2.set_yticks(np.arange(self.capital_M1.shape[0]))
        ax2.set_xlabel('Fluxgate positions')
        ax2.set_ylabel('Coils')
        #ax2.set_yticklabels(row_labels)
        #ax2.set_title('Pseudoinverse of M (*'+str(mps)+') (A/nT) ('+str(len(np.arange(self.m.shape[0])))+'coils * '+str(len(np.arange(self.m.shape[1])))+'sensors)')
        ax2.set_title('Pseudoinverse of M ('+str(len(np.arange(self.m.shape[0])))+'coils * '+str(len(np.arange(self.m.shape[1])))+'sensors)')

        ax3.set_xticks(np.arange(self.Vmat.shape[1]))
        ax3.set_yticks(np.arange(self.Vmat.shape[0]))
        ax3.set_title('V-Sqrt of eigenvalues of M*M & MM* ('+str(len(np.arange(self.m.shape[1])))+'* '+str(len(np.arange(self.m.shape[0])))+')')

        ax4.set_xticks(np.arange(self.U.shape[1]))
        ax4.set_yticks(np.arange(self.U.shape[0]))
        #ax4.set_title('U-Orthonormal eigenvectors(*'+str(mp4s)+') of MM* ('+str(len(np.arange(self.m.shape[1])))+'*'+str(len(np.arange(self.m.shape[1])))+')')
        ax4.set_title('U-Orthonormal eigenvectors of MM* ('+str(len(np.arange(self.m.shape[1])))+'*'+str(len(np.arange(self.m.shape[1])))+')')
        
        ax6.set_xticks(np.arange(self.Wt.shape[1]))
        ax6.set_yticks(np.arange(self.Wt.shape[0]))
        #ax6.set_title('W*-Orthonormal eigenvectors(*'+str(mp6s)+') of M*M ('+str(len(np.arange(self.m.shape[0])))+'*'+str(len(np.arange(self.m.shape[0])))+')')
        ax6.set_title('W*-Orthonormal eigenvectors of M*M ('+str(len(np.arange(self.m.shape[0])))+'*'+str(len(np.arange(self.m.shape[0])))+')')

        #for c in range (0,len(np.arange(self.m.shape[0]))):
	#    for s in range (0,len(np.arange(self.m.shape[1]))):
	#	ax1.text(s, c, int(self.m[c][s]), va='center', ha='center', rotation=90)
	#	ax2.text(s, c, int(self.capital_M1[c][s]*mp), va='center', ha='center', rotation=90)
	#	ax3.text(c, s, int(self.Vmat[s][c]), va='center', ha='center')

        #for c in range (0,len(np.arange(self.U.shape[0]))):
	#    for s in range (0,len(np.arange(self.U.shape[1]))):
	#	ax4.text(s, c, int(self.U[c][s]*mp4), va='center', ha='center')

        #for c in range (0,len(np.arange(self.Wt.shape[0]))):
	#    for s in range (0,len(np.arange(self.Wt.shape[1]))):
	#	ax6.text(s, c, int(self.Wt[c][s]*mp6), va='center', ha='center')
        '''

        plt.show()

        

def fiteven(x,p0,p2,p4,p6):
    return p0+p2*x**2+p4*x**4+p6*x**6

def fitodd(x,p1,p3,p5,p7):
    return p1*x+p3*x**3+p5*x**5+p7*x**7

def fitgraph(xdata,ydata,ax):
    popt,pcov=curve_fit(fiteven,xdata[abs(xdata)<.5],ydata[abs(xdata)<.5])
    print(popt)
    ax.plot(points1d,fiteven(xdata,*popt),'r--',label='$p_0$=%2.1e,$p_2$=%2.1e,$p_4$=%2.1e,$p_6$=%2.1e'%tuple(popt))