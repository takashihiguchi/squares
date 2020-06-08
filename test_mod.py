from squarespeed_mod import *

sp=scalarpotential(2,-3)
print("Sigma in spherical coordinates is %s"%sp.Sigma_spherical)
print("Sigma in cartesian coordinates is %s"%sp.Sigma)

print("Pix is %s"%sp.Pix)
print("Piy is %s"%sp.Piy)
print("Piz is %s"%sp.Piz)

a = 1.0
p0 = np.array([-a/2,-a/2,-a/2])
p1 = p0 + np.array([a,0,0])
p2 = p0 + np.array([0,a,0])
p3 = p0 + np.array([0,0,a])
points = (p0,p1,p2,p3)
#mycube=coilcube(5,5,5,points)
mycube=coilcube(3,3,3,points)


# test of sensorarray class
a = 0.5
p0 = np.array([-a/2,-a/2,-a/2])
p1 = p0 + np.array([a,0,0])
p2 = p0 + np.array([0,a,0])
p3 = p0 + np.array([0,0,a])
points = (p0,p1,p2,p3)
myarray = sensorarray(3,3,3,points)
print(myarray.sensors[2].pos)
print(myarray.numsensors)
print(myarray.sensors[myarray.numsensors-1].pos)
print(myarray.sensors[myarray.numsensors-2].pos)

print('the vector test')
print(myarray.vec_b(sp))
print(len(myarray.vec_b(sp)))


mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')
mycube.draw_coils(ax)
myarray.draw_sensors(ax)
#ax.legend()
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('z (m)')
plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')
mycube.draw_coils(ax)
myarray.draw_sensors(ax)
#ax.legend()
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('z (m)')
plt.show()


print(mycube.b(myarray.sensors[0].pos))
mycube.coil(0).set_current(1.0)
print(mycube.b(myarray.sensors[0].pos))
mycube.coil(0).set_current(0.0)


mymatrix = the_matrix(mycube,myarray)

print(mymatrix.condition)
mymatrix.show_matrices()

# Set up vector of desired fields

print(len(myarray.vec_b(sp)),myarray.vec_b(sp))
vec_i=mymatrix.Minvp.dot(myarray.vec_b(sp))
print(vec_i)

# Assign currents to coilcube

mycube.set_currents(vec_i)

# Check the field at the center of the coilcube
r=np.array([0,0,0])
print(mycube.b(r))
print(mycube.b_prime(0,0,0))

from scipy.optimize import curve_fit

fig7,(ax71)=plt.subplots(nrows=1)
fig8,(ax81)=plt.subplots(nrows=1)
fig9,(ax91)=plt.subplots(nrows=1)

# scans along each axis
points1d=np.mgrid[-1:1:101j]
bx1d_xscan,by1d_xscan,bz1d_xscan=mycube.b_prime(points1d,0.,0.)
bx1d_yscan,by1d_yscan,bz1d_yscan=mycube.b_prime(0.,points1d,0.)
bx1d_zscan,by1d_zscan,bz1d_zscan=mycube.b_prime(0.,0.,points1d)

# target field
bx1d_target_xscan=sp.fPix(points1d,0.,0.)*np.ones(np.shape(points1d))
bx1d_target_yscan=sp.fPix(0.,points1d,0.)*np.ones(np.shape(points1d))
bx1d_target_zscan=sp.fPix(0.,0.,points1d)*np.ones(np.shape(points1d))

by1d_target_xscan=sp.fPiy(points1d,0.,0.)*np.ones(np.shape(points1d))
by1d_target_yscan=sp.fPiy(0.,points1d,0.)*np.ones(np.shape(points1d))
by1d_target_zscan=sp.fPiy(0.,0.,points1d)*np.ones(np.shape(points1d))

bz1d_target_xscan=sp.fPiz(points1d,0.,0.)*np.ones(np.shape(points1d))
bz1d_target_yscan=sp.fPiz(0.,points1d,0.)*np.ones(np.shape(points1d))
bz1d_target_zscan=sp.fPiz(0.,0.,points1d)*np.ones(np.shape(points1d))

ax71.plot(points1d,bz1d_xscan,label='$B_z(x,0,0)$')
ax71.plot(points1d,bz1d_target_xscan,label='target $B_z(x,0,0)$')
ax71.plot(points1d,bz1d_yscan,label='$B_z(0,y,0)$')
ax71.plot(points1d,bz1d_target_yscan,label='target $B_z(0,y,0)$')
ax71.plot(points1d,bz1d_zscan,label='$B_z(0,0,z)$')
ax71.plot(points1d,bz1d_target_zscan,label='target $B_z(0,0,z)$')

ax81.plot(points1d,by1d_xscan,label='$B_y(x,0,0)$')
ax81.plot(points1d,by1d_target_xscan,label='target $B_y(x,0,0)$')
ax81.plot(points1d,by1d_yscan,label='$B_y(0,y,0)$')
ax81.plot(points1d,by1d_target_yscan,label='target $B_y(0,y,0)$')
ax81.plot(points1d,by1d_zscan,label='$B_y(0,0,z)$')
ax81.plot(points1d,by1d_target_zscan,label='target $B_y(0,0,z)$')

ax91.plot(points1d,bx1d_xscan,label='$B_x(x,0,0)$')
ax91.plot(points1d,bx1d_target_xscan,label='target $B_x(x,0,0)$')
ax91.plot(points1d,bx1d_yscan,label='$B_x(0,y,0)$')
ax91.plot(points1d,bx1d_target_yscan,label='target $B_x(0,y,0)$')
ax91.plot(points1d,bx1d_zscan,label='$B_x(0,0,z)$')
ax91.plot(points1d,bx1d_target_zscan,label='target $B_x(0,0,z)$')

min_field=-2.
max_field=+2.
#ax71.axis((-.5,.5,min_field,max_field))
ax71.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax71.legend()
ax81.legend()
ax91.legend()

plt.show()
