import vaex
from yade import pack, plot
import math
import random as rand
import numpy as np
import argparse
from scipy import signal as sg

#PARAMETERS

#Dimensionless numbers
inertialNumber = 0.5   #Inertial number imposed
restitCoef = 0.9		#Restitution coefficient of the particles, dimensionless
partFrictAngle = atan(0.5)	#friction angle of the particles, in radian

#Geometry
length = 20                    #Streamwise length of the periodic cell, in diameter
width = 5                      #Spanwise length of the periodic cell, in diameter
Nlayer = 10.                   #nb of layer of particles, in diameter
twoDimension = False           #If activated, consider a two dimensional case
if twoDimension==True:
	width = 1              #Spanwise length of the periodic cell, in diameter

#Nature of the forcing
zeroGravity = False	       #If activated, consider a system without gravity
slopedeg = 10               #Inclination angle of the channel slope, in deg
slope=slopedeg*2.*pi/360.     #Inclination angle of the channel slope, in radian


#Deduced parameters of the configuration
diameterPart = 1.             #Diameter of the particles
massPart = 1.                 #Diameter of the particles
partVolume = pi/6*pow(diameterPart,3.)	#Volume of particles
densPart = massPart/partVolume		#density of the particles
phiPartMax = 0.61	       #Value of the dense packing solid volume fraction, dimensionless



#CONFIGURATION

#Forcing

nbPartBoundaryPlates = length*width						#Nb of particles composing the top and lower plates
pressureImposed = nbPartBoundaryPlates*partVolume*densPart*50		#Impose the pressure as a weight of the top plate
if zeroGravity:
	pressureImposed= nbPartBoundaryPlates*partVolume*densPart
densPartBoundary = pressureImposed/(nbPartBoundaryPlates*partVolume)		#Determine the particle density of the top and lower plate
shearRate_carac = inertialNumber/(diameterPart*sqrt(densPart/pressureImposed))	#Evaluate the shear rate to impose to the boundary in order to recover the imposed inertial number
T_oscil= 4 #Oscillating period
R_asym = 2
if zeroGravity==True:
	gravityVector = Vector3(0.,0.,0.)
else:
	gravityVector = Vector3(9.81*sin(slope),0.0,-9.81*cos(slope)) #Gravity vector to consider a channel inclined with slope angle 'slope'


#Geometrical configuration, define useful quantities

height = 5*Nlayer*diameterPart	#heigth of the periodic cell, in m
groundPosition = height/4.0 #Definition of the position of the ground, in m
if zeroGravity==True:
	groundPosition = 2*diameterPart
initialTopPlatePosition = height - 2*diameterPart

#Particles contact law/material parameters

maxPressure = max(pressureImposed,densPart*phiPartMax*Nlayer*diameterPart*abs(gravityVector[2]))#Estimated max particle pressure from the static load
normalStiffness = maxPressure*diameterPart*1e4 #Evaluate the minimal normal stiffness to be in the rigid particle limit (cf Roux and Combe 2002)
youngMod = normalStiffness/diameterPart	#Young modulus of the particles from the stiffness wanted.
poissonRatio = 0.5	#poisson's ratio of the particles. Classical values, does not have much influence
O.materials.append(ViscElMat(en=restitCoef, et=0., young=youngMod, poisson=poissonRatio, density=densPart, frictionAngle=partFrictAngle, label='Mat'))  

O.materials.append(ViscElMat(en=restitCoef, et=0., young=youngMod, poisson=poissonRatio, density=densPartBoundary, frictionAngle=partFrictAngle, label='MatBoundary'))  


#CREATING PARTICLES AND ROUGHS

#Definition of the semi-periodic cell
O.periodic = True
O.cell.setBox(length,width,height)
centerYPosition = 0.
if twoDimension==True:
	O.cell.setBox(length,width*4,height)
	centerYPosition = 1.5*width

#Create a loose cloud of particle inside the cell
partCloud = pack.SpherePack()
partVolume = pi/6.*pow(diameterPart,3) #Volume of a particle
partNumber = int(Nlayer*phiPartMax*diameterPart*length*width/partVolume) #Volume of beads to obtain Nlayer layers of particles
polyDispersity = 0.
if twoDimension==True:
	polyDispersity = 0.2
partCloud.makeCloud(minCorner=(0,centerYPosition + diameterPart/2.,groundPosition+diameterPart),maxCorner=(length,centerYPosition + width - diameterPart/2., initialTopPlatePosition-diameterPart),rRelFuzz=polyDispersity, rMean=diameterPart/2.0, num = partNumber)
partCloud.toSimulation(material='Mat') #Send this packing to simulation with material Mat

#Evaluate the deposition time considering the free-fall time of the highest particle to the ground
depoTime = 0.
if zeroGravity==False:
	depoTime = 2*sqrt(Nlayer*3*diameterPart*2/abs(gravityVector[2]))

if twoDimension ==True:
	for b in O.bodies:
		b.state.blockedDOFs = 'y'


#Add the lower rough plate as a clump (particles sticked together)

lowerWall = []
for x in range(0,int(length)): #loop over the x direction, by step of a diameter
	for y in range(0,int(width)):	#loop over the y direction, by step of a diameter
		polydisp = ([-1,1][np.random.randint(0,2)]*np.random.random()/10)+1
		lowerWall.append(sphere((x*diameterPart + diameterPart/2. , centerYPosition + y*diameterPart + diameterPart/2,groundPosition - diameterPart/2.0 ),diameterPart/2.*polydisp,color=(1,1,1),fixed = True,material = 'MatBoundary'))				#Create a list of sphere element at different random height along a square lattice
#Add the list of sphere to the simulation, lowerWallIDS keeps tracks of the particles ID in the simulations
lowerWallIDS = O.bodies.appendClumped(lowerWall)

number_particles = len(O.bodies)


#Add the top rough as a clump

topWall = []
for x in range(0,int(length)):
	for y in range(0,int(width)):
		polydisp = ([-1,1][np.random.randint(0,2)]*np.random.random()/10)+1
		topWall.append(sphere((x*diameterPart + diameterPart/2.,centerYPosition + y*diameterPart + diameterPart/2., initialTopPlatePosition - diameterPart/2.0),diameterPart/2.*polydisp,color=(1,1,1),material = 'MatBoundary'))
topWallIDS = O.bodies.appendClumped(topWall)
O.bodies[topWallIDS[0]].state.blockedDOFs = 'xyXYZ'


#Block some degrees of freedom of the lower plate created
if zeroGravity==True:	#If there is no gravity, the shear is applied to both the lower and the top plate.
	O.bodies[lowerWallIDS[0]].state.blockedDOFs = 'xyXYZ'
else:
	O.bodies[lowerWallIDS[0]].state.blockedDOFs = 'xyzXYZ'


plot.addData(n = 0., time = 0.)

#get initial positions

pos = []
for part in O.bodies:
    pos.append(part.state.pos[0])

#define list for flow

n_total,V_total,Q_total,Shear_total, time_total, time_step=[0],[0],[0],[0],[0],[0]


#SIMULATION LOOP

O.engines = [
	# Reset the forces
	ForceResetter(),
	# Detect the potential contacts
	InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Wall_Aabb(),Bo1_Facet_Aabb(),Bo1_Box_Aabb()],label='contactDetection',allowBiggerThanPeriod = True),
	# Calculate the different interactions
	InteractionLoop(
   	[Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()],
   	[Ip2_ViscElMat_ViscElMat_ViscElPhys()],
   	[Law2_ScGeom_ViscElPhys_Basic()]
	,label = 'interactionLoop'),
	# Application of the deformation
	PyRunner(command='oscilShearRate(shearRate_carac,T_oscil,R_asym)',iterPeriod = 1,label = 'oscil'),
	PyRunner(command='deformationApplication()',iterPeriod = 1,label = 'deform'),
	PyRunner(command='countingParticles()', iterPeriod = 1, label = 'count'),
	PyRunner(command='addPlot()', iterPeriod = 4, label = 'plot'),
	#GlobalStiffnessTimeStepper, determine the time step
	GlobalStiffnessTimeStepper(defaultDt = 1e-4, viscEl = True,timestepSafetyCoefficient = 0.7,  label = 'GSTS'),
	# Integrate the equation and calculate the new position/velocities...
	NewtonIntegrator(damping=0., gravity=gravityVector, label='newtonIntegr'),
	PyRunner(command='save_data()', iterPeriod = 10000, label = 'export'),
	PyRunner(command='stopRunning()', iterPeriod = 1, label = 'stop')
	]

plot.plots={'time':('n','V','Q')}
plot.plot()
#plot.setLiveForceAlwaysUpdate(True)

#save the initial configuration to be able to recharge the simulation starting configuration easily

O.saveTmp()


#FUNCTION DEFINITION

#counting number of particles

def countingParticles():
	global pos, n_total, V_total, Q_total, time_total, time_step
	time_total.append(O.time)
	time_step.append(O.dt)

	n_temp = 0
	pos_temp = []
	for i in range(number_particles):
		pos_temp.append(O.bodies[i].state.pos[0])
	if pos[i] % 20 > 15 and pos_temp[i] % 20 < 5 :
		n_temp+=1
	if pos[i] % 20 < 5 and pos_temp[i] % 20 > 15 :
		n_temp-=1
	n_total.append(n_temp)
	V_total.append(np.sum(n_total[:-int(np.round(1/O.dt))]))
	Q_total.append(np.sum(n_total[:-int(np.round(1/O.dt))])/O.time)
	pos = pos_temp

def addPlot():
    plot.addData(n = np.sum(n_total[-int(np.round(1/O.dt)):]), V=V_total[-1], Q=Q_total[-1],time = O.time)
    #plot.addData(n = np.sum(n_total[-10000:]), Q_i = np.sum(n_total[:-5000]), Q_moy = np.sum(n_total[:-5000])/O.time,time = O.time)

#Function to modify shear rate speed and direction with time

def oscilShearRate(shearR,period_oscil,R):
	global shearRate,orientation,Shear_total
	t=O.time%((1+R)*0.5*period_oscil)
	if t<period_oscil/2:
		shearRate=abs(R*shearR*sin(t*2.*pi/period_oscil))
		orientation=1
	else:
		shearRate=abs(shearR*sin(((t-period_oscil/2.)*2.*pi/(R*period_oscil))))*(-1.)
		orientation=-1
	Shear_total.append(orientation*shearRate)


#Function to apply the deformation of the shear cell to the boundaries = the top and lower plates

def deformationApplication():
	VImposed = diameterPart*shearRate
	if zeroGravity ==True:
		O.bodies[topWallIDS[0]].state.vel[0] = VImposed/2.
		O.bodies[lowerWallIDS[0]].state.vel[0] = -VImposed/2.
		O.forces.addF(topWallIDS[0],Vector3(0,0,-pressureImposed/2.))
		O.forces.addF(lowerWallIDS[0],Vector3(0,0,pressureImposed/2.))
	else:
		O.bodies[topWallIDS[0]].state.vel[0] = VImposed
		O.forces.addF(topWallIDS[0],Vector3(0,0,-pressureImposed))


#Function to modify the inertial number during the simulation

def updateInertialNumber(I):
	global inertialNumber, shearRate_carac
	inertialNumber = I
	shearRate_carac = inertialNumber/(diameterPart*sqrt(densPart/pressureImposed))

def stopRunning():
	time_Q = sg.find_peaks(Q_total)[0]
	time_V = sg.find_peaks(V_total)[0]
	if len(time_Q)>2:
		peaks_Q = Q_total[time_Q]
		peaks_V = V_total[time_V]
		if abs(peaks_Q[-1]-2*peaks_Q[-2]+peaks_Q[-3])<1e-3 and abs(peaks_V[-1]-2*peaks_V[-2]+peaks_V[-3])<1e-3:
			O.stop()


def save_data():
	df = vaex.from_arrays(T=time_total, dt=time_step, R=Shear_total, n=n_total, V=V_total, Q=Q_total)
	df.export('Exp_slope{}.hdf5'.format(slopedeg))
