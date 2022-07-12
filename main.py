import random
import numpy as np
import matplotlib.pyplot
import os

## ---------------
def distance(x,y,xx,yy):
    return np.sqrt((x-xx)**2+(y-yy)**2)

## ---------------
class simulation_env():
    def __init__(self, phi, dt, w, b, l, vmax, vmin, tMax, maxAcceleration):
        self.flux = phi # pedestrian flux [ped/s]
        self.dt = dt # timestep [s]
        self.pedestriansAB = [] # list of pedestrians walking from A to B
        self.pedestriansCD = [] # list of pedestrians walking from C to D
        self.w=w # channels width [m]
        self.b=b # body radius of any pedestrian [m]
        self.l=l # channels length [m]
        self.vmax=vmax # vmax [m/s]
        self.vmin=vmin # vmin [m/s]
        self.tMax=tMax
        self.nPedestrian=0
        self.a=maxAcceleration #m/s^2
        self.outcsv=open("phi_"+str(self.flux)+"_dt_"+str(self.dt)+"_tmax_"+str(self.tMax)+"/"+"pedestrian_"+str(phi)+".csv", "w", buffering=1)
        # Initialize output file
        self.outcsv.write("Path,Individual")
        for i in range(0, int(tMax/dt)):
            self.outcsv.write(","+f'{i*dt:.{4}f}')
        self.outcsv.write("\n")

    def addPedestrian(self, inlet, t):
        new_rand=True
        while new_rand:
            ch_pos=(random.random()-0.5)*(self.w-2*self.b) # lateral position in channel
            if inlet=="A":
                new_ped_pos=(-self.l/2, ch_pos)
                pedestrians=self.pedestriansAB
                theta=0
            if inlet=="C":
                new_ped_pos=(ch_pos, -self.l/2)
                pedestrians=self.pedestriansCD
                theta=np.pi/2
            new_rand=False
            for ped_pos in pedestrians: # Check overlapping of new pedestrian with existing ones
                dist=(ped_pos.x-new_ped_pos[0])**2+(ped_pos.y-new_ped_pos[1])**2
                if dist<self.b**2:
                    new_rand=True
        vel=self.vmin+random.random()*(self.vmax-self.vmin)         #Define velocity of the new pedestrian
        pedestrians.append(
            pedestrian(new_ped_pos[0], new_ped_pos[1], vel, t, theta)
            )
        self.nPedestrian=self.nPedestrian+1
        return

    def removePedestrians(self):
        for i in self.pedestriansAB:
            if i.x>l/2:
                i.writeTraj(self.outcsv, "AB", self.dt, self.tMax)
                self.pedestriansAB.remove(i)
                self.nPedestrian=self.nPedestrian-1
        for i in self.pedestriansCD:
            if i.y>l/2:
                i.writeTraj(self.outcsv, "CD", self.dt, self.tMax)
                self.pedestriansCD.remove(i)
                self.nPedestrian=self.nPedestrian-1
        return
    
    def checkPedestrianCollision(self, pedestrianList, x_cand, y_cand, pedestrian):
        obstacle=False
        for ped_i in pedestrianList:
            if (not ped_i==pedestrian) and distance(x_cand, y_cand, ped_i.x, ped_i.y)<=2*self.b:
                obstacle=True
        return obstacle

    def checkWalls(self, pedestrian, x_cand, y_cand, x_wall, y_wall, wall):
        if wall=="y":
            return abs(y_cand)>=y_wall-self.b, abs(pedestrian.y)<=abs(y_cand)
        if wall=="x":
            return abs(x_cand)>=x_wall-self.b, abs(pedestrian.x)<=abs(x_cand)

    def choosePath(self, pedestrian, pedList, ped_in_cross, alphamin, alphamax, rotation):
        print("correcting traj of a pedestrian :" ,end='')
        nOptions1=50
        nOptions2=360
        alphaOptions=np.linspace(alphamin, alphamax, nOptions1)
        if pedestrian.alpha*pedestrian.v<=0.5:
            thetaOptions=np.linspace(0, np.pi*3/2, nOptions2)
            sortedOptions=np.dstack(
                np.unravel_index( np.argsort(#np.abs(
                np.outer(np.cos(thetaOptions), alphaOptions).ravel()#)
                ), (nOptions2, nOptions1)))[0]
        else:
            thetaOptions=np.linspace(0, 40*np.pi/180, nOptions2)
            sortedOptions=np.dstack(
                np.unravel_index( np.argsort(
                np.outer(np.cos(thetaOptions), alphaOptions).ravel()),
                (nOptions2, nOptions1)))[0]
        for iii in range(1, len(sortedOptions)+1):
            i, j = sortedOptions[-iii][0], sortedOptions[-iii][1] 
            thetaTest=thetaOptions[i]
            alpha=float(alphaOptions[j])

            obstacleVect=[]
            closerVect=[]
            for theta in [thetaTest, -thetaTest]:
                x_cand, y_cand = pedestrian.candSpost(alpha, theta+rotation, self.dt)
                obstacle = self.checkPedestrianCollision(pedList, x_cand, y_cand, pedestrian)
                if x_cand>=-self.w/2-self.b and x_cand<=self.w/2+self.b and y_cand>=-self.w/2-self.b and y_cand<=self.w/2+self.b:
                    obstacle_x = self.checkPedestrianCollision(ped_in_cross, x_cand, y_cand, pedestrian)
                    obstacle=(obstacle or obstacle_x)
                if rotation==np.pi/2:
                    wall="x"
                else:
                    wall="y"
                obstacle_w, closer=self.checkWalls(pedestrian, x_cand, y_cand, self.w/2, self.w/2, wall)
                obstacle=(obstacle or obstacle_w)
                obstacleVect.append(obstacle)
                closerVect.append(closer)

            if False in obstacleVect:
                obstacle=False
                if not closerVect[1] and not obstacleVect[1]:
                    theta=-thetaTest
                elif not closerVect[0] and not obstacleVect[0]:
                    theta=thetaTest
                elif not obstacleVect[1]:
                    theta=-thetaTest
                elif not obstacleVect[0]:
                    theta=thetaTest

            if not obstacle:
                pedestrian.progress(self.dt, alpha, theta+rotation)
                print(alpha, theta*180/np.pi)
                return
        raise StopIteration

    def movePedestrians(self):
        for listPed1, listPedX, rotation in zip([self.pedestriansCD, self.pedestriansAB], [self.pedestriansAB, self.pedestriansCD], [np.pi/2, 0]):
            ped_in_cross=[]
            for ped in listPedX:
                if ped.x>-self.w/2-self.b and ped.x<=self.w/2+self.b and ped.y>=-self.w/2-self.b and ped.y<=self.w/2+self.b:
                    ped_in_cross.append(ped)
            for pedestrian in listPed1:
                alpha_candidate=min(1, pedestrian.alpha+self.a*self.dt/pedestrian.v)
                x_cand, y_cand=pedestrian.candSpost(alpha_candidate, rotation, self.dt)
                obstacle = self.checkPedestrianCollision(listPed1, x_cand, y_cand, pedestrian)
                if x_cand>=-self.w/2-self.b and x_cand<=self.w/2+self.b and y_cand>=-self.w/2-self.b and y_cand<=self.w/2+self.b:
                    obstacle_x = self.checkPedestrianCollision(ped_in_cross, x_cand, y_cand, pedestrian)
                    obstacle=obstacle or obstacle_x
                if not obstacle:
                    pedestrian.progress(self.dt, alpha_candidate, rotation)
                else:
                    alphamin=0#max(0,pedestrian.alpha-self.a*self.dt/pedestrian.v)
                    alphamax=min(1, pedestrian.alpha+self.a*self.dt/pedestrian.v)
                    self.choosePath(pedestrian, listPed1, ped_in_cross, alphamin, alphamax, rotation)

    def getLists(self):
        xh=[]
        yh=[]
        vh=[]
        for i in self.pedestriansAB:
          xtemp, ytemp, vtemp=[], [], []
          for j in i.history:
            xtemp.append(j[0])
            ytemp.append(j[1])
            vtemp.append(j[2]/self.vmax)
          xh.append(xtemp)
          yh.append(ytemp)
          vh.append(vtemp)
        for i in self.pedestriansCD:
          xtemp, ytemp, vtemp=[], [], []
          for j in i.history:
            xtemp.append(j[0])
            ytemp.append(j[1])
            vtemp.append(j[2]/self.vmax)
          xh.append(xtemp)
          yh.append(ytemp)
          vh.append(vtemp)
        return [i.x for i in self.pedestriansAB]+[i.x for i in self.pedestriansCD], [i.y for i in self.pedestriansAB]+[i.y for i in self.pedestriansCD], [i.alpha*i.v for i in self.pedestriansAB]+[i.alpha*i.v for i in self.pedestriansCD], xh, yh, vh

## ---------------       
class pedestrian():
    def __init__(self, x, y, v, tin, theta):
        self.x = x
        self.y = y
        self.v = v
        self.alpha = 1
        self.history = []
        self.tin = tin
        self.theta=theta

    def candSpost(self, alpha_candidate, theta_candidate, dt):
        return self.x+alpha_candidate*self.v*dt*np.cos(theta_candidate), self.y+alpha_candidate*self.v*dt*np.sin(theta_candidate)

    def progress(self, dt, alpha, theta):
        self.alpha=alpha
        self.theta=theta
        self.x=self.x+self.v*np.cos(theta)*dt*self.alpha
        self.y=self.y+self.v*np.sin(theta)*dt*self.alpha
        self.history.append((self.x, self.y, self.v*self.alpha))
        return
    
    def writeTraj(self, filename, path, dt, tMax):
        filename.write(path)
        decprecision=5
        tend=len(self.history)*dt+self.tin
        for i in range(0, int(self.tin/dt)):
            filename.write(",(NaN NaN NaN)")
        for i in range(0, int((tend-self.tin)/dt)):
            filename.write(
                ",("+ f'{self.history[i][0]:.{decprecision}f}'+
                " "+f'{self.history[i][1]:.{decprecision}f}'+
                " "+f'{self.history[i][2]:.{decprecision}f}'+")")
        for t in range(int(tend/dt), int(tMax/dt)):
            filename.write(",(NaN NaN NaN)")
        filename.write("\n")
        return

## ---------------
## Initialize
w=5 
l=25 
b=0.25
vmin=1.1
vmax=1.3
dt=0.1
tMax=60
phi=2
maxAcceleration=2 #m/s^2

try:
    os.mkdir("phi_"+str(phi)+"_dt_"+str(dt)+"_tmax_"+str(tMax))
except:
    pass

sim=simulation_env(phi, dt, w, b, l, vmax, vmin, tMax, maxAcceleration)

## ---------------
# First pedestrians
t=0
tAdd=0
sim.addPedestrian("A", t)
sim.addPedestrian("C", t)

## ---------------
## Solution Loops
isave=1
while t<=tMax:
    # Add pedestrians
    if t-tAdd >= 1/phi:
        sim.addPedestrian("A", t)
        sim.addPedestrian("C", t)
        tAdd=t
    print("Time", f'{t:.{3}f}', "\tNumber of pedestrians:",  sim.nPedestrian)
    # Move forward
    sim.movePedestrians()
    # Rempve pedestrian outside domain
    sim.removePedestrians()
    t=t+dt

    # Plot timestep
    xl=[]
    yl=[]
    xl, yl, colors, xlh, ylh, colorsh=sim.getLists()
    matplotlib.pyplot.plot([-l/2, -w/2, -w/2], [-w/2, -w/2, -l/2], 'k')
    matplotlib.pyplot.plot([w/2, w/2, l/2], [-l/2, -w/2, -w/2], 'k')
    matplotlib.pyplot.plot([-l/2, -w/2, -w/2], [w/2, w/2, l/2], 'k')
    matplotlib.pyplot.plot([w/2, w/2, l/2], [l/2, w/2, w/2], 'k')
    lines=[matplotlib.pyplot.plot(xlh[i], ylh[i], lw=0.2, c='k') for i in range(len(xlh))] #color=matplotlib.pyplot.cm.jet(colorsh[i])
    matplotlib.pyplot.scatter(xl, yl, s=14, c=colors, vmin=0, vmax=sim.vmax , zorder=2)
    matplotlib.pyplot.text(-sim.l/2*1.3, -sim.l/2*1.3, "t="+f'{t:.{3}f}'+" s", fontsize=12)
    matplotlib.pyplot.colorbar(label='m/s')
    matplotlib.pyplot.savefig("phi_"+str(sim.flux)+"_dt_"+str(sim.dt)+"_tmax_"+str(sim.tMax)+"/"+"{:03d}".format(isave)+".png", dpi=500)
    isave=isave+1
    matplotlib.pyplot.clf()

# Export remaining pedestrian
for i in  sim.pedestriansAB:
    i.writeTraj(sim.outcsv, "AB", sim.dt, sim.tMax)
for i in  sim.pedestriansCD:
    i.writeTraj(sim.outcsv, "CD", sim.dt, sim.tMax)
