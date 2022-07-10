import random
import math
from gekko import GEKKO
import matplotlib.pyplot

def distance(x,y,xx,yy):
    return math.sqrt((x-xx)**2+(y-yy)**2)
def tangents(x1, y1, x2, y2, r):
    dx, dy = x1-x2, y1-y2
    dxr, dyr = -dy, dx
    d = math.sqrt(dx**2+dy**2)
    if d >= r :
        rho = r/d
        ad = rho**2
        bd = rho*math.sqrt(1-rho**2)
        T1x = x2 + ad*dx + bd*dxr
        T1y = y2 + ad*dy + bd*dyr
        T2x = x2 + ad*dx - bd*dxr
        T2y = y2 + ad*dy - bd*dyr
        if (d/r-1) < 1E-8:
            print('Point on circumference\n')
            return False, 0, 0
        else:
            return True, math.atan2(T1y-y1, T1x-x1), math.atan2(T2y-y1, T2x-x1)
    else:
        print('Point P is inside the circle. No tangent is possible...\n')


class simulation_env():
    def __init__(self, phi, dt, w, b, l, vmax, vmin, tMax, deltaThetaMax, maxAcceleration):
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
        self.dtheta=deltaThetaMax*math.pi/180/2 #rad
        self.a=maxAcceleration #m/s^2
        self.outcsv=open("pedestrian_"+str(phi)+".csv", "w", buffering=1)
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
                theta=math.pi/2
            new_rand=False
            for ped_pos in pedestrians: # Check overlapping of new pedestrian with existing ones
                dist=(ped_pos.x-new_ped_pos[0])**2+(ped_pos.y-new_ped_pos[1])**2
                if dist<self.b**2:
                    new_rand=True
        vel=self.vmin+random.random()*(self.vmax-self.vmin)         #Define velocity of the new pedestrian
        pedestrians.append(
            pedestrian(new_ped_pos[0], new_ped_pos[1], vel, theta, t)
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
    
    def movePedestrians(self):
        for pedestrian in self.pedestriansCD:
            print("x", end='', flush=True)
            #check wall
            dx=pedestrian.v*self.dt*math.sin(self.dtheta)*(pedestrian.alpha+self.a*self.dt/pedestrian.v)
            rightwall=(dx+pedestrian.x>=-self.b+self.w/2)
            leftwall=(-dx+pedestrian.x<=self.b-self.w/2)

            ped_in_candidatespace=[]
            #check other pedestrians in CD
            for ped_i in self.pedestriansCD:
                if (not ped_i==pedestrian) and distance(ped_i.x, ped_i.y, pedestrian.x, pedestrian.y)-2*self.b>=pedestrian.v*self.dt:
                    tangentsBool, th1, th2 = tangents(pedestrian.x, pedestrian.y, ped_i.x, ped_i.y, 2*self.b)
                    th1=th1-math.pi/2
                    th2=th2-math.pi/2
                    if abs(th1)<self.dtheta/2 and abs(th2)<self.dtheta/2:
                        ped_in_candidatespace.append(ped_i)
            
            if pedestrian.y>=-w and pedestrian.y<=w:
                for ped_j in self.pedestriansAB:
                    if ped_j.y>=-w and ped_j.y<=w:
                        ped_in_candidatespace.append(ped_j)
            
            if ped_in_candidatespace or rightwall or leftwall:
                m = GEKKO()
                alpha_p1,theta = [m.Var(lb=max(0,pedestrian.alpha-self.a*self.dt/pedestrian.v), ub=min(1, pedestrian.alpha+self.a*self.dt/pedestrian.v)), m.Var(lb=-self.dtheta/2, ub=self.dtheta/2)]
                alpha_p1.value = pedestrian.alpha
                theta.value = 0
                #Equations
                if leftwall:
                    m.Equation(-alpha_p1*pedestrian.v*self.dt*m.sin(theta)+pedestrian.x-2*self.b>=-self.w/2)
                if rightwall:
                    m.Equation(alpha_p1*pedestrian.v*self.dt*m.sin(theta)+pedestrian.x+2*self.b<=self.w/2)

                for ped_i in ped_in_candidatespace:
                    m.Equation((pedestrian.x+alpha_p1*pedestrian.v*self.dt*m.sin(theta)-ped_i.x)**2+
                                (pedestrian.y+alpha_p1*pedestrian.v*self.dt*m.cos(theta)-ped_i.y)**2>=4*(self.b)**2)

                #Objective
                m.Minimize(m.abs2(theta)-alpha_p1)

                #Set global options
                m.options.IMODE = 3 #steady state optimization
                #Solve simulation
                m.solve(disp=False)
                pedestrian.progress(self.dt, alpha_p1.value[0], theta.value[0]+math.pi/2)
                m.cleanup()
            else:
                pedestrian.progress(self.dt, min(1, pedestrian.alpha+self.a*self.dt/pedestrian.v), math.pi/2)

        print("\n")

    def getLists(self):
        return [i.x for i in self.pedestriansCD], [i.y for i in self.pedestriansCD]
            
class pedestrian():
    def __init__(self, x, y, v, theta, tin):
        self.x = x
        self.y = y
        self.v = v
        self.theta = theta
        self.alpha = 1
        self.history = []
        self.tin = tin

    def progress(self, dt, alpha, theta):
        self.alpha=alpha
        self.theta=theta
        self.x=self.x+self.v*math.cos(self.theta)*dt*self.alpha
        self.y=self.y+self.v*math.sin(self.theta)*dt*self.alpha
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
deltaThetaMax=20 #deg
sim=simulation_env(phi, dt, w, b, l, vmax, vmin, tMax, deltaThetaMax, maxAcceleration)

# First pedestrian
t=0
tAdd=0
#sim.addPedestrian("A", t)
sim.addPedestrian("C", t)

## ---------------
## Solution Loops
while t<=tMax:
    if t-tAdd >= 1/phi:
    #    sim.addPedestrian("A", t)
        sim.addPedestrian("C", t)
        tAdd=t

    print("Time", t, "\t\t\tNumber of pedestrians:",  sim.nPedestrian)

    # Move forward
    sim.movePedestrians()

    sim.removePedestrians()
    t=t+dt

    xl=[]
    yl=[]
    xl, yl=sim.getLists()
    matplotlib.pyplot.plot([-w/2, -w/2], [-l/2, l/2], 'k')
    matplotlib.pyplot.plot([w/2, w/2], [-l/2, l/2], 'k')
    matplotlib.pyplot.scatter(xl, yl, s=0.5)
    matplotlib.pyplot.savefig(str(t)+".png")
    matplotlib.pyplot.clf()

for i in  sim.pedestriansAB:
    i.writeTraj(sim.outcsv, "AB", sim.dt, sim.tMax)
for i in  sim.pedestriansCD:
    i.writeTraj(sim.outcsv, "CD", sim.dt, sim.tMax)