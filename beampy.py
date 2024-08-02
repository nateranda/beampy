import numpy as np
import matplotlib.pyplot as plt

class Beam:
    "Beam & support attributes"
    def __init__(self, length, ei, cantilever=None, dl=None, dr=None, sections=None, rotDelta=None):

        self.cantilever = False if cantilever is None else cantilever    # cantilever or simply-suppoted
        self.length = length                                             # length, feet
        self.dl = 0 if dl is None else dl                                # left support location, feet
        self.dr = self.length if dr is None else dr                      # right support location, feet
        self.dist = self.dr-self.dl                                      # distance between supports, feet
        self.ei = ei                                                     # e*i, lb-in^2

        self.sections = 1000 if sections is None else sections           # number of integration sections
        self.rotDelta = 0.0001 if rotDelta is None else rotDelta         # how much to change rotation delta - multiplier of ei
        self.interval = np.linspace(0, self.length, num=self.sections+1) # interval array
        self.interval_width = self.length/self.sections                  # width of each interval slice
        self.support_index_l = np.abs(self.interval-self.dl).argmin()    # index of left support
        self.support_index_r = np.abs(self.interval-self.dr).argmin()    # index of right support

        self.shear = np.zeros_like(self.interval)                        # shear array
        self.moment = np.zeros_like(self.interval)                       # moment array
        self.deflection = np.zeros_like(self.interval)                   # deflection array
        self.point_loads = []                                            # point load list
        self.dist_loads = []                                             # distributed load list

    def correct(self):
        "Corrects invilad supports"
        if self.cantilever == True:
            self.dl = 0
            self.dr = self.length
        else:
            if self.dr > self.length:
                self.dr = self.length
            if self.dl < 0:
                self.dl = 0
        
        #Recalculate dependents
        self.dist = self.dr-self.dl
        self.support_index_l = np.abs(self.interval-self.dl).argmin()
        self.support_index_r = np.abs(self.interval-self.dr).argmin()
    
    def addLoad(self, load):
        if type(load) == PointLoad:
            self.point_loads.append(load)
        elif type(load) == DistLoad:
            self.dist_loads.append(load)
        else:
            raise TypeError("Invalid load type.")
    
    def calc_sm(self):
        pshear, pmoment = point_load_calc(self)
        dshear, dmoment = dist_load_calc(self)
        self.shear = np.add(pshear, dshear)
        self.moment = np.add(pmoment, dmoment)
    
    def calc_def(self):
        rot = get_rotation(self)
        self.deflection = get_deflection(self, rot)
    
    def plot_sm(self):
        "Plots shear/moment diagram"
        plt.title("Shear/Moment Diagram")
        plt.xlabel("Length (ft)")
        plt.ylabel("Stress (lb/ft-lb)")
        plt.axhline(0, color='black', linewidth=0.5)
        plt.plot(self.interval, self.shear)
        plt.plot(self.interval, self.moment)
        plt.show()
    
    def plot_def(self):
        "Plots deflection diagram"
        plt.title("Deflection Diagram")
        plt.xlabel("Length (ft)")
        plt.ylabel("Deflection (in)")
        plt.axhline(0, color='black', linewidth=0.5)
        plt.plot(self.interval, self.deflection)
        plt.show()
        
class PointLoad:
    "Point load object"
    def __init__(self, d, m, shear=None):
        self.d = d                                    # distance, ft
        self.m = m                                    # magnitude, lb
        self.shear = True if shear is None else shear # shear or moment load

class DistLoad:
    "Distributed load object"
    def __init__(self,dl, dr, ml, mr):
        self.dl = dl # start location
        self.dr = dr # end location
        self.ml = ml # start magnitude
        self.mr = mr # end magnitude
        
        self.len = self.dr-self.dl                                               # load span
        self.mag = (self.ml+self.mr)/2*self.len                                  # load magnitude
        self.pos = self.dl + (self.ml+2*self.mr)/(3*(self.ml+self.mr))*self.len  # centroid w.r.t beam end

def point_load_calc(beam):
    "Calculates shear & moment for point loads"
    shear = np.zeros_like(beam.interval)
    moment = np.zeros_like(beam.interval)

    for load in beam.point_loads:
        # Shear loads
        if load.shear:
            for i, x in enumerate(beam.interval):
                if x >= load.d:
                    shear[i] += load.m

            if beam.cantilever:
                shear = np.add(-load.m,shear)
                moment = np.add(load.d*load.m,moment)

            if not beam.cantilever:
                vr = load.m*(load.d-beam.dl)/beam.dist
                vl = load.m - vr
                for i, x in enumerate(beam.interval):
                    if x >= beam.dl:
                        shear[i] -= vl
                    if x >= beam.dr:
                        shear[i] -= vr
        
        # Moment loads
        if not load.shear:
            for i, x in enumerate (beam.interval):
                if x >= load.d:
                    moment[i] -= load.m

            if beam.cantilever:
                moment = np.add(load.m, moment)

            if not beam.cantilever:
                vr = load.m/beam.dist
                vl = -vr
                for i, x in enumerate(beam.interval):
                    if x >= beam.dl:
                        shear[i] -= vl
                    if x >= beam.dr:
                        shear[i] -= vr
    
    # Integrate moment
    moment_sum = 0
    for i in range(1, len(beam.interval)):
        moment_sum += (shear[i]+shear[i-1])/2*beam.interval_width
        moment[i] += moment_sum
    
    return shear, moment

def dist_load_calc(beam):
    "Calculates shear & moment for distributed loads"
    shear = np.zeros_like(beam.interval)
    moment = np.zeros_like(beam.interval)

    # Shear
    for load in beam.dist_loads:
        for i, x in enumerate(beam.interval):
            if x >= load.dl and x <= load.dr:
                shear[i] += load.ml*(x-load.dl)+(x-load.dl)**2*(load.mr-load.ml)/(2*load.len)
            if x > load.dr:
                shear[i] += load.mag

        if beam.cantilever:
            shear = np.add(-load.mag, shear)
            moment = np.add(load.mag*load.pos, moment)

        if not beam.cantilever:
            vr = load.mag*(load.pos-beam.dl)/beam.dist
            vl = load.mag-vr
            for i, x in enumerate(beam.interval):
                if x >= beam.dl:
                    shear[i] -= vl
                if x >= beam.dr:
                    shear[i] -= vr
    
    # Integrate moment
    moment_sum = 0
    for i in range(1, len(beam.interval)):
        moment_sum += (shear[i]+shear[i-1])/2*beam.interval_width
        moment[i] += moment_sum

    return shear, moment

def get_deflection(beam, rot):
    rotation = np.zeros_like(beam.interval)
    deflection = np.zeros_like(beam.interval)

    rotation[0] = rot
    for i in range (1,len(rotation)):
        rotation[i] += rotation[i-1] + beam.interval_width/beam.ei*(beam.moment[i-1]+beam.moment[i])/2
    
    deflection[0] = 0
    for i in range (1,len(deflection)):
        deflection[i] += deflection[i-1] + beam.interval_width*(rotation[i-1]+rotation[i])/2
    deflection_adj = np.add(-deflection[beam.support_index_l], deflection)

    return deflection_adj

def get_rotation(beam):
    "obtains initial rotation value by guess and check"
    deflection = np.zeros_like(beam.interval)
    delta = 1/beam.ei*beam.rotDelta
    def_list = [0.,0.,0.]
    for i, rot in enumerate([-delta, 0, delta]):
        deflection = get_deflection(beam, rot)
        def_list[i] = deflection[beam.support_index_r]

    if beam.cantilever or def_list[1] == 0:
        return 0
    elif abs(def_list[0]) < abs(def_list[2]):
        direction = -1
    elif abs(def_list[0]) > abs(def_list[2]):
        direction = 1
    
    def_test = def_list[1+direction]
    def_last = def_list[1]
    rot = delta*direction

    while abs(def_test) < abs(def_last):
        def_last = def_test
        rot += delta*direction
        deflection = get_deflection(beam, rot)
        def_test = deflection[beam.support_index_r]
    return rot - delta*direction

def main():
    # Initialize a beam
    beam = Beam(length=1, ei=29000000) # Defaults to simply-supported beam
    beam.correct()

    # Add point and distributed loads
    beam.addLoad(PointLoad(d=0.5, m=-1))                # Shear load acting at midpoint
    beam.addLoad(PointLoad(shear=False, d=0.25, m=-1))  # Moment load acting at 0.25 ft
    beam.addLoad(DistLoad(dl=0, dr=1, ml=-1, mr=-1))    # Distributed constant shear load
    
    # Calculate shear, moment, and deflection
    beam.calc_sm()
    beam.calc_def()

    # Plot shear/moment diagram and deflection diagram
    beam.plot_sm()
    beam.plot_def()

if __name__ == "__main__":
    main()