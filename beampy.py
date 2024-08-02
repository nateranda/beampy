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
        self.support_index_r = np.abs(self.interval-self.dr).argmin()     # index of right support

        # Correct supports
        if self.cantilever == True:
            self.dl = 0
            self.dr = self.length
        else:
            if self.dr > self.length:
                self.dr = self.length
            if self.dl < 0:
                self.dl = 0

class PointLoad:
    "Point loads"
    def __init__(self, d, m, shear=None):
        self.d = d                                    # distance, ft
        self.m = m                                    # magnitude, lb
        self.shear = True if shear is None else shear # shear or moment load

class DistLoad:
    "Distributed loads"
    def __init__(self,dl, dr, ml, mr):
        self.dl = dl # start location
        self.dr = dr # end location
        self.ml = ml # start magnitude
        self.mr = mr # end magnitude
        
        self.len = self.dr-self.dl                                               # load span
        self.mag = (self.ml+self.mr)/2*self.len                                  # load magnitude
        self.pos = self.dl + (self.ml+2*self.mr)/(3*(self.ml+self.mr))*self.len  # centroid w.r.t beam end

def point_load_calc(beam, point_loads):
    "Calculates shear & moment for point loads"
    shear = np.zeros_like(beam.interval)
    moment = np.zeros_like(beam.interval)

    for load in point_loads:
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

def dist_load_calc(beam, dist_loads):
    "Calculates shear & moment for distributed loads"
    shear = np.zeros_like(beam.interval)
    moment = np.zeros_like(beam.interval)

    # Shear
    for load in dist_loads:
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

def load_calc(beam, point_loads, dist_loads):
    "Calculates and combines point and shear loads"
    pshear, pmoment = point_load_calc(beam, point_loads)
    dshear, dmoment = dist_load_calc(beam, dist_loads)

    shear = np.add(pshear, dshear)
    moment = np.add(pmoment, dmoment)

    return shear, moment

def get_deflection(beam, moment, rot):
    rotation = np.zeros_like(moment)
    deflection = np.zeros_like(moment)

    rotation[0] = rot
    for i in range (1,len(rotation)):
        rotation[i] += rotation[i-1] + beam.interval_width/beam.ei*(moment[i-1]+moment[i])/2
    
    deflection[0] = 0
    for i in range (1,len(deflection)):
        deflection[i] += deflection[i-1] + beam.interval_width*(rotation[i-1]+rotation[i])/2
    deflection_adj = np.add(-deflection[beam.support_index_l], deflection)

    return deflection_adj

def get_rotation(beam, moment):
    "obtains initial rotation value by guess and check"
    delta = 1/beam.ei*beam.rotDelta
    def_list = [0.,0.,0.]
    for i, rot in enumerate([-delta, 0, delta]):
        deflection = get_deflection(beam, moment, rot)
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
        deflection = get_deflection(beam, moment, rot)
        def_test = deflection[beam.support_index_r]
    return rot - delta*direction
  
def plot_sm(beam, shear, moment):
    "Plots shear/moment diagram"
    plt.title("Shear/Moment Diagram")
    plt.xlabel("Length (ft)")
    plt.ylabel("Stress (lb/ft-lb)")
    plt.axhline(0, color='black', linewidth=0.5)
    plt.plot(beam.interval, shear)
    plt.plot(beam.interval, moment)
    plt.show()

def plot_def(beam, deflection):
    "Plots shear/moment diagram"
    plt.title("Shear/Moment Diagram")
    plt.xlabel("Length (ft)")
    plt.ylabel("Deflection (in)")
    plt.axhline(0, color='black', linewidth=0.5)
    plt.plot(beam.interval, deflection)
    plt.show()

def main():
    beam = Beam(length=1,ei=290000000)

    point_loads = []
    point_loads.append(PointLoad(d=beam.length/2, m=-2))

    dist_loads = []
    #dist_loads.append(DistLoad(dl=0, dr=beam.length, ml=-1, mr=-1))

    shear, moment = load_calc(beam, point_loads, dist_loads)

    rot = get_rotation(beam, moment)
    deflection = get_deflection(beam, moment, rot)

    print(deflection[beam.support_index_r])

    #plot_sm(beam, shear, moment)
    plot_def(beam, deflection)

if __name__ == "__main__":
    main()