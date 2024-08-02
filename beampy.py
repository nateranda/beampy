import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

@dataclass
class Beam:
    "Beam & support attributes"
    cantilever: bool # cantilever or simply-suppoted
    length: float    # length, feet
    dl: float        # left support location, feet
    dr: float        # right support location, feet
    ei: float        # e*i, lb-in^2

    def dist(self) -> float:
        "Returns distance between supports"
        return self.dr-self.dl
    
    def correct(self):
        "Corrects invalid support positions"
        if self.cantilever == True:
            self.dl = 0
            self.dr = self.length
        else:
            if self.dr > self.length:
                self.dr = self.length
            if self.dl < 0:
                self.dl = 0

@dataclass
class Options:
    "Analysis options"
    sections: int    # number of integration sections
    rotDelta: float  # how much to change rotation delta - multiplier of ei
    
@dataclass
class PointLoad:
    "Point loads"
    shear: bool # shear or moment load
    d: float    # distance, ft
    m: float    # magnitude, lb

@dataclass
class DistLoad:
    "Distributed loads"
    dl: float # start location
    dr: float # end location
    ml: float # start magnitude
    mr: float # end magnitude

    def len(self) -> float:
        "Returns the span of the load"
        return self.dr-self.dl
    
    def mag(self) -> float:
        "Returns the magnitude of the load"
        return (self.ml+self.mr)/2*self.len()
    
    def pos(self) -> float:
        "Returns the centroid of the load in relation to the left end of the beam"
        return self.dl + (self.ml+2*self.mr)/(3*(self.ml+self.mr))*self.len()

def gen_interval(beam, options):
    "Generates integration interval"
    interval = np.linspace(0,beam.length,num=options.sections+1)
    return interval

def point_load_calc(beam, interval, point_loads):
    "Calculates shear & moment for point loads"
    shear = np.zeros_like(interval)
    moment = np.zeros_like(interval)

    for load in point_loads:
        # Shear loads
        if load.shear:
            for i, x in enumerate(interval):
                if x >= load.d:
                    shear[i] += load.m

            if beam.cantilever:
                shear = np.add(-load.m,shear)
                moment = np.add(load.d*load.m,moment)

            if not beam.cantilever:
                vr = load.m*(load.d-beam.dl)/beam.dist()
                vl = load.m - vr
                for i, x in enumerate(interval):
                    if x >= beam.dl:
                        shear[i] -= vl
                    if x >= beam.dr:
                        shear[i] -= vr
        
        # Moment loads
        if not load.shear:
            for i, x in enumerate (interval):
                if x >= load.d:
                    moment[i] -= load.m

            if beam.cantilever:
                moment = np.add(load.m, moment)

            if not beam.cantilever:
                vr = load.m/beam.dist()
                vl = -vr
                for i, x in enumerate(interval):
                    if x >= beam.dl:
                        shear[i] -= vl
                    if x >= beam.dr:
                        shear[i] -= vr
    
    # Integrate moment
    moment_sum = 0
    interval_width = beam.length/(len(interval)-1)
    for i in range(1, len(interval)):
        moment_sum += (shear[i]+shear[i-1])/2*interval_width
        moment[i] += moment_sum
    
    return shear, moment

def dist_load_calc(beam, interval, dist_loads):
    "Calculates shear & moment for distributed loads"
    shear = np.zeros_like(interval)
    moment = np.zeros_like(interval)

    # Shear
    for load in dist_loads:
        for i, x in enumerate(interval):
            if x >= load.dl and x <= load.dr:
                shear[i] += load.ml*(x-load.dl)+(x-load.dl)**2*(load.mr-load.ml)/(2*load.len())
            if x > load.dr:
                shear[i] += load.mag()

        if beam.cantilever:
            shear = np.add(-load.mag(), shear)
            moment = np.add(load.mag()*load.pos(), moment)

        if not beam.cantilever:
            vr = load.mag()*(load.pos()-beam.dl)/beam.dist()
            vl = load.mag()-vr
            for i, x in enumerate(interval):
                if x >= beam.dl:
                    shear[i] -= vl
                if x >= beam.dr:
                    shear[i] -= vr
    
    # Integrate moment
    moment_sum = 0
    interval_width = beam.length/(len(interval)-1)
    for i in range(1, len(interval)):
        moment_sum += (shear[i]+shear[i-1])/2*interval_width
        moment[i] += moment_sum

    return shear, moment

def load_calc(beam, interval, point_loads, dist_loads):
    "calculates and combines point and shear loads"
    pshear, pmoment = point_load_calc(beam, interval, point_loads)
    dshear, dmoment = dist_load_calc(beam, interval, dist_loads)

    shear = np.add(pshear, dshear)
    moment = np.add(pmoment, dmoment)

    return shear, moment

def get_support_indices(interval, beam):
    "gets index of each support in relation to interval"
    differences = np.abs(interval-beam.dl)
    support_index_l = differences.argmin()

    differences = np.abs(interval-beam.dr)
    support_index_r = differences.argmin()

    return [support_index_l, support_index_r]

def get_deflection(beam, moment, rot, support_indices):
    rotation = np.zeros_like(moment)
    deflection = np.zeros_like(moment)
    interval_width = beam.length/(len(moment)-1)

    rotation[0] = rot
    for i in range (1,len(rotation)):
        rotation[i] += rotation[i-1] + interval_width/beam.ei*(moment[i-1]+moment[i])/2
    
    deflection[0] = 0
    for i in range (1,len(deflection)):
        deflection[i] += deflection[i-1] + interval_width*(rotation[i-1]+rotation[i])/2
    deflection_adj = np.add(-deflection[support_indices[0]], deflection)

    return deflection_adj

def get_rotation(beam, options, moment, support_indices):
    "obtains initial rotation value by guess and check"
    delta = 1/beam.ei*options.rotDelta
    def_list = [0.,0.,0.]
    for i, rot in enumerate([-delta, 0, delta]):
        deflection = get_deflection(beam, moment, rot, support_indices)
        def_list[i] = deflection[support_indices[1]]

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
        deflection = get_deflection(beam, moment, rot, support_indices)
        def_test = deflection[support_indices[1]]
    return rot - delta*direction


    
def plot_sm(interval, shear, moment):
    "Plots shear/moment diagram"
    plt.title("Shear/Moment Diagram")
    plt.xlabel("Length (ft)")
    plt.ylabel("Stress (lb/ft-lb)")
    plt.axhline(0, color='black', linewidth=0.5)
    plt.plot(interval, shear)
    plt.plot(interval, moment)
    plt.show()

def plot_def(interval, deflection):
    "Plots shear/moment diagram"
    plt.title("Shear/Moment Diagram")
    plt.xlabel("Length (ft)")
    plt.ylabel("Deflection (in)")
    plt.axhline(0, color='black', linewidth=0.5)
    plt.plot(interval, deflection)
    plt.show()

def main():
    beam = Beam(cantilever=False, length=1, dl=0, dr=1, ei=290000000)
    beam.correct()
    options = Options(sections=1000, rotDelta=0.0001)

    point_loads = []
    point_loads.append(PointLoad(shear=False, d=beam.length/2, m=-2))

    dist_loads = []
    #dist_loads.append(DistLoad(dl=0, dr=beam.length, ml=-1, mr=-1))

    interval = gen_interval(beam, options)
    print(type(interval))

    shear, moment = load_calc(beam, interval, point_loads, dist_loads)

    support_indices = get_support_indices(interval, beam)
    rot = get_rotation(beam, options, moment, support_indices)
    deflection = get_deflection(beam, moment, rot, support_indices)
    print(deflection[support_indices[1]])

    #plot_sm(interval, shear, moment)
    plot_def(interval, deflection)

if __name__ == "__main__":
    main()