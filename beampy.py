import numpy as np
import matplotlib.pyplot as plt

class Beam:
    "Beam & support attributes"
    def __init__(self, length, ei, method=None, cantilever=None, dl=None, dr=None, sections=None, rotDelta=None):

        self.cantilever = False if cantilever is None else cantilever    # cantilever or simply-suppoted
        self.length = length                                             # length, feet
        self.dl = 0 if dl is None else dl                                # left support location, feet
        self.dr = self.length if dr is None else dr                      # right support location, feet
        self.dist = self.dr-self.dl                                      # distance between supports, feet
        self.ei = ei                                                     # e*i, lb-in^2
        self.method = "lrfd" if method is None else method               # analysis method, asd or lrfd

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
    
    def calc_sm(self, summary=None, lc=None):
        pshear, pmoment = point_load_calc(self, lc)
        dshear, dmoment = dist_load_calc(self, lc)
        self.shear = np.add(pshear, dshear)
        self.moment = np.add(pmoment, dmoment)

        if summary == None or summary == True:
            max_shear = round(np.max(self.shear), 6)
            min_shear = round(np.min(self.shear), 6)
            max_moment = round(np.max(self.moment), 6)
            min_moment = round(np.min(self.moment), 6)
            print(f"Max Shear: {max_shear} lb")
            print(f"Min Shear: {min_shear} lb")
            print(f"Max Moment: {max_moment} lb-ft")
            print(f"Min Moment: {min_moment} lb-ft")
    
    def find_lc(self):
        if self.method == "lrfd":
            load_list = lrfd
        elif self.method == "asd":
            load_list = asd
        else:
            raise Exception("Invalid method type.")
        
        max_shear_list = []
        min_shear_list = []
        max_moment_list = []
        min_moment_list = []

        for i in range(0,len(load_list)):
            pshear, pmoment = point_load_calc(self, i+1)
            dshear, dmoment = dist_load_calc(self, i+1)
            self.shear = np.add(pshear, dshear)
            self.moment = np.add(pmoment, dmoment)

            max_shear_list.append(round(np.max(self.shear), 6))
            min_shear_list.append(round(np.min(self.shear), 6))
            max_moment_list.append(round(np.max(self.moment), 6))
            min_moment_list.append(round(np.min(self.moment), 6))
        
        max_shear = max(max_shear_list)
        min_shear = min(min_shear_list)
        max_moment = max(max_moment_list)
        min_moment = min(min_moment_list)
        
        print(f"Max Shear: {max_shear} from load combination {max_shear_list.index(max_shear)+1}")
        print(f"Min Shear: {min_shear} from load combination {min_shear_list.index(min_shear)+1}")
        print(f"Max Moment: {max_moment} from load combination {max_moment_list.index(max_moment)+1}")
        print(f"Min Moment: {min_moment} from load combination {min_moment_list.index(min_moment)+1}")

    
    def calc_def(self, summary=None):
        rot = get_rotation(self)
        self.deflection = get_deflection(self, rot)

        if summary == None or summary == True:
            max_def = np.max(self.deflection)
            min_def = np.min(self.deflection)
            print(f"Max Positive Deflection: {max_def:.6e} in")
            print(f"Max Negative Deflection: {min_def:.6e} in")
    
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
    def __init__(self, d, m, shear=None, type=None):
        self.d = d                                    # distance, ft
        self.m = m                                    # magnitude, lb
        self.shear = True if shear is None else shear # shear or moment load
        self.type = type                              # ASCE load type

        # Validate type
        valid_types = [None, "Dead", "D", "Live", "L", "Roof Live", "RL", "Wind", "W", "Snow", "S", "Seismic", "E"]
        if type not in valid_types:
            raise Exception(f"Invalid load type {type}. Must be from list: {valid_types}")

class DistLoad:
    "Distributed load object"
    def __init__(self,dl, dr, ml, mr, type=None):
        self.dl = dl     # start location
        self.dr = dr     # end location
        self.ml = ml     # start magnitude
        self.mr = mr     # end magnitude
        self.type = type # ASCE load type
        
        self.len = self.dr-self.dl                                               # load span
        self.mag = (self.ml+self.mr)/2*self.len                                  # load magnitude
        self.pos = self.dl + (self.ml+2*self.mr)/(3*(self.ml+self.mr))*self.len  # centroid w.r.t beam end

        # Validate type
        valid_types = [None, "Dead", "D", "Live", "L", "Roof Live", "RL", "Snow", "S", "Rain", "R", "Wind", "W", "Seismic", "Earthquake", "E"]
        if type not in valid_types:
            raise TypeError(f"Invalid load type {type}. Must be from list: {valid_types}")

lrfd = [[1.4,0,0,0,0,0,0],
       [1.2,1.6,0.5,0.5,0.5,0,0],
       [1.2,1,1.6,1.6,1.6,0.5,0],
       [1.2,1,0.5,0.5,0.5,1,0],
       [1.2,1,0,0.2,0,0,1],
       [0.9,0,0,0,0,1,0],
       [0.9,0,0,0,0,0,1]]

asd = [[1,0,0,0,0,0,0],
       [1,1,0,0,0,0,0],
       [1,0,1,1,1,0,0],
       [1,0.75,0.75,0.75,0.75,0,0],
       [1,0,0,0,0,0.6,0.7],
       [1,0.75,0.75,0.75,0.75,0.45,0],
       [1,0.75,0,0.75,0,0,0.525],
       [0.6,0,0,0,0,0.6,0],
       [0.6,0,0,0,0,0,0.7]]

def get_mult(load, lc, method):
    if method == "lrfd":
        mult_list = lrfd[lc]
    elif method == "asd":
        mult_list = asd[lc]
    else:
        raise Exception("Invalid method type.")
    
    match load.type:
        case None:
            return 1
        case "Dead" | "D":
            return mult_list[0]
        case "Live" | "L":
            return mult_list[1]
        case "Roof Live" | "RL":
            return mult_list[2]
        case "Snow" | "S":
            return mult_list[3]
        case "Rain" | "R":
            return mult_list[4]
        case "Wind" | "W":
            return mult_list[5]
        case "Seismic" | "Earthquake" | "E":
            return mult_list[6]
        case _:
            raise Exception("Invalid load type.")

def point_load_calc(beam, lc):
    "Calculates shear & moment for point loads"
    shear = np.zeros_like(beam.interval)
    moment = np.zeros_like(beam.interval)

    for load in beam.point_loads:
        if lc == None:
            mult = 1
        else:
            mult = get_mult(load, lc-1, beam.method)
        mag = load.m * mult

        # Shear loads
        if load.shear:
            for i, x in enumerate(beam.interval):
                if x >= load.d:
                    shear[i] += mag

            if beam.cantilever:
                shear = np.add(-mag,shear)
                moment = np.add(load.d*mag,moment)

            if not beam.cantilever:
                vr = mag*(load.d-beam.dl)/beam.dist
                vl = mag - vr
                for i, x in enumerate(beam.interval):
                    if x >= beam.dl:
                        shear[i] -= vl
                    if x >= beam.dr:
                        shear[i] -= vr
        
        # Moment loads
        if not load.shear:
            for i, x in enumerate (beam.interval):
                if x >= load.d:
                    moment[i] -= mag

            if beam.cantilever:
                moment = np.add(mag, moment)

            if not beam.cantilever:
                vr = mag/beam.dist
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

def dist_load_calc(beam, lc):
    "Calculates shear & moment for distributed loads"
    shear = np.zeros_like(beam.interval)
    moment = np.zeros_like(beam.interval)

    # Shear
    for load in beam.dist_loads:
        if lc == None:
            mult = 1
        else:
            mult = get_mult(load, lc-1, beam.method)
        magl = load.ml * mult
        magr = load.mr * mult
        mag = load.mag * mult

        for i, x in enumerate(beam.interval):
            if x >= load.dl and x <= load.dr:
                shear[i] += magl*(x-load.dl)+(x-load.dl)**2*(magr-magl)/(2*load.len)
            if x > load.dr:
                shear[i] += mag

        if beam.cantilever:
            shear = np.add(-mag, shear)
            moment = np.add(mag*load.pos, moment)

        if not beam.cantilever:
            vr = mag*(load.pos-beam.dl)/beam.dist
            vl = mag-vr
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