# Beampy
An object-oriented Python package for calculating simply supported and cantilever beam stresses. It uses Numpy for array calculations and Matplotlib for plotting.

## Usage
Start by initializing a beam with the `Beam` object:
```python
from beampy import Beam

beam = Beam(length=1, ei=29000000) # Defaults to simply-supported beam
```

Then, add point and distributed loads with `addLoad`:
```python
from beampy import PointLoad, DistLoad

beam.addLoad(PointLoad(d=0.5, m=-1))                # Shear load acting at midpoint
beam.addLoad(PointLoad(shear=False, d=0.25, m=-1))  # Moment load acting at 0.25 ft
beam.addLoad(DistLoad(dl=0, dr=1, ml=-1, mr=-1))    # Distributed constant load
```

Once all loads are applied, calculate the shear and moment stresses as well as the deflection along the beam with `calc_sm` and `calc_def`:
```python
beam.calc_sm()
beam.calc_def()
```

```
Output:
Max Shear: 0.0 lb
Min Shear: -1.999 lb
Max Moment: 0.96875 lb-ft
Min Moment: -0.031001 lb-ft
Max Positive Deflection: 0.000000e+00 in
Max Negative Deflection: -2.787676e-09 in
```

Finally, plot the shear/moment diagram and deflection diagram with `plot_sm` and `plot_def`:
```python
beam.plot_sm()
beam.plot_def()
```

![plot_sm graph](https://raw.githubusercontent.com/nateranda/beampy/main/static/plot_sm.png)
![plot_def graph](https://raw.githubusercontent.com/nateranda/beampy/main/static/plot_def.png)

You can also define load types and test different ASCE load combinations. Just assign a type to each load and find critical load combinations with `find_lc`:
```python
beam.addLoad(PointLoad(d=0.5, m=-1, type="D"))                # Dead load
beam.addLoad(PointLoad(shear=False, d=0.25, m=-1, type="L"))  # Live load
beam.addLoad(DistLoad(dl=0, dr=1, ml=-1, mr=-1, type="E"))    # Seismic load

beam.find_lc()
```

```
Output:
Max Shear: 0.95 from load combination 7
Min Shear: -2.2 from load combination 2
Max Moment: 1.35 from load combination 2
Min Moment: -0.249 from load combination 2
```

## Objects

### `beampy.Beam(length float, dl float, dr float, cantilever bool, ei int, method bool, sections int, rotDelta bool)`
```
length: length of beam, ft
dl: location of left support from left end of beam, ft. Default left end of beam
dr: location of right support from left end of beam, ft. Default right end of beam
cantilever: cantilever or simply-supported beam. Default False (simply-supported)
ei: modulus of elasticity * moment of inertia, lb-in^2
method: analysis method, lrfd or asd. Default lrfd
sections: number of discrete sections to use when integrating. Default 1000
rotDelta: how much to change rotation in deflection calculation. Default 0.0001
```

### `beampy.PointLoad(d float, m float, shear bool, type str)`
```
d: distance of load from left end of beam, ft
m: magnitude of load, lb
shear: shear or moment. Default True (shear)
type: load type. Default None
```

### `beampy.DistLoad(dl float, dr float, ml float, mr float, type str)`
```
dl: distance of start from left end of beam, ft
dr: distance of end from left end of beam, ft
ml: magnitude of start, lb
mr: magnitude of end, lb
type: load type. Default None
```