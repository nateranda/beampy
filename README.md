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

Finally, plot the shear/moment diagram and deflection diagram with `plot_sm` and `plot_def`:
```python
beam.plot_sm()
beam.plot_def()
```