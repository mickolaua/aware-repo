<!--
Author: Nicolai Pankov (colinsergesen@gmail.com)
development.md (c) 2023
Desc: description
Created:  2023-03-05
Modified: 2023-03-15
!-->

# Mosaic Scanning
Given a HEALPix (sky grid) of the LVK localization region with adaptive 
resolution and UNIQ indexing scheme, implement the algorithm of scanning this 
region with a FOV of the telescope of an interest.

## Possible solution
1. Regrid the localization region in such way each pixel of the region has the 
same area $width \times height$ (i.e. convert it to const-resolution HEALPix). 
It is very intuitive to select such $width$ and $height$ corresponding to the 
telescope FOV. 

2. Find the most probable pixel and start scanning from it.

3. Continue scanning starting from most probable to least probable pixels

4. Return the list of the coordinate centers of the scanned pixels.

## FAQ
### What is HEALPix?
> HEALPix is a genuinely curvilinear partition of the sphere into exactly 
equal area quadrilaterals of varying shape. 
(c) https://healpix.jpl.nasa.gov/pdf/intro.pdf 

There are actually two types of HEALPix: constant-reolution 
HEALPix (resolution is set by $N_{side}$ parameter: 
$N{pix} = 12 \times N_{side}^2$) and new Multi-Order Coverage (MOC) HEALPix 
with adaptive resolution depending on pixel location and signal-to-noise 
level of the gravitational-wave signal.

Each HEALPix pixel has an index depending on its location and parent pixel 
information (RING and NESTED). The MOC adaptive resoltion is possible thanks 
to the UNIQ indexing scheme that combines both sky position and resolution in 
one integer $uniq = ipix + 4 N_{side}^2$. (https://emfollow.docs.ligo.org/userguide/tutorial/multiorder_skymaps.html#the-uniq-indexing-scheme). Such an indexing scheme 
gives a musch less HEALPix localization region size: a few MiB against a few 
GiB compared to "old" indexing scheme HEALPixes.

See more information on HEALPix here https://healpix.jpl.nasa.gov/pdf/intro.pdf, 
here https://healpix.sourceforge.io/, and here https://emfollow.docs.ligo.org/userguide/tutorial/.

### Where to find LVK HEALPixes?
The LVK HEALPixes are distributed via GraceDB https://gracedb.ligo.org/. You 
can search there for specific events or all public events as well as the 
latest ones.

### How to work with LVK HEALPixes (MOCs)?
Gravitational wave MOCs could be read as `astropy.table.Table` (https://docs.astropy.org/en/stable/api/astropy.table.Table.html#astropy.table.Table) tables using 
the `astropy_healpix` (https://astropy-healpix.readthedocs.io/en/latest/) 
package. Example:

```
import numpy as np
from astropy.table import Table
import astropy_healpix as ah

# Open the HEALPix file
filename = "bayestar.multiorder.fits"
healpix = Table.read(filename, format="fits")

# Information from header
order = healpix.meta["MOCORDER"]
dist_mu = healpix.meta["DISTMEAN"]
dist_std = healpix.meta["DISTSTD"]
date_obs = healpix.meta["DATE-OBS"]

# Get uniq indicies and probability density per pixel (sr^-1)
uniq = healpix["UNIQ"]
prob_dens = healpix["PROBDENSITY"]

# Get indices of HEALPix pixels
level, ipix = ah.uniq_to_level_ipix(healpix["UNIQ"])

# Get the N_side resoluition parameter
nside = ah.level_to_nside(level)

# Get probability from prob. density and pixel area
area = ah.nside_to_pixel_area(nside)
prob = (prob_dens * area).value

# Coordinates of the pixel centers
lon, lat = ah.healpix_to_lonlat(ipix, nside)
```



