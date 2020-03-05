"""                         *                          *                                            *     *

  IIIIII    lLLLL    *    lLLLL       UUU        UUU    MMMMM      MMMMM    IIIIII     NNNN     NN          AAAA
   IIII     LLLL          LLLL   *    UUU        UUU    MMMMMMM  MMMMMMM     IIII   *  NNNNN    NN        AAAaaAAA
   IIII     LLLL          LLLL        UUU *      UUU    MMM MMMMMMMM MMM     IIII      NNNNNN   NN       AAA    AAA
   IIII  *  LLLL   *      LLLL        UUU        UUU    MMM *        MMM     IIII      NNN  NNN NN     AAAAAAAAAAAAAA
   IIII     LLLl          LLLl        UUUu      uUUU    MMM          MMM     IIII      NNN   NNNNN    AAAa        aAAA
   IIII    LLLLLLLLLL    LLLLLLLLLL    UUUUUuuUUUUU     MMM     *    MMM     IIII      NNN    NNNN   aAAA    *     AAAa
  IIIIII   LLLLLLLLLLL   LLLLLLLLLLL     UUUUUUUU      mMMMm        mMMMm   IIIIII    nNNNn    NNNn  aAAA          AAAa

 **********************************************************************************************************************
 ** Illumina VERSION 3 - in Python 3.6.3                                                                             **
 ** Programmers in decreasing order of contribution  :                                                               **
 **                            Martin Aube                                                                           **
 **                            Julien-Pierre Houle                                                                   **
 **                                                                                                                  **
 **              Still having very few traces of their contributions :                                               **
 **                            Loic Franchomme-Fosse,  Mathieu Provencher, Andre Morin                               **
 **                            Alex Neron, Etienne Rousseau                                                          **
 **                            William Desroches, Maxime Girardin, Tom Neron                                         **
 **                                                                                                                  **
 **                                                                                                                  **
 **  Current version features/limitations :                                                                          **
 **                                                                                                                  **
 **    - Calculation of artificial sky radiance in a given line of sight                                             **
 **    - Calculation of the atmospheric transmittance and 1st and 2nd order of scattering                            **
 **    - Lambertian reflexion on the ground                                                                          **
 **    - Terrain slope considered (apparent surface and shadows)                                                     **
 **    - Angular photometry of a lamp is considered uniform along the azimuth                                        **
 **    - Sub-grid obstacles considered (with the mean free path of light toward ground, mean obstacle height, and    **
 **      obstacles transparency (filling factor)                                                                     **
 **    - Molecules and aerosol optics (phase function, scattering probability, aerosol absorption)                   **
 **      molecular absorption not consideres (given that we focus on the visible                                     **
 **    - Exponential concentrations vertical profile (H aerosol= 2km, H molecules= 8km  )                            **
 **    - Accounting for heterogeneity luminaires number, luminaires heights, luminaire spectrum,                     **
 **      angular photometry, obstacle properties                                                                     **
 **    - Wavelength dependant                                                                                        **
 **    - Cloud models (type and cloud base height) only the overhead clouds are considered                           **
 **    - Do not support direct observation of a source                                                               **
 **    - Direct observation of the ground not implemented                                                            **
 **    - Do not consider earth curvature (i.e. local/regional model)                                                 **
 **                                                                                                                  **
 **********************************************************************************************************************

  Copyright (C) 2019 Martin Aube PhD

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but without any warranty; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Contact: martin.aube@cegepsherbrooke.qc.ca

**********************************************************************************************************************"""


# ======================================================================
#     Import Functions and librairies
# ======================================================================

from horizon import horizon


import numpy as np
from math import pi, atan, sin


# ======================================================================
#     Variables declaration
# ======================================================================

witdh = 512                                                   # Matrix dimension in Length/width and height
nzon = 256
pi4 = 4.*pi     # pourquoi?

drefle = np.zeros((witdh, witdh))                             # Mean free path to the ground (meter)
val2d = np.zeros((witdh, witdh))                              # Temporary input array 2d
altsol = np.zeros((width,width))                              # Ground elevation (meter)
lamplu = np.zeros((width,width,nzon))                         # Source fluxes
lampal = np.zeros((width,width))                              # Height of the light sources relative to the ground (meter)
pval = np.zeros((181,nzon))                                   # Values angular photometry functions (unnormalized)
pvalno = np.zeros((181,nzon))                                 # Values angular photometry functions (normalized)
fdifa = np.zeros((181))                                       # Aerosol scattering functions (unnormalized)
fdifan = np.zeros((181))                                      # Aerosol scattering functions (normalized)
anglea = np.zeros((181))                                      # Aerosol scattering angle (degree)
inclix = np.zeros((width,width))                              # Tilt of the ground pixel along x (radian)
incliy = np.zeros((width,width))                              # Tilt of the ground pixel along y (radian)
zondif = np.zeros(3000000,3)                                  # Array for the scattering voxels positions
obsH = np.zeros((width,width)),angmin                         # Averaged height of the sub-grid obstacles
ofill = np.zeros((width,width))                               # Fill factor giving the probability to hit an obstacle when pointing in #                                                               its direction integer 0-100
ITT = np.zeros((width,width,nzon))                            # Total intensity per type of lamp
ITC = np.zeros((width,width))                                 # Total intensity per line of sight voxel
FTC = np.zeros((width,width))                                 # Fraction of the total flux at the sensor level
FCA = np.zeros((width,width))                                 # Sensor flux array
lpluto = np.zeros((width,width))                              # Total luminosity of the ground cell for all lamps

# Entier
imin = np.zeros((nzon))                                       # x and y limits containing a type of lamp
imax = np.zeros((nzon))
jmin = np.zeros((nzon))
jmax = np.zeros((nzon))

totlu = np.zeros((nzon))                                      # Total flux of a source type
flcld = np.zeros((width,width))                               # Flux crossing a low cloud

verbose = 2         # explications?                           # Very little printout=0, Many printout = 1, even more=2
diamobj = 1.        # explications?                           # Dumb value for diameter of the objective instrument of the observer.
volu = 0.
zero = 0.
un = 1.
ff = 0.
dstep = 1
ncible = 1024
stepdi = 1

if (verbose >= 1):
    print('Starting ILLUMINA computations...')


# Reading of illumina.in (input file)
# À voir (228 à 255)


siz = 10                # Explications
if (ssswit == 0):
    effdif = 0.
else:
    effdif=100000.                                             # This is apparently the minimum value to get some accuracy
    n2nd=100000

scal = 20
scalo = scal

# -----------------------------------------------------------------------------------
# Note:                                                                              |
# omemax: exclude calculations too close (<57m) this is a sustended angle of 1 deg.  |
# the calculated flux is highly sensitive to that number for a very high             |
# pixel resolution (a few 10th of meters). We assume anyway that somebody            |
# observing the sky will never lies closer than that distance to a                   |
# light fixture. This number is however somehow subjective and that means            |
# that the value of sky brightness near sources will be affected by this choice      |
# -----------------------------------------------------------------------------------

omemax = 1./((25.)**2.)
if (verbose > 1):
    print('2nd scattering grid = ', siz)
    print('2nd order scattering radius=', effdif,'m')
    print('Pixel size = ', dx,' x ', dy)
    print('Maximum radius for reflexion = ', reflsiz)

# Computing the actual AOD at the wavelength lambda      # lambda = lambd
if (verbose >= 1):
    print('500nm AOD=', taua, '500nm angstrom coeff.=', alpha)
    taua = taua*(lambd/500.)**(-1.*alpha)

"""    # Determine the Length of basenm        # Explications + traductions fortran
    lenbase = index(basenm,' ')-1
    mnaf = basenm(1:lenbase)//'_topogra.bin'                      # Determine the names of input and output files
    outfile = basenm(1:lenbase)//'.out'
    pclf = basenm(1:lenbase)//'_pcl.txt'
    pclimg = basenm(1:lenbase)//'_pcl.bin'
    pcwimg = basenm(1:lenbase)//'_pcw.bin'
    pclgp = basenm(1:lenbase)//'_pcl.gplot'
    """

# -------------------------------------------------------------------------
# Note:                                                                    |
# conversion of the geographical viewing angles toward the cartesian       |
# angle we assume that the angle in the file illumina.in                   |
# is consistent with the geographical definition                           |
# geographical: azim=0 toward north, 90 toward east, 180 toward south, etc |
# cartesian: azim=0 toward east, 90 toward north, 180 toward west, etc     |
# -------------------------------------------------------------------------

azim = 90.-azim
    if (azim < 0.):
        azim = azim+360.
    if (azim >= 360.):
        azim = azim-360.

# opening output file
# À faire!

# Check if the observation angle is above horizon
    angzen = pi/2.-angvis*pi/180.
    horizon(x, y, z, dx, dy, anga)

    if (angzen > zhoriz):                        # the line of sight is not below the horizon => we compute
        raise ValueError("""PROBLEM! You try to observe below horizon.
        No calculation will be made""")
#        write(2,*) '            Sky radiance (W/str/m**2)          '
#        write(2,2001) zero
        print('            Sky radiance (W/str/m**2)          ')
        print('                 0.0000')
#        close(2)


"""        write(2,*) 'FILE USED:'
        write(2,*) mnaf,diffil
        print*,'Wavelength (nm):',lambda,
     +       ' Aerosol optical depth:',taua
        write(2,*) 'Wavelength (nm):',lambda,
     +       ' Aerosol optical depth:',taua
        write(2,*) '2nd order scattering radius:',effdif,' m'
        print*,'2nd order scattering radius:',effdif,' m'
        write(2,*) 'Observer position (x,y,z)',x_obs,y_obs,z_o
        print*,'Observer position (x,y,z)',x_obs,y_obs,z_o
        write(2,*) 'Elevation angle:',angvis,' azim angle (clockwise fro
     +m north)',azim
        print*,'Elevation angle:',angvis,' azim angle (counterclockwise
     +from east)',azim      """


# Initialisation of the arrays and variables

if (verbose >= 1):
    print('Initializing variables...')
if (cloudt == 0):
    cloudbase=1000000000.

prmaps=1
iun=0
ideux=1
icloud=0.

for i in range(1, witdh):           # Est-ce que j'ai vrm besoin à cause de np.zeros?
    for j in range(1, witdh):
        val2d[i,j]=0.
        altsol[i,j]=0.
        obsH[i,j]=0.
        ofill[i,j]=0.
        inclix[i,j]=0.
        incliy[i,j]=0.
        lpluto[i,j]=0.
        ITC[i,j]=0.
        FTC[i,j]=0.
        FCA[i,j]=0.
        flcld[i,j]=0
        for k in range(1, nzon):
            lamplu[i,j,k]=0.
            lampal[i,j]=0.
            ITT[i,j,k]=0.

for i in range(1, 181):
    fdifa[i]=0.
    fdifan[i]=0.
    anglea[i]=0.
    for j in range(1, nzon):
        pval[i,j]=0.
        pvalno[i,j]=0.

for i in range(1, 3000000):     # Regarder si on peut set a 1 dès le début
    for j in range(1, 3):
        zondif[i, j] = 1.

    idif1=0.
    idif2=0.
    fdif2=0
    idif2p=0.
    fldir=0.
    flindi=0.
    fldiff=0.
    pdifdi=0.
    pdifin=0.
    pdifd1=0.
    pdifd2=0.
    intdir=0.
    intind=0.
    idiff2=0.
    angmin=0.
    isourc=0.
    itotty=0.
    itotci=0.
    itotrd=0.
    flcib=0.
    flrefl=0.
    irefl=0.
    irefl1=0.
    fldif1=0.
    fldif2=0.
    portio=0.
    fccld=0.
    fctcld=0.
    ometif=0.
    omefov=0.
    hh=1.

# determination of the vertical atmospheric transmittance
# tranam and tranaa are the top of atmosphere transmittance (molecules and aerosols)
transtoa(lambm, taua, pressi, tranam, tranaa)       # COMMENT LA VALEUR DES ARGUMENTS EST DÉTERMINÉ
# reading of the environment variables
# reading of the elevation file
 """call twodin(nbx,nby,mnaf,altsol)"""
# computation of the tilt of the cases along x and along y



for i in range(1, nbx):                           # Beginning of the loop over the column (longitude) of the domain.
    for j in range(1, nby):                       # Beginning of the loop over the rows (latitude) of the domain.
        if (i == 1):                              # Specific case close to the border of the domain (vertical side left).
            inclix[i,j] = atan((altsol[i+1, j]-altsol[i,j])/dx)          # Computation of the tilt along x of the surface.
        elif (i == nbx):                          # Specific case close to the border of the domain (vertical side right).
            inclix[i,j] = atan((altsol[i-1, j]-altsol[i,j])/dx)          # Computation of the tilt along x of the surface.
        else:
            #                                   VÉRIFIER QUE C'EST 2.1!
            inclix[i,j] = atan((altsol[i+1, j]-altsol[i-1,j])/(2.1*dx))  # Computation of the tilt along x of the surface.

        if (j == 1):                              # Specific case close to the border of the domain (horizontal side down).
            incliy[i,j] = atan((altsol[i, j+1]-altsol[i,j])/dy)          # Computation of the tilt along y of the surface.
        elif (j == nby):                          # Specific case close to the border of the domain (horizontal side up).
            incliy[i,j] = atan((altsol[i, j-1]-altsol[i,j])/dy))         # Computation of the tilt along y of the surface.
        else:
            incliy[i,j] = atan((altsol[i, j+1]-altsol[i,j-1])/(2.1*dy))  # Computation of the tilt along y of the surface


# Reading of the values of P(theta), height, luminosities and positions
# of the sources, obstacle height and distance
"""
        ohfile=basenm(1:lenbase)//'_obsth.bin'
        odfile=basenm(1:lenbase)//'_obstd.bin'
        alfile=basenm(1:lenbase)//'_altlp.bin'                            # setting the file name of height of the sources lumineuse.
        offile=basenm(1:lenbase)//'_obstf.bin'
        dtheta=.017453293                                                 # one degree
"""

# Reading lamp heights
"""        call twodin(nbx,nby,alfile,val2d)   """
    for i in range(1, nbx):                 # Beginning of the loop over all cells along x.
        for j in range(1, nby):             # Beginning of the loop over all cells along y.
            lampal[i,j] = val2d[i,j]        # Filling of the array for the lamps type           #   COMMENTAIRE ?

# Reading subgrid obstacles average height
"""    call twodin(nbx,nby,ohfile,val2d)      """
    for i in range(1, nbx):                 # Beginning of the loop over all cells along x.
        for j in range(1, nby):             # Beginning of the loop over all cells along y.
            obsH[i,j] = val2d[i,j]          # Filling of the array

# Reading subgrid obstacles average distance
"""    call twodin(nbx,nby,odfile,val2d)    """
    for i in range(1, nbx):                 # Beginning of the loop over all cells along x.
        for j in range(1, nby):             # Beginning of the loop over all cells along y.
            if (drefle[i,j] == 0.):         # When outside a zone, block to the size of the cell (typically 1km)
                drefle[i,j] = dx

# Reading subgrid obstacles filling factor
"""        call twodin(nbx,nby,offile,val2d)    """
    for i in range(1, nbx):                 # Beginning of the loop over all cells along x.
        for j in range(1, nby):             # Beginning of the loop over all cells along y.
            ofill[i,j]=val2d[i,j]           # Filling of the array 0-1

# Reading of the scattering parameters
"""
    Faire la section de la ligne 465 à 479
"""


# Quoi faire avec toute la section de commentaire???

for stype in range(1, ntype):               # Beginning of the loop 1 for the nzon types of sources.
    imin[stype]=nbx
    jmin[stype]=nby
    imax[stype]=1
    jmax[stype]=1
    pvalto=0.
"""       write(lampno, '(I3.3)' ) stype                                  ! support of nzon different sources (3 digits)
          pafile=basenm(1:lenbase)//'_fctem_'//lampno//'.dat'             ! setting the file name of angular photometry.
          lufile=basenm(1:lenbase)//'_lumlp_'//lampno//'.bin'             ! setting the file name of the luminosite of the cases.
"""


# Reading photometry files
#          open(UNIT=1, FILE=pafile,status='OLD')                          ! opening file pa#.dat, angular photometry.

for i in range(1, 181):
#    read(1,*) pval(i,stype)                                     ! reading of the data in the array pval.
    pvalto=pvalto+pval[i,stype]*2.*pi*sin((i-1)*dtheta)*dtheta      # Sum of the values of the  photometric function
                                                                    # (pvaleur x 2pi x sin theta x dtheta) (ou theta egale (i-1) x 1 degrees).

                                                                    # End of the loop over the 181 donnees of the fichier pa#.dat.
"""close(1)                                                        ! closing file pa#.dat, angular photometry."""
for i in range(1, 181):
    if (pvalto != 0):
        pvalno[i,stype] = pval[i,stype]/pvalto                        # Normalisation of the photometric function.



# Reading luminosity files
"""          call twodin(nbx,nby,lufile,val2d)      """


for i in range(1, nbx):                 # Beginning of the loop over all cells along x.
    for j in range(1, nby):             # Beginning of the loop over all cells along y.
        if (val2d[i,j] < 0):            # Searching of negative fluxes
            raise ValueError("Negative lamp flux!, stopping execution")

for i in range(1, nbx):
    for j in range(1, nby):             # Searching of the smallest rectangle containing the zone
        if (val2d[i,j] != 0):           # of non-null luminosity to speedup the calculation
            if (i-1 < imin[stype]):
                imin[stype]=i-2
            if (imin[stype] < 1):       # Est ce que ce if rentre dans le if au dessus?
                imin[stype] = 1

""" VÉRIFIER GOTO FONCTIONNEMENT """
# COMPLÉTER CE BLOC DE LA LIGNE 642 À 698



# Some preliminary tasks
    dy=dx
    omefov=0.00000001                         # Solid angle of the spectrometer slit on the sky. Here we only need a small value
    z_obs=z_o+altsol[x_obs,y_obs]             # z_obs = the local observer elevation plus the height of observation above ground (z_o)
    rx_obs = x_obs*dx
    ry_obs = y_obs*dy
    if (z_obs == 0.):
        z_obs = 0.001
    largx=dx*nbx                              # Computation of the Width along x of the case.
    largy=dy*nby                              # Computation of the Width along x of the case.

"""        write(2,*) 'Width of the domain [NS](m):',largx,'#cases:',nbx
        write(2,*) 'Width of the domain [EO](m):',largy,'#cases:',nby
        write(2,*) 'Size of a cell (m):',dx,' X ',dy
        write(2,*) 'latitu center:',latitu
"""
