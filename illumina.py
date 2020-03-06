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

from horizon import horizons
from transTOA import transtoa
from twodin import twodin
from anglezenithal import anglezenithal
from angleazimutal import angleazimutal
from transmitm import transmitm
from transmita import transmita
from angle3points import angle3points
from diffusion import diffusion
from cloudreflectance import cloudreflectance
from anglesolide import anglesolide
from zone_diffusion import zone_diffusion
from twodout import twodout


import numpy as np
from math import sqrt, pi, sin, cos, tan, atan

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
    largy=dy*nby                              # Computation of the Width along y of the case.

"""        write(2,*) 'Width of the domain [NS](m):',largx,'#cases:',nbx
        write(2,*) 'Width of the domain [EO](m):',largy,'#cases:',nby
        write(2,*) 'Size of a cell (m):',dx,' X ',dy
        write(2,*) 'latitu center:',latitu
"""



# FAIRE LE BLOC DE LA LIGNE 716 à 801


# Beginning of the loop over the types of light sources

for stype in range(1, ntype):               # Beginning of the loop over the source types.
    if (totlu[stype] != 0):                 # Check if there are any flux in that source type otherwise skip this lamp
        if (verbose >= 1):
            print("Turning on lamps", stype)
"""     if (verbose.ge.1) write(2,*) ' Turning on lamps' stype      """
        itotty=0.
        for x_s in range(1, nbx):           # Initialisation of the contribution of a source types to
            for y_s in range(1, nby):       # the intensity toward the sensor by a line of sight voxel.
                ITT[x_s, y_s, stype]=0.     "Vérifier que j'ai besoin de faire ça"

        for x_s in range(imin[stype], imax[stype])           # Beginning of the loop over the column (longitude the) of the domain.
            for y_s in range(jmin[stype], jmax[stype])       # Beginning of the loop over the rows (latitud) of the domain.
                intdir=0.
                itotind=0.
                itodif=0.
                itotrd=0.
                isourc=0.
                rx_s = x_s*dx
                ry_s = y_s*dy
                if (lamplu[x_s, y_s, stype] != 0):                 # If the luminosite of the case is null, ignore this case.
                    z_s = (altsol[x_s, y_s] + lampal[x_s, y_s])    # Definition of the position (metre) vertical of the source.

# *********************************************************************************************************
# * computation of the direct intensity toward the observer by a line of sight voxel from the source      *
# *********************************************************************************************************

                    dirck=0                             # Initialisation of the verification of the position of the source
                    if ((rx_s == rx_c) and (ry_s == ry_c) and (z_s == z_c)):   # If the position of the source and the line of sight
                        dirck =1                                                # voxel are the same then ...
                        if (verbose == 2):
                            if (verbose >= 1):      # Erreur? Pk faire ça?
                                print("Source = line of sight")
                    # End of the case positions x and y source and line of sight voxel identical.


                    if (dirck != 1):                  # The source is not at the line of sight voxel position
# computation of the zenithal angle between the source and the line of sight
# computation of the horizon for the resolved shadows direct              ! horizon resolution is 1 degree
                        distd = sqrt((rx_c-rx_s)**2. + (ry_c-ry_s)**2. + (z_c-z_s)**2.)
                        dho = sqrt((rx_c-rx_s)**2. + (ry_c-ry_s)**2.)
                        # Computation of the zenithal angle between the source and the line of sight voxel.
                        anglezenithal(rx_s, ry_s, z_s, rx_c, ry_c, z_c, angzen)
                        # Computation of the angle azimutal direct line of sight-source
                        angleazimutal(rx_s,ry_s,rx_c,ry_c,angzen)
                        if (angzen > pi/4.):               # 45deg. it is unlikely to have a 1km high mountain less than 1
                            horizon(x_s,y_s,z_s,dx,dy,altsol,angazi,zhoriz,dh)
                            if (dh <= dho):
                                if (angzen < zhoriz):  # Shadow the path line of sight-source not below the horizon => we compute
                                    hh=1.
                                else:
                                    hh=0.
                            else
                                hh=1.
                        else
                            hh=1.

# Sub-grid obstacles
                        angmin = pi/2.-atan((altsol[x_s,y_s] + obsH[x_s,y_s]-z_s) / drefle[x_s,y_s])
                        if (angzen < angmin):             # Condition sub-grid obstacles direct.
                                ff=0.
                        else:
                            ff=ofill[x_s,y_s]

# Computation of the transmittance between the source and the line of sight
                        transmitm(angzen,z_s,z_c,distd,transm,tranam)
                        transmita(angzen,z_s,z_c,distd,transa,tranaa)

# Computation of the solid angle of the line of sight voxel seen from the source
                        omega = 1./distd**2.
                        if (omega > omemax):
                            omega=0.
                        anglez = round(180.*angzen/pi)+1
                        P_dir = pvalno[anglez,stype]

# Computation of the flux direct reaching the line of sight voxel
                        fldir=lamplu(x_s,y_s,stype)*P_dir*omega*transm*transa*(1.-ff)*hh   # Correction for obstacle filling factor

# Computation of the scattering probability of the direct light
# Distance pour traverser la cellule unitaire parfaitement orientée
                        if (omega != 0.):
                            angle3points(rx_s,ry_s,z_s,rx_c,ry_c,z_c,rx_obs,ry_obs,z_obs,angdif)    # Scattering angle.
                            diffusion(angdif,tranam,tranaa,secdif,un,fdifan,pdifdi,z_c) # Scattering probability of the direct light.
                        else:
                            pdifdi=0.

# Computation of the source contribution to the direct intensity toward the sensor by a line of sight voxel
                        intdir = fldir*pdifdi

# Contribution of the cloud reflexion of the light coming directly from the source
                        if (cloudt. != 0)                      #Line of sight voxel = cloud
                            if (cloudbase-z_c <= iz*scal):
                                anglezenithal(rx_c,ry_c,z_c,rx_obs,ry_obs,z_obs,azcl1)  # Zenith angle from cloud to observer
                                anglezenithal(rx_c,ry_c,z_c,rx_s,ry_s,z_s,azcl2)    # Zenith angle from source to cloud
                                doc2=(rx_c-rx_obs)**2.+(ry_c-ry_obs)**2.+(z_c-z_obs)**2.
                                dsc2=(rx_s-rx_c)**2.+(ry_s-ry_c)**2.+(z_s-z_c)**2.
                                cloudreflectance(angzen,cloudt,rcloud)      # Cloud intensity from direct illum
                                icloud=icloud+fldir/omega*rcloud*doc2*omefov*abs(cos(azcl2)/cos(azcl1))/dsc2/pi

                    else:
                        intdir=0.
# End of the case Position Source is not equal to the line of sight voxel position
# End of the computation of the direct intensity


# **********************************************************************************************************************
# * Computation of the scattered light toward the observer by a line of sight voxel lighted by the ground reflexion    *
# **********************************************************************************************************************

# etablissement of the conditions ands boucles
                    itotind=0.           # Initialisation of the reflected intensity of the source
                    itotrd=0.
                    boxx=round(reflsiz/dx)           # Number of column to consider left/right of the source for the reflexion.
                    boxy=round(reflsiz/dy)           # Number of column to consider up/down of the source for the reflexion.
                    xsrmi = x_s-boxx
                    if (xsrmi < 1):
                        xsrmi=1
                        xsrma=x_s+boxx          # Vérifier indentation
                    if (xsrma > nbx):
                        xsrma=nbx
                        ysrmi=y_s-boxy          # Vérifier indentation
                    if (ysrmi < 1):
                        ysrmi=1
                        ysrma=y_s+boxy          # Vérifier indentation
                    if (ysrma > nby):
                        ysrma=nby

                    # Vérifier indentation

                    for x_sr in range(xsrmi, xsrma):        # Beginning of the loop over the column (longitude) reflecting.
                        rx_sr=real(x_sr)*dx

                        for y_sr in range(ysrmi, ysrma):    # Beginning of the loop over the rows (latitu) reflecting.
                            ry_sr = y_sr*dy
                            irefl=0.
                            z_sr = altsol[x_sr,y_srz]
                            if((x_sr > nbx) or (x_sr < 1) or (y_sr > nby) or (y_sr < 1)):
                                if (verbose == 2):
                                    print('Ground cell out of borders')
                            else:
                                if((x_s == x_sr) and (y_s == y_sr) == (z_s == z_sr)):
                                    if (verbose == 2):
                                        print('Source pos = Ground cell')
                                else:
                                    # if haut is negative, the ground cell is lighted from below
                                    haut=-(rx_s-rx_sr)*tan(inclix[x_sr,y_sr])-(ry_s-ry_sr)*tan(incliy[x_sr,y_sr])+z_s-z_sr

# Computation of the zenithal angle between the source and the surface reflectance
                                    if (haut > 0.):              # Condition: the ground cell is lighted from above
                                        anglezenithal(rx_s,ry_s,z_s,rx_sr,ry_sr,z_sr,angzen)
                                        # end of the case "observer at the same latitu/longitude than the source".
                                        # End mais on reste dans le même if?

# Computation of the transmittance between the source and the ground surface
                                        distd=sqrt((rx_s-rx_sr)**2.+(ry_s-ry_sr)**2.+(z_s-z_sr)**2.)
                                        transmitm(angzen,z_s,z_sr,distd,transm,tranam)
                                        transmita(angzen,z_s,z_sr,distd,transa,tranaa)



# Computation of the solid angle of the reflecting cell seen from the source
                                        xc=dble[x_sr]*dble[dx]            # Position in meters of the observer voxel (longitude).
                                        yc=dble[y_sr]*dble[dy]            # Position in meters of the observer voxel (latitu).
                                        zc=dble[z_sr]                     # Position in meters of the observer voxel (altitude).
                                        xn=dble[x_s]*dble[dx]             # Position in meters of the source (longitude).
                                        yn=dble[y_s]*dble[dy]             # Position in meters of the source (latitu).
                                        zn=dble[z_s]                      # Position in meters of the source (altitude).
                                        epsilx=inclix[x_sr,y_sr]          # tilt along x of the ground reflectance
                                        epsily=incliy[x_sr,y_sr]          # tilt along x of the ground reflectance
                                        # use a sub-grid surface when the reflectance radius is smaller than the cell size
                                        if (dx > reflsiz):
                                            if ((x_sr == x_s) and (y_sr == y_s)):
                                                dxp=reflsiz
                                            else:
                                                dxp=dx
                                        else:
                                            dxp=dx

                                        if (dy > reflsiz):
                                            if ((x_sr == x_s) and (y_sr == y_s)):
                                                dyp=reflsiz
                                            else:
                                                dyp=dy
                                        else:
                                            dyp=dy

                                        r1x=xc-dxp/2.-xn            # computation of the composante along x of the first vector.
                                        r1y=yc+dyp/2.-yn            # computation of the composante along y of the first vector.
                                        r1z=zc-tan(epsilx)*dxp/2.+tan(epsily)*(dyp)/2.-zn  # computation of the composante en z of the first vector.
                                        r2x=xc+(dxp)/2.-xn            # computation of the composante along x of the second vector.
                                        r2y=yc+(dyp)/2.-yn            # computation of the composante along y of the second vector.
                                        r2z=zc+tan((epsilx))*(dxp)/2.+tan((epsily))*(dyp)/2.-zn    # computation of the composante en z of the second vector.
                                        r3x=xc-(dxp)/2.-xn            # computation of the composante along x of the third vector.
                                        r3y=yc-(dyp)/2.-yn            # computation of the composante along y of the third vector.
                                        r3z=zc-tan((epsilx))*(dxp)/2.-tan((epsily))*(dyp)/2.-zn   #  computation of the composante en z of the third vector.
                                        r4x=xc+(dxp)/2.-xn            # computation of the composante along x of the fourth vector.
                                        r4y=yc-(dyp)/2.-yn            # computation of the composante along y of the fourth vector.
                                        r4z=zc+tan((epsilx))*(dxp)/2.-tan((epsily))*(dyp)/2.-zn   # computation of the composante en z of the fourth vector.

                                        # Call of the routine anglesolide to compute the angle solide.
                                        anglesolide(omidiega,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y, r3z,r4x,r4y,r4z)

                                        # VÉRIFIER INDENTATION ???
                                        if (omega < 0.):
                                            raise ValueError("ERROR: Solid angle of the reflecting surface < 0.")

# Estimation of the half of the underlying angle of the solid angle       ! this angle servira a obtenir un meilleur isime (moyenne) of
# P_dir for le cas of grans solid angles the , pvalno varie significativement sur +- ouvang.
                                        ouvang = sqrt(omega/pi)             # Angle in radian.
                                        ouvang = ouvang*180./pi             # Angle in degrees.

# computation of the photometric function of the light fixture toward the reflection surface
#=======================================================================

                                        anglez=round(180.*angzen/pi)
                                        if (anglez < 0):
                                            anglez=-anglez
                                        if (anglez > 180):
                                            anglez=360-anglez
                                        # vérifier indentation
                                        anglez=anglez+1     # Transform the angle in integer degree into the position in the array.
                                                            # average +- ouvang

                                        naz=0
                                        nbang=0.
                                        P_indir=0.

                                        for na i range(-round(ouvang), round(ouvang)):
                                            naz=anglez+na
                                            if (naz < 0) naz=-naz
                                          if (naz > 181):     # symetric function
                                             naz=362-naz
                                          if (naz == 0):
                                              naz=1
                                          P_indir = P_indir+pvalno[naz,stype]*abs(sin(pi*(naz)/180.))/2.
                                          nbang = nbang+1.*abs(sin(pi*(naz)/180.))/2.

                                        P_indir = P_indir/nbang

# Computation of the flux reaching the reflecting surface
                                        flrefl = lamplu[x_s,y_s,stype]*P_indir*omega*transm*transa
# Computation of the reflected intensity leaving the ground surface
                                        irefl1 = flrefl*srefl/pi      # The factor 1/pi comes from the normalisation of the fonction

# Rendu à la ligne 1074
