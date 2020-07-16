"""=======================================================================
  Routine zone diffusion 2010
  Determine les cellules se trouvant dans la zone de diffusion (effet)
  des cellules (x1,y1,z1) et (x2,y2,z2)
  Retourne la matrice des cellules diffusantes (diffusion) ainsi que le nombre
  de cellules diffusantes (ncellule)
------------------------------------------------------------------------
    Copyright (C) 2010  Martin Aube
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
=======================================================================
Statut:  demander explications Martin
======================================================================="""

import numpy as np
from math import sqrt, pi

def zone_diffusion(effet, siz):       # erreur dstep pas utilisé + vérifier ncell
# car inintialisé à zéro

    zondif = np.zeros((3000000, 3))     # pk une matrice de cette grandeur ?

    neffet = round(effet/siz)   # arrondir à combien ?
    dmin = effet
    stepdi = 1

    # limits of the calculations loops
    imin = -neffet
    imax = neffet
    jmin = -neffet
    jmax = neffet
    kmin = -neffet
    kmax = neffet
    ncell = 0   # 10 dans le .f?
    keep = 0


    for i in range(imin, imax):
        x0 = i*siz
        for j in range(jmin, jmax):   # pk ERREUR?
            y0 = j*siz
            for k in range(1, kmax):
                z0 = k*siz

                d = sqrt((x0)**2.+(y0)**2.+(z0)**2.)

                if (d <= dmin):
                    keep += 1
                    if (keep == stepdi):

                        keep = 0
                        ncell = ncell+1
                        zondif[ncell,0] = x0
                        zondif[ncell,1] = y0
                        zondif[ncell,2] = z0

    return zondif   # retourner juste stepdi?

print(zone_diffusion(404, 32))
