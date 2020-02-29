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



======================================================================="""

import numpy as np
from math import sqrt, pi

def zone_diffusion(x1, y1, z1, x2, y2, z2, effet, alts, cloudbase, siz):       # erreur step pas utilisé + vérifier ncell

    zondif = np.zeros((3000000, 3))     # pk une matrice de cette grandeur ?

    neffet = round(effet/siz)   # arrondir à combien ?
    dmin = math.sqrt((x1-x2)**2.+(y1-y2)**2.+(z1-z2)**2.)

#    find an approximate value to stepdi
#    stepdi = 90000000

    stepdi=round((dmin+effet)* pi/siz) * neffet/n2nd*neffet/2      # importé fromnumeric pour round?

    if (stepdi = 0):
        stepdi=1

    step = round(stepdi**(1./3.))   # pourquoi stepdi real ?
    if (step <= 0):
        step=1

    # limits of the calculations loops
    x_1 = round(x1/siz)
    y_1 = round(y1/siz)
    x_2 = round(x2/siz)
    y_2 = round(y2/siz)
    z_2 = round(z2/siz)+1

    if (x_1 < x_2):
        imin = x_1 - neffet
        imax = x_2 + neffet
    else:
        imin = x_2 - neffet
        imax = x_1 + neffet

    if (y_1 < y_2):
        jmin = y_1 - neffet
        jmax = y_2 + neffet
    else:
        jmi = y_2 - neffet
        jmax = y_1 + neffet

    kmax=z_2+neffet
    if (z2 > cloudbase):            # que représente cloudbase? --> pas déclaré en haut
        kmax = round(cloudbase/siz)

    ncell = 0       # 10 dans le .f?
    keep = 0


    for i in range(imin, imax, step):
        for j in range(jmin, jmax, step):
            for k in range(1, kmax, step):
                x0 = i*siz
                y0 = j*six      # pourquoi écrire real? --> redéfinition des int ?
                z0 = k*siz
                d1 = sqrt((x1-x0)**2.+(y1-y0)**2.+(z1-z0)**2.)
                d2 = sqrt((x2-x0)**2.+(y2-y0)**2.+(z2-z0)**2.)
                d = d1+d2

                if ((z0. > alts) and (z0 < 35000.)):
                    if (d < dmin+2. *effet):
                        ncell += 1
                        zondif(ncell,1) = x0      # vérifier syntaxe numpy
                        zondif(ncell,2) = y0
                        zondif(ncell,3) = z0
    stepdi=step**3

return zondif, stepdi,
