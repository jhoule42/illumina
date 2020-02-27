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
import math

def zone_diffusion(x1, y1, z1, x2, y2, z2):

    zondif = np.zeros((3000000, 3))     # pk une matrice de cette grandeur ?

    neffet = round(effet/siz)   # arrondir à combien ?
    dmin = math.sqrt((x1-x2)**2.+(y1-y2)**2.+(z1-z2)**2.)

#       find an approximate value to stepdi
#       stepdi=nint((dmin+effet)*3.14159/siz)*neffet/n2nd*neffet/2
    stepdi = 90000000
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
    if (z2 > cloudbase):            # que représente cloudbase --> pas déclaré en haut
        kmax = round(cloudbase/siz)

    ncell = 0
    keep = 0

    
