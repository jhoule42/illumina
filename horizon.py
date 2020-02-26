"""======================================================================
  Routine horizon (Martin Aube 2017)


  Determine si la lumiere est bloquee par l'horizon

  pour utilisation avec Illumina
-------------------------------------------------------------------------

    Copyright (C) 2019  Martin Aube

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
-----------------------------------------------------------------------""""


import math
import numpy as np

def horizon(x, y, z, dx, dy, anga):

# Variables:

    width = 512
    altsol = np.zeros((width, width))       # matrice de 512x512 --> vérifier quoi à l'intérieur ?
    angaz1 = anga                           # angaz1 = (pi*anga)/180.  ?
    ix = (math.cos(angaz1))
    iy = (math.sin(angaz1))
    dist = width * math.sqrt(1. + math.tan(angaz1)**2.)
    scalef = dx/3.
    posx = x * dx       # dx -> taille du pixel en x
    posy = y * dy       # dy -> taille du pixel en y
    zhoriz = math.pi


    while (((posx <= (width * dx) and (posx > dx)) and ((posy <= (width * dy)) and (posy > dx)):

        posx = posx + ix * scalef
        posy = posy + iy * scalef
        nx = round(posx/dx)
        ny = round(posy/dy)

        if (altsol(nx, ny) > z):    # la matrice altsol n'a pas déja été défini ???

            zout = (math.pi/2) - (math.atan((altsol(nx, ny)-z)/sqrt(dx**2 * ((nx-x))**2 + dy**2 *(ny-y)**2)))
            d = sqrt(dx**2 * (nx-x)**2 + dy**2 *(ny-y)**2)

        else:
            zout = (math.pi/2) - 0.5 * math.pi/180
            d = (width)*dx                          # bug for zhoriz=pi, anyway in the real world pi is almost impossible

        if (zout < zhoriz):
            zhoriz = zout

return zhoriz, d