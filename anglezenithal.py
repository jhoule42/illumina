"""=======================================================================
  Routine anglezenithal (Andre Morin 2004)
  debugge par Martin Aube 2004 (cette routine ne calculait pas du tout
  l'angle zenithal
  Determine l'angle zenithal entre les points (x1,y1,z1) et (x2,y2,z2)
  Retourne l'angle angzen en radians.
   Angle zenithal == angle a partir du haut
  pour utilisation avec Illumina
-----------------------------------------------------------------------
    Copyright (C) 2009  Martin Aube
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
** Fonctionnel
======================================================================="""


import math

def anglezenithal(x1, y1, z1, x2, y2, z2):

    hdist = math.sqrt((x2-x1)**2. + (y2-y1)**2.)    # pythagore

    if (z2 - z1 != 0):
        angzen = math.atan(hdist/abs(z2-z1))

        if (z2 - z1 < 0):
            angzen = math.pi - angzen   # Rotation de pi pour etre dans le domaine

    else:
        angzen = math.pi / 2.

    if ((angzen < 0) or (angzen > math.pi)):    # 0 pointe vers le haut et pi vers le bas
        print("ERREUR angzen2=", angzen)

    return angzen


print(anglezenithal(43,32,31,48,54,64))
