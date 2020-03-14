"""=======================================================================
  Routine angleazimutal (Martin Aube 2010)


  Determine l'angle azimutal entre les points (x1,y1,z1) et (x2,y2,z2)
  Retourne l'angle anglezen en radians

  pour utilisation avec Illumina
-----------------------------------------------------------------------

    Copyright (C) 2010  Martin Aube

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but without any warranry; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: martin.aube@cegepsherbrooke.qc.ca

=======================================================================

** Fonctionnel --> comparer à .f

======================================================================="""

from math import atan, pi

# rotation selon le plan surface autour de y
def angleazimutal(x1, y1, x2, y2):

    if (x2 - x1 != 0.):
        angazi = abs(atan((y2-y1)/(x2-x1)))   

    if ((x2 -x1 == 0.) and (y2 -y1 == 0.)):
        angazi = 0.

    else:
        if (x2 - x1 > 0.):
            if (y2 - y1 < 0.):
                angazi = 2.* (pi-angazi)

        elif (x2 - x1 < 0.):
            if (y2 - y1 < 0.):
                angazi = angazi + pi
            else:
                angazi = pi - angazi

        else:
            if (y2 > y1):
                angazi = pi/2
            if (y2 < y1):
                angazi = 3.* pi/2

        if ((angazi < 0.) or (angazi > 2.*pi)):  # doit être entre 0 et 2pi
            print("'ERREUR angazi = ", angazi, x1, y2, x2, y2)

    return angazi

print(angleazimutal(10, 12, 20, 21))
