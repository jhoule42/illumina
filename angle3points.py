"""-----------------------------------------------------------------------
=======================================================================
  Routine angle3points (Andre Morin 2004)
  debuggee et modifiee par Martin Aube 2004
  Determine l'angle entre 3 points (x1,y1,z1), (x2,y2,z2) et (x3,y3,z3)
  dont le sommet est au point 2
  Retourne l'angle en radians
  pour utilisation avec Illumina
-----------------------------------------------------------------------
    Copyright (C) 2009  Martin Aube
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but without amy warranty; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    Contact: martin.aube@cegepsherbrooke.qc.ca
=======================================================================
Statut : Fonctionnel (petite marge erreur)
======================================================================="""

from math import sqrt, pi, acos

def angle3points(x1,y1,z1,x2,y2,z2,x3,y3,z3):

    xu=x2-x1         # Voici les composantes du vecteur u.
    yu=y2-y1
    zu=z2-z1

    xv=x3-x2         # Voici les composantes du vecteur v.
    yv=y3-y2
    zv=z3-z2


    if (xv == 0.) and (yv == 0) and (zv == 0):
        raise ValueError("Erreur vecteur de sortie nul")
        print(x1,y1,z1,x2,y2,z2,x3,y3,z3)

    if (xu == 0) and (yu == 0) and (zu == 0):
        raise ValueError("Erreur vecteur d'entree nul")
        print(x1,y1,z1,x2,y2,z2,x3,y3,z3)

    argume = ((xu*xv) + (yu*yv) + (zu*zv)) / (sqrt(xu**2.+ yu**2.+ zu**2.)* sqrt(xv**2.+yv**2.+zv**2.))

    if (argume >= 1):
        an3pts = 0.

    elif (argume <= -1):
        an3pts = pi

    else:
        an3pts = acos(argume)

    if (an3pts < 0):
        ValueError("Error: an3pts < 0")

    return an3pts

print(angle3points(2.0, 4.0, 3.0, 2.0, 3.0, 16.0, 9.0, 1.0, 7.0))
