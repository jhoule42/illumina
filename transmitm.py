"""=======================================================================
 Routine transmitm (Andre Morin, Alex Neron, Etienne Rousseau 2004)
 debuggee par Martin Aube 2004
 Determine la transmittance des molecules atmospheriques sur un parcours
 entre les cellules (x_i,y_i,z_i) et (x_f,y_f,z_f)
 Recoit des longueurs d'ondes en nanometre et les transforme en microns.
 La pressi doit etre en KPa
 Retourne la transmittance transm

  *** J'ai valide le calcul zenith tout atm avec modtran et
      cela concorde M. Aubé mars 2010

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

** Demander explications Martin

======================================================================="""

import math

def transmitm(angz,z_i,z_f,distd,transm,tranam):

    if (z_i > z_f):         # Vérifier code fortran, car certaine variables ne semble pas être utilisé
        z2 = z_1
        z1 = z_f

    else:
        z1 = z_i
        z2 = z_f

    if (z1 != z2):                                                                # vérifier signe - + dans .f
        transm = math.exp((math.log(tranam) / abs(math.cos(angz)))*(math.exp(-1. * z1/8000.) - exp(-1. *z2/8000.)))

    else:
        transm = math.exp((math.log(tranam)) * math.exp(-1. * z1/8000.) * distd)

    if ((transm < 0.) or (transm > 1)):
        print("ERREUR avec transm", transm, tranam, z_f, z_i, distd, angz, airm)     # vérifier print*

# return    vérifier quoi retourner --> pas clair
