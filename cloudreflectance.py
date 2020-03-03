"""================================================================================
    Copyright (C) 2015 Martin Aube

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

 **      http://cegepsherbrooke.qc.ca/~aubema/index.php/Prof/IllumEn?action=download&upname=intensite_lumineuse.pdf  **

================================================================================

Statut : Fonctionnel --> comparer à .f

================================================================================"""

#--------------------------------------
# Types de nuages:                     |
#1 : thin Cirrus/Cirrostratus          |
#2 : thick Thick Cirrus/Cirrostratus   |
#3 : Altostratus/Altocumulus           |
#4 : Stratotcumuls & Stratus           |
#5 : Cumulus ou Cumulonimbus           |
#--------------------------------------

#   fitted parameters for the cloud reflectance as a function of the incident zenith angle
#   rho(z)=a0+a1*cos z + a2 * cos^2 z + a3 * cos^3 z according to Shapiro 1982 Table 10

import numpy as np
from math import cos

def cloudreflectance(angzen, cloud):

    if (cloud < 0 or cloud > 4):
        raise ValueError("Error cloud type not define")

# (type de nuages, coefficiant de la formule)
    rhocld = np.zeros((5, 4))

    rhocld[0,0]=0.25674
    rhocld[0,1]=-0.18077
    rhocld[0,2]=-0.21961
    rhocld[0,3]=0.252724
    rhocld[1,0]=0.60540
    rhocld[1,1]=-0.55142
    rhocld[1,2]=-0.23389
    rhocld[1,3]=0.43648
    rhocld[2,0]=0.66152
    rhocld[2,1]=-0.14863
    rhocld[2,2]=-0.08193
    rhocld[2,3]=0.13442
    rhocld[3,0]=0.71214
    rhocld[3,1]=-0.15033
    rhocld[3,2]=0.00696
    rhocld[3,3]=0.03904
    rhocld[4,0]=0.67072
    rhocld[4,1]=-0.13805
    rhocld[4,2]=-0.10895
    rhocld[4,3]=0.09460

    rcloud = rhocld[cloud,0] + rhocld[cloud,1]*cos(angzen) + rhocld[cloud,2]*(cos(angzen))**2.+rhocld[cloud,3]*(cos(angzen))**3.
    print("rcloud:", rcloud, "angzen:", angzen, "cos(angzen)", cos(angzen))


print(cloudreflectance(1.52, 0))
