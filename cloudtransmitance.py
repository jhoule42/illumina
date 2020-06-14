
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
Statut : Fonctionnel (petite marge d'erreur)
================================================================================"""

#--------------------------------------
# Types de nuages:                     |
#1 : thin Cirrus/Cirrostratus          |
#2 : thick Thick Cirrus/Cirrostratus   |
#3 : Altostratus/Altocumulus           |
#4 : Stratotcumuls & Stratus           |
#5 : Cumulus ou Cumulonimbus           |
#--------------------------------------

#   Fitted parameters for the cloud reflectance as a function of the incident zenith angle
#   rho(z)=a0+a1*cos z + a2 * cos^2 z + a3 * cos^3 z according to Shapiro 1982 Table 11

import numpy as np
from math import cos

def cloudreflectance(angzen, cloud):

    if (cloud < 0 or cloud > 4):
        raise ValueError("Error cloud type not define")

# (type de nuages, coefficiant de la formule)
    thocld = np.zeros((5, 4))

    thocld[0,0] = 0.63547
    thocld[0,1] = 0.35229
    thocld[0,2] = 0.08709
    thocld[0,3] = 0.22902
    thocld[1,0] = 0.26458
    thocld[1,1] = 0.66829
    thocld[1,2] = 0.24228
    thocld[1,3] = 0.49357
    thocld[2,0] = 0.19085
    thocld[2,1] = 0.32817
    thocld[2,2] = 0.08613
    thocld[2,3] = 0.08197
    thocld[3,0] = 0.13610
    thocld[3,1] = 0.29964
    thocld[3,2] = 0.14041
    thocld[3,3] = 0.00952
    thocld[4,0] = 0.17960
    thocld[4,1] = 0.34855
    thocld[4,2] = 0.14875
    thocld[4,3] = 0.01962


    tcloud = thocld[cloud,0] + thocld[cloud,1]*cos(angzen) + thocld[cloud,2]*(cos(angzen))**2.+thocld[cloud,3]*(cos(angzen))**3.
    print("rcloud:", tcloud, "angzen:", angzen, "cos(angzen)", cos(angzen))

print(cloudreflectance(1.15, 4))
