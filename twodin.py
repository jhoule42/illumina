"""===================================================================
Convertir les données sont forme binaire

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
-----------------------------------------------------------------------

Statut: demander questions de programmation Martin

-----------------------------------------------------------------------"""

import numpy as np

def twodin(nbx, nby, filename):

    width = 512
    bindata = np.zeros((width, width))

    with open(filename, 'r') as rfile:

        if (nbx > width) or (nby > width):
            raise ValueError("""You try to use a domain larger than the maximum allowed.
            Please restrict it to no more that 512 x 512. Your domain size is:"""
            nbx,'x',nby)

        for j in range(nby, 1, -1):
            for i in range(1, nbx):



# quoi retourner?
