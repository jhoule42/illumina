 """===================================================================
 Calculation of the total vertical transmittance of the atmosphere
 (aerosols and molecules separately

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

Statut: À vérifier

-----------------------------------------------------------------------""""


def transtoa(lambm, taua, pressi):

    m = (pressi/101.3)  # masse d'air

    #  transmittance tiree de Kneizys et al. (1980)
    # trop de parenthèses dans verion O.G
    tranam = exp(-1. * m/(((lambm/1000.)**4.) * (115.6406-(lambm/1000.)**2.)))  # molecules
    tranaa=exp(-1.*taua)  # aerosols

return tranam, tranaa
