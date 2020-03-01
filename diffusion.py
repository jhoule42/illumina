"""------------------------------------------------------------------------

=======================================================================
  Routine diffusion (Alex Neron 2004) (Modifie par Martin Aube, Valerie Houle, Philippe Robert-Staehler)

  Determine la probabilite de diffusion de la lumiere par unite d'angle solide
  dans la direction angdif. Les parametres de diffusion sont donnes
  par secdif (rapport de la section efficace de diffusion sur la section
  efficace d'extinction totale), fonc_anorm (fonction de diffusion
  normalisees des aerosols)
  Retourne la probabilite de diffusion pdif

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
======================================================================="""


import math
from math import pi


def diffusion(angdif, tranam, tranaa, un, secdif, fonc_a, altit):

#--------------------------------------------------------
#  Calcul et normalisation des fonctions de diffusion
#--------------------------------------------------------

    if (angdif < 0.):
        angdif = -angdif   # c'est * -1 ?
    if (angdif - pi > 0.00001):
        angdif = pi
    angdeg = ((angdif*180.)/pi)
    rang = angdeg + 1

#=======================================================================
#        Calcul de la fonction d'emission de la source vers la cible
#=======================================================================

    fonc_ae = fonc_a(rang)       # fonc_ae est une liste de 181 éléments ?
    fctmol = 0.75 * (1.+((math.cos()(angdif))**2.)) / (4.*pi)

#--------------------------------------------------------------------
#  Calcul des probabilites de diffusion par unite d'angle solide
#-----------------------------------------------------------------------

# Les fonctions utilisees ici sont deja normalisees
    prob_a = (1.-exp(log(tranaa) * exp(-1.*altit/2000.) * un/2000.))*secdif * fonc_ae
    prob_m=(1.-exp(log(tranam)*exp(-1.*altit/8000.)*un/8000.))* fctmol

    # Fonc_ae normalisee dans le MAIN, fctmol dans la routine (voir la division par 4 pi)


    pdif = prob_a+prob_m                  # Ce calcul est approximatif et bon seulement si 1-transa et  1-transm sont tres petits.


    if (prob_a > 1.):
        raise ValueError("a > 1")

    if (prob_a < 0.):
        raise ValueError("a < 0")

    if (prob_m > 1.):
        raise ValueError("m > 1")

    if (prob_m < 0.):
        raise ValueError("m < 0")       # caractères étranges dans .f

    if (pdif > 1.):
        raise ValueError("pdif > 1")

    if (pdif < 0.):
        raise ValueError("pdif < 0")

    # vérifier le commentaire dans .f

    return pdif

print(diffusion(1, 0.32, 0.12, 1,32, 5, 2))
