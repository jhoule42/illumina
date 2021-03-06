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
=======================================================================
Statut: Demander questions a Martin --> fonc_a semble tjrs vide
======================================================================="""


import numpy as np
from math import pi, cos, exp, log

def diffusion(angdif, tranam, tranaa, un, secdif, fonc_a, altit):

    """ Calcul et normalisation des fonctions de diffusion"""


    if (angdif < 0.):
        angdif = -angdif
    if (angdif - pi > 0.00001):
        angdif = pi
    angdeg = ((angdif*180.)/pi)
    rang = int(angdeg + 1)

#=======================================================================
#        Calcul de la fonction d'emission de la source vers la cible
#=======================================================================

    fonc_a = np.zeros((180))
    print("fonc_a", fonc_a)
    fonc_ae = fonc_a[rang]       # fonc_ae est une liste de 181 elements ou de size(rang)
    print("fonc_ae", fonc_ae)
    fctmol = 0.75 * (1.+((cos((angdif))**2.)) / (4.*pi))

#---------------------------------------------------------------------
#  Calcul des probabilites de diffusion par unite d'angle solide
#---------------------------------------------------------------------

    # Les fonctions utilisees ici sont deja normalisees
    prob_a = (1.-exp(log(tranaa) * exp(-1.*altit/2000.) * un/2000.))*secdif * fonc_ae
    print("prob_a", prob_a)
    prob_m = (1.-exp(log(tranam) * exp(-1.*altit/8000.) * un/8000.))* fctmol
    print("prob_m", prob_m)

    # Fonc_ae normalisee dans le MAIN, fctmol dans la routine (voir la division par 4 pi)
    pdif = prob_a+prob_m                  # Ce calcul est approximatif et bon seulement si 1-transa et  1-transm sont tres petits.



    if (prob_a > 1.):
        raise ValueError("a > 1")

    if (prob_a < 0.):
        raise ValueError("a < 0")

    if (prob_m > 1.):
        raise ValueError("m > 1")

    if (prob_m < 0.):
        raise ValueError("m < 0")

    if (pdif > 1.):
        raise ValueError("pdif > 1", pdif,prob_a,prob_m)

    if (pdif < 0.):
        raise ValueError("pdif < 0", pdif,prob_a,prob_m)

    return pdif

print(diffusion(2.0, 0.62, 0.22, 1.0 , 0.05, 0.65, 242.0))
