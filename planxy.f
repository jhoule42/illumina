c-----------------------------------------------------------------------
c
c=======================================================================
c  Routine planxy (Martin Aube 2010)
c
c
c  pour utilisation avec Illumina 
c-----------------------------------------------------------------------
c   
c    Copyright (C) 2010  Martin Aube
c
c    This program is free software: you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation, either version 3 of the License, or
c    (at your option) any later version.
c
c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.
c
c    Contact: martin.aube@cegepsherbrooke.qc.ca
c
c
      subroutine planxy(dx1,dy1,xc,xn,yc,yn,zc,zn,
     +r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
      real*8 r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z              ! Composantes des vecteurs utilises dans la routine angle solide.     
      real*8 xc,yc,zc,xn,yn,zn                                            ! Position (metre) des elements (arrivee, depart) pour le calcul       
      real dx1,dy1
        
                        r1x=xc-dble(dx1)/2.-xn                             ! Calcul de la composante en x du premier vecteur.
                        r1y=yc+dble(dy1)/2.-yn                             ! Calcul de la composante en y du premier vecteur.
                        r1z=zc-zn                                         ! Calcul de la composante en z du premier vecteur.
                        r2x=xc+dble(dx1)/2.-xn                             ! Calcul de la composante en x du deuxieme vecteur.
                        r2y=yc+dble(dy1)/2.-yn                             ! Calcul de la composante en y du deuxieme vecteur.
                        r2z=zc-zn                                         ! Calcul de la composante en z du deuxieme vecteur.
                        r3x=xc-dble(dx1)/2.-xn                             ! Calcul de la composante en x du troisieme vecteur.
                        r3y=yc-dble(dy1)/2.-yn                             ! Calcul de la composante en y du troisieme vecteur.
                        r3z=zc-zn                                         ! Calcul de la composante en z du troisieme vecteur.
                        r4x=xc+dble(dx1)/2.-xn                             ! Calcul de la composante en x du quatrieme vecteur.
                        r4y=yc-dble(dy1)/2.-yn                             ! Calcul de la composante en y du quatrieme vecteur.
                        r4z=zc-zn                                         ! Calcul de la composante en z du quatrieme vecteur.
      return
      end

