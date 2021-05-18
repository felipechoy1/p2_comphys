!    2021-05-17
!    electro.f90
!    Felipe Ixcamparic (felipechoy1@gmail.com)

!    Este es un programa que resuelve el problema de 3 cuerpos para 4 particulas cargadas 
!    además calcula el momentum angular además de graficar las trayectorias de las particulas.

!    Codificación del texto: UTF8
!    Compiladores probados: GNU Fortran (WSl ubuntu Linux 20.04) 7.5.0
!    Instrucciones de compilación: no requiere nada mas
!    gfortran -Wall -pedantic  -c -o electro.o electro.f90
!    gfortran -o electro.x electro.o 
!
!
!    !    Copyright (C) 2021 
!    Este código es de autoria totalmente original de 
!    Bryant Morazán(bryant.morazan@gmail.com) quien lo presentó en el foro perteneciente
!    al curso de Física Computacional en el año 2021.
!    el permitió el uso para realizar las simulaciones como 
!    pruebas de rendimiento de este mismo código.
!
!    This program is free software: you can redistribute it and/or
!    modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of
!    the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see
!    <http://www.gnu.org/licenses/>.
!


!-----------------------------------------------------------------------------------------------------
!----------------------------------------Inicia el programa--------------------------------------------
program electro
!-----------------------Variables-------------------------------
    implicit none
    !-----Eje x----- 
    real(16), dimension(4) :: x0, x !Posiciones 
    real(16), dimension(4) :: vx0, vx !Velocidades
    real(16), dimension(4) :: ax !Aceleraciones
    !-----Eje y-----
    real(16), dimension(4) :: y0, y !Posiciones
    real(16), dimension(4) :: vy0, vy !Velocidades 
    real(16), dimension(4) :: ay !Aceleraciones
    !---Escalares---
    real(16), dimension(4) :: q, m !Masas y cargas
    real(16) :: mom_ang
    !---Distancia---
    real(16), dimension(6) :: r !Distancia
    !--Constantes---
    real(16) :: G = 6.674E-11
    real(16) :: k = 8.9874E+9
    !----Tiempo----
    real(16) :: t !Intervalo de tiempo
    !--Iteradores--
    integer :: iteraciones, i  
!---------------------------------------------------------------
!--------------------------Datos--------------------------------
!------------Masas[kg]-------------
m(1) = 9.11E-31 
m(2) = 9.11E-31 
m(3) = 9.11E-31 
m(4) = 9.11E-31
!------------Cargas[C]-------------
q(1) = 1.602E-19
q(2) = 1.602E-19
q(3) = 1.602E-19
q(4) = 1.602E-19
!--------Posicion x[m]-----------
x0(1) = 0.0 
x0(2) = 5 
x0(3) = 0.0
x0(4) = 5 
!---------------------------------
!--------Posicion y[km]-----------
y0(1) = 0.0 
y0(2) = 0.0 
y0(3) = 1
y0(4) = 1
!---------------------------------
!--------Velocidad x[km/s]--------
vx0(1) = 0.0
vx0(2) = 0.0
vx0(3) = 0.0
vx0(4) = 0.0
!---------------------------------
!--------Velocidad y[km/s]--------
vy0(1) = 0.0
vy0(2) = 0.0
vy0(3) = 0.0
vy(4) = 0.0
!---------------------------------
!-------Tiempo e iteraciones------
t = 1
iteraciones = 3e6
!---------------------------------
!--------------------------------------------Bloque principal---------------------------------------------------
!------Arhivos de resultados-------
open(1, file='movimiento.txt')
open(2, file='momentum_angular.txt')
!----------------------------------
do i=0, iteraciones
!-----------------------------Distancias------------------------------------- 
    r(1) = ( (x0(1) - x0(2))**2 + (y0(1) - y0(2))**2 )**0.5 !Entre 1 y 2
    r(2) = ( (x0(1) - x0(3))**2 + (y0(1) - y0(3))**2 )**0.5 !Entre 1 y 3
    r(3) = ( (x0(1) - x0(4))**2 + (y0(1) - y0(4))**2 )**0.5 !Entre 1 y 4
    r(4) = ( (x0(2) - x0(3))**2 + (y0(2) - y0(3))**2 )**0.5 !Entre 2 y 3
    r(5) = ( (x0(2) - x0(4))**2 + (y0(2) - y0(4))**2 )**0.5 !Entre 2 y 4
    r(6) = ( (x0(3) - x0(4))**2 + (y0(3) - y0(4))**2 )**0.5 !Entre 3 y 4
!----------------------------------------------------------------------------
!----------------------------Aceleraciones-----------------------------------
!-------------Eje x-----------------
    ax(1) = (x0(1) - x0(2))*( (k*q(1)*q(2))/m(1) - G*m(2) )/r(1)**3 + (x0(1) - x0(3))*( (k*q(1)*q(3))/m(1)&
    & - G*m(3) )/r(2)**3 + (x0(1) - x0(4))*( (k*q(1)*q(4))/m(1) - G*m(4) )/r(3)**3
    ax(2) = (x0(2) - x0(1))*( (k*q(2)*q(1))/m(2) - G*m(1) )/r(1)**3 + (x0(2) - x0(3))*( (k*q(2)*q(3))/m(2)&
    & - G*m(3) )/r(4)**3 + (x0(2) - x0(4))*( (k*q(2)*q(4))/m(2) - G*m(4) )/r(5)**3
    ax(3) = (x0(3) - x0(1))*( (k*q(3)*q(1))/m(3) - G*m(1) )/r(2)**3 + (x0(3) - x0(2))*( (k*q(3)*q(2))/m(3)&
    & - G*m(2) )/r(4)**3 + (x0(3) - x0(4))*( (k*q(3)*q(4))/m(3) - G*m(4) )/r(6)**3
    ax(4) = (x0(4) - x0(1))*( (k*q(4)*q(1))/m(4) - G*m(1) )/r(3)**3 + (x0(4) - x0(2))*( (k*q(4)*q(2))/m(4)&
    & - G*m(2) )/r(5)**3 + (x0(4) - x0(3))*( (k*q(4)*q(3))/m(4) - G*m(3) )/r(6)**3 
!-----------------------------------
!-------------Eje y-----------------
    ay(1) = (y0(1) - y0(2))*( (k*q(1)*q(2))/m(1) - G*m(2) )/r(1)**3 + (y0(1) - y0(3))*( (k*q(1)*q(3))/m(1)&
    & - G*m(3) )/r(2)**3 + (y0(1) - y0(4))*( (k*q(1)*q(4))/m(1) - G*m(4) )/r(3)**3
    ay(2) = (y0(2) - y0(1))*( (k*q(2)*q(1))/m(2) - G*m(1) )/r(1)**3 + (y0(2) - y0(3))*( (k*q(2)*q(3))/m(2)&
    & - G*m(3) )/r(4)**3 + (y0(2) - y0(4))*( (k*q(2)*q(4))/m(2) - G*m(4) )/r(5)**3
    ay(3) = (y0(3) - y0(1))*( (k*q(3)*q(1))/m(3) - G*m(1) )/r(2)**3 + (y0(3) - y0(2))*( (k*q(3)*q(2))/m(3)&
    & - G*m(2) )/r(4)**3 + (y0(3) - y0(4))*( (k*q(3)*q(4))/m(3) - G*m(4) )/r(6)**3
    ay(4) = (y0(4) - y0(1))*( (k*q(4)*q(1))/m(4) - G*m(1) )/r(3)**3 + (y0(4) - y0(2))*( (k*q(4)*q(2))/m(4)&
    & - G*m(2) )/r(5)**3 + (y0(4) - y0(3))*( (k*q(4)*q(3))/m(4) - G*m(3) )/r(6)**3
!-----------------------------------
!-----------------------------------------------------------------------------
!---------------------------Velocidades---------------------------------------
!-----------Eje x-------------------
    vx = vx0 + t*ax
!-----------Eje y-------------------
    vy = vy0 + t*ay
!-----------------------------------------------------------------------------
!------------------------------Posiciones-------------------------------------
!-------------Eje x---------------
    x = x0 + t*vx0 + (t**2)*ax/2
!-------------Eje y---------------
    y = y0 + t*vy0 + (t**2)*ay/2
!-----------------------------------------------------------------------------
!-------------------------Reasignacion de valores-----------------------------
!-----------Eje x--------------
    x0 = x
    vx0 = vx
!-----------Eje y--------------
    y0 = y
    vy0 = vy
!-----------------------------------------------------------------------------
!---------------------------Momentum angular----------------------------------
     mom_ang = m(1)*( x0(1)*vy0(1) - y(1)*vx0(1) ) + m(2)*( x0(2)*vy0(2) - y(2)*vx0(2) ) +&
     & m(3)*( x0(3)*vy0(3) - y(3)*vx0(3) )+ m(4)*( x0(4)*vy0(4) - y(4)*vx0(4) )
!-----------------------------------------------------------------------------
!-------------------------Escritura de resultados-----------------------------
!-------------Trayectorias-----------
    if (mod(i,10000) == 0) then
    write(1,*) x(1),';',y(1),';',x(2),';',y(2),';', x(3),';',y(3),';', x(4),';',y(4)
    end if 
!----------Momentum angular----------
    write(2,*) mom_ang
!-------------------------------------
end do
!-------------------------------------------------------------------------------------------------------
close(1)
close(2)
!-------------------------------------------------------------------------------------------------------
end program electro
!-------------------------------------------------------------------------------------------------------