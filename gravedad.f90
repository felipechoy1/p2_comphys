!    2021-05-17
!    gravedad.f90
!    Felipe Ixcamparic (felipechoy1@gmail.com)

!    Este programa es una simulación del problema de los
!    n cuerpos con la Tierra, Luna y sol. 

!    Codificación del texto: UTF8
!    Compiladores probados: GNU Fortran (WSl ubuntu Linux 20.04) 7.5.0
!    Instrucciones de compilación: no requiere nada mas
!    gfortran -Wall -pedantic  -c -o gravedad.o gravedad.f90
!    gfortran -o gravedad.x gravedad.o 

!    Copyright (C) 2021 
!    Este código es de autoria totalmente original de 
!    Bryant Morazán(bryant.morazan@gmail.com) quien lo presentó en el foro perteneciente
!    al curso de Física Computacional en el año 2021.
!    el permitió el uso para realizar las simulaciones como 
!    pruebas de rendimiento de este mismo código.

!    Felipe Ixcamparic
!    felipechoy1@gmail.com

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



!----------------------------------------Inicia el programa--------------------------------------------
program gravedad
    !-----------------------Variables-------------------------------
        implicit none
        !-----Eje x----- 
        real(16), dimension(3) :: x0, x !Posiciones 
        real(16), dimension(3) :: vx0, vx !Velocidades
        real(16), dimension(3) :: ax !Aceleraciones
        !-----Eje y-----
        real(16), dimension(3) :: y0, y !Posiciones
        real(16), dimension(3) :: vy0, vy !Velocidades 
        real(16), dimension(3) :: ay !Aceleraciones
        !---Escalares---
        real(16), dimension(3) :: m !Masas
        real(16) :: mom_ang
        !---Distancia---
        real(16), dimension(3) :: r !Distancia
        !--Constantes---
        real(16) :: G = 6.674E-20
        !----Tiempo----
        real(16) :: t !Intervalo de tiempo
        !--Iteradores--
        integer :: iteraciones, i  
    !---------------------------------------------------------------
    !--------------------------Datos--------------------------------
    !------------Masas[kg]-------------
    m(1) = 1.989*10.0**30.0 !Sol
    m(2) = 5.972*10.0**24.0 !Tierra
    m(3) = 7.349*10.0**22.0 !Luna
    !---------------------------------
    !--------Posicion x[km]-----------
    x0(1) = 0.0 !Sol
    x0(2) = 152098290 !Tierra
    x0(3) = 152504990 !Luna
    !---------------------------------
    !--------Posicion y[km]-----------
    y0(1) = 0.0 !Sol
    y0(2) = 0.0 !Tierra
    y0(3) = 0.0 !Luna
    !---------------------------------
    !--------Velocidad x[km/s]--------
    vx0(1) = 0.0 !Sol
    vx0(2) = 0.0 !Tierra
    vx0(3) = 0.0 !Luna
    !---------------------------------
    !--------Velocidad y[km/s]--------
    vy0(1) = 0 !Sol
    vy0(2) = 28.76 !Tierra
    vy0(3) = 29.72 !Luna
    !---------------------------------
    !-------Tiempo e iteraciones------
    t = 20
    iteraciones = 3000000
    !---------------------------------
    !--------------------------------------------Bloque principal---------------------------------------------------
    !------Arhivos de resultados-------
    open(1, file='tierra.txt')
    open(4, file='tierra2.txt')
    open(5, file='tierra3.txt')
    open(6, file='tierra4.txt')
    open(2, file='tierra-luna.txt')
    open(3, file='momentum_angular.txt')
    !----------------------------------
    do i=0, iteraciones
    !-----------------------------Distancias------------------------------------- 
        r(1) = ( (x0(1) - x0(2))**2 + (y0(1) - y0(2))**2 )**0.5 !Entre 1 y 2
        r(2) = ( (x0(1) - x0(3))**2 + (y0(1) - y0(3))**2 )**0.5 !Entre 1 y 3
        r(3) = ( (x0(2) - x0(3))**2 + (y0(2) - y0(3))**2 )**0.5 !Entre 2 y 3
    !----------------------------------------------------------------------------
    !----------------------------Aceleraciones-----------------------------------
    !-------------Eje x-----------------
        ax(1) = -G*( (x0(1) - x0(2))*m(2)/r(1)**3 + (x0(1) - x0(3))*m(3)/r(2)**3 ) 
        ax(2) = -G*( (x0(2) - x0(1))*m(1)/r(1)**3 + (x0(2) - x0(3))*m(3)/r(3)**3 ) 
        ax(3) = -G*( (x0(3) - x0(1))*m(1)/r(2)**3 + (x0(3) - x0(2))*m(2)/r(3)**3 ) 
    !-----------------------------------
    !-------------Eje y-----------------
        ay(1) = -G*( (y0(1) - y0(2))*m(2)/r(1)**3 + (y0(1) - y0(3))*m(3)/r(2)**3 ) 
        ay(2) = -G*( (y0(2) - y0(1))*m(1)/r(1)**3 + (y0(2) - y0(3))*m(3)/r(3)**3 ) 
        ay(3) = -G*( (y0(3) - y0(1))*m(1)/r(2)**3 + (y0(3) - y0(2))*m(2)/r(3)**3 )
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
         mom_ang = m(1)*( x0(1)*vy0(1) - y(1)*vx0(1) ) + m(2)*( x0(2)*vy0(2) - y(2)*vx0(2) ) + m(3)*( x0(3)*vy0(3) - y(3)*vx0(3) )
    !-----------------------------------------------------------------------------
    !-------------------------Escritura de resultados-----------------------------
    !---------Orbita terrestre-----------
        write(1,*) -x(2), y(2)
        write(4,*) x(2), y(2)
        write(5,*) -x(2), -y(2)
        write(6,*) x(2), -y(2)
    !-----Orbita lunar y terrestre-------
        if (mod(i,10000) == 0) then
        write(2,*) x0(1),';', y0(1),';', x0(2),';', y0(2),';', x0(3) + 70*( x0(3) - x0(2) ) ,';', y0(3) + 70*( y0(3) - y0(2) )
        end if 
    !----------Momentum angular----------
        write(3,*) mom_ang
    !-------------------------------------
    end do
    !-------------------------------------------------------------------------------------------------------
    close(1)
    close(2)
    !-------------------------------------------------------------------------------------------------------
    end program gravedad
    !-------------------------------------------------------------------------------------------------------