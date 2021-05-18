!    2021-05-17
!    pendulo_simple.f90
!    Felipe Ixcamparic (felipechoy1@gmail.com)

!    Este programa es una simulación del movimiento
!    de un péndulo simple utilizando el método de runge-Kutta
!    para dar con su posición angular como con su velocidad  
!

!    Codificación del texto: UTF8
!    Compiladores probados: GNU Fortran (WSl ubuntu Linux 20.04) 7.5.0
!    Instrucciones de compilación: no requiere nada mas
!    gfortran -Wall -pedantic  -c -o pendulo_simple.o pendulo_simple.f90
!    gfortran -o pendulo_simple.x pendulo_simple.o 
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


!Constantes otorgadas, incluyendo el step h
MODULE ctes
    IMPLICIT NONE
    REAL, PARAMETER    :: l=1 !longitud de la cuerda
    REAL, PARAMETER    :: g=9.8 !g
    REAL, PARAMETER    :: m=1  !masa
    REAL, PARAMETER    :: theta10=0.2 !Posición angular inicial
    REAL, PARAMETER    :: theta20=0.0 !Velocidad angular inicial
    !Constantes de muestreo
    REAL, PARAMETER    :: h=0.001 !step size
    REAL, PARAMETER    :: tfinal=5 !Tiempo final
    INTEGER,PARAMETER  :: N=CEILING(tfinal/h) !Número de iteraciones
END MODULE ctes



PROGRAM runge
    USE ctes
    IMPLICIT NONE
    REAL::t
    !Vectores de almacenamiento de datos
    REAL, dimension (N+1) :: o
    REAL, dimension (N+1) :: p
    !Constantes del método de Runge-Kutta
    REAL::k1theta1,k1theta2,k2theta1,k2theta2,k3theta1,k3theta2,k4theta1,k4theta2
    !Se llaman a las funciones de las derivadas
    REAL,EXTERNAL::theta1prim,theta2prim
    INTEGER::i !iteardor

    !Se almacenan las condiciones iniciales
    o(1)=theta10
    p(1)=theta20
    t=0
    !Se abre un archivo de texto
     OPEN(13,file="pendulo_simple.txt")
        DO i=1,N
            
            !Implementación del método de Runge-Kutta por definición
            k1theta1=theta1prim(t,           o(i)                ,p(i)         )
            k1theta2=theta2prim(t,           o(i)                ,p(i)          )
            k2theta1=theta1prim(t+ (h/2),    o(i)+(h/2)*(k1theta1)    ,p(i)+(h/2)*k1theta2 )
            k2theta2=theta2prim(t+ (h/2),    o(i)+(h/2)*(k1theta1)    ,p(i)+(h/2)*k1theta2 )
            k3theta1=theta1prim(t+ (h/2),    o(i)+(h/2)*(k2theta1)    ,p(i)+(h/2)*k2theta2 )
            k3theta2=theta2prim(t+ (h/2),    o(i)+(h/2)*(k2theta1)    ,p(i)+(h/2)*k2theta2 )
            k4theta1=theta1prim(t+ h,        o(i)+(h)*(k3theta1)      ,p(i)+(h)*k3theta2 )
            k4theta2=theta2prim(t+ h,        o(i)+(h)*(k3theta1)      ,p(i)+(h)*k3theta2 )
            !Se imprime la iteración 
            WRITE(13,*) t,";",o(i),";",p(i)

            !Se guarda la próxima iteración en los vectores
            o(i+1)=o(i)+(h/6)*(k1theta1+2*k2theta1+2*k3theta1+k4theta1)
            p(i+1)=p(i)+(h/6)*(k1theta2+2*k2theta2+2*k3theta2+k4theta2)
            !Se aumenta el tiempo de muestreo 
            t=t+h
            
        END DO
    CLOSE(13)   
   
   
    
END PROGRAM runge

!Función para theta1prima
REAL FUNCTION theta1prim(t,theta1,theta2)
USE ctes
REAL::t,theta1,theta2
theta1prim =  theta2
END FUNCTION theta1prim

!Función para theta2 prima 
REAL FUNCTION theta2prim(t,theta1,theta2)
USE ctes
REAL ::t,theta1,theta2
theta2prim = -(g/l)*sin(theta1)
END FUNCTION theta2prim