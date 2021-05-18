!    2021-05-17
!    pendsimeuler.f90
!    Felipe Ixcamparic (felipechoy1@gmail.com)

!    Este programa es una simulación del movimiento
!    de un péndulo simple utilizando el método de Euler
!    para dar con su posición angular como con su velocidad  
!

!    Codificación del texto: UTF8
!    Compiladores probados: GNU Fortran (WSl ubuntu Linux 20.04) 7.5.0
!    Instrucciones de compilación: no requiere nada mas
!    gfortran -Wall -pedantic  -c -o pendsimeuler.o pendsimeuler.f90
!    gfortran -o pendsimeuler.x pendsimeuler.o 
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


PROGRAM psimple
    USE ctes
    IMPLICIT NONE
    !Vectores para almacenar posición y velocidad
    REAL, dimension (N+1) :: o
    REAL, dimension (N+1) :: p

    !Se llaman a las funciones
    REAL,EXTERNAL::theta1prim,theta2prim
    INTEGER::i
    REAL::t=0
    !Condiciones iniciales
    o(1)=theta10
    p(1)=theta20


    OPEN(69,file="simple_euler.txt")    
    DO i=1,N

        !Se anota la i-esima componente en el archivo de texto
        WRITE(69,*) t,";",o(i),";",p(i)
        
        !Se implementa el método de EULER para encontrar el paso siguiente
        o(i+1)=o(i)+h*theta1prim(p(i))
        p(i+1)=p(i)+h*theta2prim(o(i),p(i))
        !Se incrimenta el valor de t en h.
        t=t+h

    END DO
    CLOSE(69)

END PROGRAM psimple





!Funciones despejadas del sistema de Ecuaciones diferenciales

REAL FUNCTION theta1prim(theta2)
USE ctes
REAL,INTENT(IN)::theta2
theta1prim =  theta2
END FUNCTION theta1prim

!Función para theta2 prima 
REAL FUNCTION theta2prim(theta1,theta2)
USE ctes
REAL,INTENT(IN) ::theta1,theta2
theta2prim = -(g/l)*sin(theta1)
END FUNCTION theta2prim