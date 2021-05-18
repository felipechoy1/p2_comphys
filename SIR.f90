!    2021-05-17
!    SIR.f90
!    Felipe Ixcamparic (felipechoy1@gmail.com)

!    Este programa hace una implementación del
!    modelo SIR para panemias utilizando el método
!    de Runge-Kutta.    
!

!    Codificación del texto: UTF8
!    Compiladores probados: GNU Fortran (SUSE Linux) 7.5.0
!    Instrucciones de compilación: no requiere nada mas
!    gfortran -Wall -pedantic -std=f95 -c -o SIR.o SIR.f90
!    gfortran -o SIR.x SIR.o 


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
MODULE constantes
IMPLICIT NONE 
   
   REAL, PARAMETER :: beta=0.5
   REAL, PARAMETER :: gamma=0.1111111
   REAL, PARAMETER :: mu=0.001167
   REAL, PARAMETER :: A= 1449401
   REAL, PARAMETER :: N= 1449401
   REAL, PARAMETER :: S0 = 1446093 / N
   REAL, PARAMETER :: I0 = 1885 / N
   REAL, PARAMETER :: R0 = 1423 /N
   REAL, PARAMETER::  h=0.01
   REAL, PARAMETER::  tfinal=20
   INTEGER,PARAMETER  :: eNe=CEILING(tfinal/h)
END MODULE constantes

PROGRAM test
    !Se llama al módulo de constantes
    USE constantes
    IMPLICIT NONE
    !Se llaman a las funciones de SIR
    REAL,EXTERNAL::S,I,R
    !Iteradores
    INTEGER::w
    !Vectores para almacenar cada valor de SIR
    REAL, dimension (eNe+1) :: o
    REAL, dimension (eNe+1) :: p
    REAL, dimension (eNe+1) :: q
        
    !Condiciones iniciales
    o(1)=S0
    p(1)=I0
    q(1)=R0

    OPEN(13,file="SIR.txt")

    !Si se desea plotear en GNUPLOT
    ! OPEN(14,file="I.txt")
    ! OPEN(15,file="R.txt")
    
    DO w=2,eNe
        !Se calculan los S_j+1 utilizando los valores enteriores
        o(w)=S(o(w-1),p(w-1),q(w-1))
        p(w)=I(o(w-1),p(w-1),q(w-1))
        q(w)=R(o(w-1),p(w-1),q(w-1))

    END DO
                                                                    
    DO w=1,eNe
        !Se anotan los valores obtenidos en un archivo de texto
        WRITE(13,*) w*h,";",o(w),";",p(w),";",q(w)

        !Si se desea realizar en GNUPLOT
        ! WRITE(13,*) w*h," ",o(w)
        ! WRITE(14,*) w*h," ",p(w)
        ! WRITE(15,*) w*h," ",q(w)
    END DO
    CLOSE(13)

    ! CLOSE(14)
    ! CLOSE(15)
    


END PROGRAM test


!Funciones para el incremento en cada uno de los contenedores


REAL FUNCTION S(x,y,z)
    USE constantes
    REAL::x,y,z
    REAL,EXTERNAL::k1,k2,k3,k4
    S=x+ (h/6)*(k1(x,y,z)+2*k2(x,y,z)+2*k3(x,y,z)+k4(x,y,z))
END FUNCTION


REAL FUNCTION I(x,y,z)
    USE constantes
    REAL::x,y,z
    REAL,EXTERNAL:: l1,l2,l3,l4
    I = y + (h/6)*(l1(x,y,z)+2*l2(x,y,z)+2*l3(x,y,z)+l4(x,y,z))
END FUNCTION


REAL FUNCTION R(x,y,z)
    USE constantes
    REAL::x,y,z
    REAL,EXTERNAL:: m1,m2,m3,m4
    R = z + (h/6)*(m1(x,y,z)+2*m2(x,y,z)+2*m3(x,y,z)+m4(x,y,z))
END FUNCTION



!Coeficientes k's, l's y m's.

REAL FUNCTION k1(x,y,z)
    USE constantes
    REAL::x,y,z
    k1= A-mu*x-(beta*x*y/n)
END FUNCTION k1

REAL FUNCTION l1(x,y,z)
    USE constantes
    REAL::x,y,z
    l1= (beta*x*y/n)-(mu+gamma)*y
END FUNCTION l1

REAL FUNCTION m1(x,y,z)
    USE constantes
    REAL::x,y,z
    m1 = gamma*y-mu*z
END FUNCTION m1

REAL FUNCTION k2(x,y,z)
    USE constantes
    REAL::x,y,z
    REAL,EXTERNAL::k1,l1
    k2= A-mu*(x+(h*k1(x,y,z)/2))-((beta/n)*(y+(h*l1(x,y,z)/2))*(x+(h*k1(x,y,z)/2)) )
END FUNCTION k2

REAL FUNCTION l2(x,y,z)
    USE constantes
    REAL::x,y,z
    REAL,EXTERNAL::k1,l1
    l2=(beta/n)*(y+(h*l1(x,y,z)/2))*(x+(h*k1(x,y,z)/2))-((mu+gamma)*(y+(h*l1(x,y,z)/2)))
END FUNCTION l2

REAL FUNCTION m2(x,y,z)
    USE constantes
    REAL::x,y,z
    REAL,EXTERNAL::l1,m1
    m2 = gamma*(y+(h*l1(x,y,z)/2))-mu*(z+(h*m1(x,y,z)/2))
END FUNCTION m2

REAL FUNCTION k3(x,y,z)
    USE constantes
    REAL::x,y,z
    REAL, EXTERNAl::k2,l2
    k3= A- mu*(x+(h*k2(x,y,z)/2))-(beta/n)*(y+(h*l2(x,y,z)/2))*(x+(h*k2(x,y,z)/2))
END FUNCTION k3

REAL FUNCTION l3(x,y,z)
    USE constantes
    REAL::x,y,z
    REAL, EXTERNAL::l2,k2
    l3 = (beta/n)*(y+(h*l2(x,y,z)/2))*(x+(h*k2(x,y,z)/2))-(mu+gamma)*(y+(h*l2(x,y,z)/2))
END FUNCTION l3

REAL FUNCTION m3(x,y,z)
    USE constantes
    REAL::x,y,z
    REAL, EXTERNAL::l2,m2
    m3 = gamma*(y+(h*l2(x,y,z)/2))-mu*(z+(h*m2(x,y,z)/2))
END FUNCTION m3

REAL FUNCTION k4(x,y,z)
    USE constantes
    REAL::x,y,z
    REAL,EXTERNAL::k3,l3
    k4= A-mu*(x+h*k3(x,y,z))-(beta/n)*(y+h*l3(x,y,z))*(x+(h*k3(x,y,z)))
END FUNCTION k4

REAL FUNCTION l4(x,y,z)
    USE constantes
    REAL::x,y,z
    REAL, EXTERNAL::l3,k3
    l4= (beta/n)*(y+(h*l3(x,y,z)))*(x+(h*k3(x,y,z)))-(mu+gamma)*(y+h*l3(x,y,z))
END FUNCTION l4

REAL FUNCTION m4(x,y,z)
    USE constantes
    REAL::x,y,z
    REAL, EXTERNAL::l3,m3
    m4= gamma*(y+h*l3(x,y,z))-mu*(z+h*m3(x,y,z))
END FUNCTION m4