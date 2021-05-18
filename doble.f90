!    2021-05-17
!    doble.f90
!    Felipe Ixcamparic (felipechoy1@gmail.com)

!    Este programa es una simulación del movimiento
!    de un péndulo doble implementando el método de 
!    Runge-Kutta de cuarto orden (RK4). 
!

!    Codificación del texto: UTF8
!    Compiladores probados: GNU Fortran (WSl ubuntu Linux 20.04) 7.5.0
!    Instrucciones de compilación: no requiere nada mas
!    gfortran -Wall -pedantic  -c -o doble.o doble.f90
!    gfortran -o doble.x doble.o 
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
    !Constantes/propiedades de los péndulos
    REAL,PARAMETER      ::l1=10.00
    REAL,PARAMETER      ::l2=15.00
    REAL,PARAMETER      ::m1=15.00
    REAL,PARAMETER      ::m2=3.00
    REAL,PARAMETER      ::g=9.8

    !Tasa de muestreo
    REAL, PARAMETER     ::h=0.01
    REAL, PARAMETER     ::tfinal=50.00
    INTEGER,PARAMETER   ::N=CEILING(tfinal/h)


    !CONDICIONES INICIALES
    REAL, PARAMETER     ::phi10=1.2
    REAL, PARAMETER     ::phi20=0.47
    REAL, PARAMETER     ::omega10=0
    REAL, PARAMETER     ::omega20=0
END MODULE ctes

PROGRAM pdoble

    USE ctes
    IMPLICIT NONE
    !Arreglos para almcenar valores
    REAL, DIMENSION (N+1)::o
    REAL, DIMENSION (N+1)::p
    REAL, DIMENSION (N+1)::q
    REAL, DIMENSION (N+1)::r
    !Coeficientes para el método de runge-kutta
    REAL:: k1phi1,k1phi2,k1omega1,k1omega2,k2phi1,k2phi2,k2omega1,k2omega2
    REAL:: k3phi1,k3phi2,k3omega1,k3omega2, k4phi1,k4phi2,k4omega1,k4omega2

    !Se llaman a las funciones 
    REAL, EXTERNAL::phi1prim,phi2prim,omega1prim,omega2prim
    !Condiciones iniciales
   
    REAL::t=0
    INTEGER::i=0
    !Se inicializa las condiciones iniciales 
    
    o(1)=phi10
    p(1)=phi20
    q(1)=omega10
    r(1)=omega20
    
    

    OPEN(7,file="pendulo_doble.txt")
     
        
        DO i=1,N
            

            !k1's descritos por el método 
            k1phi1 = phi1prim(q(i))
            k1phi2 = phi2prim(r(i))
            k1omega1 = omega1prim(o(i),p(i),r(i))
            k1omega2 = omega2prim(o(i),p(i),q(i),r(i))
            

            !k2's descritos por el método
            k2phi1      = phi1prim(q(i)+(h/2)*k1omega1)
            k2phi2      = phi2prim(r(i)+(h/2)*k1omega2)
            k2omega1    = omega1prim(o(i)+(h/2)*k1phi1,p(i)+(h/2)*k1phi2, r(i)+(h/2)*k1omega2 )
            k2omega2    = omega2prim(o(i)+(h/2)*k1phi1,p(i)+(h/2)*k1phi2, q(i)+(h/2)*k1omega1 ,r(i)+(h/2)*k1omega2 )


            !k3's descritos por el método 
            k3phi1      = phi1prim(q(i)+(h/2)*k2omega1)
            k3phi2      = phi2prim(r(i)+(h/2)*k2omega2)
            k3omega1    = omega1prim(o(i)+(h/2)*k2phi1, p(i)+(h/2)*k2phi2 , r(i)+(h/2)*k2omega2)
            k3omega2    = omega2prim(o(i)+(h/2)*k2phi1, p(i)+(h/2)*k2phi2 , q(i)+(h/2)*k2omega1 , r(i)+(h/2)*k2omega2)



            !k4's descritos por el método
            k4phi1      =phi1prim(q(i)+h*k3omega1)
            k4phi2      =phi2prim(r(i)+h*k3omega2)
            k4omega1    =omega1prim(o(i)+h*k3phi1,p(i)+h*k3phi2,r(i)+h*k3omega2)
            k4omega2    =omega2prim(o(i)+h*k3phi1,p(i)+h*k3phi2,q(i)+h*k3omega1,r(i)+h*k3omega2)
            
            !Se escriba la i-ésima componente en el archivo de texto 
            WRITE(7,*) t,";",o(i),";",p(i),";",q(i),";",r(i)

            !Se encuentra el siguiente valor para ángulos y velocidades angulares
            !descrito por el método de runge kutta por definición
            o(i+1)=o(i)+(h/6)*(k1phi1+2*k2phi1+2*k3phi1+k4phi1)
            p(i+1)=p(i)+(h/6)*(k1phi2+2*k2phi2+2*k3phi2+k4phi2)
            q(i+1)=q(i)+(h/6)*(k1omega1+2*k2omega1+2*k3omega1+k4omega1)
            r(i+1)=r(i)+(h/6)*(k1omega2+2*k2omega2+2*k3omega2+k4omega2)
            t=t+h

        

        END DO
        !Se cierra el archivo de texto
     CLOSE(7)     
   



END PROGRAM pdoble



!FUNCIONES A UTILIZAR descritas por el sistema de ecuaciones diferenciales

REAL FUNCTION phi1prim(omega1)
    USE ctes
    IMPLICIT NONE
    REAL,INTENT(IN)::omega1
    phi1prim=omega1
END FUNCTION phi1prim


REAL FUNCTION phi2prim(omega2)
    USE ctes
    IMPLICIT NONE
    REAL,INTENT(IN)::omega2
    phi2prim=omega2
END FUNCTION phi2prim


REAL FUNCTION omega1prim(phi1,phi2,omega2)
    USE ctes
    IMPLICIT NONE
    REAL,INTENT(IN)::phi1,phi2,omega2
    omega1prim= (-g*(2*m1+m2)*sin(phi1)+m2*g*sin(phi1-2*phi2)-2*sin(phi1-phi2)*m2*(omega2**2*l2+omega2**2*&
    &l1*cos(phi1-phi2))) / (l1*2*(m1+m2-m2*cos(2*phi1-2*phi2)))
END FUNCTION omega1prim



REAL FUNCTION omega2prim(phi1,phi2,omega1,omega2)
    USE ctes
    IMPLICIT NONE
    REAL,INTENT(IN)::phi1,phi2,omega1,omega2
    omega2prim  =(2*sin(phi1-phi2)*(omega1**2*l1*(m1+m2)+g*(m1+m2)*cos(phi1)+omega2**2*l2*m2*cos(phi1-phi2))) &
    &/(l2*(2*m1+m2-m2*cos(2*phi1-2*phi2)) )
END FUNCTION omega2prim