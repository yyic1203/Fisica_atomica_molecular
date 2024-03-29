program HEHF
  implicit none
  INTEGER :: i, j, k, l, iter
  INTEGER, PARAMETER :: n= 2 !!Número atómico para iterar
  DOUBLE PRECISION, PARAMETER :: Z=2.0D0, re=-2.904D0 !! Número atómico
  DOUBLE PRECISION, dimension(n) :: zeta !! Coeficiente de Slater
  DOUBLE PRECISION, dimension(n,n) :: P, S, H, G, F, C, oldP, X_dag, X, CP, FP, tem !!Matrices de 2 dimensiones
	DOUBLE PRECISION, dimension(n,n,n,n) :: T !!Matriz de 4 dimensiones
	DOUBLE PRECISION :: su, mul, E, Delta, crit, er!! Valores
 	
  zeta(1)=1.45363D0 !!Parametro de slater 1
  zeta(2)=2.91093D0 !!Parametro de slater 2

  DO i= 1, n
    DO j= 1, n
      P(i,j)=0.0D0 !!Densidad de probabilidad
    END DO
  END DO

  DO i= 1, n
    DO j= 1, n
      S(i,j)=(8.0D0*(zeta(i)*zeta(j))**(1.5D0))/(zeta(i)+zeta(j))**(3.0D0) !!Matriz de ortogonalidad
    END DO
  END DO

  DO i= 1, n
    DO j= 1, n
      H(i,j)= (4.0d0 * (zeta(i) * zeta(j))**(2.5d0) / (zeta(i) + zeta(j))**(3.0d0)) - &
              (4.0d0 * Z * (zeta(i) * zeta(j))**(1.5d0) / (zeta(i) + zeta(j))**(2.0d0)) !!Matriz de energía de 1 electrón
    END DO
  END DO
  
  DO i= 1, n
    DO j= 1, n
      DO k= 1, n
      	DO l= 1, n
      		su= zeta(i)+zeta(j)+zeta(k)+zeta(l)
      		mul= zeta(i)*zeta(j)*zeta(k)*zeta(l)
      		T(i,j,k,l)= 16.0D0*((mul)**(1.5D0))*((2/((zeta(i)+zeta(k))**(3.0D0)*(zeta(j)+zeta(l))**(2.0D0))) - & 
          (2/((zeta(i)+zeta(k))**(2.0D0)*(su)**(3.0D0)))- & 
          (2/((zeta(i)+zeta(k))**(3.0D0)*(su)**(2.0D0)))) !!Matriz de energía de 2 electrón
      	END DO
      END DO
    END DO
  END DO
  
  call find_X_and_X_dagger(S, X, X_dag, n) !!Llama la subrutina para calcular la matriz de transformación
  
  Open(12, file="energías", status='unknown')
  write(12,*) '*********************************************************'
  write(12,15)
  15 format (4x, 'Iteración', 2x, 'E_total u.a', 2x, 'E_real u.a', 2x, 'Error %',/)
  iter=0
  crit= 1.0D-15 !!Criterio de estabilidad
  1 iter=iter+1
  G=0.0D0
  DO i= 1, n
    DO j= 1, n
      DO k= 1, n
      	DO l= 1, n
      		G(i,j) = G(i,j) + P(k,l)*(T(i,k,j,l) - 0.5D0*T(i,k,l,j))!!Aporte interacción electrón-electrón
      	END DO
      END DO
    END DO
  END DO
  F=0.0D0
  DO i= 1, n
    DO j= 1, n
      F(i,j)= H(i,j) + G(i,j) !!Matriz de Fock
    END DO
  END DO
  
  E=0.0D0
  DO i= 1, n
    DO j= 1, n
      E = E + 0.5D0*P(i,j)*(H(i,j)+F(i,j)) !!Energía del sistema
    END DO
  END DO
	 
	er=(abs(E+2.904D0)/2.904D0)*100 !!Error relativo porcentual
	
	tem=matmul(X_dag,F)
	FP=matmul(tem,X) !!Transformación de F a la base ortonormal
	
	call autovector(FP, CP, n) !!Calculo de los autovectores de F en el sistema primado 
	
	C=0.0D0
	C=matmul(X,transpose(CP)) !!Nuevos coeficientes de peso
	
	oldP=P
	
	DO i= 1, n
    DO j= 1, n
      DO k= 1, n
      	P(i,j)=(2.0D0*C(k,i)*C(j,k)) !!Nueva matriz de densidad poblacional 
      END DO
    END DO
  END DO
  
  Delta=0.0D0
  DO i= 1, n
    DO j= 1, n
      Delta = Delta + ABS(oldP(i,j)-P(i,j)) !!Diferencia poblacional
    END DO
  END DO
 
  write (12,18) iter, E, re, er
  18 format(4x,I3,8x,F9.6,3x,F9.6,3x,F9.6)
  
  
  if (Delta .gt. crit) then
		go to 1 !!Retorna el ciclo
	end if
	return
end program HEHF

subroutine autovector(FP, CP, n) !!Subrutina para calcula autovectores
   	use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    !!Parametros para llamar la función dsyevd
    integer, intent(in) :: n 
    real(dp), intent(in) :: FP(n, n)
    real(dp), intent(out) :: CP(n, n)
    real(dp), dimension(:), allocatable :: work
    real(dp), dimension(n) :: eigenvalues
    integer, dimension(:), allocatable :: iwork
    integer :: info, lwork, liwork
    integer :: i

    ! Calculamos los tamaños de los arrays de trabajo
    lwork = -1
    liwork = -1
    allocate(work(1))
    allocate(iwork(1))
		
		 ! Llamada de trabajo para dsyevd para calcular el tamaño óptimo de lwork y liwork
    call dsyevd('V', 'U', n, FP, n, eigenvalues, work, lwork, iwork, liwork, info)

    ! Asignamos los tamaños de los arrays de trabajo
    lwork = int(work(1))
    liwork = iwork(1)
    deallocate(work)
    deallocate(iwork)
    allocate(work(lwork))
    allocate(iwork(liwork))
		
    ! Llamamos a dsyevd para calcular los valores propios y los vectores propios de S
    call dsyevd('V', 'U', n, FP, n, eigenvalues, work, lwork, iwork, liwork, info)

    if (info /= 0) then
        print *, 'Error en la diagonalizacion: Info =', info
        stop
    end if



    ! Formamos la matriz C en la base ortonormal
    CP = FP

end subroutine autovector

subroutine find_X_and_X_dagger(S, X, X_dag, n) !!Subrutina para ortonormalizar la base
   	use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    !!Parametros para usar el comando dsyevd
    integer, intent(in) :: n
    real(dp), intent(in) :: S(n, n)
    real(dp), intent(out) :: X(n, n), X_dag(n, n)
    real(dp), dimension(:), allocatable :: work
    real(dp), dimension(n) :: eigenvalues
    real(dp), dimension(n, n) :: D, SP
    integer, dimension(:), allocatable :: iwork
    integer :: info, lwork, liwork
    integer :: i,j

    ! Calculamos los tamaños de los arrays de trabajo
    lwork = -1
    liwork = -1
    allocate(work(1))
    allocate(iwork(1))
		
		 ! Llamada de trabajo para dsyevd para calcular el tamaño óptimo de lwork y liwork
    call dsyevd('V', 'U', n, S, n, eigenvalues, work, lwork, iwork, liwork, info)

    ! Asignamos los tamaños de los arrays de trabajo
    lwork = int(work(1))
    liwork = iwork(1)
    deallocate(work)
    deallocate(iwork)
    allocate(work(lwork))
    allocate(iwork(liwork))
		
    ! Llamamos a dsyevd para calcular los valores propios y los vectores propios de S
    call dsyevd('V', 'U', n, S, n, eigenvalues, work, lwork, iwork, liwork, info)

    if (info /= 0) then
        print *, 'Error en la diagonalizacion: Info =', info
        stop
    end if

    ! Inicializamos la matriz D como cero y asignamos la inversa de la raíz cuadrada
    ! de los valores propios en su diagonal
    D = 0.0_dp
    do i = 1, n
        D(i, i) = 1.0_dp / sqrt(eigenvalues(i))
    end do

    ! Formamos la matriz X multiplicando los vectores propios de S por D
    
    X = matmul(S, D)
		X = transpose (X)
    ! X_dagger es simplemente la transpuesta de X
    X_dag = transpose(X)
end subroutine find_X_and_X_dagger

