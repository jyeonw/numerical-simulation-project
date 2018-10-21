!!! Shelly CLEMENTE et Jieyeon WOO

MODULE mymodule

IMPLICIT NONE

REAL :: eps = 1.E-10

CONTAINS

!Affichage d'une matrice M de dimension NxN
SUBROUTINE affichage(M,N)
REAL, DIMENSION(:,:), INTENT(IN) :: M   !Matrice M
INTEGER, INTENT(IN) :: N                !Dimension de la matrice M
INTEGER :: i, j
PRINT*, "La matrice est:"
DO i = 1, N
    DO j = 1, N
        WRITE (*,"(F10.4)", ADVANCE='no') M(i,j)    !Affichage des valeurs dans la matrice en 10bits avec 4 chiffres décimales
    END DO
WRITE (*,*)
END DO
END SUBROUTINE affichage

!Méthode de la Descente
SUBROUTINE descente (A,B,X)
REAL, DIMENSION(:,:), INTENT(IN) :: A               !Matrice A
REAL, DIMENSION(:), INTENT(IN) :: B                 !Vecteur B
REAL, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: X   !Vecteur inconnu X
INTEGER :: n, i, j, ok
REAL :: s

n=SIZE(A,1)
ALLOCATE(X(n), STAT=ok)
IF (ok/=0) STOP "ERREUR : descente"

DO i = 1, n
	s=0.
	DO j=1, i-1
		s=s+A(i,j)*X(j)
	end DO
	X(i)= (B(i)-s)/A(i,i)
end DO
end SUBROUTINE descente

!Méthode de la Remontée
SUBROUTINE remontee (A,B,X)
REAL, DIMENSION(:,:), INTENT(IN) :: A               !Matrice A
REAL, DIMENSION(:),INTENT(IN) :: B                  !Vecteur B
REAL, DIMENSION(:), INTENT(OUT) :: X                !Vecteur inconnu X
INTEGER :: n, i, j, ok
REAL :: s

n=SIZE(A,1)

DO i = n,1,-1
	s=0.
	DO j=n, i+1, -1
		s=s+A(i,j)*X(j)
	end DO
	X(i)= (B(i)-s)/A(i,i)
end DO
end SUBROUTINE remontee

!Decomposition L,U
SUBROUTINE decomp (A,L,U)
REAL, DIMENSION(:,:), INTENT(IN) :: A               !Matrice A
REAL, DIMENSION(:,:), INTENT(OUT) :: L, U           !Matrices de décomposition L et U
INTEGER :: n, i, j, k, ok
REAL :: s

n=SIZE(A,1)

DO j = 1, n
	L(j,j) = 1.
	DO i = 1, j
		s=0.
		DO k = 1, i-1
			s = s + L(i,k)*U(k,j)
		end DO
		U(i,j) = A(i,j) - s
	end DO
	DO i = j+1, n
		s=0.
		DO k = 1, j-1
			s = s + L(i,k)*U(k,j)
		end DO
		L(i,j) = (A(i,j) - s) / U(j,j)
	end DO
end DO
end SUBROUTINE decomp

!Generation de la Matrice M (symétrique)
SUBROUTINE generation_M(M,N)
REAL, DIMENSION(:,:), INTENT(OUT) :: M              !Matrice M
INTEGER, INTENT(IN) :: N                            !Dimension de la matrice M
INTEGER :: i, j, ok

DO i = 1, N
    DO j = 1, N
        IF (i == j) THEN
            M(i,j) = 2
        ELSE IF ((i==j+1) .OR. (i==j-1)) THEN
            M(i,j) = -1
        ELSE
            M(i,j) = 0
        END IF
    END DO
END DO
!CALL affichage(M,N)
END SUBROUTINE generation_M

!Méthode de Décomposition de Cholesky (M=LLt)
SUBROUTINE cholesky(N,M,L,Lt)
REAL, DIMENSION(:,:), INTENT(IN) :: M               !Matrice M
REAL, DIMENSION(:,:), INTENT(OUT) :: L, Lt          !Matrices de décomposition L et Lt
INTEGER, INTENT(IN) :: N                            !Dimension de la matrice
INTEGER :: i, j, k, ok

DO j = 1, N
    DO k = 1, j-1
        L(j,j) = L(j,j) + L(j,k)*L(j,k)
    END DO
    L(j,j) = SQRT(M(j,j) - L(j,j))
    DO i = j+1, N
        DO k = 1, j-1
            L(i,j) = L(i,j) + L(i,k)*L(j,k)
        END DO
        L(i,j) = (M(i,j) - L(i,j))/L(j,j)
    END DO
END DO
Lt = TRANSPOSE(L)
END SUBROUTINE cholesky

!Methode de la puissance inverse avec Cholesky
SUBROUTINE puis_iter_inv_cholesky(M, vp, lambda)

REAL, DIMENSION(:,:), INTENT(IN) :: M               !Matrice M
REAL, DIMENSION(:), ALLOCATABLE,INTENT(OUT) :: vp   !Vecteur propre
REAL, INTENT(OUT) :: lambda                         !Valeur propre
REAL, DIMENSION(:,:), ALLOCATABLE :: L, Lt          !Matrices de décomposition L et Lt
REAL, DIMENSION(:), ALLOCATABLE :: u, y             !Vecteurs vecteur intermédiaire
REAL :: oldlambda, oldmu, mu                        !Valeurs intermédiaire
INTEGER :: i=0, n, ok

n=SIZE(M,1)
ALLOCATE(L(n,n), Lt(n,n), u(n), y(n), STAT=ok)
IF (ok/=0) STOP "ERREUR : puis_iter_inv_cholesky"

CALL cholesky(n,M,L,Lt)
CALL RANDOM_NUMBER(u)                               !Nombre aléatoire de u
vp = u/ SQRT(SUM(u*u))
mu = -10.
oldmu = mu + 0.1

DO WHILE (ABS(mu-oldmu) > eps)                      !Condition de sortie
    i = i + 1
    vp = u / SQRT(SUM(u*u))
    CALL descente(L,vp,y)
    CALL remontee(Lt,y,u)
    oldmu= mu
    oldlambda = lambda
    mu = DOT_PRODUCT (vp, u)
END DO

vp = u / SQRT(SUM(u*u))
lambda= 1./mu

DEALLOCATE(u,y)
END SUBROUTINE puis_iter_inv_cholesky

!Methode SOR avec un pas d'iteration
SUBROUTINE iter_SOR(omega,M,r,e)
REAL, INTENT(IN) :: omega                           !Paramètre de relaxation
REAL, DIMENSION(:,:), INTENT(IN) :: M               !Matrice M
REAL, DIMENSION(:), INTENT(IN) :: r                 !Second membre des méthodes itératives
REAL, DIMENSION(:), INTENT(OUT) :: e                !Système inconnue
REAL :: s
INTEGER :: i, j, n

n = SIZE(M,1)

DO i = 1, n
    s=0.
    DO j = 1, i-1
        s = s + M(i,j)*e(j)
    END DO
    e(i)=omega*(r(i)-s)/M(i,i)
END DO
end SUBROUTINE iter_SOR

!Methode de la puissance inverse avec SOR
SUBROUTINE puis_iter_inv_SOR(M, vp, lambda)

REAL, DIMENSION(:,:), INTENT(IN) :: M               !Matrice M
REAL, DIMENSION(:), INTENT(OUT) :: vp               !Vecteur propre
REAL, INTENT(OUT) :: lambda                         !Valeur propre
REAL, DIMENSION(:), ALLOCATABLE :: u, y, r, e       !Vecteurs vecteur intermédiaire, Second membre et Système inconnue
REAL :: mu, oldmu, cond, omega                      !Valeurs intermédiaire, Paramètre de relaxation
INTEGER :: i = 0, n, k, ok

n = SIZE(M,1)

IF(n<=10) then              !Nombre d'itérations le plus petit pour omega = 10 : 1.5->6 iter, 100 : 1.99->13 iter, 1000 : 1.99->122 iter
    omega = 1.5
ELSEIF(n<=1000) then
    omega = 1.99
ELSE
    omega = 1.
END IF

ALLOCATE (u(n),y(n),r(n),e(n), STAT=ok)
IF (ok/=0) STOP "ERREUR : puis_iter_inv_SOR"

CALL RANDOM_NUMBER(u)                               !Nombre aléatoire de u
vp = u / SQRT(SUM(u*u))
mu = -10.
oldmu = mu + 0.1
k=1

DO WHILE (ABS(mu - oldmu) > eps)                    !Condition de sortie de 1ère boucle
    i = i + 1
    vp = u / SQRT(SUM(u*u))
    cond = eps*2
    DO WHILE (cond > eps)                           !Condition de sortie de 2ème boucle
        r = vp - MATMUL(M,u)
        CALL iter_SOR(omega,M,r,e)
        u = u + e
        cond = MAXVAL(ABS(r))
    END DO
    oldmu = mu
    mu = DOT_PRODUCT (vp, u)
    k = k+1
!print*, (N+1)/sqrt(mu)
END DO
!print*, k

vp = u / SQRT(SUM(u*u))
lambda= 1./mu

DEALLOCATE (u, y)
END SUBROUTINE puis_iter_inv_SOR

!Méthode de la puissance iteree
SUBROUTINE puis_iter(M, vp, lambda)
REAL, DIMENSION(:,:), INTENT(IN) :: M               !Matrice M
REAL, DIMENSION(:), INTENT(OUT) :: vp               !Vecteur propre
REAL, INTENT(OUT) :: lambda                         !Valeur propre
REAL, DIMENSION(:), ALLOCATABLE :: u                !Vecteur intermédiaire
REAL :: anclambda                                   !Valeur intermédiaire
INTEGER :: i = 0, n

n = SIZE(M,1)
ALLOCATE (u(n))
CALL RANDOM_NUMBER(u)                               !Nombre aléatoire de u
vp = u / SQRT(SUM(u*u))
u = MATMUL(M, vp)
lambda = DOT_PRODUCT (vp, u)
anclambda = lambda + 2*eps
DO WHILE (abs(anclambda-lambda) > eps)              !Condition de sortie
    i = i + 1
    vp = u / SQRT(SUM(u*u))
    u = MATMUL(M, vp)
    anclambda = lambda
    lambda = DOT_PRODUCT (vp, u)
END DO
vp = u / SQRT(SUM(u*u))

DEALLOCATE (u)
End subroutine puis_iter

!Méthode de Déflation
SUBROUTINE deflation(M, vp, lambda)
REAL, DIMENSION(:,:), INTENT(INOUT) :: M            !Matrice M
REAL, DIMENSION(:), INTENT(IN) :: vp                !Vecteur propre
REAL, INTENT(IN) :: lambda                          !Valeur propre
REAL, DIMENSION(:), ALLOCATABLE :: vpt              !Vecteur propre transposé
REAL :: lambdat, fact                               !Valeur propre transposée et Valeur intermédiaire
INTEGER :: i, j, n

n = SIZE(M, 1)
ALLOCATE(vpt(n))
!On cherche le vecteur propre dominant vpt de M transposée
CALL puis_iter (TRANSPOSE(M), vpt, lambdat)
DO WHILE (abs(lambdat - lambda) < eps)             !Condition de sortie

fact = lambda/dot_product(vp, vpt)

DO i = 1, n
    DO j = 1, n
        M(i,j)=M(i,j)-fact*vp(i)*vpt(j)
    END DO
END DO
END DO

DEALLOCATE(vpt)
END SUBROUTINE deflation

END MODULE mymodule
