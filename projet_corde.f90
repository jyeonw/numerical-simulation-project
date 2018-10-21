!!! Shelly CLEMENTE et Jieyeon WOO

PROGRAM projet_corde

    USE mymodule
    IMPLICIT NONE

    !Les variables pour compter le temps d'execution
    integer(kind = 4) clock_count
    integer(kind = 4) clock_count1
    integer(kind = 4) clock_count2
    integer(kind = 4) clock_max
    integer(kind = 4) clock_rate

    REAL, DIMENSION(:,:), ALLOCATABLE :: M, Mchg, Stockvp, MC, L, U                            !Matrice fixe M, Matrice changeable Mchg, Matrice de stockage des valeurs propres et Matrice des Moindres Carrés
    REAL, DIMENSION(:), ALLOCATABLE :: vpCHO, vpSOR, vpDEF, Stocklambda, x, y, b, inter, Res   !Vecteurs propres et Vecteur de stockage des lambdas, Vecteurs des Moindres Carrés
    REAL :: lambdaCHO, lambdaSOR, lambdaDEF, h, s                                              !Valeurs propres, Pas d'espace constant h, valeur intermédiaire
    INTEGER :: N, i, j, ok1, ok2, choix
!Ex1 :
    REAL, PARAMETER :: sigma=167.51, rho=7200., long=0.88                                      !Tension d'intensité(N), Masse volumique(kgm-3) et Longuer(m)
    REAL, PARAMETER :: pi = ACOS(-1.)                                                          !Valeur numérique de Pi

    h = long/N

    !Choix de programme
    print*, "Veuillez choisir : "
    print*, "   1 : Methode de la puissance inverse avec Cholesky pour la pulsation propre la plus petite donnee. (Q.A.2)"
    print*, "   2 : Methode de la puissance inverse avec S.O.R. pour la pulsation propre la plus petite donnee. (Q.A.3)"
    print*, "   3 : Methode de deflation avec puissance iteree pour les 5 modes propres les plus bas. (Q.A.4)"
    print*, "       avec l'interpolation avec la techinque des moindres carrées pour le 5eme vecteur propres. (Q.A.5)"

    read*, choix

    !Saisir la dimension de la matrice M
    print*, "Veuillez entrer le nombre de divisions."
    read*, N

    ALLOCATE (M(N,N), Mchg(N,N), Stockvp(N,N), vpCHO(N), vpSOR(N), vpDEF(N), Stocklambda(N), stat=ok1)
    ALLOCATE (x(n), y(N), b(N), inter(N), MC(N,N), L(N,N), U(N,N), Res(N), stat=ok2)
    IF (ok1/=0) STOP "Erreur : main"
    IF (ok2/=0) STOP "Erreur : main"

    !Generation de la matrice M
    CALL generation_M(M,N)
    !M = M/(h**2)

    !Diriger ver le programme choisit
    SELECT CASE(choix)

!Ex2 :
    !Methode de la puissance inverse avec Cholesky - 1er Choix
    CASE(1)
        print*, 'Methode de la puissance inverse avec Cholesky : '
        !On commence à compter le temps d'exécution
        CALL SYSTEM_CLOCK(clock_count1, clock_rate, clock_max)

        CALL puis_iter_inv_cholesky(M,vpCHO,lambdaCHO)
        print*,'vpCHO = [',vpCHO(:),']'                                     !Vecteurs propres les plus petits de Cholesky
        print*,'lambdaCHO = ',lambdaCHO                                     !Valeur propre la plus petite de Cholesky
        print*, 'pulCHO = ', (N+1)*SQRT(lambdaCHO)/long*SQRT(sigma/rho)     !Pulsation propre la plus petite de Cholesky
        print*, ""

        !On finit de compter le temps d'exécution
        CALL system_clock(clock_count2,clock_rate,clock_max)
        !On affiche le temps d'exécution
        write(*,*)'Elapsed real time = ',real(clock_count2-clock_count1)/real(clock_rate)
        print*, ""

        !Calcul de la valeur analythique de la pulsation la plus petite
        print*, 'Valeur Analythique de pul0 : ', pi/long*SQRT(sigma/rho)


!Ex3 :
    !Methode de la puissance inverse avec SOR - 2ème choix
    CASE(2)
        print*, 'Methode de la puissance inverse avec S.O.R. : '
        !On commence à compter le temps d'exécution
        CALL SYSTEM_CLOCK(clock_count1, clock_rate, clock_max)

        CALL puis_iter_inv_SOR(M,vpSOR,lambdaSOR)
        print*,'vpSOR = [',vpSOR(:),']'                                     !Vecteurs propres les plus petits de S.O.R.
        print*,'lambdaSOR = ',lambdaSOR                                     !Valeur propre la plus petite de S.O.R.
        print*, 'pulSOR = ', (N+1)*SQRT(lambdaSOR)/long*SQRT(sigma/rho)     !Pulsation propre la plus petite de S.O.R.
        print*, ""

        !On finit de compter le temps d'exécution
        CALL system_clock(clock_count2,clock_rate,clock_max)
        !On affiche le temps d'exécution
        write(*,*)'Elapsed real time = ',real(clock_count2-clock_count1)/real(clock_rate)
        print*, ""

        !Calcul de la valeur analythique de la pulsation la plus petite
        print*, 'Valeur Analythique de pul0 : ', pi/long*SQRT(sigma/rho)


!Ex4 :
    !Deflation pour calculer 5 valeurs propres (lambda) et 5e vecteur propre (vp)

    !Puissance Itérée :
    CASE(3)
        print*, 'Methode de deflation avec la puissance iteree : '
        !On commence a compter le temps d'execution
        CALL SYSTEM_CLOCK(clock_count1, clock_rate, clock_max)

        Mchg = M
        !Determination des 1er vecteur propre(vp) et valeur propre (lambda)
        i = 1
        CALL puis_iter(Mchg, vpDEF, lambdaDEF)
        Stockvp(:,i) = vpDEF                    !Stockage des vecteurs propres
        Stocklambda(i) = lambdaDEF              !Stockage des valeurs propres associees

        !Determination des vecteurs propres et valeurs propres restants
        DO i=2,n
        CALL deflation (Mchg, vpDEF, lambdaDEF)
        CALL puis_iter (Mchg, vpDEF, lambdaDEF)
        Stockvp(:,i) = vpDEF                    !Stockage des vecteurs propres
        Stocklambda(i) = lambdaDEF              !Stockage des valeurs propres associees
        END DO


        DO i=1,5
            print*, 'vp',i,'= [', Stockvp(:,i),']'
            print*, 'lambda',i, '=', Stocklambda(i)
        END DO

        print*, 'pul5 : ',(N+1)/long*SQRT(sigma/rho)*lambdaDEF


        print*, ""
        CALL system_clock(clock_count2,clock_rate,clock_max)
        !On affiche le temps d'execution
        write(*,*)'Elapsed real time = ',real(clock_count2-clock_count1)/real(clock_rate)
        print*, ""

        print*, 'Valeur Analythique de pul5 : ', 5*pi/long*SQRT(sigma/rho)
        print*, ""

!Ex5 :
    !Interpolation du cinquième vecteur propre en utilisant la méthode des moindres carrés de degrée 3
        s = 0.
        y = Stockvp(:,5)

        DO i=0,N-1
            x(i+1) = i*h
        END DO

        print*, "x = ", x
        print*, ""

        MC = reshape((/(N+1)*1.,SUM(x),SUM(x**2),SUM(x),SUM(x**2),SUM(x**3),SUM(x**2),SUM(x**3),SUM(x**4)/),(/3,3/))

        print*,"Moindres Carres :"
        CALL affichage(MC,3)
        b = (/SUM(y),SUM(y*x),SUM(y*(x**2))/)
        print*, ""

        print*, "b = ", b
        print*, ""

        CALL decomp(MC,L,U)
        CALL descente(L,b,inter)
        CALL remontee(U,inter,Res)

        PRINT*,"Res =",Res

    CASE DEFAULT
        print*, "Veuillez choisir un nombre entre 1 et 3!"

    END SELECT

DEALLOCATE (M, Mchg, vpCHO, vpSOR, vpDEF, Stockvp, Stocklambda, x, y, b, inter, MC, L, U, Res)
END PROGRAM projet_corde
