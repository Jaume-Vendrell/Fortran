program answer
    real(4) terror1
    terro1=0;
    do n=4, 200
        call Error_Schemes(n,terror1)
        print*, terror1
    enddo
    
end program

subroutine Error_Schemes (n,terror1)
    implicit none
    integer, intent(in) :: n
    real(4), intent(out) :: terror1
    real(8), dimension(n) :: b,x,d,a,c,w,anal,error
    integer::m,i,scheme !! Integers used to go through loops
    real(8)::aE,aW,Dif,F,bL,bR,gamma,Ax,Velocity,Density,s,sp1,sp2,zero,Peclet,Length,terror
    zero=0
    terror1=0
    terror=0
    s=n
    Length=1 !!Length of the rod
    gamma=0.1 !!Difussion coefficient
    Velocity=30!!Velocity of the fluid in m/s
    Density=1 !!Density of the fluid in kg/m^3
    Ax=Length/s !! Width of the grid cell
    Dif=gamma/Ax   
    F=Density*Velocity !!
    bL=100 ! Value of Phi at the left boundary
    bR=20 ! Value of Phi at the right boundary
    Peclet=Density*Velocity*Ax/gamma
       if (F>=0) then
                aE=Dif 
                aW=Dif+F
            else if (F<=0) then
                aE=Dif-F
                aW=Dif
        endif
            sp1=(2*Dif)+F
            sp2=2*Dif
        
    do i=1 , n !Initialize vector space
                if (i==1) then
                    w(i)=Ax/2
                else if (i>1) then
                    w(i)=Ax*i-Ax/2                           
                end if
        enddo
    do i=1 , n !Analytical Solution
            anal(i)=100+(exp(Density*Velocity*w(i)/gamma)-1)/(exp(Density*Velocity*Length/gamma)-1)*(bR-bL)
    enddo   
   
    !Initialize Matrix vectors
    
    do i=1, n !!Initialize vector a
        if (i==1) then
         a(i)=0
        else if (i/=1) then
            a(i)=(-aW)
        end if
    enddo
    
    do i=1 , n !Initialize vector b
        if (i==1) then
         b(i)=0+aE+(sp1)
        else if (i/=1 .OR. i/=n) then
            b(i)=(aW+aE)
        if (i==n) then
            b(i)=aW+0+(sp2)
        end if
        end if
    enddo
    
     do i=1, n !Initialize vector c
        if (i/=n) then
            c(i)=(-aE)
        else if (i==n) then
            c(i)=0
        end if
    enddo
    
    do i=1, n !Initialize vector d
        if (i==1) then
            d(i)=(sp1)*bL
        else if (i>=1 .AND. i<=n) then
            d(i)=0
        if (i==n) then
            d(i)=(sp2)*bR
        end if
        end if
    enddo
        
    call solve_tridiag(a,b,c,d,x,n)
    
    do m=1, n
            error(m)=ABS(((x(m)-anal(m))/anal(m)))
    enddo

    do m=1, n
            terror=terror+error(m)
    enddo
    terror1=terror*100/s
end subroutine Error_Schemes

subroutine solve_tridiag(a, b, c, d, x, n)
    implicit none
    !	 a - sub-diagonal (means it is the diagonal below the main diagonal)
    !	 b - the main diagonal
    !	 c - sup-diagonal (means it is the diagonal above the main diagonal)
    !	 d - right part
    !	 x - the answer
    !	 n - number of equations

    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: a, b, c, d
    real(8), dimension(n), intent(out) :: x
    real(8), dimension(n) :: cp, dp
    real(8) :: m
    integer i

    ! initialize c-prime and d-prime
    cp(1) = c(1)/b(1)
    dp(1) = d(1)/b(1)
    ! solve for vectors c-prime and d-prime
    do i = 2, n
        m = b(i) - cp(i - 1) * a(i)
        cp(i) = c(i)/m
        dp(i) = (d(i) - dp(i - 1) * a(i))/m
    enddo
    ! initialize x
    x(n) = dp(n)
    ! solve for x from the vectors c-prime and d-prime
    do i = n - 1, 1, -1
        x(i) = dp(i) - cp(i) * x(i + 1)
    end do

end subroutine solve_tridiag
   