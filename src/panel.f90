module panel_m
    use math_m
    implicit none
    
    type panel_t
        integer :: npts
        
        !Options
        integer :: PointsPlacement  ! 1 = Standard Method Presented by Phillips 
                                    ! 2 = Cosine Clustered on Panel using Panel normal
                                    ! 3 = Cosine Clustered on Surface using Panel Normal
                                    ! 4 = Cosine Clustered on Surface using Surface normal
        integer :: OpenTE           ! 1 Use the standard open trailing edge NACA 4-Digit Airfoil
                                    ! 0 Use the modified closed trailing edge NACA 4-Digit Airfoil
        integer :: readfile         ! 1 = read geometry from a file
        character(LEN=50):: readfilename
        
        character(LEN=2) :: load   != 'UL' if uniform load is specified
        real :: CLd     ! Design lift coefficient
        
        real :: alpha   !angle of attack in radians
        real :: Vinf(2) !Freestream velocity vector
        real :: Vmag    !Magnitude of freestream velocity
!        real :: pi      !pi defined in setup subroutine
        
        CHARACTER(LEN=4) :: naca

        real :: LE(2)   !Leading Edge position
        real :: TE(2)   !Trailing Edge position
        real :: xmc     !maximum camber x-location
        real :: ymc     !maximum camber y-percent
        real :: tm      !maximum thickness
        
        real,allocatable,dimension(:,:) :: Points
        real,allocatable,dimension(:,:) :: Control
        real,allocatable,dimension(:,:) :: Normals
        real,allocatable,dimension(:,:) :: Camber
        real,allocatable,dimension(:) :: Gamma
        real,allocatable,dimension(:,:) :: Wake
    end type panel_t
    
contains

!-----------------------------------------------------------------------------------------------------------
    subroutine panel_allocate(t)
        implicit none
        type(panel_t) :: t

        allocate(t%Points(t%npts,2))
        allocate(t%Control(t%npts,2))
        allocate(t%Normals(t%npts,2))
        allocate(t%Camber(t%npts,2))
        allocate(t%Gamma(t%npts))
        
    end subroutine panel_allocate

!-----------------------------------------------------------------------------------------------------------
    subroutine panel_cleanup(t)
        implicit none
        type(panel_t) :: t
        deallocate(t%Points)
        deallocate(t%Control)
        deallocate(t%Normals)
        deallocate(t%Camber)
        deallocate(t%Gamma)
    end subroutine panel_cleanup

!-----------------------------------------------------------------------------------------------------------
    subroutine panel_create_from_file(t)
        implicit none
        type(panel_t) :: t
        t%readfile = 1
        call panel_read_airfoil(t)
        call panel_setup(t)
    end subroutine panel_create_from_file

!-----------------------------------------------------------------------------------------------------------
    subroutine panel_create_from_data(t,datasize,rawdata)
        implicit none
        type(panel_t) :: t
        integer :: datasize
        real :: rawdata(datasize,2)
        t%readfile = 0

        t%npts = datasize
        call panel_allocate(t)
        t%Points(:,:) = rawdata(:,:)

        call panel_setup(t)

    end subroutine panel_create_from_data

!-----------------------------------------------------------------------------------------------------------
    subroutine panel_create_from_naca(t)
        implicit none
        type(panel_t) :: t
        real :: dtheta,x,xc,yc
        integer :: i
        
        t%readfile = 0

        if(t%load.eq.'UL') then
            t%xmc = 0.5
            t%ymc = -t%CLd/(4.0*pi)*log(0.5)
            read(t%naca(3:4),*) t%tm;   t%tm = t%tm/100.0
        else
            read(t%naca(1:1),*) t%ymc; t%ymc = t%ymc/100.0
            read(t%naca(2:2),*) t%xmc; t%xmc = t%xmc/10.0
            read(t%naca(3:4),*) t%tm;   t%tm = t%tm/100.0
        end if

        call panel_allocate(t)
        ! Points along surface
        dtheta = 2.0*pi/real(t%npts-1)
        do i=1,t%npts/2
            x = 0.5*(1.0-cos((real(i)-0.5)*dtheta))
            CALL panel_geometry(t,x,xc,yc,t%Points(t%npts/2+i,1),t%Points(t%npts/2+i,2),t%Points(t%npts/2+1-i,1),&
                                &t%Points(t%npts/2+1-i,2))
        end do

        call panel_setup(t)

    end subroutine panel_create_from_naca


!-----------------------------------------------------------------------------------------------------------
    subroutine panel_setup(t)
        implicit none
        type(panel_t) :: t
        integer :: i,ierror
        real :: dtheta,x,xc,yc,xu,yu,xl,yl,li,D1,D2,percent,dummy
        100 FORMAT (1X, 1000ES25.16)
        
        t%LE(:) = [0.0,0.0]
        t%TE(:) = 0.5*(t%Points(1,:) + t%Points(t%npts,:))

        if(t%PointsPlacement.eq.1) then !Standard Method Presented by Phillips - automatically used if t%readfile.eq.1
            do i=1,t%npts-1
                t%Control(i,1) = 0.5*(t%Points(i,1) + t%Points(i+1,1))
                t%Control(i,2) = 0.5*(t%Points(i,2) + t%Points(i+1,2))
                li = sqrt((t%Points(i,1)-t%Points(i+1,1))**2 + (t%Points(i,2)-t%Points(i+1,2))**2)
                t%Normals(i,1) = - (t%Points(i+1,2)-t%Points(i,2))/li
                t%Normals(i,2) =   (t%Points(i+1,1)-t%Points(i,1))/li
            end do
        else if(t%PointsPlacement.eq.2) then !Cosine Clustered on Panel using Panel normal
            dtheta = pi/(real((t%npts-1))/2.0)
            t%Control(t%npts/2,:) = 0.0
            do i=1,t%npts/2-1
                x = 0.5*(1.0-cos((real(i))*dtheta))
                CALL panel_geometry(t,x,xc,yc,t%Control(t%npts/2+i,1),t%Control(t%npts/2+i,2),t%Control(t%npts/2-i,1),&
                                    &t%Control(t%npts/2-i,2))
            end do
            do i=1,t%npts-1
                D1 = panel_distance(t%Control(i,:),t%Points(i,:))
                D2 = panel_distance(t%Control(i,:),t%Points(i+1,:))
                percent = D1/(D1+D2)
                t%Control(i,:) = t%Points(i,:) + percent*(t%Points(i+1,:)-t%Points(i,:))
                li = sqrt((t%Points(i,1)-t%Points(i+1,1))**2 + (t%Points(i,2)-t%Points(i+1,2))**2)
                t%Normals(i,1) = - (t%Points(i+1,2)-t%Points(i,2))/li
                t%Normals(i,2) =   (t%Points(i+1,1)-t%Points(i,1))/li
            end do
        else if(t%PointsPlacement.eq.3) then !Cosine Clustered on Surface using Panel Normal
            dtheta = pi/(real((t%npts-1))/2.0)
            t%Control(t%npts/2,:) = 0.0
            do i=1,t%npts/2-1
                x = 0.5*(1.0-cos((real(i))*dtheta))
                CALL panel_geometry(t,x,xc,yc,t%Control(t%npts/2+i,1),t%Control(t%npts/2+i,2),t%Control(t%npts/2-i,1),&
                                    &t%Control(t%npts/2-i,2))
            end do
            do i=1,t%npts-1
                li = sqrt((t%Points(i,1)-t%Points(i+1,1))**2 + (t%Points(i,2)-t%Points(i+1,2))**2)
                t%Normals(i,1) = - (t%Points(i+1,2)-t%Points(i,2))/li
                t%Normals(i,2) =   (t%Points(i+1,1)-t%Points(i,1))/li
            end do
        else if(t%PointsPlacement.eq.4) then !Cosine Clustered on Surface using Surface normal
            dtheta = pi/(real((t%npts-1))/2.0)
            t%Control(t%npts/2,:) = 0.0
            CALL panel_surface_normal(t,t%LE(1),t%Normals(t%npts/2,:),t%Normals(t%npts/2,:))
            do i=1,t%npts/2-1
                x = 0.5*(1.0-cos((real(i))*dtheta))
                CALL panel_geometry(t,x,xc,yc,t%Control(t%npts/2+i,1),t%Control(t%npts/2+i,2),t%Control(t%npts/2-i,1),&
                                    &t%Control(t%npts/2-i,2))
                CALL panel_surface_normal(t,x,t%Normals(t%npts/2+i,:),t%Normals(t%npts/2-i,:))
            end do
        end if
        do i=1,t%npts
            CALL panel_geometry(t,t%Points(i,1),t%Camber(i,1),t%Camber(i,2),dummy,dummy,dummy,dummy)
        end do

        !write airfoil to airfoil.txt
!        open(unit = 10, File = 'geom.txt', status="replace", action = "write", iostat = ierror)
!        write(10,*) '   xpoint                   ypoint                   xcamber                  ycamber               &
!                    &   xcontrol                 ycontrol                 xnorm                    ynorm'
!        do i=1,t%npts-1
!            write(10,100) t%Points(i,:),t%Camber(i,:),t%Control(i,:),t%Normals(i,:)
!        end do
!        write(10,100) t%Points(t%npts,:),t%Camber(t%npts,:)
!        close(10)
    end subroutine panel_setup

!-----------------------------------------------------------------------------------------------------------
    SUBROUTINE panel_write_airfoil(t)
        implicit none
        type(panel_t) :: t
        integer :: i,ierror
        CHARACTER(LEN=100)::fn
    100 FORMAT (1X, 1000ES25.16)
        write(fn,*) t%npts
        fn = trim(adjustl(t%naca))//'_'//trim(adjustl(fn))//'_pts.txt'
        
        open(unit = 10, File = fn, status="replace", action = "write", iostat = ierror)
        write(10,*) t%npts
        do i=1,t%npts
            write(10,100) t%Points(i,:)
        end do
        close(10)
        write(*,*) 'Airfoil written to file: ',fn
    end subroutine panel_write_airfoil

!-----------------------------------------------------------------------------------------------------------
    SUBROUTINE panel_read_airfoil(t)
        implicit none
        type(panel_t) :: t
        integer :: i,ierror
        CHARACTER(LEN=100)::fn
        
        open(unit = 10, File = t%readfilename, action = "read", iostat = ierror)
        read(10,*) t%npts
        call panel_allocate(t)
        do i=1,t%npts
            read(10,*) t%Points(i,:)
        end do
        close(10)
    end subroutine panel_read_airfoil
    
!-----------------------------------------------------------------------------------------------------------
    SUBROUTINE panel_geometry(t,x,xc,yc,xu,yu,xl,yl) !--------- NACA 4-Digit
        IMPLICIT NONE
        type(panel_t) :: t
        REAL :: x,xc,yc,xu,yu,xl,yl,dydx,thick,theta
        
        xc = x
        if(x.eq.0.0) then
            yc = 0.0; xu = 0.0; yu = 0.0; xl = 0.0; yl = 0.0
        else
            if(t%load.eq.'UL') then !Uniform Load Camber Line
                if(x.eq.1.0) then
                    yc = 0.0
                    dydx = -t%CLd/(4.0*pi)*(log(1.0-1.0e-15)-log(1.0e-15))
                else
                    yc = -t%CLd/(4.0*pi)*((1.0-x)*log(1.0-x) + x*log(x))
                    dydx = -t%CLd/(4.0*pi)*(log(x)-log(1.0-x))
                end if
            else    !Standard NACA 4-digit camber line
                if(x < t%xmc) then
                    yc = t%ymc*(2.0*(x/t%xmc) - (x/t%xmc)**2)
                    dydx = 2.0*t%ymc/t%xmc*(1.0 - x/t%xmc)
                else
                    yc = t%ymc*(2.0*(1.0-x)/(1.0-t%xmc) - ((1.0-x)/(1.0-t%xmc))**2)
                    dydx = -2.0*t%ymc/(1.0-t%xmc)*(1.0 - (1.0-x)/(1.0-t%xmc))
                end if
            end if
            
            if(t%OpenTE .eq. 1) then
                thick = t%tm*(2.969*sqrt(x) - 1.260*x - 3.516*x**2 + 2.843*x**3 - 1.015*x**4) !Original NACA
            else
                thick = t%tm*(2.969*sqrt(x) - 1.260*x - 3.523*x**2 + 2.836*x**3 - 1.022*x**4) !Modified to close TE
            end if
            xu = x  - thick/(2.0*sqrt(1.0 + dydx**2))*dydx
            yu = yc + thick/(2.0*sqrt(1.0 + dydx**2))
            xl = x  + thick/(2.0*sqrt(1.0 + dydx**2))*dydx
            yl = yc - thick/(2.0*sqrt(1.0 + dydx**2))
        end if
        
    END SUBROUTINE panel_geometry

!-----------------------------------------------------------------------------------------------------------
    SUBROUTINE panel_velocity(t,x,y,Vx,Vy,Cp)
    !   This subroutine returns the axial and normal components of velocity and the pressure coefficient (Vx,Vy,Cp)
    !   at the Cartesian coordinates (x,y)
        IMPLICIT NONE
        type(panel_t) :: t
        REAL :: x,y,Vx,Vy,Cp,V(2),P(2,2)
        INTEGER :: i

        V(:) = t%Vinf(:)
        do i=1,t%npts-1
            CALL panel_P_mat(t%Points(i,1),t%Points(i,2),t%Points(i+1,1),t%Points(i+1,2),x,y,P)
            V(:) = V(:) + matmul(P,t%Gamma(i:i+1))
        end do
        Vx = V(1); Vy = V(2)
        Cp = 1.0-(Vx**2 + Vy**2)/t%Vmag**2
    END SUBROUTINE panel_velocity

!-----------------------------------------------------------------------------------------------------------
    SUBROUTINE panel_coefficients(t,Mach,CL,CMle,CMc4)
        implicit none
        type(panel_t) :: t
        integer :: i
        real :: Mach, CL, CMle,CMc4,li,x1,x2,y1,y2,G1,G2
        
        CL = 0.0
        CMle = 0.0
        do i=1,t%npts-1
            x1 = t%Points(i,1); x2 = t%Points(i+1,1)
            y1 = t%Points(i,2); y2 = t%Points(i+1,2)
            G1 = t%Gamma(i)
            G2 = t%Gamma(i+1)
            li = sqrt((x2-x1)**2 + (y2-y1)**2)
            CL = CL + li*(G1+G2)/t%Vmag
            CMle = CMle + li/t%Vmag*((2.0*x1*G1 + x1*G2 + x2*G1 + 2.0*x2*G2)*cos(t%alpha) &
                               & + (2.0*y1*G1 + y1*G2 + y2*G1 + 2.0*y2*G2)*sin(t%alpha))
        end do
        CMle = -1.0/3.0*CMle
        CMc4 = CMle + 0.25*CL*cos(t%alpha)
        
        if(Mach < 1.0) then !Modify for subsonic flight
            CL = CL/sqrt(1.0-Mach**2)
            CmLE = CmLE/sqrt(1.0-Mach**2)
            Cmc4 = Cmc4/sqrt(1.0-Mach**2)
        end if

    END SUBROUTINE panel_coefficients

!-----------------------------------------------------------------------------------------------------------
    SUBROUTINE panel_solve(t)
        IMPLICIT NONE
        type(panel_t) :: t
        integer :: i,j,ierror
        real :: x,xc,yc,dtheta
        real :: P(2,2) !,Amat(t%npts,t%npts),B(t%npts),Ainv(t%npts,t%npts)
        
        real,allocatable,dimension(:,:) :: Amat,Ainv
        real,allocatable,dimension(:) :: Bvec
        
        allocate(Amat(t%npts,t%npts))
        allocate(Ainv(t%npts,t%npts))
        allocate(Bvec(t%npts))
        
        t%Vinf(1) = cos(t%alpha)
        t%Vinf(2) = sin(t%alpha)
        t%Vmag    = 1.0

        Amat(:,:) = 0.0; Bvec(:) = 0.0; Ainv(:,:) = 0.0
        do i=1,t%npts-1
            do j=1,t%npts-1
                CALL panel_P_mat(t%Points(j,1),t%Points(j,2),t%Points(j+1,1),&
                                &t%Points(j+1,2),t%Control(i,1),t%Control(i,2),P)
                Amat(i,j) =   Amat(i,j)   + t%Normals(i,2)*P(2,1) + t%Normals(i,1)*P(1,1)
                Amat(i,j+1) = Amat(i,j+1) + t%Normals(i,2)*P(2,2) + t%Normals(i,1)*P(1,2)
            end do
            Bvec(i) = -dot_product(t%Vinf(:),t%Normals(i,:))
        end do
        Amat(t%npts,1) = 1.0; Amat(t%npts,t%npts) = 1.0; Bvec(t%npts) = 0.0

!    100 FORMAT (1X, 1000ES25.16)
!    open(unit = 20, File = 'A.txt', status="replace", action = "write", iostat = ierror)
!    do i=1,t%npts
!        write(20,100) Amat(i,:)
!    end do
!    close(20)
!    open(unit = 20, File = 'B.txt', status="replace", action = "write", iostat = ierror)
!    do i=1,t%npts
!        write(20,100) Bvec(i)
!    end do
!    close(20)

        CALL math_matinv(t%npts,Amat,Ainv)
        t%Gamma = matmul(Ainv,Bvec)
!        call math_AXB_LUD(t%npts,Amat,Bvec,t%Gamma) ! Either works fine.
        deallocate(Amat)
        deallocate(Ainv)
        deallocate(Bvec)
    END SUBROUTINE panel_solve

!-----------------------------------------------------------------------------------------------------------
    SUBROUTINE panel_P_mat(x1,y1,x2,y2,x,y,P)
        IMPLICIT NONE
        real :: x1,y1,x2,y2,x,y,P(2,2),li,A(2,2),B(2,2),vec1(2),vec2(2)
        real :: dx,dy,phi,psi,eta,xi
        li = sqrt((x2-x1)**2 + (y2-y1)**2)
        dx = x2-x1; dy = y2-y1
        A(1,1) =  dx; A(1,2) = dy
        A(2,1) = -dy; A(2,2) = dx
        vec1(1) = x-x1; vec1(2) = y-y1
        vec2 = 1.0/li*matmul(A,vec1)
        xi = vec2(1); eta = vec2(2)
        phi = atan2(eta*li,eta**2 + xi**2 - xi*li)
        psi = 0.5*log((xi**2 + eta**2)/((xi-li)**2 + eta**2))
        A(1,1) = dx; A(1,2) = -dy
        A(2,1) = dy; A(2,2) =  dx
!        B(1,1) = -eta*phi + (li-xi)*psi + li;     B(1,2) = eta*phi + xi*psi - li   ! Source Panels
!        B(2,1) = (li-xi)*phi + eta*psi;           B(2,2) = xi*phi - eta*psi        ! Source Panels
        B(1,1) = (li-xi)*phi + eta*psi;     B(1,2) = xi*phi - eta*psi               ! Vortex Panels
        B(2,1) = eta*phi - (li-xi)*psi - li; B(2,2) = -eta*phi - xi*psi + li        ! Vortex Panels

        P = 1.0/(2.0*pi*li**2)*matmul(A,B)
            
    END SUBROUTINE panel_P_mat

!--------------------------------- General Functions and Subroutines
REAL FUNCTION panel_distance(P1,P2)
    implicit none
    real :: P1(2),P2(2)
    panel_distance = sqrt((P2(1)-P1(1))**2 + (P2(2)-P1(2))**2)
END FUNCTION panel_distance

!-----------------------------------------------------------------------------------------------------------
SUBROUTINE panel_aseq(t,astart,aend,da)
    implicit none
    type(panel_t) :: t
    integer :: i,ierror
    real :: alfa,astart,aend,da,CL,CMle,CMc4
    100 FORMAT (1X, 1000ES25.16)
    
    open(unit = 10, File = 'aseq.txt', status="replace", action = "write", iostat = ierror)
    write(10,*) '  alpha (deg)         CL                  CM_le               CM_c/4'
    alfa = astart
    do while (alfa <= aend)
        t%alpha = pi/180.0*alfa
        CALL panel_solve(t)
        CALL panel_coefficients(t,0.0,CL,CMle,CMc4)
        write(10,100) alfa,CL,CMle,CMc4
        alfa = alfa + da
    end do
    close(10)
END SUBROUTINE panel_aseq

!-----------------------------------------------------------------------------------------------------------
SUBROUTINE panel_plot_pressure(t)
    implicit none
    type(panel_t) :: t
    real :: dtheta,x,xc,yc,xu,yu,xl,yl,Vx,Vy,Cp,offset
    integer :: i,ierror
    100 FORMAT (1X, 1000ES25.16)
    
    dtheta = 2.0*pi/real(t%npts-1)
    offset = 1.0e-3
    open(unit = 20, File = 'pressure.txt', status="replace", action = "write", iostat = ierror)
    write(20,*) '   x                        y                        Vx                       Vy                       Cp'
    do i=1,t%npts-1
        CALL panel_velocity(t,t%Control(i,1)+offset*t%Normals(i,1),t%Control(i,2)+offset*t%Normals(i,2),Vx,Vy,Cp)
        write(20,100) t%Control(i,:),Vx,Vy,Cp
    end do
    close(20)

    open(unit = 20, File = 'cp.plt', status="replace", action = "write", iostat = ierror)
    write(20,*) 'set encoding iso_8859_1'
    write(20,*) 'set terminal postscript eps enhanced monochrome solid "Times-Roman" 20'
    write(20,*) 'set output "cp.eps"'
    write(20,*) 'set multiplot'
    write(20,*) 'set size 1,0.75'
    write(20,*) 'set origin 0.0,0.0'
    write(20,*) 'set format x "%4.1f"'
    write(20,*) 'set format y "%4.1f"'
    write(20,*) 'set xlabel "{/Times-Italic x/c}"'
    write(20,*) 'set ylabel "{/Times-Italic C_p}"'
    write(20,*) 'set bmargin 3'
    write(20,*) 'set tmargin 0'
    write(20,*) 'set xrange [-0.05:1.05]'
    write(20,*) 'set yrange [1.0:-1.1]'
    write(20,*) 'plot "pressure.txt" u 1:5 noti w lines'
    write(20,*) 'set size 1,.204851752'
    write(20,*) 'set origin 0.0,0.75'
    write(20,*) 'set format x ""'
    write(20,*) 'set xlabel ""'
    write(20,*) 'set ylabel "{/Times-Italic y/c}"'
    write(20,*) 'set ytic -0.1,0.1,0.1'
    write(20,*) 'set bmargin 0'
    write(20,*) 'set tmargin 0'
    write(20,*) 'set xrange [-0.05:1.05]'
    write(20,*) 'set yrange [-0.1:0.1]'
    write(20,*) 'plot "pressure.txt" u 1:2 noti w lines'
    write(20,*) 'unset multiplot'
    close(20)

END SUBROUTINE panel_plot_pressure

!-----------------------------------------------------------------------------------------------------------
SUBROUTINE panel_plot_streamlines(t,xstart,x_lower_limit,x_upper_limit,ds,nlines,delta)
    IMPLICIT NONE
    type(panel_t) :: t
    REAL :: xstart,x_lower_limit,x_upper_limit,ds,delta
    REAL :: xend,yend,x,Vu,Vl,fx,f2,xnew,m,xstag,ystag,dx,Vu1,Vl1,Vu2,Vl2,Vx,Vy,Cp,x2
    REAL :: xc,yc,xu,yu,xl,yl,ystart,nu(2),nl(2),norm(2),dtheta
    INTEGER :: nlines,i,ierror
    CHARACTER(LEN=5)::search
    CHARACTER(LEN=10)::fn
    100 FORMAT (1X, 1000ES25.16)

    !Plot geometry
    dtheta = 2.0*pi/real(t%npts-1)
    open(unit = 10, File = 'geom_old.txt', status="replace", action = "write", iostat = ierror)
    write(10,*) '  x                   x_camber            y_camber            x_upper &
                &            y_upper             x_lower             y_lower'
    do i=1,t%npts/2
        x = 0.5*(1.0-cos((real(i)-0.5)*dtheta))
        CALL panel_geometry(t,x,xc,yc,xu,yu,xl,yl)
        write(10,100) x,xc,yc,xu,yu,xl,yl
    end do
    close(10)

    dx = 1.0e-1
    !Find forward stagnation streamline
    x = t%LE(1)
    if(t%readfile.eq.1) then
        xstag = -1.0e-14; ystag = 0.0
    else
        CALL panel_surface_tangential_velocity(t,x,Vu,Vl)
        fx = Vu
        if(abs(fx)<1.0e-12) then !Stagnation at leading edge
            write(*,*) 'Forward stagnation point at leading edge'
            CALL panel_geometry(t,x,xstag,ystag,xu,yu,xl,yl)
            CALL panel_surface_normal(t,x,nl,norm)
        else
            if(fx > 0.0) search = 'lower'
            if(fx < 0.0) search = 'upper'
            write(*,*) 'searching ',search,' surface for forward stagnation point.'
            x2 = x + dx
            do while(abs(fx) > 1.0e-8)
                CALL panel_surface_tangential_velocity(t,x2,Vu,Vl)
                if(search == 'upper') f2 = Vu
                if(search == 'lower') f2 = Vl
                m = (f2-fx)/(x2-x)
                xnew = x - fx/m
                if(xnew.ne.xnew) then
                    write(*,*) 'Note: Leading stagnation point not found!'
                    exit
                end if
                if(xnew < 0.0) xnew = 0.0
                x = x2; x2 = xnew; fx = f2
            end do
            if(search == 'upper') then
                CALL panel_geometry(t,x,xc,yc,xstag,ystag,xl,yl)
                CALL panel_surface_normal(t,x,norm,nl)
            end if
            if(search == 'lower') then
                CALL panel_geometry(t,x,xc,yc,xu,yu,xstag,ystag)
                CALL panel_surface_normal(t,x,nl,norm)
            end if
        end if
        xstag = xstag + 0.01*norm(1); ystag = ystag + 0.01*norm(2)
    end if
    CALL panel_streamline(t,0,xstag,ystag,x_lower_limit,x_upper_limit,-ds,xend,yend)
    ystart = yend
    
    !aft stagnation point at trailing edge
    x = t%TE(1)
    write(*,*) 'Aft stagnation point at trailing edge'
    if(t%readfile.eq.1) then
        xstag = t%Points(1,1); ystag = t%Points(1,2)
        norm(:) = 0.5*(t%Normals(1,:) + t%Normals(t%npts-1,:))
        call panel_normalize(norm)
    else
        CALL panel_geometry(t,x,xstag,ystag,xu,yu,xl,yl)
        CALL panel_surface_normal(t,x,nl,norm)
    end if
    write(*,*) 'Trailing edge normal: ',norm(:)
    xstag = xstag + 0.01*norm(1); ystag = ystag + 0.01*norm(2)
    CALL panel_streamline(t,2*nlines+1,xstag,ystag,x_lower_limit,x_upper_limit,ds,xend,yend)

    !Plot the rest of the streamlines
    write(*,*) 'Calculating Streamlines'
    do i=1,nlines
        CALL panel_streamline(t,i,xstart,ystart+real(i)*delta/cos(t%alpha),x_lower_limit,x_upper_limit,ds,xend,yend)
        CALL panel_streamline(t,nlines+i,xstart,ystart-real(i)*delta/cos(t%alpha),x_lower_limit,x_upper_limit,ds,xend,yend)
    end do
    
    !write gnuplot file
    open(unit = 10, File = 'stream.plt', status="replace", action = "write", iostat = ierror)
    write(10,*) 'set encoding iso_8859_1'
    write(10,*) 'set terminal postscript eps enhanced monochrome solid "Times-Roman" 20'
    write(10,*) 'set output "stream.eps"'
    write(10,*) 'set size square'
    write(10,*) 'set format x "%4.1f"'
    write(10,*) 'set format y "%4.1f"'
    write(10,*) 'set xlabel "{/Times-Italic x}"'
    write(10,*) 'set label 1 "{/Times-Italic y}" at graph -.115, graph .5'
    write(10,*) 'set xrange [-0.2:1.2]'
    write(10,*) 'set yrange [-0.7:0.7]'
    write(10,*) 'plot "geom.txt" u 1:2 noti w l,\'
    write(10,*) '     "geom.txt" u 3:4 noti w l,\'
    write(10,*) '     "geom.txt" u 1:2 noti w p pt 6 ps 0.2,\'
    do i=0,2*nlines
        write(fn,*) i
        fn = trim(adjustl(fn))//'.txt"'
        write(10,*) '     "',fn,' u 1:2 noti w l,\'
    end do
    write(fn,*) 2*nlines+1
    fn = trim(adjustl(fn))//'.txt"'
    write(10,*) '     "',fn,' u 1:2 noti w l'
    close(10)
END SUBROUTINE panel_plot_streamlines

!-----------------------------------------------------------------------------------------------------------
SUBROUTINE panel_streamline(t,n,xstart,ystart,x_lower_limit,x_upper_limit,ds,xend,yend)
!   This subroutine plots a single streamline from the starting point (xstart,ystart) to the end point where x=xlimit
!   and returns the end point as (xend,yend).
!   The step size, ds, is measured along the streamline and may be either positive or negative.
!   The solution is based on a fourth-order Runge-Kutta numerical integration.
    IMPLICIT NONE
    type(panel_t) :: t
    INTEGER :: n,ierror
    REAL :: xstart,ystart,x_lower_limit,x_upper_limit,ds,xend,yend
    REAL :: pos1(2),pos2(2),s
    CHARACTER(LEN=100)::fn
    
    100 FORMAT (1X, 1000ES30.16)
    write(fn,*) n
    fn = trim(adjustl(fn))//'.txt'
    open(unit = 10, File = fn, status="replace", action = "write", iostat = ierror)
    write(10,*) '   x              y'
    s = 0.0
    pos1(1) = xstart; pos1(2) = ystart
    write(10,100) pos1
    do while((pos1(1)>=x_lower_limit).and.(pos1(1)<=x_upper_limit))
        CALL panel_rnkta4(t,2,s,pos1,ds,pos2)
        write(10,100) pos2
        pos1=pos2
        s = s + ds
    end do
    close(10)
    xend = pos1(1); yend = pos1(2)

END SUBROUTINE panel_streamline

!-----------------------------------------------------------------------------------------------------------
integer function panel_make_wake_streamline(t,x_limit,wakeLength)
    IMPLICIT NONE
    type(panel_t) :: t
    real :: x_limit,wakeLength,x,xstart,ystart,norm(2),ds,pos1(2),pos2(2),s,percent
    integer :: i,nWakePoints

    write(*,*) 'Computing wake streamline from trailing edge...'
    ds = 0.003*wakeLength
    !aft stagnation point at trailing edge
    x = t%TE(1)
    xstart = t%Points(1,1); ystart = t%Points(1,2)
    norm(:) = 0.5*(t%Normals(1,:) + t%Normals(t%npts-1,:))
    call panel_normalize(norm)
!    write(*,*) '    Trailing edge normal: ',norm(:)
    
    !Run once to find how long
    nWakePoints = 2
    s = 0.0
    pos1(1) = xstart + ds*norm(1); pos1(2) = ystart + ds*norm(2)
!    write(*,*) pos1
    do while(pos1(1)<x_limit)
        nWakePoints = nWakePoints + 1
        CALL panel_rnkta4(t,2,s,pos1,ds,pos2)
!        write(*,*) pos2
        pos1=pos2
        s = s + ds
    end do

    allocate(t%Wake(nWakePoints,2))
!    write(*,*) 'Second Run'
    !Run again to store data in t%Wake
    s = 0.0
    pos1(1) = xstart + ds*norm(1); pos1(2) = ystart + ds*norm(2)
    t%Wake(1,:) = [xstart,ystart]
    t%Wake(2,:) = pos1(:)
    do i=3,nWakePoints
        CALL panel_rnkta4(t,2,s,pos1,ds,pos2)
        t%Wake(i,:) = pos2(:)
        pos1=pos2
        s = s + ds
!        write(*,*) t%Wake(i,:)
    end do
    
    !Cut off end to be exactly x_limit
    i = nWakePoints
    percent = (x_limit - t%Wake(i-1,1)) / (t%Wake(i,1) - t%Wake(i-1,1))
    t%Wake(nWakePoints,:) = percent*(t%Wake(i,:) - t%Wake(i-1,:)) + t%Wake(i-1,:)
    panel_make_wake_streamline = nWakePoints
end function panel_make_wake_streamline

!-----------------------------------------------------------------------------------------------------------
SUBROUTINE panel_surface_normal(t,x,nu,nl)
    IMPLICIT NONE
    type(panel_t) :: t
    REAL :: x,nu(2),nl(2)
    REAL :: dx,xc,yc,xu1,yu1,xl1,yl1,xu2,yu2,xl2,yl2,vec(2)

    dx = 1.0e-12
    if(abs(x-t%LE(1)) < 1.0e-10) then !leading edge
        CALL panel_geometry(t,x+dx,xc,yc,xu1,yu1,xl1,yl1)
        vec(1) = xu1-xl1; vec(2) = yu1-yl1
        CALL panel_cross_product_2D(vec,nu)
        CALL panel_normalize(nu)
        nu = - nu
        nl = nu
    else if(abs(x-t%TE(1)) < 1.0e-10) then !trailing edge
        CALL panel_geometry(t,x-dx,xc,yc,xu1,yu1,xl1,yl1)
        vec(1) = xu1-xl1; vec(2) = yu1-yl1
        CALL panel_cross_product_2D(vec,nu)
        CALL panel_normalize(nu)
        nl = nu
    else
        CALL panel_geometry(t,x-dx,xc,yc,xu1,yu1,xl1,yl1)
        CALL panel_geometry(t,x+dx,xc,yc,xu2,yu2,xl2,yl2)
        vec(1) = xu2-xu1; vec(2) = yu2-yu1
        CALL panel_cross_product_2D(vec,nu)
        CALL panel_normalize(nu)
        nu = -nu
        vec(1) = xl2-xl1; vec(2) = yl2-yl1
        CALL panel_cross_product_2D(vec,nl)
        CALL panel_normalize(nl)
    end if
END SUBROUTINE panel_surface_normal

!-----------------------------------------------------------------------------------------------------------
SUBROUTINE panel_surface_tangent(t,x,tu,tl)
    IMPLICIT NONE
    type(panel_t) :: t
    REAL :: x,tu(2),tl(2)
    REAL :: dx,xc,yc,xu1,yu1,xl1,yl1,xu2,yu2,xl2,yl2

    dx = 1.0e-12
    if(abs(x-t%LE(1)) < 1.0e-10) then !leading edge
        CALL panel_geometry(t,x+dx,xc,yc,xu1,yu1,xl1,yl1)
        tu(1) = xu1-xl1; tu(2) = yu1-yl1
        CALL panel_normalize(tu)
        tl = tu
    else if(abs(x-t%TE(1)) < 1.0e-10) then !trailing edge
        CALL panel_geometry(t,x-dx,xc,yc,xu1,yu1,xl1,yl1)
        tu(1) = xl1-xu1; tu(2) = yl1-yu1
        CALL panel_normalize(tu)
        tl = tu
    else
        CALL panel_geometry(t,x-dx,xc,yc,xu1,yu1,xl1,yl1)
        CALL panel_geometry(t,x+dx,xc,yc,xu2,yu2,xl2,yl2)
        tu(1) = xu2-xu1; tu(2) = yu2-yu1
        CALL panel_normalize(tu)
        tl(1) = xl1-xl2; tl(2) = yl1-yl2
        CALL panel_normalize(tl)
    end if
END SUBROUTINE panel_surface_tangent

!-----------------------------------------------------------------------------------------------------------
SUBROUTINE panel_surface_tangential_velocity(t,x,Vu,Vl)
    IMPLICIT NONE
    type(panel_t) :: t
    REAL :: x,Vu,Vl,xc,yc,xu,yu,xl,yl,tu(2),tl(2),vec(2),Vx,Vy,Cp
    
    CALL panel_geometry(t,x,xc,yc,xu,yu,xl,yl)
    CALL panel_surface_tangent(t,x,tu,tl)

    CALL panel_velocity(t,xu,yu,Vx,Vy,Cp)
    vec(1) = Vx; vec(2) = Vy
    Vu = dot_product(tu,vec)

    CALL panel_velocity(t,xl,yl,Vx,Vy,Cp)
    vec(1) = Vx; vec(2) = Vy
    Vl = dot_product(tl,vec)
END SUBROUTINE panel_surface_tangential_velocity

!-----------------------------------------------------------------------------------------------------------
SUBROUTINE panel_normalize(v)
    IMPLICIT NONE
    REAL :: v(2),mag
    mag = sqrt(v(1)**2 + v(2)**2)
    v = v/mag
END SUBROUTINE panel_normalize

!-----------------------------------------------------------------------------------------------------------
SUBROUTINE panel_cross_product_2D(a,b)
    IMPLICIT NONE
    REAL :: a(2),b(2)
    b(1) =  a(2)
    b(2) = -a(1)
END SUBROUTINE panel_cross_product_2D


!-----------------------------------------------------------------------------------------------------------
      subroutine panel_rnkta4(t,n,t0,y0,dt,y)
      implicit none
!     This working precision subroutine computes a value for the n component
!     vector y(t0+dt) from a known value of the vector y(t0)=y0.  The solution
!     is based on a fourth order Runge-Kutta solution to the system of n
!     differential equations,
!
!                dy(i)/dt = f(i,a,t,y)              i = 1,2,3,...,n
!
!     where a is a coefficient array passed to the function f.  The working
!     precision function subprogram f(i,a,t,y) must be provided by the user.
!
      type(panel_t) :: t
      integer n,j,i
      real t0,y0(n),dt,y(n),c(4) !,dy(n,4)
      real,allocatable,dimension(:,:) :: dy
      allocate(dy(n,4))
      
      c(1) = 1.0/6.0
      c(2) = 1.0/3.0
      c(3) = c(2)
      c(4) = c(1)
      do j=1,n
         dy(j,1)=panel_fy(t,j,t0,y0)*dt
         y(j)=y0(j)+dy(j,1)/2.
      end do
      do j=1,n
         dy(j,2)=panel_fy(t,j,t0+dt/2.,y)*dt
      end do
      do j=1,n
         y(j)=y0(j)+dy(j,2)/2.
      end do
      do j=1,n
         dy(j,3)=panel_fy(t,j,t0+dt/2.,y)*dt
      end do
      do j=1,n
         y(j)=y0(j)+dy(j,3)
      end do
      do j=1,n
         dy(j,4)=panel_fy(t,j,t0+dt,y)*dt
      end do
      do j=1,n
         y(j)=y0(j)
         do i=1,4
            y(j)=y(j)+c(i)*dy(j,i)
         end do
      end do
      deallocate(dy)
      return
      end subroutine panel_rnkta4

!-----------------------------------------------------------------------------------------------------------
REAL FUNCTION panel_fy(t,i,s,pos)
    IMPLICIT NONE
    type(panel_t) :: t
    INTEGER :: i
    REAL :: s,pos(2),Vx,Vy,Cp,V
    CALL panel_velocity(t,pos(1),pos(2),Vx,Vy,Cp)
    V = sqrt(Vx**2 + Vy**2)
    if(i==1) panel_fy = Vx/V
    if(i==2) panel_fy = Vy/V
    
END FUNCTION panel_fy

end module panel_m
