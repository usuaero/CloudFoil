module view_m
    use grid_m
    implicit none
    
contains

!-----------------------------------------------------------------------------------------------------------
subroutine view_plotmtv(t,flag)
    type(grid_t) :: t
    character(100) :: filename,command
    integer :: i,j,flag(5)
    integer :: ierror = 0

    filename = trim(adjustl(t%master_filename))//'_PlotData.txt'
    open(unit = 10, File = filename, action = 'write', iostat = ierror)

    write(10,*) ' $ DATA=CURVE3D'
    write(10,*) ' % toplabel   = "Planform"'
    write(10,*) ' % xlabel     = "x"'
    write(10,*) ' % ylabel     = "y"'
    write(10,*) ' % zlabel     = "z"'
    write(10,*) ' % grid       = False'
!    write(10,*) ' % axislabel  = False'
    write(10,*) ' % equalscale = True'
    write(10,*) ' % axisscale  = False'
!    write(10,*) ' % axisguides = False'
    write(10,*) ' % eyepos.x   = 0.50'
    write(10,*) ' % eyepos.y   = 0.75'
    write(10,*) ' % eyepos.z   = 0.25'
!    write(10,*) ' % xmin   = -0.60'
!    write(10,*) ' % xmax   =  0.60'
!    write(10,*) ' % ymin   = -0.70'
!    write(10,*) ' % ymax   =  0.2'
    write(10,*) ' % dlinecolor = 6'
    write(10,*) ' '

    if(flag(1) .eq. 1) then !show organic grid
        write(*,*) 'Displaying organic grid in yellow.'
        do j=1,t%nj
            do i=1,t%ni-1
                write(10,*) ' % linecolor = 1' !Quarter Chord
                write(10,*) t%organicGrid(i,j,:),0.0
                write(10,*) t%organicGrid(i+1,j,:),0.0
                write(10,*) ' '
            end do
        end do

        do j=1,t%nj-1
            do i=1,t%ni
                write(10,*) ' % linecolor = 1' !Quarter Chord
                write(10,*) t%organicGrid(i,j,:),0.0
                write(10,*) t%organicGrid(i,j+1,:),0.0
                write(10,*) ' '
            end do
        end do
    end if

    if(flag(2) .eq. 1) then !show algebraic grid
        write(*,*) 'Displaying algebraic grid in light blue.'
        do j=1,t%nj
            do i=1,t%ni-1
                write(10,*) ' % linecolor = 2' !Quarter Chord
                write(10,*) t%algebraicGrid(i,j,:),0.0
                write(10,*) t%algebraicGrid(i+1,j,:),0.0
                write(10,*) ' '
            end do
        end do

        do j=1,t%nj-1
            do i=1,t%ni
                write(10,*) ' % linecolor = 2' !Quarter Chord
                write(10,*) t%algebraicGrid(i,j,:),0.0
                write(10,*) t%algebraicGrid(i,j+1,:),0.0
                write(10,*) ' '
            end do
        end do
    end if

    if(flag(3) .eq. 1) then !show master
        write(*,*) 'Displaying master grid in green.'
        do j=1,t%nj
            do i=1,t%ni-1
                write(10,*) ' % linecolor = 3' !Quarter Chord
                write(10,*) t%masterGrid(i,j,:),0.0
                write(10,*) t%masterGrid(i+1,j,:),0.0
                write(10,*) ' '
            end do
        end do

        do j=1,t%nj-1
            do i=1,t%ni
                write(10,*) ' % linecolor = 3' !Quarter Chord
                write(10,*) t%masterGrid(i,j,:),0.0
                write(10,*) t%masterGrid(i,j+1,:),0.0
                write(10,*) ' '
            end do
        end do
    end if

    if(flag(4) .eq. 1) then !show medium
        write(*,*) 'Displaying medium grid in red.'
        do j=1,(t%nj-1)/2+1
            do i=1,(t%ni-1)/2
                write(10,*) ' % linecolor = 4' !Quarter Chord
                write(10,*) t%mediumGrid(i,j,:),0.0
                write(10,*) t%mediumGrid(i+1,j,:),0.0
                write(10,*) ' '
            end do
        end do

        do j=1,(t%nj-1)/2
            do i=1,(t%ni-1)/2+1
                write(10,*) ' % linecolor = 4' !Quarter Chord
                write(10,*) t%mediumGrid(i,j,:),0.0
                write(10,*) t%mediumGrid(i,j+1,:),0.0
                write(10,*) ' '
            end do
        end do
    end if

    if(flag(5) .eq. 1) then !show coarse
        write(*,*) 'Displaying coarse grid in blue.'
        do j=1,(t%nj-1)/4+1
            do i=1,(t%ni-1)/4
                write(10,*) ' % linecolor = 5' !Quarter Chord
                write(10,*) t%coarseGrid(i,j,:),0.0
                write(10,*) t%coarseGrid(i+1,j,:),0.0
                write(10,*) ' '
            end do
        end do

        do j=1,(t%nj-1)/4
            do i=1,(t%ni-1)/4+1
                write(10,*) ' % linecolor = 5' !Quarter Chord
                write(10,*) t%coarseGrid(i,j,:),0.0
                write(10,*) t%coarseGrid(i,j+1,:),0.0
                write(10,*) ' '
            end do
        end do
    end if
    
    close(10)
!    command = '/sw/bin/Plotmtv -fg white -bg black '//trim(filename)//' > plotmtv_console.out'
!    call system(command)
    command = 'Plotmtv -fg white -bg black '//trim(filename)//' > plotmtv_console.out'
    call system(command)
end subroutine view_plotmtv

!-----------------------------------------------------------------------------------------------------------
end module view_m
