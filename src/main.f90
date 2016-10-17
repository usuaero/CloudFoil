!gfortran -fdefault-real-8 math.f90 dataset.f90 airfoil.f90 section.f90 wing.f90 plane.f90 special_functions.f90 view.f90 main.f90 -o AF3D.out
program main
    use grid_m
    use view_m
    use su2_m
    implicit none
    type(grid_t) :: mygrid
    type(json_value),pointer :: json_this, json_run, json_sub
    character(100) :: filename,run_command
    integer :: i,k,nRunCommands,viewGrid(5)

    real :: time1,time2
    
    call cpu_time(time1)
    write(*,*) '-----------------------------------------------'
    write(*,*) '|               CloudFoil 2.0                 |'
    write(*,*) '|                                             |'
    write(*,*) '|           (c) Doug Hunsaker, 2016           |'
    write(*,*) '|                                             |'
    write(*,*) '|          This software comes with           |'
    write(*,*) '| ABSOLUTELY NO WARRANTY EXPRESSED OR IMPLIED |'
    write(*,*) '|                                             |'
    write(*,*) '|        Submit bug reports online at:        |'
    write(*,*) '|               aero.go.usu.edu               |'
    write(*,*) '-----------------------------------------------'

    call get_command_argument(1,filename)
    mygrid%master_filename = filename
    
    !These are all done no matter what
    call grid_set_defaults(mygrid)
    call grid_load_json(mygrid)
    call grid_set_conditions(mygrid)
    call grid_create_airfoil_surface(mygrid)

    !Run Commands
    call mygrid%json%get('run', json_this)
    nRunCommands = json_value_count(json_this)

    do i=1,nRunCommands
        call json_value_get(json_this,i,json_run)
        run_command = trim(json_run%name)

        write(*,*)
        write(*,*) 'Running command : ',run_command

        select case (run_command)
            case('surface')
                call grid_save_surface_geometry(mygrid)
            case('characterize')
                call grid_characterize_airfoil(mygrid)
            case('mesh')
                call grid_create_wake_surface(mygrid)
                call grid_create_inner_grids(mygrid)
                call grid_merge(mygrid)
                call grid_coarsen(mygrid)
            case ('view')
                viewGrid = 0
                do k=1,json_value_count(json_run)
                    call json_value_get(json_run,k,json_sub)
                    if(trim(json_sub%name) .eq. 'organic')   viewGrid(1) = 1;
                    if(trim(json_sub%name) .eq. 'algebraic') viewGrid(2) = 1;
                    if(trim(json_sub%name) .eq. 'master')    viewGrid(3) = 1;
                    if(trim(json_sub%name) .eq. 'medium')    viewGrid(4) = 1;
                    if(trim(json_sub%name) .eq. 'coarse')    viewGrid(5) = 1;
                end do
                call view_plotmtv(mygrid,viewGrid)
            case ('saveSU2')
                do k=1,json_value_count(json_run)
                    call json_value_get(json_run,k,json_sub)
                    call grid_save_SU2(mygrid,trim(json_sub%name));
                end do
            case ('saveVTK')
                do k=1,json_value_count(json_run)
                    call json_value_get(json_run,k,json_sub)
                    call grid_save_VTK(mygrid,trim(json_sub%name));
                end do
            case ('saveWeb')
                do k=1,json_value_count(json_run)
                    call json_value_get(json_run,k,json_sub)
                    call grid_save_web(mygrid,trim(json_sub%name));
                end do
            case('setupSU2')
                call su2_setup(mygrid)
            case('command')
                do k=1,json_value_count(json_run)
                    call json_value_get(json_run,k,json_sub)
                    write(*,*) trim(json_sub%name)
                    call system(trim(json_sub%name));
                end do
            case default
                write(*,*) 'Command not found!'
        end select
    end do


    call grid_deallocate(mygrid)

    call cpu_time(time2)
    write(*,*) 'CPU time total (sec): ',time2-time1

end program main
