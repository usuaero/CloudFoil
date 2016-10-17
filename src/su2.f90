module su2_m
    use grid_m
    implicit none
    

contains


!-----------------------------------------------------------------------------------------------------------
subroutine su2_setup(t)
    type(grid_t) :: t
    character(100) :: filename,turbModel,gradType,linearSolver,formatOut,MGlevel,MGcycle,convType,spacialOrder
    character(100) :: convTurb,spacialOrderTurb
    character(len=:),allocatable :: cval
    integer :: ios,maxIter,residualReduction,residualMinimum
    real :: CFL,MGdamp

    110 FORMAT (2ES25.16, I10)
    120 FORMAT (A, I10)

    filename = trim(t%master_filename)//'.cfg' ! use this for Windows Compilation
!    filename = 'SU2.cfg'                        ! use this for hpc compilation
    write(*,*) 'Saving SU2 Config to file: ',trim(filename)

    call t%json%get('settings.turbModel',    cval);   call json_check(); turbModel = trim(cval)
    call t%json%get('settings.gradType',    cval);   call json_check(); gradType = trim(cval)
    call t%json%get('settings.CFL',           CFL);   call json_check()
    call t%json%get('settings.linearSolver',    cval);   call json_check(); linearSolver = trim(cval)

    call t%json%get('settings.MGlevel',    cval);   call json_check(); MGlevel = trim(cval)
    call t%json%get('settings.MGcycle',    cval);   call json_check(); MGcycle = trim(cval)
    call t%json%get('settings.MGdamp',           MGdamp);   call json_check()

    call t%json%get('settings.convType',    cval);   call json_check(); convType = trim(cval)
    call t%json%get('settings.spacialOrder',    cval);   call json_check(); spacialOrder = trim(cval)

    call t%json%get('settings.convTurb',    cval);   call json_check(); convTurb = trim(cval)
    call t%json%get('settings.spacialOrderTurb',    cval);   call json_check(); spacialOrderTurb = trim(cval)

    call t%json%get('settings.residualReduction',  residualReduction);   call json_check()
    call t%json%get('settings.residualMinimum',  residualMinimum);   call json_check()
    call t%json%get('settings.maxIter',   maxIter);   call json_check()

    call t%json%get('settings.formatOut',    cval);   call json_check(); formatOut = trim(cval)

    open(unit = 100, File = trim(filename), action = "write", iostat = ios)
    write(100,'(a)') '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(100,'(a)') '% SU2 configuration file                                                       %'
    write(100,'(a)') '% Case description: 2D Airfoil (compressible)                                  %'
    write(100,'(a)') '% Author: CloudFoil / Utah State University AeroLab : aero.go.usu.edu          %'
    write(100,'(a)') '% Case Name: '//trim(t%tagName)
    write(100,'(a)') '% Case UUID: '//trim(t%tagUUID)//'                              %'
    write(100,'(a)') '% Case Date: '//trim(t%tagDate)//'                                          %'
    write(100,'(a)') '% File Version 4.3.0 "Cardinal"                                                %'
    write(100,'(a)') '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(100,*)
    write(100,*)
    write(100,'(a)') '% ------------- DIRECT, ADJOINT, AND LINEARIZED PROBLEM DEFINITION ------------%'
    if(trim(turbModel) .eq. 'EULER') then
        write(100,'(a)') 'PHYSICAL_PROBLEM= EULER'
    else
        write(100,'(a)') 'PHYSICAL_PROBLEM= NAVIER_STOKES'
        write(100,'(a)') 'KIND_TURB_MODEL= '//trim(turbModel)
    end if
    write(100,'(a)') 'MATH_PROBLEM= DIRECT'
    write(100,'(a)') 'RESTART_SOL= NO'
    write(100,'(a)') 'LOW_MEMORY_OUTPUT= NO'
    write(100,'(a)') 'REGIME_TYPE= COMPRESSIBLE'
    write(100,'(a)') 'SYSTEM_MEASUREMENTS= SI'
    write(100,*)
    write(100,'(a)') '% -------------------- COMPRESSIBLE FREE-STREAM DEFINITION --------------------%'
    write(100,'(a ES25.16)') 'MACH_NUMBER= ',t%MachNumber
    write(100,'(a ES25.16)') 'AoA= ',t%alpha*180.0/pi
    write(100,'(a)') 'INIT_OPTION= REYNOLDS'
    write(100,'(a)') 'FREESTREAM_OPTION= TEMPERATURE_FS'
    write(100,'(a ES25.16)') 'FREESTREAM_PRESSURE= ',t%pressure
    write(100,'(a ES25.16)') 'FREESTREAM_TEMPERATURE= ',t%temperature
    write(100,'(a ES25.16)') 'REYNOLDS_NUMBER= ',t%ReynoldsNumber
    write(100,'(a ES25.16)') 'REYNOLDS_LENGTH= ',t%afLength
    write(100,*)
    write(100,'(a)') '% ---------------------- REFERENCE VALUE DEFINITION ---------------------------%'
    write(100,'(a ES25.16)') 'REF_ORIGIN_MOMENT_X = ',0.25*t%afLength
    write(100,'(a ES25.16)') 'REF_ORIGIN_MOMENT_Y = ',0.0*t%afLength
    write(100,'(a ES25.16)') 'REF_ORIGIN_MOMENT_Z = ',0.0*t%afLength
    write(100,'(a ES25.16)') 'REF_LENGTH_MOMENT= ',t%afLength
    write(100,'(a ES25.16)') 'REF_AREA= ',t%afLength
    write(100,'(a)') 'REF_DIMENSIONALIZATION= FREESTREAM_PRESS_EQ_ONE'
    write(100,*)
    write(100,'(a)') '% ---- IDEAL GAS, POLYTROPIC, VAN DER WAALS AND PENG ROBINSON CONSTANTS -------%'
    write(100,'(a)') 'FLUID_MODEL= STANDARD_AIR'
    write(100,'(a)') 'GAMMA_VALUE= 1.4'
    write(100,'(a)') 'GAS_CONSTANT= 287.058'
    write(100,'(a)') 'CRITICAL_TEMPERATURE= 131.00'
    write(100,'(a)') 'CRITICAL_PRESSURE= 3588550.0'
    write(100,'(a)') 'CRITICAL_DENSITY= 263.0'
    write(100,'(a)') 'ACENTRIC_FACTOR= 0.035'
    write(100,*)
    write(100,'(a)') '% --------------------------- VISCOSITY MODEL ---------------------------------%'
    write(100,'(a)') 'VISCOSITY_MODEL= SUTHERLAND'
    write(100,'(a)') 'MU_CONSTANT= 1.716E-5'
    write(100,'(a)') 'MU_REF= 1.716E-5'
    write(100,'(a)') 'MU_T_REF= 273.15'
    write(100,'(a)') 'SUTHERLAND_CONSTANT= 110.4'
    write(100,*)
    write(100,'(a)') '% -------------------- BOUNDARY CONDITION DEFINITION --------------------------%'
    if(trim(turbModel) .eq. 'EULER') then
        write(100,'(a)') 'MARKER_EULER= ( airfoil )'
    else
        write(100,'(a)') 'MARKER_HEATFLUX= ( airfoil, 0.0 )'
    end if
    write(100,'(a)') 'MARKER_FAR= ( farfield, outlet )'
    write(100,'(a)') 'MARKER_PLOTTING= ( airfoil )'
    write(100,'(a)') 'MARKER_MONITORING= ( airfoil )'
    write(100,*)
    write(100,'(a)') '% ------------- COMMON PARAMETERS DEFINING THE NUMERICAL METHOD ---------------%'
    write(100,'(a)') 'NUM_METHOD_GRAD= '//trim(gradType)
    write(100,'(a ES25.16)') 'CFL_NUMBER= ',CFL
    write(100,'(a)') 'MAX_DELTA_TIME= 1E10'
    write(100,'(a)') 'CFL_ADAPT= NO'
    write(100,'(a)') 'CFL_ADAPT_PARAM= ( 1.5, 0.5, 1.0, 100.0 )'
    write(100,'(a)') 'RK_ALPHA_COEFF= ( 0.66667, 0.66667, 1.000000 )'
    write(100,*)
    write(100,'(a)') '% ----------------------- SLOPE LIMITER DEFINITION ----------------------------%'
    write(100,'(a ES25.16)') 'REF_ELEM_LENGTH= ',0.1*t%afLength
    write(100,'(a)') 'LIMITER_COEFF= 0.3'
    write(100,'(a)') 'LIMITER_ITER= 99999'
    write(100,'(a)') 'SHARP_EDGES_COEFF= 3'
    write(100,*)
    write(100,'(a)') '% ------------------------ LINEAR SOLVER DEFINITION ---------------------------%'
    write(100,'(a)') 'LINEAR_SOLVER= '//trim(linearSolver)
    write(100,'(a)') 'LINEAR_SOLVER_PREC= LU_SGS'
    write(100,'(a)') 'LINEAR_SOLVER_ERROR= 1E-6'
    write(100,'(a)') 'LINEAR_SOLVER_ITER= 5'
    write(100,*)
    write(100,'(a)') '% -------------------------- MULTIGRID PARAMETERS -----------------------------%'
    write(100,'(a)') 'MGLEVEL= '//trim(MGlevel)
    write(100,'(a)') 'MGCYCLE= '//trim(MGcycle)
    write(100,'(a)') 'MG_PRE_SMOOTH= ( 1, 2, 3, 3 )'
    write(100,'(a)') 'MG_POST_SMOOTH= ( 0, 0, 0, 0 )'
    write(100,'(a)') 'MG_CORRECTION_SMOOTH= ( 0, 0, 0, 0 )'
    write(100,'(a ES25.16)') 'MG_DAMP_RESTRICTION= ',MGdamp
    write(100,'(a ES25.16)') 'MG_DAMP_PROLONGATION= ',MGdamp
    write(100,*)
    write(100,'(a)') '% -------------------- FLOW NUMERICAL METHOD DEFINITION -----------------------%'
    write(100,'(a)') 'CONV_NUM_METHOD_FLOW= '//trim(convType)
    write(100,'(a)') 'SPATIAL_ORDER_FLOW= '//trim(spacialOrder)
    write(100,'(a)') 'SLOPE_LIMITER_FLOW= VENKATAKRISHNAN'
    write(100,'(a)') 'AD_COEFF_FLOW= ( 0.15, 0.5, 0.02 )'
    write(100,'(a)') 'TIME_DISCRE_FLOW= EULER_IMPLICIT'
    write(100,*)
    write(100,'(a)') '% -------------------- TURBULENT NUMERICAL METHOD DEFINITION ------------------%'
    write(100,'(a)') 'CONV_NUM_METHOD_TURB= '//trim(convTurb)
    write(100,'(a)') 'SPATIAL_ORDER_TURB= '//trim(spacialOrderTurb)
    write(100,'(a)') 'SLOPE_LIMITER_TURB= VENKATAKRISHNAN'
    write(100,'(a)') 'TIME_DISCRE_TURB= EULER_IMPLICIT'
    write(100,'(a)') 'CFL_REDUCTION_TURB= 1.0'
    write(100,*)
    write(100,'(a)') '% --------------------------- CONVERGENCE PARAMETERS --------------------------%'
    write(100,'(a I10)') 'EXT_ITER= ',maxIter
    write(100,'(a)') 'CONV_CRITERIA= RESIDUAL'
    write(100,'(a I3)') 'RESIDUAL_REDUCTION= ',residualReduction
    write(100,'(a I3)') 'RESIDUAL_MINVAL= ',residualMinimum
    write(100,'(a)') 'STARTCONV_ITER= 10'
    write(100,'(a)') 'CAUCHY_ELEMS= 100'
    write(100,'(a)') 'CAUCHY_EPS= 1E-6'
    write(100,'(a)') 'CAUCHY_FUNC_FLOW= DRAG'
    write(100,*)
    write(100,'(a)') '% ------------------------- INPUT/OUTPUT INFORMATION --------------------------%'
    write(100,'(a)') 'MESH_FILENAME= mesh.su2'
    write(100,'(a)') 'MESH_FORMAT= SU2'
    write(100,'(a)') 'MESH_OUT_FILENAME= mesh_out.su2'
    write(100,'(a)') 'SOLUTION_FLOW_FILENAME= solution_flow.dat'
    write(100,'(a)') 'SOLUTION_ADJ_FILENAME= solution_adj.dat'
    write(100,'(a)') 'OUTPUT_FORMAT= '//trim(formatOut)
    write(100,'(a)') 'CONV_FILENAME= history'
    write(100,'(a)') 'BREAKDOWN_FILENAME= forces_breakdown.dat'
    write(100,'(a)') 'RESTART_FLOW_FILENAME= restart_flow.dat'
    write(100,'(a)') 'RESTART_ADJ_FILENAME= restart_adj.dat'
    write(100,'(a)') 'VOLUME_FLOW_FILENAME= flow'
    write(100,'(a)') 'VOLUME_ADJ_FILENAME= adjoint'
    write(100,'(a)') 'GRAD_OBJFUNC_FILENAME= of_grad.dat'
    write(100,'(a)') 'SURFACE_FLOW_FILENAME= surface'
    write(100,'(a)') 'SURFACE_ADJ_FILENAME= surface_adjoint'
    write(100,'(a)') 'WRT_SOL_FREQ= 10000'
    write(100,'(a)') 'WRT_CON_FREQ= 10'

    close(100)

end subroutine su2_setup

end module su2_m
