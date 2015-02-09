!*TRANSPORT CODE SOLVER INTERFACE FOR AHOTN DGFEM AND SCTSTEP
!> @brief Solves the neutron transport equation for the specified solver with the passed in input information
!> @detail Standard usage of these solveres should be through the pyne spatial_solver interface.
!> 
!Code originally by Sebastian Schunert, modified by Josh Howland.
!
!> @param titlein desired title for solver run 
!> @param solver_in spatial solver class
!> @param solver_type_in spatial solver specific type
!> @param spatial_order_in lambda order value
!> @param angular_quadrature_order_in quadrature order
!> @param angular_quadrature_type_in quadrature type
!> @param qdtypin type of quadrature scheme the code uses.
!> @param nodes_x_in number of nodes in x direction
!> @param nodes_y_in number of nodes in y direction
!> @param nodes_z_in number of nodes in z direction
!> @param num_groups_in number of groups of ?????????????
!> @param num_materials_in number of materials
!> @param x_cell_widths_in x node size
!> @param y_cell_widths_in y node size
!> @param z_cell_widths_in z node size
!> @param x_boundry_condition_1_in x south boundry conditions
!> @param x_boundry_condition_2_in x east boundry conditions
!> @param y_boundry_condition_1_in y north boundry conditions
!> @param y_boundry_condition_2_in y west boundry conditions
!> @param z_boundry_condition_1_in z upper boundry conditions
!> @param z_boundry_condition_2_in z lower boundry conditions
!> @param material_id_in material file in
!> @param quadrature_file
!> @param xs_file_in
!> @param source_input_file_in
!> @param bc_input_filein
!> @param flux_output_filein
!> @param convergence_criterion_in
!> @param max_iterations_in (itmxin)
!> @param moments_converged_in
!> @param converge_tolerence_in
!> @param max_mom_printed_in
!> @param moment_sum_flag_in
!> @param mom_at_a_pt_flag_in
!> @param quad_flux_print_flag_in

!OUTPUT NOT PARAM
!> @param fluxout: ouput array of final solution
!> @param error_code: error code if solver failed
!> @param tsolve_out: system time when solver began
!> @param ttosolve_out: system time when solver terminated
!> @param tend_out: solver runtime

SUBROUTINE main(qdfile, xsfile, srcfile, mtfile,inflow_file,phi_file, titlein,&
 solver_in, solver_type_in, spatial_order_in,&
 angular_quadrature_order_in, qdtypin, nodes_x_in, nodes_y_in, nodes_z_in,&
 num_groups_in, num_materials_in, x_cell_widths_in, y_cell_widths_in,&
 z_cell_widths_in, x_boundry_condition_1_in, x_boundry_condition_2_in,&
 y_boundry_condition_1_in, y_boundry_condition_2_in, z_boundry_condition_1_in,&
 z_boundry_condition_2_in, material_id_in, quadrature_file, xs_file_in,&
 source_input_file_in, bc_input_filein, flux_output_filein, &
 convergence_criterion_in, itmxin, moments_converged_in, converge_tolerence_in, &
 max_mom_printed_in, moment_sum_flag_in,&
 mom_at_a_pt_flag_in, quad_flux_print_flag_in,fluxout,error_code_out,&
 tsolve_out, ttosolve_out, tend_out) 

!-------------------------------------------------------------
!
!    Read the input parameters from input files
!
!    Comments below demonstrate order of the reading
!
!    Dependency: 
!           angle   = gets the angular quadrature data
!           readmt  = reads the material map from file
!           readxs  = reads the cross sections
!           readsrc = reads the source distribution
!           check   = input check on all the values
!
!    This code can by dynamically allocated. It uses a module to hold all 
!    input variables: invar
!
!
!    Solver types =  "AHOTN", "DGFEM", and "SCTSTEP"
!                    - AHOTN solvers: "LL" "LN" and "NEFD"
!                    - DGFEM solvers: "LD" "DENSE" and "LAGRANGE"
!         - SCTSTEP solvers: SCTSTEP (only one)
!
! Some problem size specifications that are passed in:
!   lambda => LAMDBA, the AHOT spatial order
!   qdord => Angular quadrature order
!   qdtyp => Angular quadrature type = 0/1/2 = TWOTRAN/EQN/Read-in
!   nx    => Number of 'x' cells
!   ny    => Number of 'y' cells
!   nz    => Number of 'z' cells
!   ng    => Number of groups
!   nm    => Number of materials
! NOTE: we should probably define or describe all of the input varialbes using 
! doxygen syntax so that the api will be clear.
!-------------------------------------------------------------

USE invar
USE solvar
USE timevar
USE precision_module
IMPLICIT NONE
  
INTEGER :: i, j, k, n
! File Names
CHARACTER(30), INTENT(OUT) :: qdfile, xsfile, srcfile, mtfile,inflow_file,&
                              phi_file
LOGICAL :: ex1, ex2, ex3, ex4
REAL*8 :: wtsum

CHARACTER(80), INTENT(IN) :: titlein
CHARACTER(30), INTENT(IN) :: solver_in, solver_type_in
INTEGER, INTENT(IN) :: spatial_order_in, angular_quadrature_order_in,&
 qdtypin, nodes_x_in, nodes_y_in, nodes_z_in, num_groups_in, num_materials_in

REAL*8, INTENT(IN), DIMENSION(:) :: x_cell_widths_in, y_cell_widths_in, z_cell_widths_in

INTEGER, INTENT(IN) :: x_boundry_condition_1_in, x_boundry_condition_2_in,&
 y_boundry_condition_1_in, y_boundry_condition_2_in, z_boundry_condition_1_in,&
 z_boundry_condition_2_in 

! Cell materials
INTEGER, INTENT(IN), DIMENSION(:,:,:) :: material_id_in

CHARACTER(30), INTENT(IN) :: quadrature_file, xs_file_in, source_input_file_in,&
 bc_input_filein, flux_output_filein

! Iteration Controls
REAL*8, INTENT(IN) :: convergence_criterion_in, converge_tolerence_in
INTEGER, INTENT(IN) :: itmxin, moments_converged_in

! Editing data
INTEGER, INTENT(IN) :: max_mom_printed_in, moment_sum_flag_in, mom_at_a_pt_flag_in,&
 quad_flux_print_flag_in

REAL*8, INTENT(OUT), DIMENSION(nodes_x_in,num_groups_in*nodes_y_in,num_groups_in*nodes_z_in) :: fluxout
! Works for all solvers!

INTEGER, INTENT(OUT) :: error_code_out
REAL*8, INTENT(OUT) :: tsolve_out, ttosolve_out, tend_out

! Set error codes to 0 initmoments_convergedy
error_code = 0
error_code_out = error_code

! Set all of the input values
title = titlein
solver = solver_in
solvertype = solver_type_in
lambda = spatial_order_in
qdord = angular_quadrature_order_in
qdtyp = qdtypin
nx = nodes_x_in
ny = nodes_y_in
nz = nodes_z_in
ng = num_groups_in
nm = num_materials_in
dx = x_cell_widths_in
dy = y_cell_widths_in
dz = z_cell_widths_in
xsbc = x_boundry_condition_1_in
xebc = x_boundry_condition_2_in
ysbc = y_boundry_condition_1_in
yebc = y_boundry_condition_2_in
zsbc = z_boundry_condition_1_in
zebc = z_boundry_condition_2_in
mat = material_id_in
inflow_file = bc_input_filein
phi_file = flux_output_filein
convergence_criterion = convergence_criterion_in
converge_tolerence = converge_tolerence_in
itmx = itmxin
moments_converged = moments_converged_in
momp = max_mom_printed_in
momsum = moment_sum_flag_in
mompt = mom_at_a_pt_flag_in
qdflx = quad_flux_print_flag_in

CALL version

IF (solver == "DGFEM") THEN
    IF (solvertype == "LD") THEN
        lambda=1
        WRITE(8,*) "DGFEM LN SOLVER" 
    ELSE IF (solvertype == "DENSE") THEN
        WRITE(8,*) "DGFEM DENSE SOLVER" 
    ELSE IF (solvertype == "LAGRANGE") THEN
        WRITE(8,*) "DGFEM LAGRANGE SOLVER"
    END IF
ELSE IF (solver == "AHOTN") THEN
    IF (solvertype == "LN" .or. solvertype == "LL") THEN
        IF (lambda .ne. 1) then
            WRITE(8,*) "ERROR: Lambda must be equal to one." 
            error_code_out = 1001
            RETURN
            !STOP
        END IF
    END IF
    IF (solvertype == "LN") THEN
        WRITE(8,*) "AHOTN LN SOLVER" 
    ELSE IF (solvertype == "LL") THEN
        WRITE(8,*) "AHOTN LL SOLVER" 
    ELSE IF (solvertype == "NEFD") THEN
        WRITE(8,*) "AHOTN NEFD SOLVER" 
    END IF
ELSE IF (solver == "SCTSTEP") THEN
    WRITE(8,*) "SCT STEP SOLVER" 
    lambda = 0
END IF

! Check that the order given is greater than zero and is even
IF (qdord <= 0) THEN
    WRITE(8,'(/,3x,A)') "ERROR: Illegal value for qdord. Must be greater than zero."
    error_code_out = 1002
    RETURN
    !STOP
ELSE IF (MOD(qdord,2) /= 0) THEN
    WRITE(8,'(/,3x,A)') "ERROR: Illegal value for the quadrature order. Even #s only."
    error_code_out = 1003
    return
    !STOP
END IF

!INQUIRE(FILE = xs_file_in, EXIST = ex1)
!INQUIRE(FILE = source_input_file_in, EXIST = ex2)
!IF (ex1 .eqv. .FALSE. .OR. ex2 .eqv. .FALSE.) THEN
!   WRITE(8,'(/,3x,A)') "ERROR: File does not exist for reading."
!   STOP
!END IF

! Set up the extra needed info from the read input
apo = (qdord*(qdord+2))/8
IF (solver == "AHOTN") THEN
    order = lambda+1
    ordsq = order**2
    ordcb = order**3
ELSE IF (solver == "DGFEM") THEN
    IF (solvertype == "LD") THEN
        dofpc = 4
    ELSE IF (solvertype == "DENSE") THEN
        dofpc = (lambda+3)*(lambda+2)*(lambda+1)/6
    ELSE IF (solvertype == "LAGRANGE") THEN
        order = lambda+1
        ordsq = order**2
        ordcb = order**3
    END IF
ELSE IF (solver=="SCTSTEP") THEN
  order = lambda+1
  ordsq = order**2
  ordcb = order**3
END IF

! Angular quadrature
ALLOCATE(ang(apo,3), w(apo))
ang = 0.0d0
   write(8,*) "angle initial: ", ang
IF (qdtyp == 2) THEN
  INQUIRE(FILE=quadrature_file, EXIST=ex3)
  IF (qdfile == '        ' .OR. ex3 .eqv. .FALSE.) THEN
    WRITE(8,'(/,3x,A)') "ERROR: illegal entry for the qdfile name."
    error_code_out = 1004
    RETURN
    !STOP
   END IF
   OPEN(UNIT=10, FILE=quadrature_file)
   READ(10,*)
   READ(10,*) (ang(n,1),ang(n,2),w(n),n=1,apo)
   write(8,*) ang
   ! Renormalize all the weights
   wtsum = SUM(w)
   DO n = 1, apo
     w(n) = w(n) * 0.125/wtsum
   END DO
ELSE
  CALL angle
END IF

IF (qdtyp == 2) CLOSE(UNIT=10)
! Call for the input check
CALL check

IF (error_code /= 0) THEN
  error_code_out = error_code
  RETURN
END IF

! Setting orpc value for sweep.
IF (solver == "DGFEM") THEN
    IF (solvertype == "LD" .or. solvertype == "DENSE") THEN
        orpc = dofpc
    ELSE IF (solvertype == "LAGRANGE") THEN
        orpc = ordcb
    END IF
END IF

CALL readxs(xs_file_in)
CALL readsrc(source_input_file_in)

IF (xsbc .eq. 2) THEN
    IF (solver == "AHOTN") THEN
        CALL read_inflow_ahotn(inflow_file)
    ELSE IF (solver == "DGFEM") THEN
        CALL read_inflow_dgfem(inflow_file)
    ELSE IF (solver == "SCTSTEP") THEN
        CALL read_inflow_sct_step(inflow_file)
    END IF
END IF

!CALL echo
WRITE(8,*) "Solver type: ", solvertype
CALL solve
CALL output
fluxout = flux_out

! Time the end of the job
CALL CPU_TIME(tend)

tsolve_out = tsolve
ttosolve_out = ttosolve
tend_out = tend

!Cleanup previously found in old main file
IF( allocated(flux_out)) deallocate(flux_out)
IF( allocated(ang)) deallocate(ang)
IF( allocated(w)) deallocate(w)
IF( allocated(sigt)) deallocate(sigt)
IF( allocated(sigs)) deallocate(sigs)
IF( allocated(s)) deallocate(s)
IF( allocated(frbc)) deallocate(frbc)
IF( allocated(babc)) deallocate(babc)
IF( allocated(lebc)) deallocate(lebc)
IF( allocated(ribc)) deallocate(ribc)
IF( allocated(bobc)) deallocate(bobc)
IF( allocated(tobc)) deallocate(tobc)
IF( allocated(tfrbc)) deallocate(tfrbc)
IF( allocated(tbabc)) deallocate(tbabc)
IF( allocated(tlebc)) deallocate(tlebc)
IF( allocated(tribc)) deallocate(tribc)
IF( allocated(tbobc)) deallocate(tbobc)
IF( allocated(ttobc)) deallocate(ttobc)

IF( allocated(ssum)) deallocate(ssum)
IF( allocated(dx)) deallocate(dx)
IF( allocated(dy)) deallocate(dy)
IF( allocated(dz)) deallocate(dz)
IF( allocated(mat)) deallocate(mat)
IF( allocated(f)) deallocate(f)
IF( allocated(e)) deallocate(e)
IF( allocated(f_ahot_l)) deallocate(f_ahot_l)
IF( allocated(e_ahot_l)) deallocate(e_ahot_l)
IF( allocated(cnvf)) deallocate(cnvf)
IF( allocated(amat)) deallocate(amat)
IF( allocated(bmat)) deallocate(bmat)
IF( allocated(gmat)) deallocate(gmat)
IF( allocated(jmat)) deallocate(jmat)
IF( allocated(gaa)) deallocate(gaa)
IF( allocated(gaxy)) deallocate(gaxy)
IF( allocated(gaxz)) deallocate(gaxz)
IF( allocated(gayz)) deallocate(gayz)
IF( allocated(gxya)) deallocate(gxya)
IF( allocated(gxyxy)) deallocate(gxyxy)
IF( allocated(gxyxz)) deallocate(gxyxz)
IF( allocated(gxyyz)) deallocate(gxyyz)
IF( allocated(gxza)) deallocate(gxza)
IF( allocated(gxzxy)) deallocate(gxzxy)
IF( allocated(gxzxz)) deallocate(gxzxz)
IF( allocated(gxzyz)) deallocate(gxzyz)
IF( allocated(gyza)) deallocate(gyza)
IF( allocated(gyzxy)) deallocate(gyzxy)
IF( allocated(gyzxz)) deallocate(gyzxz)
IF( allocated(gyzyz)) deallocate(gyzyz)
IF( allocated(xmat)) deallocate(xmat)
IF( allocated(ymat)) deallocate(ymat)
IF( allocated(zmat)) deallocate(zmat)
IF( allocated(phisum)) deallocate(phisum)
IF( allocated(refl_left)) deallocate(refl_left)
IF( allocated(refl_right)) deallocate(refl_right)
IF( allocated(refl_front)) deallocate(refl_front)
IF( allocated(refl_back)) deallocate(refl_back)
IF( allocated(refl_top)) deallocate(refl_top)
IF( allocated(refl_bottom)) deallocate(refl_bottom)

RETURN 

END SUBROUTINE main
