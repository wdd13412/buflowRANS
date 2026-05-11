!1.dataStructures.jl:数据结构   2.MAIN SOLVER:主求解函数(subroutine solve) +  应变率张量，湍流粘性，分子粘性
!3. BOUNDARY CONDITIONS:边界条件(wall通过无量纲距离进行了壁面修正)
!4.CONSTITUTIVE RELATIONS:状态变量与通量推导公式   5.JST SCHEME:通量残差   
!6. NEW RANS-SPECIFIC FUNCTIONS + SAME FUNCTIONS(就是之前的mesh.jl):网格 + 更新湍流场，计算wall距离，湍流梯度，粘性通量，湍流应力张量，源项
!7.MODIFIED NUMERICS:梯度,面值插值,面通量到单元格通量 + turbulence source terms to residuals
!8.TIME INTEGRATION:更新解   8.SECTION 9: OUTPUT:输出文件 
!10.Vector Functions:向量计算函数
!11.EXAMPLE INITIALIZATION:初始化

!===============================================================================
! RANS k-omega TURBULENCE MODEL ADDITION


module BuFlowModule	
    implicit none
    
    !===========================================================================
    ! RANS k-omega: Added turbulence model parameters
    !===========================================================================
    real(kind=8), parameter :: SIGMA_K = 0.5d0      ! k diffusion coefficient
    real(kind=8), parameter :: SIGMA_OMEGA = 0.5d0  ! omega diffusion coefficient
    real(kind=8), parameter :: BETA_STAR = 0.09d0   ! Turbulent production coefficient
    real(kind=8), parameter :: BETA = 0.075d0       ! omega destruction coefficient
    real(kind=8), parameter :: ALPHA = 5.0d0/9.0d0  ! Production coefficient
    real(kind=8), parameter :: KAPPA = 0.41d0       ! von Karman constant
    real(kind=8), parameter :: CAPPA = 0.09d0       ! k-omega model constant
    real(kind=8), parameter :: SMALL_NUM = 1.0d-20  ! Small number for divisions
    real(kind=8), parameter :: Y_PLUS_LAMINAR = 11.63d0  ! y+ transition
    
    !===========================================================================
    ! Solver type switch
    !===========================================================================
    integer, parameter :: SOLVER_DENSITY_BASED = 1
    integer, parameter :: SOLVER_SIMPLE = 2
    integer, parameter :: SOLVER_TYPE = SOLVER_SIMPLE
    
    !===========================================================================
    ! Turkel low-Mach preconditioning (for density-based solver)
    !===========================================================================
    logical, parameter :: LOW_MACH_PRECOND = .true.
    real(kind=8), parameter :: M_CUTOFF = 0.01d0
    
    !===========================================================================
    ! SIMPLE solver parameters
    !===========================================================================
    real(kind=8), parameter :: SIMPLE_ALPHA_U = 0.02d0
    real(kind=8), parameter :: SIMPLE_ALPHA_P = 0.01d0
    real(kind=8), parameter :: SIMPLE_ALPHA_K = 0.5d0
    real(kind=8), parameter :: SIMPLE_ALPHA_OMEGA = 0.5d0
    integer, parameter :: SIMPLE_MAX_ITER = 3000
    integer, parameter :: SIMPLE_MAX_INNER_U = 10
    integer, parameter :: SIMPLE_MAX_INNER_P = 1000
    real(kind=8), parameter :: SIMPLE_TOLERANCE = 1.0d-6
    real(kind=8), parameter :: SIMPLE_PCG_TOLERANCE = 1.0d-10
    real(kind=8), parameter :: SIMPLEC_PSEUDO_CFL = 0.05d0
    real(kind=8), parameter :: SIMPLEC_D_DIAG_FLOOR = 0.5d0
    real(kind=8), parameter :: SIMPLEC_RC_DAMPING = 0.3d0
    real(kind=8), parameter :: SIMPLEC_MAX_PRESSURE_STEP_FACTOR = 0.001d0
    real(kind=8), parameter :: SIMPLE_MAX_SPEED = 17.0d0
    real(kind=8), parameter :: SIMPLE_MAX_PRESSURE_RANGE = 200.0d0
    integer, parameter :: SIMPLE_OUTPUT_INTERVAL = 50
    logical, parameter :: DEBUG_MESH_VERBOSE = .false.
    !===========================================================================
    
	type Celll
		integer, allocatable :: faceIndices(:)
		integer, allocatable :: pointIndices(:)
	end type Celll
	
	type SolutionState
		!--- Original inviscid variables ---
		! real(kind=8), allocatable :: cellState(:,:)      ! Original: [ρ, ρu, ρv, ρw, ρE] (5 vars)
		! real(kind=8), allocatable :: cellFluxes(:,:)     ! Original: 15 components
		! real(kind=8), allocatable :: cellPrimitives(:,:) ! Original: [P, T, u, v, w] (5 vars)
		! real(kind=8), allocatable :: fluxResiduals(:,:)  ! Original: 5 vars
		! real(kind=8), allocatable :: faceFluxes(:,:)     ! Original: 15 components
		
		!--- RANS k-omega: Modified to include turbulence variables ---
		real(kind=8), allocatable :: cellState(:,:)      ! RANS: [ρ, ρu, ρv, ρw, ρE, ρk, ρω] (7 vars)
		real(kind=8), allocatable :: cellFluxes(:,:)     ! RANS: 21 components (7 vars × 3 directions)
		real(kind=8), allocatable :: cellPrimitives(:,:) ! RANS: [P, T, u, v, w, k, ω] (7 vars)
		real(kind=8), allocatable :: fluxResiduals(:,:)  ! RANS: 7 vars
		real(kind=8), allocatable :: faceFluxes(:,:)     ! RANS: 21 components
		
		!--- RANS k-omega: Additional turbulence fields ---
		real(kind=8), allocatable :: cellMuT(:)          ! Turbulent viscosity (eddy viscosity)
		real(kind=8), allocatable :: cellMuL(:)          ! Laminar viscosity
		real(kind=8), allocatable :: cellYplus(:)        ! Wall distance in wall units
		real(kind=8), allocatable :: cellWallDist(:)     ! Distance to nearest wall
		real(kind=8), allocatable :: cellStrainRate(:)   ! Strain rate magnitude
		real(kind=8), allocatable :: cellVorticity(:)    ! Vorticity magnitude
		real(kind=8), allocatable :: cellProduction(:)   ! Turbulent production term
	end type SolutionState
	
	type Meshh
		integer(kind=8), allocatable :: cells(:,:)
		real(kind=8), allocatable :: cVols(:)
		real(kind=8), allocatable :: cCenters(:,:)
		real(kind=8), allocatable :: cellSizes(:,:)		
		integer(kind=8), allocatable :: faces(:,:)
		real(kind=8), allocatable :: fAVecs(:,:)
		real(kind=8), allocatable :: fCenters(:,:)
		integer(kind=8), allocatable :: boundaryFaces(:,:)
		character(len=100), allocatable :: boundaryNames(:)
		
		!--- RANS k-omega: Added wall distance field ---
		real(kind=8), allocatable :: wallDistance(:)     ! Distance from cell center to nearest wall
	end type Meshh
	
	type Fluidd
		real(kind=8) :: Cp
		real(kind=8) :: R
		real(kind=8) :: gammaa
		
		!--- RANS k-omega: Added laminar transport properties ---
		real(kind=8) :: mu_laminar    ! Laminar dynamic viscosity
		real(kind=8) :: Pr            ! Prandtl number
		real(kind=8) :: Pr_turb       ! Turbulent Prandtl number
	end type Fluidd
	
	type SolverStatus
		real(kind=8) :: currentTime
		integer(kind=8) :: nTimeSteps
		real(kind=8) :: nextOutputTime
		real(kind=8) :: endTime
	end type SolverStatus
	
	type FaceType
		integer, allocatable :: points(:)
	end type FaceType
	
	type FaceArray
		type(FaceType), allocatable :: faces(:)
	end type FaceArray
	
	type BoundaryCondition
		integer :: type
		real(kind=8), allocatable :: params(:)
	end type BoundaryCondition
	
	integer, parameter :: wallBoundary = 1
	integer, parameter :: emptyBoundary = 2
	integer, parameter :: InletBoundary = 3
	integer, parameter :: OutletBoundary = 4
	
	type CelllArray
		type(Celll), allocatable :: cells(:)
	end type CelllArray
	
	type RestrictResult
		logical :: needTruncate
		real(kind=8) :: actualDt
	end type RestrictResult
	
	type SIMPLEState
		integer :: nCells, nFaces, nBoundaries, nBdryFaces
		real(kind=8), allocatable :: p(:)
		real(kind=8), allocatable :: u(:), v(:), w(:)
		real(kind=8), allocatable :: k(:), omega(:)
		real(kind=8), allocatable :: mu_t(:), mu_l(:)
		real(kind=8), allocatable :: p_prime(:)
		real(kind=8), allocatable :: aP_u(:)
		real(kind=8), allocatable :: rho_field(:)
		real(kind=8), allocatable :: phi_f(:)  ! 面质量通量 (kg/s)
		real(kind=8) :: rho_ref
		real(kind=8) :: p_ref
	end type SIMPLEState
	
	interface dot
		module procedure dot_vec_vec
		module procedure dot_vec_mat
	end interface dot
	
	type MeshData
		real(kind=8), allocatable :: points(:,:)
		type(FaceArray) :: faces
		integer, allocatable :: owner(:), neighbour(:)
		character(len=100), allocatable :: boundaryNames(:)
		integer, allocatable :: boundaryNumFaces(:), boundaryStartFaces(:)
	end type MeshData
	
contains

!===============================================================================
! SECTION 2: MAIN SOLVER - RANS k-omega modifications
!===============================================================================

	
	!--- RANS k-omega: Modified to include k and omega initialization ---
	function initializeUniformSolution3D(mesh, P, T, Ux, Uy, Uz, k_init, omega_init)
		implicit none
		type(Meshh), intent(in) :: mesh
		! real(kind=8), intent(in) :: P, T, Ux, Uy, Uz  ! Original
		real(kind=8), intent(in) :: P, T, Ux, Uy, Uz, k_init, omega_init  ! RANS
		real(kind=8), allocatable :: initializeUniformSolution3D(:,:)
		integer(kind=8) :: nCells, meshInfo(4) 
		integer :: nVars, c, nDims
		real(kind=8), allocatable :: initialValues(:)
		    !--- VALIDATION: Check input parameters ---
		print *, "=== initializeUniformSolution3D: Validating inputs ==="
		print *, "  P     = ", P
		print *, "  T     = ", T
		print *, "  Ux    = ", Ux
		print *, "  Uy    = ", Uy
		print *, "  Uz    = ", Uz
		print *, "  k     = ", k_init
		print *, "  omega = ", omega_init
    
		if (P <= 0.0d0 .or. P /= P) then
			print *, "ERROR: Invalid pressure P = ", P
			stop
		end if
		if (T <= 0.0d0 .or. T /= T) then
			print *, "ERROR: Invalid temperature T = ", T
			stop
		end if
		if (k_init < 0.0d0 .or. k_init /= k_init) then
			print *, "ERROR: Invalid k_init = ", k_init
			stop
		end if
		if (omega_init <= 0.0d0 .or. omega_init /= omega_init) then
			print *, "ERROR: Invalid omega_init = ", omega_init
			stop
		end if
    !--- END VALIDATION ---
		
		nDims = 3
		meshInfo = unstructuredMeshInfo(mesh)
		nCells = meshInfo(1)
		! nVars = 2 + nDims  ! Original: P, T, Ux, Uy, Uz (5 vars)
		nVars = 2 + nDims + 2  ! RANS: P, T, Ux, Uy, Uz, k, omega (7 vars)
		print *, "  nCells = ", nCells
		print *, "  nVars  = ", nVars
		
		allocate(initializeUniformSolution3D(nCells, nVars))
		allocate(initialValues(nVars))
		
		! initialValues = [P, T, Ux, Uy, Uz]  ! Original
		initialValues = [P, T, Ux, Uy, Uz, k_init, omega_init]  ! RANS
		initialValues = initialValues(1:nVars)
		print *, "  initialValues = ", initialValues
		
		do c = 1, nCells
			initializeUniformSolution3D(c, :) = initialValues
		end do
		print *, "=== Initialization complete ==="
		print *, "  Output array size: ", size(initializeUniformSolution3D, 1), "x", size(initializeUniformSolution3D, 2)
		print *, "  First row: ", initializeUniformSolution3D(1, :)
		
	end function initializeUniformSolution3D

	!===========================================================================
	! RANS k-omega: New function to compute laminar viscosity (Sutherland's law)
	!===========================================================================
	function computeLaminarViscosity(T) result(mu)
		implicit none
		real(kind=8), intent(in) :: T
		real(kind=8) :: mu
		real(kind=8), parameter :: mu_ref = 1.7894d-5  ! Reference viscosity at T_ref
		real(kind=8), parameter :: T_ref = 288.15d0    ! Reference temperature (K)
		real(kind=8), parameter :: S = 110.4d0         ! Sutherland constant (K)
		
		! Sutherland's law
		mu = mu_ref * (T / T_ref)**1.5d0 * (T_ref + S) / (T + S)
	end function computeLaminarViscosity

	!===========================================================================
	! RANS k-omega: New function to compute turbulent viscosity
	!===========================================================================
	function computeTurbulentViscosity(rho, k, omega) result(mu_t)
		implicit none
		real(kind=8), intent(in) :: rho, k, omega
		real(kind=8) :: mu_t
		
		! Prevent division by zero
		if (omega < SMALL_NUM) then
			mu_t = 0.0d0
		else
			! mu_t = rho * k / omega
			mu_t = rho * k / max(omega, SMALL_NUM)
		end if
		
		! Limit turbulent viscosity (optional stability measure)
		mu_t = min(mu_t, 1.0d5)
	end function computeTurbulentViscosity

	!===========================================================================
	! RANS k-omega: New function to compute strain rate tensor magnitude
	!===========================================================================
	function computeStrainRateMagnitude(gradU) result(S_mag)
		implicit none
		real(kind=8), intent(in) :: gradU(3,3)  ! Velocity gradient tensor
		real(kind=8) :: S_mag
		real(kind=8) :: S(3,3)  ! Strain rate tensor
		integer :: i, j
		
		! S_ij = 0.5 * (du_i/dx_j + du_j/dx_i)
		do i = 1, 3
			do j = 1, 3
				S(i,j) = 0.5d0 * (gradU(i,j) + gradU(j,i))
			end do
		end do
		
		! S_mag = sqrt(2 * S_ij * S_ij)
		S_mag = 0.0d0
		do i = 1, 3
			do j = 1, 3
				S_mag = S_mag + S(i,j) * S(i,j)
			end do
		end do
		S_mag = sqrt(2.0d0 * S_mag)
	end function computeStrainRateMagnitude

	!===========================================================================
	! RANS k-omega: Modified CFL computation (accounts for viscous terms)
	!===========================================================================
	function CFL(mesh, sln, fluid, dt) result(CFLL)
		implicit none
		real(kind=8), allocatable :: CFLL(:), cellRhoT(:,:)
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(in) :: sln
		type(Fluidd), intent(in) :: fluid
		real(kind=8), intent(in), optional :: dt
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4)
		integer(kind=8) :: f, ownerCell, neighbourCell, stat
		real(kind=8), allocatable :: faceRhoT(:,:)
		real(kind=8) :: faceRho, faceT, flux, a, dt_val
		real(kind=8), allocatable :: faceVel(:), positionn(:)
		! RANS k-omega: Additional variables for viscous CFL
		real(kind=8) :: mu_eff, nu_eff, visc_cfl
		
		dt_val = 1.0d0
		if(present(dt)) dt_val = dt
		
		meshInfo = unstructuredMeshInfo(mesh)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		allocate(CFLL(nCells))
		CFLL = 0.0d0
		
		allocate(faceRhoT(nFaces, 2))		
		allocate(cellRhoT(nCells, 2))
		cellRhoT(:, 1) = sln%cellState(:, 1)
		cellRhoT(:, 2) = sln%cellPrimitives(:, 2)
		call linInterp_3D_RANS(mesh, cellRhoT, faceRhoT)
		
		allocate(faceVel(3), positionn(3), stat=stat)
		if (stat /= 0) stop "Error: 分配 faceVel/positionn 失败"
		
		do f = 1, nFaces
			ownerCell = mesh%faces(f,1)
			neighbourCell = mesh%faces(f,2)
			
			faceRho = faceRhoT(f, 1)
			faceT = faceRhoT(f, 2)
			
			if(neighbourCell == -1) then
				faceRho = sln%cellState(ownerCell, 1)
				faceT = sln%cellPrimitives(ownerCell, 2)
			end if
			
			faceVel = sln%faceFluxes(f, 1:3) / faceRho
			flux = abs(dot_product(faceVel, mesh%fAVecs(f,:))) * dt_val
			
			if(faceT <= 0.0d0) then
				positionn = mesh%fCenters(f,:)
			end if
			
			a = sqrt(fluid%gammaa * fluid%R * faceT)
			if (LOW_MACH_PRECOND) then
				flux = flux + mag(mesh%fAVecs(f,:)) * &
					preconditioned_spectral_radius(mag(faceVel), a, fluid%gammaa) * dt_val - &
					abs(dot_product(faceVel, mesh%fAVecs(f,:))) * dt_val
			else
				flux = flux + mag(mesh%fAVecs(f,:)) * a * dt_val
			end if
			
			!--- RANS k-omega: Add viscous contribution to CFL ---
			! mu_eff = mu_laminar + mu_turbulent
			mu_eff = sln%cellMuL(ownerCell) + sln%cellMuT(ownerCell)
			nu_eff = mu_eff / faceRho
			visc_cfl = nu_eff * mag(mesh%fAVecs(f,:))**2 / mesh%cVols(ownerCell) * dt_val
			flux = flux + visc_cfl
			!---
			
			CFLL(ownerCell) = CFLL(ownerCell) + flux
			if(neighbourCell > -1) then
				CFLL(neighbourCell) = CFLL(neighbourCell) + flux
			end if
		end do
		
		deallocate(faceVel, positionn)
		CFLL = CFLL / (2.0d0 * mesh%cVols)
		deallocate(faceRhoT, cellRhoT)
	end function CFL

	!===========================================================================
	! RANS k-omega: Modified populateSolution 
	!===========================================================================
	function populateSolution(cellPrimitives, nCells, nFaces, fluid, nDims) result(sln)
		implicit none
		real(kind=8), intent(in) :: cellPrimitives(:,:)
		integer(kind=8), intent(in) :: nCells, nFaces
		type(Fluidd), intent(in) :: fluid
		integer, intent(in), optional :: nDims
		type(SolutionState) :: sln
		integer :: nConservedVars, nFluxes, nDims_val
		integer :: i 
		
		nDims_val = 3
		if (present(nDims)) nDims_val = nDims
		
		! nConservedVars = 2 + nDims_val  ! Original: 5 vars
		nConservedVars = 2 + nDims_val + 2  ! RANS: 7 vars (added k and omega)
		nFluxes = nConservedVars * nDims_val  ! 21 components
		
		! Allocate arrays
		allocate(sln%cellPrimitives(nCells, 7))     ! RANS: 7 primitives
		allocate(sln%cellState(nCells, nConservedVars))
		allocate(sln%cellFluxes(nCells, nFluxes))
		allocate(sln%fluxResiduals(nCells, nConservedVars))
		allocate(sln%faceFluxes(nFaces, nFluxes))
		
		!--- RANS k-omega: Allocate turbulence fields ---
		allocate(sln%cellMuT(nCells))
		allocate(sln%cellMuL(nCells))
		allocate(sln%cellYplus(nCells))
		allocate(sln%cellWallDist(nCells))
		allocate(sln%cellStrainRate(nCells))
		allocate(sln%cellVorticity(nCells))
		allocate(sln%cellProduction(nCells))
		!--- CRITICAL CHECK: Verify input array ---
		print *, "=== populateSolution: Checking input array ==="
		print *, "  cellPrimitives size: ", size(cellPrimitives, 1), "x", size(cellPrimitives, 2)

		if (size(cellPrimitives, 2) /= 7) then
			print *, "ERROR: cellPrimitives must have 7 columns!"
			stop
		end if

		! Check for NaN in input
		if (any(cellPrimitives /= cellPrimitives)) then
			print *, "ERROR: NaN detected in input cellPrimitives!"
			print *, "  Searching for NaN locations..."
			do i = 1, min(10, size(cellPrimitives, 1))
				if (any(cellPrimitives(i,:) /= cellPrimitives(i,:))) then
					print *, "  NaN found in row ", i, ": ", cellPrimitives(i, :)
				end if
			end do
			stop
		end if

		print *, "  Input validation passed"
		
		
		! Initialize
		sln%cellPrimitives = cellPrimitives
		
		!--- Validate input primitives ---
		print *, "DEBUG populateSolution: Validating input primitives"
		print *, "  Pressure range:", minval(cellPrimitives(:,1)), maxval(cellPrimitives(:,1))
		print *, "  Temperature range:", minval(cellPrimitives(:,2)), maxval(cellPrimitives(:,2))
		print *, "  k range:", minval(cellPrimitives(:,6)), maxval(cellPrimitives(:,6))
		print *, "  omega range:", minval(cellPrimitives(:,7)), maxval(cellPrimitives(:,7))

		sln%cellState      = encodePrimitives3D_RANS(cellPrimitives, fluid)  ! Modified
		!--- Validate encoded state ---
		print *, "DEBUG populateSolution: Validating encoded state"
		print *, "  Density range:", minval(sln%cellState(:,1)), maxval(sln%cellState(:,1))
		
		sln%cellFluxes     = 0.0d0
		sln%fluxResiduals  = 0.0d0
		sln%faceFluxes     = 0.0d0
		
		!--- RANS k-omega: Initialize turbulence fields ---
		sln%cellMuT = 0.0d0
		sln%cellMuL = 0.0d0
		sln%cellYplus = 0.0d0
		sln%cellWallDist = 0.0d0
		sln%cellStrainRate = 0.0d0
		sln%cellVorticity = 0.0d0
		sln%cellProduction = 0.0d0
		print *, "DEBUG populateSolution: Calling decodeSolution_3D_RANS"
		
		call decodeSolution_3D_RANS(sln, fluid)  ! Modified
		
		print *, "DEBUG populateSolution: After decode"
		print *, "  Pressure range:", minval(sln%cellPrimitives(:,1)), maxval(sln%cellPrimitives(:,1))
		print *, "  Velocity U range:", minval(sln%cellPrimitives(:,3)), maxval(sln%cellPrimitives(:,3))
	end function populateSolution

	function restrictTimeStep(status_val, desiredDt) result(res)
		implicit none
		type(SolverStatus), intent(in) :: status_val
		real(kind=8), intent(in) :: desiredDt
		type(RestrictResult) :: res
		real(kind=8) :: maxStep
		
		maxStep = min(status_val%endTime - status_val%currentTime, &
		              status_val%nextOutputTime - status_val%currentTime)
		
		if (desiredDt > maxStep) then
			res%needTruncate = .true.
			res%actualDt = maxStep
		else
			res%needTruncate = .false.
			res%actualDt = desiredDt
		end if
	end function restrictTimeStep

	function adjustTimeStep_LTS(targetCFL, dt, status_val) result(output)
		implicit none
		real(kind=8), intent(in) :: targetCFL
		real(kind=8), intent(inout) :: dt(:)
		type(SolverStatus), intent(in) :: status_val
		real(kind=8) :: output(3)
		type(RestrictResult) :: res
		real(kind=8) :: CFLL
		
		if (status_val%nTimeSteps < 100) then
			CFLL = (status_val%nTimeSteps + 1) * targetCFL / 100.0d0
		else
			CFLL = targetCFL
		end if
		
		res = restrictTimeStep(status_val, CFLL)
		CFLL = res%actualDt
		dt(1) = CFLL
		
		!output = [merge(1.0d0, 0.0d0, res%needTruncate), &
		!	  dt(1),&
		!	  CFLL]
		output(1) = merge(1.0d0, 0.0d0, res%needTruncate)
		output(2) = dt(1)
		output(3) = CFLL
	end function adjustTimeStep_LTS

	subroutine advanceStatus(status_val, dt, CFLL, timeIntegrationFn, silent)
		implicit none
		type(SolverStatus), intent(inout) :: status_val
		real(kind=8), intent(in) :: dt(:), CFLL
		character(len=*), intent(in) :: timeIntegrationFn
		logical, intent(in) :: silent
		
		status_val%nTimeSteps = status_val%nTimeSteps + 1
		status_val%currentTime = status_val%currentTime + CFLL
		
		if (.not. silent) then
			write(*,'(A,I5,A,F9.4,A,F9.4)') 'Timestep: ', status_val%nTimeSteps, &
			    ', simTime: ', status_val%currentTime, ', Max CFLL: ', CFLL
		end if
	end subroutine advanceStatus

	function fluiddd() result(fluid)
		type(Fluidd) :: fluid
		fluid%Cp = 1005.0d0
		fluid%R = 287.05d0
		fluid%gammaa = 1.4d0
		
		!--- RANS k-omega: Initialize transport properties ---
		fluid%mu_laminar = 1.7894d-5  ! Air at 15°C (kg/(m·s))
		fluid%Pr = 0.72d0              ! Prandtl number for air
		fluid%Pr_turb = 0.9d0          ! Turbulent Prandtl number
	end function fluiddd

	!===========================================================================
	! RANS k-omega: Modified main solve routine
	!===========================================================================
	subroutine solve(mesh, meshPath, cellPrimitives, boundaryConditions)
		implicit none
		type(Meshh), intent(inout) :: mesh
		character(len=*), intent(in) :: meshPath
		real(kind=8), intent(inout) :: cellPrimitives(:,:)
		type(BoundaryCondition), intent(in) :: boundaryConditions(:)
		
		integer :: dbg_nshow, dbg_i
		real(kind=8) :: dbg_min, dbg_max		
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4), tsResult(3)
		type(SolutionState), allocatable :: sln
		real(kind=8), allocatable :: dt(:)
		type(SolverStatus) :: status_val
		real(kind=8), allocatable :: CFLvec(:)
		real(kind=8) :: CFLL, LTSResult(3)
		logical :: writeOutputThisIteration, silent, restart, createRestartFile, createVTKOutput
		character(len=32) :: timeIntegrationFn, fluxFunction
		real(kind=8) :: initDt, endTime, outputInterval, targetCFL
		type(Fluidd) :: fluid
		character(len=256) :: restartFile
		integer :: i, j
		
		! Hardcoded parameters
		timeIntegrationFn = 'LTSEuler'
		fluxFunction = 'unstructured_JSTFlux_RANS'  ! Modified for RANS
		initDt = 0.0000001
		endTime = 350
		outputInterval = 25
		targetCFL = 0.5d0
		fluid = fluiddd()
		silent = .false.
		restart = .false.
		createRestartFile = .true.
		createVTKOutput = .true.
		restartFile = 'FvCFDRestart.txt'
		
		if(.not. silent) print*, 'Initializing RANS k-omega Simulation'
		
		meshInfo = unstructuredMeshInfo(mesh)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		
		!--- RANS k-omega: Compute wall distances ---
		if(.not. silent) print*, 'Computing wall distances...'
		call computeWallDistances(mesh)
		if(.not. silent) print*, 'Wall distances computed successfully'
		
		sln = populateSolution(cellPrimitives, nCells, nFaces, fluid, 3)
		
		if(timeIntegrationFn == 'LTSEuler') then
			allocate(dt(nCells))
			dt = 0.0d0
		else
			dt = initDt
		end if
		
		status_val = SolverStatus(0, 0, outputInterval, endTime)
		allocate(CFLvec(nCells))
		CFLvec = 0.0d0
		
		if(.not. silent) print*, 'Starting RANS iterations'
		
		do while(status_val%currentTime < status_val%endTime)
			LTSResult = adjustTimeStep_LTS(targetCFL, dt, status_val)
			writeOutputThisIteration = (ltsResult(1) == 1.0d0)
			dt(1) = ltsResult(2)
			CFLL = ltsResult(3)
			
			!--- RANS k-omega: Update turbulence fields before flux computation ---
			if (status_val%nTimeSteps == 0) then
				print *, "DEBUG: First call to updateTurbulenceFields"
				print *, "  sln%cellState(:,1) range:", minval(sln%cellState(:,1)), maxval(sln%cellState(:,1))
				print *, "  sln%cellPrimitives(:,1) range:", minval(sln%cellPrimitives(:,1)), maxval(sln%cellPrimitives(:,1))
			end if
			
			!--- RANS k-omega: Update turbulence fields before flux computation ---
			call updateTurbulenceFields(mesh, sln, fluid)
			
			if (status_val%nTimeSteps == 0) then
				print *, "DEBUG: After updateTurbulenceFields"
				print *, "  cellMuL range:", minval(sln%cellMuL), maxval(sln%cellMuL)
				print *, "  cellMuT range:", minval(sln%cellMuT), maxval(sln%cellMuT)
			end if
			
			!--- Modified time integration for RANS ---
			call LTSEuler_RANS(mesh, sln, boundaryConditions, fluid, dt)
			
			if (status_val%nTimeSteps == 0) then
				print *, "DEBUG: After first LTSEuler_RANS"
				print *, "  sln%cellState(:,1) range:", minval(sln%cellState(:,1)), maxval(sln%cellState(:,1))
				print *, "  sln%cellPrimitives(:,1) range:", minval(sln%cellPrimitives(:,1)), maxval(sln%cellPrimitives(:,1))
			end if
			
			call advanceStatus(status_val, dt, CFLL, timeIntegrationFn, silent)
			
			if(writeOutputThisIteration) then
				status_val%nextOutputTime = status_val%nextOutputTime + outputInterval
			end if
		end do
		
		! Output
		print *, new_line('a')//"=== RANS Solution Complete ==="
		print *, "Total cells: ", nCells
		print *, "Solution dimensions: ", size(sln%cellPrimitives,1), "×", size(sln%cellPrimitives,2)
		print *, "Format: Cell | P | T | Ux | Uy | Uz | k | omega"
		
		do i = 1, min(15, nCells)
			write(*, '(I5, 7F12.6)') i, &
			    sln%cellPrimitives(i,1), sln%cellPrimitives(i,2), &
			    sln%cellPrimitives(i,3), sln%cellPrimitives(i,4), sln%cellPrimitives(i,5), &
			    sln%cellPrimitives(i,6), sln%cellPrimitives(i,7)
		end do

if ( allocated(sln%cellPrimitives) ) then
  dbg_nshow = min(5, size(sln%cellPrimitives,1))
  write(*,*) 'DEBUG_A: sln%cellPrimitives size =', size(sln%cellPrimitives)
  dbg_min = minval(sln%cellPrimitives(:,1))
  dbg_max = maxval(sln%cellPrimitives(:,1))
  write(*,'(A,2F12.6)') 'DEBUG_A: Pressure min/max = ', dbg_min, dbg_max
  dbg_min = minval(sln%cellPrimitives(:,3))
  dbg_max = maxval(sln%cellPrimitives(:,3))
  write(*,'(A,2F12.6)') 'DEBUG_A: Velocity_x min/max = ', dbg_min, dbg_max
  write(*,*) 'DEBUG_A: sample sln%cellPrimitives rows 1..', dbg_nshow
  do dbg_i = 1, dbg_nshow
     write(*,'(I6,7F12.6)') dbg_i, sln%cellPrimitives(dbg_i,1), sln%cellPrimitives(dbg_i,2), &
          sln%cellPrimitives(dbg_i,3), sln%cellPrimitives(dbg_i,4), sln%cellPrimitives(dbg_i,5), &
          sln%cellPrimitives(dbg_i,6), sln%cellPrimitives(dbg_i,7)
  end do
else
  write(*,*) 'DEBUG_A: sln%cellPrimitives is NOT allocated'
end if

if ( allocated(sln%cellState) ) then
  dbg_nshow = min(5, size(sln%cellState,1))
  write(*,*) 'DEBUG_A: sln%cellState size =', size(sln%cellState)
  dbg_min = minval(sln%cellState(:,1))
  dbg_max = maxval(sln%cellState(:,1))
  write(*,'(A,2F12.6)') 'DEBUG_A: cellState(rho) min/max = ', dbg_min, dbg_max
  write(*,*) 'DEBUG_A: sample sln%cellState rows 1..', dbg_nshow
  do dbg_i = 1, dbg_nshow
     write(*,'(I6,7F12.6)') dbg_i, sln%cellState(dbg_i,1), sln%cellState(dbg_i,2), sln%cellState(dbg_i,3), &
          sln%cellState(dbg_i,4), sln%cellState(dbg_i,5), sln%cellState(dbg_i,6), sln%cellState(dbg_i,7)
  end do
else
  write(*,*) 'DEBUG_A: sln%cellState is NOT allocated'
end if
! --- end DEBUG BLOCK A ---

		
		call writeOutput_RANS(sln, restartFile, meshPath, createRestartFile, createVTKOutput)

		
		deallocate(dt, CFLvec)
	end subroutine solve

!===============================================================================
! SECTION 3: BOUNDARY CONDITIONS - RANS k-omega modifications
!===============================================================================

	!===========================================================================
	! RANS k-omega: Modified inlet boundary (includes k and omega)
	!===========================================================================
	subroutine updateInletBoundary_RANS(mesh, sln, boundaryNumber, inletConditions, fluid)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(inout) :: sln
		integer, intent(in) :: boundaryNumber
		! real(kind=8), intent(in) :: inletConditions(4)  ! Original: [Pt, Tt, nx, ny]
		real(kind=8), intent(in) :: inletConditions(6)  ! RANS: [Pt, Tt, nx, ny, k_inlet, omega_inlet]
		type(Fluidd), intent(in) :: fluid
		
		integer :: face, ownerCell, stat
		integer(kind=8), allocatable :: currentBoundary(:)
		real(kind=8) :: Pt, Tt, nx, ny, nz, gammaa, R
		real(kind=8) :: adjustedVelocity, machNum, e
		! real(kind=8) :: primitives(5), state(5), velUnitVector(3), boundaryFluxes(15)  ! Original
		real(kind=8) :: primitives(7), state(7), velUnitVector(3), boundaryFluxes(21)  ! RANS
		real(kind=8) :: k_inlet, omega_inlet
		integer(kind=8) :: nCells, nFaces, nFacePerBdry, face_idx
		
		nFaces = size(mesh%faces, 1)
		nCells = size(sln%cellState, 1)
		gammaa = fluid%gammaa
		R = fluid%R
		
		Pt = inletConditions(1)
		Tt = inletConditions(2)
		nx = inletConditions(3)
		ny = inletConditions(4)
		k_inlet = inletConditions(5)      ! RANS: Inlet turbulent kinetic energy
		omega_inlet = inletConditions(6)  ! RANS: Inlet specific dissipation rate
		nz = 0.0d0
		
		velUnitVector = [nx, ny, nz]
		
		nFacePerBdry = size(mesh%boundaryFaces, 2)
		allocate(currentBoundary(nFacePerBdry), stat=stat)
		currentBoundary = mesh%boundaryFaces(boundaryNumber, :)
		
		do face = 1, nFacePerBdry
			face_idx = currentBoundary(face)
			if (face_idx == 0) exit
			ownerCell = mesh%faces(face_idx, 1)
			
			! Compute static conditions
			adjustedVelocity = dot_product(velUnitVector, sln%cellPrimitives(ownerCell, 3:5))
			primitives(2) = Tt - (gammaa - 1.0d0) / (2.0d0 * gammaa * R) * adjustedVelocity**2
			machNum = abs(adjustedVelocity) / sqrt(gammaa * R * primitives(2))
			primitives(1) = Pt * (1.0d0 + (gammaa - 1.0d0) / 2.0d0 * machNum**2)**(-gammaa / (gammaa - 1.0d0))
			primitives(3:5) = adjustedVelocity * velUnitVector
			
			!--- RANS k-omega: Set inlet turbulence values ---
			primitives(6) = k_inlet
			primitives(7) = omega_inlet
			
			! Compute conservative variables
			state(1) = idealGasRho(primitives(2), primitives(1), R)
			state(2:4) = primitives(3:5) * state(1)
			e = calPerfectEnergy(primitives(2), fluid)
			state(5) = state(1) * (e + 0.5d0 * sum(primitives(3:5)**2))
			
			!--- RANS k-omega: Conservative turbulence variables ---
			state(6) = state(1) * primitives(6)  ! rho*k
			state(7) = state(1) * primitives(7)  ! rho*omega
			
			boundaryFluxes = 0.0d0
			! call calculateFluxes3D(boundaryFluxes, primitives, state)  ! Original
			call calculateFluxes3D_RANS(boundaryFluxes, primitives, state, fluid, sln, ownerCell)  ! RANS
			
			sln%faceFluxes(face_idx, :) = boundaryFluxes
		end do
		
		deallocate(currentBoundary)
	end subroutine updateInletBoundary_RANS

	!===========================================================================
	! RANS k-omega: Modified outlet boundary
	!===========================================================================
	subroutine updateOutletBoundary_RANS(mesh, sln, boundaryNumber, outletPressure, fluid)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(inout) :: sln
		integer, intent(in) :: boundaryNumber
		real(kind=8), intent(in) :: outletPressure(:)
		type(Fluidd), intent(in) :: fluid
		
		! integer :: nFluxes, face, ownerCell, flux, d, stat  ! Original
		integer :: nFluxes, face, ownerCell, flux, d, stat, v  ! RANS
		integer(kind=8), allocatable :: currentBoundary(:)
		real(kind=8) :: outletPressuree, origP
		integer(kind=8) :: nCells, nFaces, nFacePerBdry, face_idx
		
		nFaces = size(mesh%faces, 1)
		nCells = size(sln%cellState, 1)
		nFluxes = size(sln%cellFluxes, 2)
		nFacePerBdry = size(mesh%boundaryFaces, 2)
		
		allocate(currentBoundary(nFacePerBdry), stat=stat)
		if (stat /= 0) stop "Error: Allocation currentBoundary failed"
		currentBoundary = mesh%boundaryFaces(boundaryNumber, :)
		outletPressuree = outletPressure(1)
		
		do face = 1, nFacePerBdry
			face_idx = currentBoundary(face)
			if (face_idx == 0) exit
			
			ownerCell = mesh%faces(face_idx, 1)
			
			
			do flux = 1, nFluxes
				sln%faceFluxes(face_idx, flux) = sln%cellFluxes(ownerCell, flux)
			end do
			
			! Correct pressure in momentum fluxes
			origP = sln%cellPrimitives(ownerCell, 1)
			
			
			do d = 1, 3
				flux = 4 * d  ! Indices: 4, 8, 12
				sln%faceFluxes(face_idx, flux) = sln%faceFluxes(face_idx, flux) + &
				                                 (outletPressuree - origP)
			end do
			
			
			do d = 1, 3
				flux = 12 + d  ! Indices: 13, 14, 15
				sln%faceFluxes(face_idx, flux) = sln%faceFluxes(face_idx, flux) + &
				                                 sln%cellPrimitives(ownerCell, 2 + d) * (outletPressuree - origP)
			end do
			
			
		end do
		
		deallocate(currentBoundary)
	end subroutine updateOutletBoundary_RANS

	!===========================================================================
	! RANS k-omega: Modified wall boundary (with wall functions)
	!===========================================================================
	subroutine updateWallBoundary_RANS(mesh, sln, boundaryNumber, dummy1, dummy2)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(inout) :: sln
		integer, intent(in) :: boundaryNumber
		real(kind=8), intent(in), optional :: dummy1(:)
		type(Fluidd), intent(in) :: dummy2
		
		integer :: face, stat, i
		integer(kind=8), allocatable :: currentBoundary(:)
		real(kind=8) :: faceP
		integer(kind=8) :: nCells, nFaces, nFacePerBdry, face_idx, ownerCell
		integer :: realFaceCount
		
		!--- RANS k-omega: Wall boundary variables ---
		real(kind=8) :: y_plus, u_tau, tau_w, k_wall, omega_wall
		real(kind=8) :: rho_wall, mu_wall, nu_wall, u_mag
		real(kind=8) :: wall_dist, C_mu
		
		C_mu = BETA_STAR  ! = 0.09
		
		nFaces = size(mesh%faces, 1)
		nCells = size(sln%cellState, 1)
		nFacePerBdry = size(mesh%boundaryFaces, 2)
		
		allocate(currentBoundary(nFacePerBdry), stat=stat)
		if (stat /= 0) stop "Error: Allocation currentBoundary failed"
		currentBoundary = mesh%boundaryFaces(boundaryNumber, :)
		
		! Count real faces
		realFaceCount = 0
		do i = 1, nFacePerBdry
			if (currentBoundary(i) == 0) exit
			realFaceCount = realFaceCount + 1
		end do
		
		do face = 1, realFaceCount
			face_idx = currentBoundary(face)
			
			if (face_idx < 1 .or. face_idx > nFaces) cycle
			
			ownerCell = max(mesh%faces(face_idx, 1), mesh%faces(face_idx, 2))
			
			!--- Original wall boundary: No-slip, impermeability ---
			faceP = sln%cellPrimitives(ownerCell, 1)
			
			! sln%faceFluxes(face_idx, 4) = faceP    ! x momentum  ! Original
			! sln%faceFluxes(face_idx, 8) = faceP    ! y momentum  ! Original
			! sln%faceFluxes(face_idx, 12) = faceP   ! z momentum  ! Original
			! sln%faceFluxes(face_idx, 1:3) = 0.0d0  ! mass flux   ! Original
			! sln%faceFluxes(face_idx, 13:15) = 0.0d0 ! energy flux ! Original
			
			!--- RANS k-omega: Set momentum fluxes (include pressure and wall shear) ---
			sln%faceFluxes(face_idx, 4) = faceP   ! x momentum
			sln%faceFluxes(face_idx, 8) = faceP   ! y momentum
			sln%faceFluxes(face_idx, 12) = faceP  ! z momentum
			sln%faceFluxes(face_idx, 1:3) = 0.0d0   ! mass flux = 0
			sln%faceFluxes(face_idx, 13:15) = 0.0d0 ! energy flux
			
			!--- RANS k-omega: Wall boundary conditions for turbulence ---
			! Get wall properties
			rho_wall = sln%cellState(ownerCell, 1)
			mu_wall = sln%cellMuL(ownerCell)
			nu_wall = mu_wall / rho_wall
			wall_dist = sln%cellWallDist(ownerCell)
			
			! Velocity magnitude in adjacent cell
			u_mag = sqrt(sln%cellPrimitives(ownerCell, 3)**2 + &
			             sln%cellPrimitives(ownerCell, 4)**2 + &
			             sln%cellPrimitives(ownerCell, 5)**2)
			
			! Compute y+
			u_tau = sqrt(C_mu) * sqrt(sln%cellPrimitives(ownerCell, 6))  ! u_tau = sqrt(C_mu * k)
			if (u_tau < SMALL_NUM) u_tau = SMALL_NUM
			y_plus = rho_wall * u_tau * wall_dist / mu_wall
			sln%cellYplus(ownerCell) = y_plus
			
			!--- Wall function for k ---
			! k at wall = 0 (Dirichlet), but we use wall function for first cell
			if (y_plus > Y_PLUS_LAMINAR) then
				k_wall = u_tau**2 / sqrt(C_mu)  ! Log-layer assumption
			else
				k_wall = 0.0d0  ! Laminar sublayer
			end if
			
			!--- Wall function for omega ---
			! omega_wall = u_tau^2 / (nu * sqrt(C_mu))  for y+ > 2.5
			! omega_wall = 6*nu / (beta * y^2)          for y+ <= 2.5
			if (y_plus > 2.5d0) then
				omega_wall = u_tau / (sqrt(C_mu) * KAPPA * wall_dist)
			else
				omega_wall = 6.0d0 * nu_wall / (BETA * wall_dist**2)
			end if
			
			! Prevent unrealistic values
			omega_wall = max(omega_wall, 1.0d-10)
			omega_wall = min(omega_wall, 1.0d10)
			
			!--- Set turbulence fluxes at wall ---
			! k flux (indices 16, 17, 18): zero flux
			sln%faceFluxes(face_idx, 16:18) = 0.0d0
			
		
			sln%faceFluxes(face_idx, 19:21) = 0.0d0
			
			! Store wall omega in adjacent cell (for production/destruction terms)
			sln%cellPrimitives(ownerCell, 7) = omega_wall
			sln%cellState(ownerCell, 7) = rho_wall * omega_wall
		end do
		
		deallocate(currentBoundary)
	end subroutine updateWallBoundary_RANS

	!===========================================================================
	! RANS k-omega: Empty boundary (unchanged)
	!===========================================================================
	subroutine updateEmptyBoundary(mesh, sln, boundaryNumber, dummy1, dummy2)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(inout) :: sln
		integer, intent(in) :: boundaryNumber
		real(kind=8), intent(in), optional :: dummy1(:)
		type(Fluidd), intent(in) :: dummy2
		! Empty boundary does nothing
	end subroutine updateEmptyBoundary

	!===========================================================================
	! RANS k-omega: Modified unified boundary condition dispatcher
	!===========================================================================
	subroutine updateBoundaryConditions_RANS(mesh, sln, boundaryConditions, nBoundaries, fluid)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(inout) :: sln
		type(BoundaryCondition), intent(in) :: boundaryConditions(:)
		integer(kind=8), intent(in) :: nBoundaries
		type(Fluidd), intent(in) :: fluid
		
		integer :: boundaryNumber
		integer(kind=8) :: nCells, nFaces, nVars, dim1, dim2
		
		nCells = size(sln%cellState, 1)
		nFaces = size(mesh%faces, 1)
		nVars = size(sln%cellState, 2)
		
		do boundaryNumber = 1, nBoundaries
			select case (boundaryConditions(boundaryNumber)%type)
			case (wallBoundary)
				call updateWallBoundary_RANS(mesh, sln, boundaryNumber, &
				                             boundaryConditions(boundaryNumber)%params, fluid)
			case (emptyBoundary)
				call updateEmptyBoundary(mesh, sln, boundaryNumber, &
				                         boundaryConditions(boundaryNumber)%params, fluid)
			case (InletBoundary)
				call updateInletBoundary_RANS(mesh, sln, boundaryNumber, &
				                              boundaryConditions(boundaryNumber)%params, fluid)
			case (OutletBoundary)
				call updateOutletBoundary_RANS(mesh, sln, boundaryNumber, &
				                               boundaryConditions(boundaryNumber)%params, fluid)
			case default
				stop "Error: Unknown boundary type"
			end select
		end do
	end subroutine updateBoundaryConditions_RANS

	!===============================================================================
	! SECTION 4: CONSTITUTIVE RELATIONS - RANS k-omega modifications
	!===============================================================================

		
	function idealGasRho(T, P, R)
		implicit none
		real(kind=8), intent(in) :: T, P
		real(kind=8), intent(in), optional :: R
		real(kind=8) :: idealGasRho
		real(kind=8) :: R_val
    
		R_val = 287.05d0
		if(present(R)) R_val = R
    
    !--- SAFETY CHECKS ---
		if (T <= 0.0d0 .or. T /= T) then
			print *, "ERROR in idealGasRho: Invalid temperature T = ", T
			idealGasRho = 0.0d0
			return
		end if
		if (P <= 0.0d0 .or. P /= P) then
			print *, "ERROR in idealGasRho: Invalid pressure P = ", P
			idealGasRho = 0.0d0
			return
		end if
		if (R_val <= 0.0d0) then
			print *, "ERROR in idealGasRho: Invalid gas constant R = ", R_val
			idealGasRho = 0.0d0
			return
		end if
    !
    
    idealGasRho = P/(R_val*T)
    

		if (idealGasRho /= idealGasRho .or. idealGasRho <= 0.0d0) then
			print *, "ERROR in idealGasRho: Computed invalid density"
			print *, "  Input: P=", P, " T=", T, " R=", R_val
			print *, "  Output: rho=", idealGasRho
			idealGasRho = SMALL_NUM
		end if
    
	end function idealGasRho

	function idealGasP(rho, T, R)
		implicit none
		real(kind=8), intent(in) :: rho, T
		real(kind=8), intent(in), optional :: R
		real(kind=8) :: idealGasP
		real(kind=8) :: R_val
		
		R_val = 287.05d0
		if(present(R)) R_val = R
		
		idealGasP = rho*R_val*T
	end function idealGasP



	function calPerfectEnergy(T, fluid)
		implicit none
		real(kind=8), intent(in) :: T
		type(Fluidd), intent(in) :: fluid
		real(kind=8) :: calPerfectEnergy
		
		calPerfectEnergy = T * (fluid%Cp - fluid%R)
	end function calPerfectEnergy

	function calPerfectT(e, fluid)
		implicit none
		real(kind=8), intent(in) :: e
		type(Fluidd), intent(in) :: fluid
		real(kind=8) :: calPerfectT
		
		calPerfectT = e / (fluid%Cp - fluid%R)
	end function calPerfectT

	!===========================================================================
	! RANS k-omega: Modified decodePrimitives3D (now handles 7 variables)
	!===========================================================================
	
subroutine decodePrimitives3D_RANS(primitives, cellState, fluid)
    implicit none
    real(kind=8), intent(inout) :: primitives(7)  ! RANS: [P, T, Ux, Uy, Uz, k, omega]
    real(kind=8), intent(in)    :: cellState(7)   ! RANS: [rho, rho*u, rho*v, rho*w, rho*E, rho*k, rho*omega]
    type(Fluidd), intent(in)    :: fluid
    
    real(kind=8) :: e, vel_mag, rho, vel_mag_sq
    
  
    real(kind=8), parameter :: MIN_TEMPERATURE = 100.0d0   ! Minimum physical temperature (K)
    real(kind=8), parameter :: MAX_TEMPERATURE = 5000.0d0  ! Maximum physical temperature (K)
    real(kind=8), parameter :: MIN_PRESSURE = 1000.0d0     ! Minimum physical pressure (Pa)
    real(kind=8), parameter :: MAX_VELOCITY = 1000.0d0     ! Maximum physical velocity (m/s)
   
    
    rho = cellState(1)
    
  
    if (rho < SMALL_NUM .or. rho /= rho) then  ! Check for zero or NaN
        print *, "WARNING: Invalid density detected, rho = ", rho
        print *, "         cellState = ", cellState
        rho = max(rho, SMALL_NUM)  ! Use minimum safe value
    end if
  
    ! Velocity components
    primitives(3) = cellState(2) / rho  ! Ux
    primitives(4) = cellState(3) / rho  ! Uy
    primitives(5) = cellState(4) / rho  ! Uz
    
    
    vel_mag_sq = primitives(3)**2 + primitives(4)**2 + primitives(5)**2
    vel_mag = sqrt(vel_mag_sq)
    
    if (vel_mag > MAX_VELOCITY) then
        print *, "WARNING: Excessive velocity detected, |V| = ", vel_mag
        print *, "         Clipping to MAX_VELOCITY = ", MAX_VELOCITY
        ! Scale down velocity components
        primitives(3) = primitives(3) * MAX_VELOCITY / vel_mag
        primitives(4) = primitives(4) * MAX_VELOCITY / vel_mag
        primitives(5) = primitives(5) * MAX_VELOCITY / vel_mag
        vel_mag = MAX_VELOCITY
        vel_mag_sq = vel_mag**2
    end if
   
    
   
    e = (cellState(5) / rho) - 0.5d0 * vel_mag_sq
    
   
    if (e < 0.0d0 .or. e /= e) then
        ! Energy became negative or NaN
        print *, "WARNING: Invalid energy detected in decodePrimitives"
        print *, "  Cell energy e = ", e
        print *, "  Total energy = ", cellState(5) / rho
        print *, "  Kinetic energy = ", 0.5d0 * vel_mag_sq
        print *, "  Velocity mag = ", vel_mag
        
   
        e = (fluid%Cp - fluid%R) * MIN_TEMPERATURE
        print *, "  Corrected to e = ", e, " (T = ", MIN_TEMPERATURE, "K)"
    end if
   
    ! Compute temperature from internal energy
    primitives(2) = calPerfectT(e, fluid)  ! T
    
    !--- ADD TEMPERATURE BOUNDS CHECK (CRITICAL FIX) ---
    if (primitives(2) < MIN_TEMPERATURE .or. primitives(2) /= primitives(2)) then
        print *, "WARNING: Invalid temperature detected, T = ", primitives(2)
        print *, "  Clipping to MIN_TEMPERATURE = ", MIN_TEMPERATURE, "K"
        primitives(2) = MIN_TEMPERATURE
    else if (primitives(2) > MAX_TEMPERATURE) then
        print *, "WARNING: Excessive temperature detected, T = ", primitives(2)
        print *, "  Clipping to MAX_TEMPERATURE = ", MAX_TEMPERATURE, "K"
        primitives(2) = MAX_TEMPERATURE
    end if
    !--- END TEMPERATURE BOUNDS CHECK ---
    
    ! Pressure from equation of state
    primitives(1) = idealGasP(rho, primitives(2), fluid%R)  ! P
    
    !--- ADD PRESSURE BOUNDS CHECK ---
    if (primitives(1) < MIN_PRESSURE .or. primitives(1) /= primitives(1)) then
        print *, "WARNING: Invalid pressure detected, P = ", primitives(1)
        print *, "  Clipping to MIN_PRESSURE = ", MIN_PRESSURE, "Pa"
        primitives(1) = MIN_PRESSURE
    end if
    !--- END PRESSURE BOUNDS CHECK ---
    
    !--- RANS k-omega: Decode turbulence variables ---
    primitives(6) = cellState(6) / rho  ! k
    primitives(7) = cellState(7) / rho  ! omega
    
    !--- ADD TURBULENCE BOUNDS CHECK ---
    ! Ensure positive turbulence quantities with reasonable upper bounds
    primitives(6) = max(primitives(6), SMALL_NUM)
    primitives(6) = min(primitives(6), 0.1d0 * vel_mag_sq)  ! k < 10% of kinetic energy
    
    primitives(7) = max(primitives(7), SMALL_NUM)
    primitives(7) = min(primitives(7), 1.0d10)  ! Cap omega at reasonable value
    !--- END TURBULENCE BOUNDS CHECK ---
    
end subroutine decodePrimitives3D_RANS

	!===========================================================================
	! RANS k-omega: Modified encodePrimitives3D (now handles 7 variables)
	!===========================================================================
	! function encodePrimitives3D(cellPrimitives, fluid)  ! Original
	function encodePrimitives3D_RANS(cellPrimitives, fluid)
    implicit none
    real(kind=8), intent(in) :: cellPrimitives(:,:)
    type(Fluidd), intent(in) :: fluid
    real(kind=8), allocatable :: encodePrimitives3D_RANS(:,:)
    integer :: nCells, c
    real(kind=8) :: e, rho, P, T, vel_mag_sq
    logical :: error_found
    
    nCells = size(cellPrimitives, 1)
    
    print *, "=== encodePrimitives3D_RANS: Starting ==="
    print *, "  nCells = ", nCells
    print *, "  Input dimensions: ", size(cellPrimitives, 1), "x", size(cellPrimitives, 2)
    
    ! Validate input array
    if (size(cellPrimitives, 2) /= 7) then
        print *, "ERROR: cellPrimitives must have 7 columns, got ", size(cellPrimitives, 2)
        stop
    end if
    
    allocate(encodePrimitives3D_RANS(nCells, 7))
    error_found = .false.
    
    do c = 1, nCells
        ! Extract values for debugging
        P = cellPrimitives(c, 1)
        T = cellPrimitives(c, 2)
        
        ! Check for NaN in input
        if (P /= P .or. T /= T) then
            print *, "ERROR in encodePrimitives: NaN detected in input at cell ", c
            print *, "  P = ", P
            print *, "  T = ", T
            print *, "  Full row = ", cellPrimitives(c, :)
            error_found = .true.
            cycle
        end if
        
        ! Check for invalid values
        if (P <= 0.0d0 .or. T <= 0.0d0) then
            print *, "ERROR in encodePrimitives: Invalid P or T at cell ", c
            print *, "  P = ", P, " (must be > 0)"
            print *, "  T = ", T, " (must be > 0)"
            error_found = .true.
            cycle
        end if
        
        ! Density
        rho = idealGasRho(T, P, fluid%R)
        
        ! Check computed density
        if (rho /= rho .or. rho <= 0.0d0) then
            print *, "ERROR: idealGasRho returned invalid density at cell ", c
            print *, "  Input: P=", P, " T=", T, " R=", fluid%R
            print *, "  Output: rho=", rho
            error_found = .true.
            rho = max(rho, SMALL_NUM)
        end if
        
        encodePrimitives3D_RANS(c,1) = rho
        
        ! Momentum
        encodePrimitives3D_RANS(c,2:4) = cellPrimitives(c,3:5) * rho
        
        ! Total energy
        e = calPerfectEnergy(T, fluid)
        vel_mag_sq = cellPrimitives(c,3)**2 + cellPrimitives(c,4)**2 + cellPrimitives(c,5)**2
        encodePrimitives3D_RANS(c,5) = rho * (e + 0.5d0 * vel_mag_sq)
        
        ! Turbulence conservative variables
        encodePrimitives3D_RANS(c,6) = rho * cellPrimitives(c,6)  ! rho*k
        encodePrimitives3D_RANS(c,7) = rho * cellPrimitives(c,7)  ! rho*omega
        
        ! Debug first cell
        if (c == 1) then
            print *, "  Cell 1 encoding:"
            print *, "    P, T = ", P, T
            print *, "    rho = ", rho
            print *, "    e = ", e
            print *, "    vel_mag_sq = ", vel_mag_sq
            print *, "    Encoded state = ", encodePrimitives3D_RANS(c, :)
        end if
    end do
    
    if (error_found) then
        print *, "FATAL: Errors found during encoding, stopping"
        stop
    end if
    
    print *, "=== encodePrimitives3D_RANS: Complete ==="
    print *, "  Output rho range: ", minval(encodePrimitives3D_RANS(:,1)), maxval(encodePrimitives3D_RANS(:,1))
    
	end function encodePrimitives3D_RANS

	!===========================================================================
	! RANS k-omega: Modified calculateFluxes3D (includes viscous and turbulent fluxes)
	!===========================================================================
	! subroutine calculateFluxes3D(out_fluxes, prim, state)  ! Original
	subroutine calculateFluxes3D_RANS(out_fluxes, prim, state, fluid, sln, cellIdx)
		implicit none
		! real(kind=8), intent(out) :: out_fluxes(15)  ! Original
		! real(kind=8), intent(in) :: prim(5), state(5)  ! Original
		real(kind=8), intent(out) :: out_fluxes(21)  ! RANS: 7 vars × 3 directions
		real(kind=8), intent(in) :: prim(7), state(7)  ! RANS: includes k, omega
		type(Fluidd), intent(in) :: fluid
		type(SolutionState), intent(in) :: sln
		integer, intent(in) :: cellIdx
		
		real(kind=8) :: mu_eff, k_eff, rho
		
		rho = state(1)
		
		!--- Original inviscid fluxes ---
		! Mass fluxes
		out_fluxes(1) = state(2)  ! rho*Ux
		out_fluxes(2) = state(3)  ! rho*Uy
		out_fluxes(3) = state(4)  ! rho*Uz
		
		! Momentum fluxes (inviscid part)
		out_fluxes(4) = state(2)*prim(3) + prim(1)   ! x-momentum x-flux
		out_fluxes(7) = state(2)*prim(4)              ! y-momentum x-flux
		out_fluxes(10) = state(2)*prim(5)             ! z-momentum x-flux
		out_fluxes(5) = out_fluxes(7)                 ! x-momentum y-flux (symmetric)
		out_fluxes(8) = state(3)*prim(4) + prim(1)   ! y-momentum y-flux
		out_fluxes(11) = state(3)*prim(5)             ! z-momentum y-flux
		out_fluxes(6) = out_fluxes(10)                ! x-momentum z-flux (symmetric)
		out_fluxes(9) = out_fluxes(11)                ! y-momentum z-flux (symmetric)
		out_fluxes(12) = state(4)*prim(5) + prim(1)  ! z-momentum z-flux
		
		! Energy fluxes (inviscid part)
		out_fluxes(13) = prim(3)*state(5) + prim(1)*prim(3)
		out_fluxes(14) = prim(4)*state(5) + prim(1)*prim(4)
		out_fluxes(15) = prim(5)*state(5) + prim(1)*prim(5)
		
		!--- RANS k-omega: Turbulence fluxes ---
		! k fluxes (convective part)
		out_fluxes(16) = state(6) * prim(3)  ! rho*k*Ux
		out_fluxes(17) = state(6) * prim(4)  ! rho*k*Uy
		out_fluxes(18) = state(6) * prim(5)  ! rho*k*Uz
		
		! omega fluxes (convective part)
		out_fluxes(19) = state(7) * prim(3)  ! rho*omega*Ux
		out_fluxes(20) = state(7) * prim(4)  ! rho*omega*Uy
		out_fluxes(21) = state(7) * prim(5)  ! rho*omega*Uz
		
				
	end subroutine calculateFluxes3D_RANS

!===============================================================================
! SECTION 5: JST SCHEME - RANS k-omega modifications
!===============================================================================

	!===========================================================================
	! RANS k-omega: Modified JST artificial dissipation (same structure, more variables)
	!===========================================================================
	function unstructured_JSTEps_RANS(mesh, sln, fluid, k2, k4, c4) result(eps)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(in) :: sln
		type(Fluidd), intent(in) :: fluid
		real(kind=8), intent(in), optional :: k2, k4, c4
		real(kind=8), allocatable :: eps(:,:)
		
		
		
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4)
		real(kind=8), allocatable :: P(:), gradP(:,:), P_matrix(:,:), temp_grad(:,:,:)
		real(kind=8), allocatable :: sj(:), sjCount(:), rj(:), rjsjF(:,:)
		real(kind=8), allocatable :: eps2(:), eps4(:)
		integer(kind=8) :: f, ownerCell, neighbourCell, c
		real(kind=8) :: d(3), oP, nP, farOwnerP, farNeighbourP, epsilonn
		real(kind=8) :: k2_val, k4_val, c4_val
		
		k2_val = 0.5d0
		k4_val = 1.0d0/32.0d0
		c4_val = 1.0d0
		if(present(k2)) k2_val = k2
		if(present(k4)) k4_val = k4
		if(present(c4)) c4_val = c4
		
		meshInfo = unstructuredMeshInfo(mesh)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		
		
		allocate(P(nCells))
		P = sln%cellPrimitives(:,1)
		
		allocate(P_matrix(nCells, 1))
		P_matrix(:, 1) = P
		allocate(temp_grad(nCells, 1, 3))
		temp_grad = greenGaussGrad_RANS(mesh, P_matrix, .false.)
		allocate(gradP(nCells, 3))
		! Density-based JST pressure sensor follows the original x-gradient-only form.
		gradP = 0.0d0
		gradP(:,1) = temp_grad(:,1,1)
		
		allocate(sj(nCells), sjCount(nCells))
		sj = 0.0d0
		sjCount = 0.0d0
		epsilonn = 1.0d-10
		
		
		do f = 1, nFaces - nBdryFaces
			ownerCell = mesh%faces(f, 1)
			neighbourCell = mesh%faces(f, 2)
			
			d = mesh%cCenters(neighbourCell, :) - mesh%cCenters(ownerCell, :)
			
			oP = P(ownerCell)
			nP = P(neighbourCell)
			
			farOwnerP = nP - 2.0d0 * dot_product(d, gradP(ownerCell, :))
			farNeighbourP = oP + 2.0d0 * dot_product(d, gradP(neighbourCell, :))
			
			sj(ownerCell) = sj(ownerCell) + (abs(nP - 2.0d0*oP + farOwnerP) / &
			    max(abs(nP - oP) + abs(oP - farOwnerP), epsilonn))**2
			sjCount(ownerCell) = sjCount(ownerCell) + 1
			
			sj(neighbourCell) = sj(neighbourCell) + (abs(oP - 2.0d0*nP + farNeighbourP) / &
			    max(abs(farNeighbourP - nP) + abs(nP - oP), epsilonn))**2
			sjCount(neighbourCell) = sjCount(neighbourCell) + 1
		end do
		
		
		allocate(rj(nCells))
		do c = 1, nCells
			rj(c) = preconditioned_spectral_radius( &
				mag(sln%cellPrimitives(c,3:5)), &
				sqrt(fluid%gammaa * fluid%R * sln%cellPrimitives(c,2)), &
				fluid%gammaa)
			if(sjCount(c) > 0) then
				sj(c) = sj(c) / sjCount(c)
			end if
		end do
		
		allocate(rjsjF(nFaces, 2))
		rjsjF = maxInterp(mesh, sj, rj)
		
		! Compute dissipation coefficients
		allocate(eps2(nFaces), eps4(nFaces))
		eps2 = 0.0d0
		eps4 = 0.0d0
		
		do f = 1, nFaces - nBdryFaces
			eps2(f) = k2_val * rjsjF(f,2) * rjsjF(f,1)
			eps4(f) = max(0.0d0, k4_val * rjsjF(f,1) - c4_val * eps2(f))
		end do
		
		allocate(eps(nFaces, 2))
		eps(:,1) = eps2
		eps(:,2) = eps4
		
		deallocate(P, P_matrix, temp_grad, gradP)
		deallocate(sj, sjCount, rj, rjsjF, eps2, eps4)
	end function unstructured_JSTEps_RANS

	!===========================================================================
	! RANS k-omega: Modified JST flux computation (with viscous terms)
	!===========================================================================
	function unstructured_JSTFlux_RANS(mesh, sln, boundaryConditions, fluid) result(unstructured_JSTFluxx)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(inout) :: sln
		type(BoundaryCondition), intent(in) :: boundaryConditions(:)
		type(Fluidd), intent(in) :: fluid
		real(kind=8), allocatable :: unstructured_JSTFluxx(:,:)
		
		! integer(kind=8) :: nVars  ! Original: 5
		integer(kind=8) :: nVars  ! RANS: 7
		real(kind=8) :: d(3), grad_v(3)
		real(kind=8), allocatable :: fDeltas(:,:), fDGrads(:,:,:)
		real(kind=8), allocatable :: eps2(:), eps4(:), diffusionFlux(:), unitFA(:)
		real(kind=8), allocatable :: fD(:), farOwnerfD(:), farNeighbourfD(:), eps(:,:)
		real(kind=8), allocatable :: gradU(:,:,:), gradK(:,:,:), gradOmega(:,:,:)
		real(kind=8), allocatable :: viscousFlux(:), turbulentDiffFlux(:)
		integer(kind=8) :: f, ownerCell, neighbourCell, v, i1, i2
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4), i, j
		real(kind=8) :: mu_eff, rho_face, mu_L_face, mu_T_face, sigma_k_eff, sigma_omega_eff
		real(kind=8) :: grad_k_n, grad_omega_n, face_area
		
		meshInfo = unstructuredMeshInfo(mesh)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		nVars = size(sln%cellState, 2)  ! 7 for RANS
		
		!--- 1. Apply boundary conditions ---
		call updateBoundaryConditions_RANS(mesh, sln, boundaryConditions, nBoundaries, fluid)
		
		!--- 2. Compute central difference fluxes ---
		call linInterp_3D_RANS(mesh, sln%cellFluxes, sln%faceFluxes)
		
		!--- 3. Compute JST artificial dissipation ---
		allocate(fDeltas(nFaces, nVars))
		fDeltas = faceDeltas_RANS(mesh, sln)
		
		allocate(fDGrads(nCells, nVars, 3))
		fDGrads = greenGaussGrad_RANS(mesh, fDeltas, .false.)
		
		allocate(eps2(nFaces), eps4(nFaces))
		eps = unstructured_JSTEps_RANS(mesh, sln, fluid)
		eps2 = eps(:,1)
		eps4 = eps(:,2)
		deallocate(eps)
		
		allocate(diffusionFlux(nVars), unitFA(3), fD(nVars))
		allocate(farOwnerfD(nVars), farNeighbourfD(nVars))
		allocate(viscousFlux(nVars), turbulentDiffFlux(nVars))
		
		!--- 4. Compute velocity gradients for viscous terms ---
		allocate(gradU(nCells, 3, 3))  ! gradU(cell, velocity_component, space_direction)
		call computeVelocityGradients(mesh, sln, gradU)
		
		!--- 5. Compute turbulence gradients ---
		allocate(gradK(nCells, 1, 3))
		allocate(gradOmega(nCells, 1, 3))
		call computeTurbulenceGradients(mesh, sln, gradK, gradOmega)
		
		!--- 6. Loop over interior faces ---
		do f = 1, nFaces-nBdryFaces
			ownerCell = mesh%faces(f, 1)
			neighbourCell = mesh%faces(f, 2)
			d = mesh%cCenters(neighbourCell, :) - mesh%cCenters(ownerCell, :)
			
			!--- JST artificial dissipation (same as original) ---
			fD = fDeltas(f,:)
			do v = 1, nVars
				grad_v = fDGrads(ownerCell, v, :)
				farOwnerfD(v) = fD(v) - dot_product(d, grad_v)
				farNeighbourfD(v) = fD(v) + dot_product(d, grad_v)
			end do
			
			diffusionFlux = eps2(f) * fD - eps4(f) * (farNeighbourfD - 2*fD + farOwnerfD)
			
			unitFA = normalize(mesh%fAVecs(f,:))
			face_area = mag(mesh%fAVecs(f,:))
			
			!--- Add JST dissipation ---
			do v = 1, nVars
				i1 = (v-1)*3 + 1
				i2 = i1 + 2
				sln%faceFluxes(f,i1:i2) = sln%faceFluxes(f,i1:i2) - (diffusionFlux(v) * unitFA)
			end do
			
			!--- RANS k-omega: Add viscous fluxes ---
			call computeViscousFluxes_RANS(mesh, sln, f, ownerCell, neighbourCell, &
			                               gradU, fluid, viscousFlux)
			
			
			do v = 1, 3  ! 3 momentum equations
				i1 = (v)*3 + 1  ! Start index for momentum v
				i2 = i1 + 2
				sln%faceFluxes(f,i1:i2) = sln%faceFluxes(f,i1:i2) + &
				                          viscousFlux(v+1) * unitFA * face_area
			end do
			
			! Energy: index 13-15
			sln%faceFluxes(f,13:15) = sln%faceFluxes(f,13:15) + &
			                          viscousFlux(5) * unitFA * face_area
			
			!--- RANS k-omega: Add turbulent diffusion ---
			! Interpolate properties to face
			rho_face = 0.5d0 * (sln%cellState(ownerCell,1) + sln%cellState(neighbourCell,1))
			mu_L_face = 0.5d0 * (sln%cellMuL(ownerCell) + sln%cellMuL(neighbourCell))
			mu_T_face = 0.5d0 * (sln%cellMuT(ownerCell) + sln%cellMuT(neighbourCell))
			
			! Effective diffusivity for k
			sigma_k_eff = mu_L_face + mu_T_face / SIGMA_K
			
			! Effective diffusivity for omega
			sigma_omega_eff = mu_L_face + mu_T_face / SIGMA_OMEGA
			
			! Gradient normal to face
			grad_k_n = dot_product(0.5d0*(gradK(ownerCell,1,:) + gradK(neighbourCell,1,:)), unitFA)
			grad_omega_n = dot_product(0.5d0*(gradOmega(ownerCell,1,:) + gradOmega(neighbourCell,1,:)), unitFA)
			
			! Turbulent diffusion fluxes
			! k: indices 16-18
			sln%faceFluxes(f,16:18) = sln%faceFluxes(f,16:18) - &
			                          sigma_k_eff * grad_k_n * unitFA * face_area
			
			! omega: indices 19-21
			sln%faceFluxes(f,19:21) = sln%faceFluxes(f,19:21) - &
			                          sigma_omega_eff * grad_omega_n * unitFA * face_area
		end do
		
		!--- 7. Integrate fluxes to get cell residuals ---
		allocate(unstructured_JSTFluxx(nCells, nVars))
		unstructured_JSTFluxx = integrateFluxes_unstructured3D_RANS(mesh, sln)
		
		!--- Cleanup ---
		deallocate(fDeltas, fDGrads, eps2, eps4, diffusionFlux, unitFA, fD)
		deallocate(farOwnerfD, farNeighbourfD, viscousFlux, turbulentDiffFlux)
		deallocate(gradU, gradK, gradOmega)
	end function unstructured_JSTFlux_RANS

!===============================================================================
! SECTION 6: NEW RANS-SPECIFIC FUNCTIONS
!===============================================================================

	!===========================================================================
	! RANS k-omega: Compute wall distances
	!===========================================================================
subroutine computeWallDistances(mesh)
    implicit none
    type(Meshh), intent(inout) :: mesh
    integer(kind=8) :: nCells, nFaces, nBoundaries, meshInfo(4)
    integer :: c, b, f, face_idx
    real(kind=8) :: min_dist, dist
    real(kind=8) :: cell_center(3), wall_center(3)
    integer :: total_wall_faces, boundary_face_count
    integer :: nFacesInBoundary
    
    meshInfo = unstructuredMeshInfo(mesh)
    nCells = meshInfo(1)
    nFaces = meshInfo(2)
    nBoundaries = meshInfo(3)
    
    print *, "DEBUG: Computing wall distances for ", nCells, " cells"
    print *, "DEBUG: Number of boundaries: ", nBoundaries
    
    ! Allocate wall distance array
    if (.not. allocated(mesh%wallDistance)) then
        allocate(mesh%wallDistance(nCells))
    end if
    
    ! Initialize to large value
    mesh%wallDistance = 1.0d10
    total_wall_faces = 0
    
    ! First pass: count total wall faces across all boundaries
    do b = 1, nBoundaries
        boundary_face_count = 0
        do f = 1, size(mesh%boundaryFaces, 2)
            if (mesh%boundaryFaces(b, f) == 0) exit
            boundary_face_count = boundary_face_count + 1
        end do
        print *, "  Boundary ", b, " has ", boundary_face_count, " faces"
        total_wall_faces = total_wall_faces + boundary_face_count
    end do
    
    print *, "DEBUG: Total boundary faces to process: ", total_wall_faces
    
    if (total_wall_faces == 0) then
        print *, "WARNING: No boundary faces found! Using default distance"
        mesh%wallDistance = 0.01d0
        return
    end if
    
    ! Second pass: compute actual distances
    do c = 1, nCells
        cell_center = mesh%cCenters(c,:)
        min_dist = 1.0d10
        
        ! Loop over all boundaries
        do b = 1, nBoundaries
            ! Loop over faces in this boundary
            nFacesInBoundary = 0
            do f = 1, size(mesh%boundaryFaces, 2)
                face_idx = mesh%boundaryFaces(b, f)
                if (face_idx == 0) exit  ! No more faces in this boundary
                if (face_idx < 1 .or. face_idx > nFaces) cycle  ! Invalid face index
                
                nFacesInBoundary = nFacesInBoundary + 1
                
                ! Get wall face center
                wall_center = mesh%fCenters(face_idx,:)
                
                ! Compute distance
                dist = mag(cell_center - wall_center)
                
                ! Update minimum distance
                if (dist < min_dist) then
                    min_dist = dist
                end if
            end do
        end do
        
        mesh%wallDistance(c) = min_dist
        
        ! Debug first few cells
        if (c <= 3) then
            print *, "  Cell ", c, " wall distance: ", mesh%wallDistance(c)
        end if
    end do
    
    print *, "Wall distances computed. Min:", minval(mesh%wallDistance), &
             " Max:", maxval(mesh%wallDistance)
    print *, "Total boundary faces found:", total_wall_faces
    
end subroutine computeWallDistances

	!===========================================================================
	! RANS k-omega: Update turbulence fields
	!===========================================================================
	subroutine updateTurbulenceFields(mesh, sln, fluid)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(inout) :: sln
		type(Fluidd), intent(in) :: fluid
		integer :: c, nCells
		real(kind=8) :: rho, k, omega, T
		
		nCells = size(sln%cellState, 1)
		  !--- SAFETY CHECK: Ensure wall distances are computed ---
		if (.not. allocated(mesh%wallDistance)) then
			print *, "ERROR: Wall distances not computed before updateTurbulenceFields!"
			stop
		end if
		if (size(mesh%wallDistance) /= nCells) then
			print *, "ERROR: Wall distance array size mismatch!"
			print *, "  Expected:", nCells, " Got:", size(mesh%wallDistance)
			stop
		end if
    !--- END SAFETY CHECK ---
		
		do c = 1, nCells
			rho = sln%cellState(c, 1)
			k = sln%cellPrimitives(c, 6)
			omega = sln%cellPrimitives(c, 7)
			T = sln%cellPrimitives(c, 2)
			
			     !--- Additional safety checks ---
			if (rho < SMALL_NUM .or. rho /= rho) then
				print *, "WARNING: Invalid rho in cell", c, " rho=", rho
				rho = max(rho, SMALL_NUM)
				sln%cellState(c, 1) = rho
			end if
			if (T < SMALL_NUM .or. T /= T) then
				print *, "WARNING: Invalid T in cell", c, " T=", T
				T = 300.0d0  ! Default temperature
				sln%cellPrimitives(c, 2) = T
			end if
        !--- End additional checks ---
			
			! Ensure positive values
			k = max(k, SMALL_NUM)
			omega = max(omega, SMALL_NUM)
			
			! Update laminar viscosity
			sln%cellMuL(c) = computeLaminarViscosity(T)
			
			! Update turbulent viscosity
			sln%cellMuT(c) = computeTurbulentViscosity(rho, k, omega)
			
			! Update wall distance (already computed in mesh)
			sln%cellWallDist(c) = mesh%wallDistance(c)
		end do
	end subroutine updateTurbulenceFields

	!===========================================================================
	! RANS k-omega: Compute velocity gradients
	!===========================================================================
	subroutine computeVelocityGradients(mesh, sln, gradU)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(in) :: sln
		real(kind=8), intent(out) :: gradU(:,:,:)  ! (nCells, 3 components, 3 directions)
		real(kind=8), allocatable :: velocity(:,:), temp_grad(:,:,:)
		integer :: nCells, i
		
		nCells = size(sln%cellPrimitives, 1)
		
		! Extract velocity components
		allocate(velocity(nCells, 3))
		velocity = sln%cellPrimitives(:, 3:5)  ! [Ux, Uy, Uz]
		
		! Compute gradients using Green-Gauss
		allocate(temp_grad(nCells, 3, 3))
		temp_grad = greenGaussGrad_RANS(mesh, velocity, .false.)
		
		! Copy to output array
		gradU = temp_grad
		
		deallocate(velocity, temp_grad)
	end subroutine computeVelocityGradients

	!===========================================================================
	! RANS k-omega: Compute turbulence gradients
	!===========================================================================
	subroutine computeTurbulenceGradients(mesh, sln, gradK, gradOmega)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(in) :: sln
		real(kind=8), intent(out) :: gradK(:,:,:), gradOmega(:,:,:)
		real(kind=8), allocatable :: k_matrix(:,:), omega_matrix(:,:)
		integer :: nCells
		
		nCells = size(sln%cellPrimitives, 1)
		
		! Extract k
		allocate(k_matrix(nCells, 1))
		k_matrix(:,1) = sln%cellPrimitives(:, 6)
		gradK = greenGaussGrad_RANS(mesh, k_matrix, .false.)
		deallocate(k_matrix)
		
		! Extract omega
		allocate(omega_matrix(nCells, 1))
		omega_matrix(:,1) = sln%cellPrimitives(:, 7)
		gradOmega = greenGaussGrad_RANS(mesh, omega_matrix, .false.)
		deallocate(omega_matrix)
	end subroutine computeTurbulenceGradients

	!===========================================================================
	! RANS k-omega: Compute viscous fluxes
	!===========================================================================
subroutine computeViscousFluxes_RANS(mesh, sln, faceIdx, owner, neighbour, &
                                     gradU, fluid, viscousFlux)
    implicit none
    type(Meshh), intent(in) :: mesh
    type(SolutionState), intent(in) :: sln
    integer(kind=8), intent(in) :: faceIdx, owner, neighbour
    real(kind=8), intent(in) :: gradU(:,:,:)
    type(Fluidd), intent(in) :: fluid
    real(kind=8), intent(out) :: viscousFlux(:)
    
    real(kind=8) :: mu_eff, tau(3,3), face_normal(3), area
    real(kind=8) :: grad_U_face(3,3), T_owner, T_neighbour, T_face
    real(kind=8) :: grad_T(3), k_thermal, q_dot_n
    real(kind=8) :: u_face(3), tau_dot_u
    integer :: i, j
    
    ! Average gradients to face
    do i = 1, 3
        do j = 1, 3
            grad_U_face(i,j) = 0.5d0 * (gradU(owner,i,j) + gradU(neighbour,i,j))
        end do
    end do
    
    ! Effective viscosity at face
    mu_eff = 0.5d0 * (sln%cellMuL(owner) + sln%cellMuT(owner) + &
                      sln%cellMuL(neighbour) + sln%cellMuT(neighbour))
    
    ! Compute stress tensor
    call computeStressTensor(grad_U_face, mu_eff, tau)
    
    ! Face normal and area
    face_normal = normalize(mesh%fAVecs(faceIdx,:))
    area = mag(mesh%fAVecs(faceIdx,:))
    
    ! Initialize
    viscousFlux = 0.0d0
    
    ! Viscous momentum flux: tau · n
    viscousFlux(2) = dot_product(tau(1,:), face_normal)  ! x-momentum
    viscousFlux(3) = dot_product(tau(2,:), face_normal)  ! y-momentum
    viscousFlux(4) = dot_product(tau(3,:), face_normal)  ! z-momentum
    
   
    ! Face velocity (interpolated)
    u_face(1) = 0.5d0 * (sln%cellPrimitives(owner,3) + sln%cellPrimitives(neighbour,3))
    u_face(2) = 0.5d0 * (sln%cellPrimitives(owner,4) + sln%cellPrimitives(neighbour,4))
    u_face(3) = 0.5d0 * (sln%cellPrimitives(owner,5) + sln%cellPrimitives(neighbour,5))
    
    ! Work done by viscous stresses: tau · u
    tau_dot_u = 0.0d0
    do i = 1, 3
        tau_dot_u = tau_dot_u + dot_product(tau(i,:), face_normal) * u_face(i)
    end do
    
    ! Temperature gradient (simple approximation)
    T_owner = sln%cellPrimitives(owner, 2)
    T_neighbour = sln%cellPrimitives(neighbour, 2)
    T_face = 0.5d0 * (T_owner + T_neighbour)
    
    ! Temperature gradient normal to face (simple finite difference)
    grad_T = face_normal * (T_neighbour - T_owner) / &
             mag(mesh%cCenters(neighbour,:) - mesh%cCenters(owner,:))
    
    ! Thermal conductivity (k = mu * Cp / Pr)
    k_thermal = mu_eff * fluid%Cp / fluid%Pr
    
    ! Heat flux: q = -k * grad(T)
    q_dot_n = -k_thermal * dot_product(grad_T, face_normal)
    
    ! Viscous energy flux = tau·u - q·n
    viscousFlux(5) = tau_dot_u - q_dot_n
    
end subroutine computeViscousFluxes_RANS

	!===========================================================================
	! RANS k-omega: Compute stress tensor
	!===========================================================================
	subroutine computeStressTensor(gradU, mu_eff, tau)
		implicit none
		real(kind=8), intent(in) :: gradU(3,3), mu_eff
		real(kind=8), intent(out) :: tau(3,3)
		real(kind=8) :: div_U
		integer :: i, j
		
		! Divergence of velocity
		div_U = gradU(1,1) + gradU(2,2) + gradU(3,3)
		
		! Stress tensor: tau_ij = mu * (du_i/dx_j + du_j/dx_i) - 2/3 * mu * div(U) * delta_ij
		do i = 1, 3
			do j = 1, 3
				tau(i,j) = mu_eff * (gradU(i,j) + gradU(j,i))
				if (i == j) then
					tau(i,j) = tau(i,j) - (2.0d0/3.0d0) * mu_eff * div_U
				end if
			end do
		end do
	end subroutine computeStressTensor

	!===========================================================================
	! RANS k-omega: Compute source terms for k and omega equations
	!===========================================================================
subroutine computeTurbulenceSourceTerms(mesh, sln, fluid, sources)
    implicit none
    type(Meshh), intent(in) :: mesh
    type(SolutionState), intent(inout) :: sln
    type(Fluidd), intent(in) :: fluid
    real(kind=8), intent(out) :: sources(:,:)  ! (nCells, 2) for k and omega
    
    integer :: c, nCells
    real(kind=8) :: rho, k, omega, mu_t, S_mag
    real(kind=8) :: P_k, D_k, P_omega, D_omega
    real(kind=8), allocatable :: gradU(:,:,:)
    
    !--- ADD THESE PARAMETERS ---
    real(kind=8), parameter :: RELAX_K = 0.5d0       ! Under-relaxation for k source
    real(kind=8), parameter :: RELAX_OMEGA = 0.5d0   ! Under-relaxation for omega source
    real(kind=8), parameter :: PROD_LIMIT = 10.0d0   ! Limit production/dissipation ratio
    !--- END NEW PARAMETERS ---
    
    nCells = size(sln%cellState, 1)
    allocate(gradU(nCells, 3, 3))
    call computeVelocityGradients(mesh, sln, gradU)
    
    do c = 1, nCells
        rho = sln%cellState(c, 1)
        k = sln%cellPrimitives(c, 6)
        omega = sln%cellPrimitives(c, 7)
        mu_t = sln%cellMuT(c)
        
        ! Strain rate magnitude
        S_mag = computeStrainRateMagnitude(gradU(c,:,:))
        sln%cellStrainRate(c) = S_mag
        
        !--- k equation source terms ---
        ! Production: P_k = mu_t * S^2
        P_k = mu_t * S_mag**2
        
        ! Dissipation: D_k = beta_star * rho * omega * k
        D_k = BETA_STAR * rho * omega * k
        
        
        if (D_k > SMALL_NUM) then
            P_k = min(P_k, PROD_LIMIT * D_k)
        end if
        !--- END LIMITER ---
        
        !--- APPLY UNDER-RELAXATION TO k SOURCE ---
        sources(c, 1) = RELAX_K * (P_k - D_k)
        !--- OLD CODE WAS: sources(c, 1) = P_k - D_k
        
        sln%cellProduction(c) = P_k
        
        !--- omega equation source terms ---
        ! Production: P_omega = alpha * omega/k * P_k
        if (k > SMALL_NUM) then
            P_omega = ALPHA * (omega / k) * P_k
        else
            P_omega = 0.0d0
        end if
        
        ! Dissipation: D_omega = beta * rho * omega^2
        D_omega = BETA * rho * omega**2
        
        !--- ADD PRODUCTION LIMITER FOR OMEGA ---
        if (D_omega > SMALL_NUM) then
            P_omega = min(P_omega, PROD_LIMIT * D_omega)
        end if
        !--- END LIMITER ---
        
        !--- APPLY UNDER-RELAXATION TO omega SOURCE ---
        sources(c, 2) = RELAX_OMEGA * (P_omega - D_omega)
        !--- OLD CODE WAS: sources(c, 2) = P_omega - D_omega
        
    end do
    
    deallocate(gradU)
end subroutine computeTurbulenceSourceTerms
	!===============================================================================
	! SAME FUNCTIONS OF SECTION 6
	!===============================================================================
function triangleCentroid(points)
		implicit none
		real(kind=8), intent(in) :: points(:,:)
		real(kind=8) :: triangleCentroid(3)
		integer :: nPts, pt
		
		triangleCentroid = [0.0d0, 0.0d0, 0.0d0]
		nPts = size(points, 1)
		
		do pt = 1, nPts
			triangleCentroid = triangleCentroid + points(pt,:)   ! 累加每个点坐标
		end do
		
		triangleCentroid = triangleCentroid / nPts   ! 求平均，即几何中心
	end function triangleCentroid

	! 为几何中心函数设置一个别名
	function geometricCenter(points)
		implicit none
		real(kind=8), intent(in) :: points(:,:)
		real(kind=8) :: geometricCenter(3)
		geometricCenter = triangleCentroid(points)
	end function geometricCenter

	! 计算三角形面积向量!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function triangleArea(points)
		implicit none
		real(kind=8), intent(in) :: points(:,:)
		real(kind=8) :: triangleArea(3)
		real(kind=8) :: side1(3), side2(3)
		
		side1 = points(2,:) - points(1,:)
		side2 = points(3,:) - points(1,:)
		triangleArea = cross(side1, side2) / 2.0d0   ! 使用向量叉乘除以2
	end function triangleArea

	! 计算多边形面元的面积向量与几何中心!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine faceAreaCentroid(points, first5_refs, fAVec, centroid)
	!subroutine faceAreaCentroid(points, fAVec, centroid)
		implicit none
		real(kind=8), intent(in) :: points(:,:)
		real(kind=8), intent(out) :: fAVec(3), centroid(3)
		real(kind=8) :: gC(3), subTriPts(3,3), triCentroid(3), subFAVec(3)
		integer :: nPts, i
		logical :: print_output  ! 标记是否需要输出
		real(kind=8), intent(in) :: first5_refs(5,3)   ! 前5个面的参考顶点坐标!!!!!!!!!!!!!!
		integer :: current_face_id, f ! 当前面的索引（1-5）!!!!!!!!!!!!!!!!!!!
		logical :: is_in_first5!!!!!!!!!
		! 判断是否为前5个面（通过第一个顶点坐标匹配）!!!!!!!!
		is_in_first5 = .false.
		current_face_id = 0
		if (DEBUG_MESH_VERBOSE) then
		do f = 1, 5
		    ! 允许微小数值误差（1e-6）
		    if (abs(points(1,1)-first5_refs(f,1)) < 1d-6 .and. &
		        abs(points(1,2)-first5_refs(f,2)) < 1d-6 .and. &
		        abs(points(1,3)-first5_refs(f,3)) < 1d-6) then
		        is_in_first5 = .true.
		        current_face_id = f
		        exit
		    end if
		end do!!!!!!!!!!!!!!!!
		end if
		gC = geometricCenter(points)       ! 得到整个面的几何中心
		nPts = size(points, 1)        
		fAVec = [0.0d0, 0.0d0, 0.0d0]      ! 初始化面面积向量
		centroid = [0.0d0, 0.0d0, 0.0d0]   ! 初始化面中心        
	 	! 仅当处理第1个面时才输出
		print_output = (f == 1)!!!!!!!!!!!!!!!
		! 将多边形面分割成多个由几何中心和两个相邻顶点组成的三角形!!!!!!!!!!!!!!!!!!!!
		do i = 1, nPts
			if (i < nPts) then
			    subTriPts(1,:) = gC
			    subTriPts(2,:) = points(i,:)
			    subTriPts(3,:) = points(i+1,:)
			else
			    subTriPts(1,:) = gC
			    subTriPts(2,:) = points(i,:)
			    subTriPts(3,:) = points(1,:)  ! 闭合面
			end if
            !subTriPts没问题
			triCentroid = triangleCentroid(subTriPts)    ! 子三角形中心
			subFAVec = triangleArea(subTriPts)           ! 子三角形面积向量
			! 前5个面输出subFAVec
		    if (is_in_first5) then
		        print *, "=== 第", current_face_id, "个面的子三角形 ", i, " ==="
		        print *, "  subFAVec (x,y,z): (", subFAVec(1), ",", subFAVec(2), ",", subFAVec(3), ")"
		        print *, "  面积大小: ", mag(subFAVec)  ! 面积大小（标量）
		    end if
			fAVec = fAVec + subFAVec
			centroid = centroid + triCentroid * mag(subFAVec)  ! 质心加权（面积加权）
		end do

		centroid = centroid / mag(fAVec)   ! 总质心为加权平均
		if (is_in_first5) then
		        print *, "centroid=", centroid 
        end if
	end subroutine faceAreaCentroid

	! 原调用处:cVols[c], cCenters[c] = cellVolCentroid(pts, cell_fAVecs, fCs),fCs为一个单元中各个面的几何中心,cell_fAVecs为这个包对应面的法向量
	! 计算单元体积和中心位置!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine cellVolCentroid(points, fAVecs, faceCentroids, vol, centroid)
		implicit none
		! 假设输入数组维度：
		! points: (nPoints, 3)  单元所有顶点坐标
		! fAVecs: (nFaces, 3)   每个面的面积向量（每行一个面）
		! faceCentroids: (nFaces, 3) 每个面的中心（每行一个面）
		real(kind=8), intent(in) :: points(:,:), fAVecs(:,:), faceCentroids(:,:)
		real(kind=8), intent(out) :: vol, centroid(3)
		real(kind=8) :: gC(3), cellCenterVec(3), subPyrVol, subPyrCentroid(3)
		integer :: nFaces, f
		
		! 计算单元几何中心（确保与Julia的geometricCenter逻辑一致）
		gC = geometricCenter(points)
		
		! 关键修正：正确获取面数（若fAVecs是(nFaces,3)，则用size(fAVecs,1)）
		! 若fAVecs是(3,nFaces)，则改为 nFaces = size(fAVecs, 2)
		nFaces = size(fAVecs, 1)  

		vol = 0.0d0
		centroid = [0.0d0, 0.0d0, 0.0d0]

		do f = 1, nFaces
		    ! 计算面中心到单元几何中心的向量（修正数组访问）
		    cellCenterVec = faceCentroids(f,:) - gC  ! 若faceCentroids是(3,nFaces)，则用faceCentroids(:,f)
		    
		    ! 计算子金字塔体积（点积+绝对值）
		    subPyrVol = abs( dot_product(fAVecs(f,:), cellCenterVec) ) / 3.0d0  
		    ! 改用dot_product确保点积计算正确（与Julia的sum(.+)一致）
		    
		    ! 计算子金字塔中心
		    subPyrCentroid = 0.75d0 * faceCentroids(f,:) + 0.25d0 * gC
		    
		    ! 累加体积和中心
		    vol = vol + subPyrVol
		    centroid = centroid + subPyrCentroid * subPyrVol
		end do

		centroid = centroid / vol
	end subroutine cellVolCentroid

	! 计算每个面从其所属单元中心指向面中心的向量!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function cellCentroidToFaceVec(faceCentroids, cellCentroids)
		implicit none
		real(kind=8), intent(in) :: faceCentroids(:,:), cellCentroids(3)
		real(kind=8) :: cellCentroidToFaceVec(size(faceCentroids,1), 3)
		integer :: nFaces, f
		
		nFaces = size(faceCentroids, 1)

		do f = 1, nFaces
			cellCentroidToFaceVec(f,:) = faceCentroids(f,:) - cellCentroids  ! 注意：cellCentroids 应是单个中心点？
		end do
	end function cellCentroidToFaceVec

	! ######################### 工具函数 ###########################

	! 返回网格的基本信息，包括单元数、面数、边界数量、边界面数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function unstructuredMeshInfo(mesh) result(info)
		implicit none
		type(Meshh), intent(in) :: mesh  ! 输入网格对象（与Julia的Mesh对应）
		integer(kind=8) :: info(4)       ! 返回数组：[nCells, nFaces, nBoundaries, nBdryFaces]
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces
		integer(kind=8) :: bdry, countt, i         ! 循环变量（边界索引）
		
		! 计算单元总数（网格单元的行数）
		nCells = size(mesh%cells, 1)
		! 计算面总数（网格面的行数）
		nFaces = size(mesh%faces, 1)
		! 计算边界总数（边界面的行数）
		nBoundaries = size(mesh%boundaryFaces, 1)		
		! 统计边界面总数（累加每个边界的面数量）
		! 正确计算边界面数（忽略填充的0）
		nBdryFaces = 0
		do bdry = 1, nBoundaries
		    ! 统计当前边界的非零元素个数
		    !i = 0
		    !do while (i < size(mesh%boundaryFaces, 2) .and. mesh%boundaryFaces(bdry, i+1) /= 0)
		    !    i = i + 1
		    !end do
   		    i = 0
    		    do
        		if (i+1 > size(mesh%boundaryFaces, 2)) exit
       			if (mesh%boundaryFaces(bdry, i+1) == 0) exit
        		i = i + 1
   		    end do
		    nBdryFaces = nBdryFaces + i
		end do
		! 赋值返回结果（保持与Julia相同的顺序）
		info = [nCells, nFaces, nBoundaries, nBdryFaces]
	end function unstructuredMeshInfo
	
    ! 核心函数：读取并处理OpenFOAM网格
    function OpenFOAMMesh(polyMeshPath) result(mesh)
        character(len=*), intent(in) :: polyMeshPath
        type(Meshh) :: mesh        
        type(MeshData) :: tempMesh
        integer :: nCells, nFaces, nBoundaries, nPoints
        integer :: f, c, b, i, j, k, pt_id, startF, endF, faceCount, max_faces_per_cell, j_max
        integer, allocatable :: cellFaceCount(:)
        real(kind=8), allocatable :: facePts(:,:)
        real(kind=8) :: maxCoords(3), minCoords(3)
        real(kind=8), allocatable :: cellPts(:,:), fCs(:,:), cell_fAVecs(:,:)
        real(kind=8), allocatable :: first5_refs(:,:)  ! 形状(5,3)!!!!!!!!!!!!
        integer, allocatable :: f_list(:) !!!!!!!!!!!!!!!!!!!!!!!!!!
        integer :: new_size
		real(kind=8), allocatable :: temp_cellPts(:,:) 
		
		allocate(first5_refs(5,3))
		first5_refs = 0.0d0 !!!!!!!!!!
        print*, 'Reading mesh: ', polyMeshPath
        tempMesh = readOpenFOAMMesh(polyMeshPath)    
        nPoints = size(tempMesh%points, 1)
        nFaces = size(tempMesh%owner)
        nCells = maxval(tempMesh%owner)
        nBoundaries = size(tempMesh%boundaryNames)     
        allocate(mesh%boundaryNames(nBoundaries))
        mesh%boundaryNames = tempMesh%boundaryNames
        allocate(mesh%faces(nFaces, 2))
        allocate(mesh%fAVecs(nFaces, 3), mesh%fCenters(nFaces, 3))
        
        allocate(cellFaceCount(nCells))
        cellFaceCount = 0
        do f = 1, nFaces
            cellFaceCount(tempMesh%owner(f)) = cellFaceCount(tempMesh%owner(f)) + 1
            if (allocated(tempMesh%neighbour) .and. f <= size(tempMesh%neighbour)) then
                cellFaceCount(tempMesh%neighbour(f)) = cellFaceCount(tempMesh%neighbour(f)) + 1
            end if
        end do
        max_faces_per_cell = maxval(cellFaceCount)
        allocate(mesh%cells(nCells, max_faces_per_cell))
        mesh%cells = 0
        
        allocate(mesh%cVols(nCells), mesh%cCenters(nCells, 3), mesh%cellSizes(nCells, 3))
        ! 计算每个面的面积向量和几何中心
        do f = 1, nFaces
            allocate(facePts(size(tempMesh%faces%faces(f)%points), 3))
            do i = 1, size(tempMesh%faces%faces(f)%points)
                pt_id = tempMesh%faces%faces(f)%points(i)
                facePts(i,:) = tempMesh%points(pt_id,:)
            end do
            ! 保存前5个面的第一个顶点坐标作为参考!!!!!!
			if (f <= 5) then
				first5_refs(f,:) = facePts(1,:)
			end if!!!!!!!!!!!!!!!!
			
			! 调用面计算子程序，传递前5个面的参考特征
			call faceAreaCentroid(facePts, first5_refs, mesh%fAVecs(f,:), mesh%fCenters(f,:))!!!!!!!
            !call faceAreaCentroid(facePts, mesh%fAVecs(f,:), mesh%fCenters(f,:))
            deallocate(facePts)
        end do
        deallocate(first5_refs)!!!!!!!!!!!!!
        
        cellFaceCount = 0
        do f = 1, nFaces
            c = tempMesh%owner(f)
            cellFaceCount(c) = cellFaceCount(c) + 1
            mesh%cells(c, cellFaceCount(c)) = f
            mesh%faces(f, 1) = c
            mesh%faces(f, 2) = -1
            
            if (allocated(tempMesh%neighbour) .and. f <= size(tempMesh%neighbour)) then
                c = tempMesh%neighbour(f)
                cellFaceCount(c) = cellFaceCount(c) + 1
                mesh%cells(c, cellFaceCount(c)) = f
                mesh%faces(f, 2) = c
            end if
        end do
        deallocate(cellFaceCount)
        ! 对每个单元，利用其所有面的信息计算体积和中心
        do c = 1, nCells
            faceCount = 0
            do i = 1, size(mesh%cells, 2)
                if (mesh%cells(c, i) == 0) exit
                faceCount = faceCount + 1
            end do
            
            allocate(fCs(faceCount, 3), cell_fAVecs(faceCount, 3))
            allocate(cellPts(0, 3))
            do i = 1, faceCount
                f = mesh%cells(c, i)
                fCs(i,:) = mesh%fCenters(f,:)
                cell_fAVecs(i,:) = mesh%fAVecs(f,:)
                
                do j = 1, size(tempMesh%faces%faces(f)%points)
                    pt_id = tempMesh%faces%faces(f)%points(j)
                    if (.not. isInArray(cellPts, tempMesh%points(pt_id,:))) then
						! 1. 计算新尺寸（原行数 + 1）
						new_size = size(cellPts, 1) + 1						
						! 2. 分配临时数组（新尺寸）
						allocate(temp_cellPts(new_size, 3))						
						! 3. 复制原有顶点（若存在）
						if (size(cellPts, 1) > 0) then
							temp_cellPts(1:new_size-1, :) = cellPts  ! 前n行复制原cellPts
						end if						
						! 4. 添加新顶点到最后一行（关键：按行存储）
						temp_cellPts(new_size, :) = tempMesh%points(pt_id,:)						
						! 5. 替换原cellPts
						call move_alloc(temp_cellPts, cellPts)
					end if
                end do
            end do
            ! 仅输出前两个单元的信息
			if (DEBUG_MESH_VERBOSE .and. c <= 2) then
			    allocate(f_list(faceCount))
				do i = 1, faceCount
					f_list(i) = mesh%cells(c, i)
				end do
			    ! 输出单元包含的面索引（前2个）
				print *, "包含的面索引（前2个）: ["
				do i = 1, min(6, faceCount)
					if (i > 1) then
						print "(A)", ", "  ! 元素间加逗号
					end if
					print "(I0)", f_list(i)  ! 逐个输出整数
				end do
				print *, "]"  ! 闭合括号
				! 对每个面，输出其包含的顶点ID（前3个）
				do i = 1, min(2, faceCount)
					f = f_list(i)  ! 面索引
					j_max = size(tempMesh%faces%faces(f)%points)
					print "(A, I0, A)", "  面", f, "的顶点ID（前3个）: ["  ! 先输出前缀
					! 逐个输出顶点ID（前3个）
					do j = 1, min(4, j_max)
						if (j > 1) then
							write(*, "(A)", advance="no") ", "  ! 元素间加逗号（不换行）
						end if
						write(*, "(I0)", advance="no") tempMesh%faces%faces(f)%points(j)  ! 输出单个顶点ID（不换行）
					end do
					print *, "]"  ! 闭合括号（换行）
				end do
				deallocate(f_list)
				print *, " "
				print *, "=== 单元 ", c, " 输入检验 ==="
				! 前4个顶点坐标
				print *, "全部顶点坐标："
				do k = 1, size(cellPts, 1)  ! 循环所有顶点（取消min限制）
					print "(A, I0, A, F0.6, A, F0.6, A, F0.6, A)", "  顶点", k, ": (", &
						  cellPts(k,1), ", ", cellPts(k,2), ", ", cellPts(k,3), ")"
				end do
				! 前4个面中心
				print *, "前4个面中心："
				do k = 1, min(6, faceCount)
				    print "(A, I0, A, F0.6, A, F0.6, A, F0.6, A)", "  面", k, ": (", &
				          fCs(k,1), ", ", fCs(k,2), ", ", fCs(k,3), ")"
				end do
				! 前4个面面积向量
				print *, "前4个面面积向量："
				do k = 1, min(6, faceCount)
				    print "(A, I0, A, F0.6, A, F0.6, A, F0.6, A)", "  面", k, "面积向量: (", &
				          cell_fAVecs(k,1), ", ", cell_fAVecs(k,2), ", ", cell_fAVecs(k,3), ")"
				end do
			end if
            call cellVolCentroid(cellPts, cell_fAVecs, fCs, mesh%cVols(c), mesh%cCenters(c,:))
            deallocate(cellPts, fCs, cell_fAVecs)
        end do
        ! 生成 boundaryFaces 数组
		! 恢复原始分配方式
		! 生成 boundaryFaces 数组时，先全部初始化为0
		allocate(mesh%boundaryFaces(nBoundaries, maxval(tempMesh%boundaryNumFaces)))
		mesh%boundaryFaces = 0  ! 关键：将所有元素初始化为0

		do b = 1, nBoundaries
			startF = tempMesh%boundaryStartFaces(b)
			endF = startF + tempMesh%boundaryNumFaces(b) - 1
			do f = startF, endF
				mesh%boundaryFaces(b, f - startF + 1) = f  ! 填充有效面索引（无效位置保持0）
			end do
		end do
        ! 计算每个单元在 x、y、z 方向上的尺寸（用边界包围盒法）
        do c = 1, nCells
            !maxCoords = -huge(1.0d0)
            !minCoords = huge(1.0d0)
	    maxCoords = -1.0d30
            minCoords = 1.0d30
            do i = 1, size(mesh%cells, 2)
                f = mesh%cells(c, i)
                if (f == 0) exit
                do j = 1, size(tempMesh%faces%faces(f)%points)
                    pt_id = tempMesh%faces%faces(f)%points(j)
                    do k = 1, 3
                        maxCoords(k) = max(maxCoords(k), tempMesh%points(pt_id, k))
                        minCoords(k) = min(minCoords(k), tempMesh%points(pt_id, k))
                    end do
                end do
            end do
            mesh%cellSizes(c,:) = maxCoords - minCoords
        end do
        
        call deallocateMeshData(tempMesh)
    end function OpenFOAMMesh

    ! 读取OpenFOAM网格主函数
    function readOpenFOAMMesh(polyMeshPath) result(mesh)
        character(len=*), intent(in) :: polyMeshPath
        type(MeshData) :: mesh
        character(len=256) :: pointsFilePath, facesFilePath, ownerFilePath, neighbourFilePath, boundaryFilePath
        
        pointsFilePath = trim(polyMeshPath) // "/points"
        facesFilePath = trim(polyMeshPath) // "/faces"
        ownerFilePath = trim(polyMeshPath) // "/owner"
        neighbourFilePath = trim(polyMeshPath) // "/neighbour"
        boundaryFilePath = trim(polyMeshPath) // "/boundary"
        
        mesh%points = readOFPointsFile(pointsFilePath)
        mesh%faces = readOFFacesFile(facesFilePath)
        mesh%owner = readOFOwnerFile(ownerFilePath)
        mesh%neighbour = readOFNeighbourFile(neighbourFilePath)
        call readOFBoundaryFile(boundaryFilePath, mesh%boundaryNames, mesh%boundaryNumFaces, mesh%boundaryStartFaces)
    end function readOpenFOAMMesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 主函数：从OpenFOAM网格路径读取并处理单元点索引
    subroutine OpenFOAMMesh_findCellPts(polyMeshPath, pointLocations, cells)
		character(len=*), intent(in) :: polyMeshPath
		real(kind=8), allocatable, intent(out) :: pointLocations(:,:)
		type(CelllArray), intent(out) :: cells
		
		! 主程序变量（包含pointIndicesByFace，避免全局变量）
		type(MeshData) :: meshDataa  ! 网格数据
		integer :: nCells, nFaces, nBoundaries
		integer :: i, j, fi, nFacesCell, quadFaceCount
		integer :: maxFacePts, npts_face
		! 核心：点索引数据（对应Julia的pointIndicesByFace，在主程序内定义）
		integer, allocatable :: pointIndicesByFace(:,:)  ! (面索引, 点序号)
		integer, allocatable :: faceNpts(:)
		integer :: ierr, len_common
		integer, allocatable :: commonPoints(:)		! 1. 读取网格数据
		meshDataa = readOpenFOAMMesh(polyMeshPath)
		
		! 2. 提取pointIndicesByFace（关键：在主程序内定义，供内置函数使用）
		! 假设从meshDataa中提取点索引（根据实际数据结构调整）
		nFaces = size(meshDataa%faces%faces)
		maxFacePts = 0
		do i = 1, nFaces
			npts_face = size(meshDataa%faces%faces(i)%points)
			if (npts_face > maxFacePts) maxFacePts = npts_face
		end do
		if (maxFacePts <= 0) then
			print *, "Error: 无效面顶点数（maxFacePts<=0）"
			stop
		end if
		allocate(pointIndicesByFace(nFaces, maxFacePts))
		allocate(faceNpts(nFaces))
		pointIndicesByFace = 0
		faceNpts = 0
		do i = 1, nFaces
			npts_face = size(meshDataa%faces%faces(i)%points)
			faceNpts(i) = npts_face
			if (npts_face > 0) then
				pointIndicesByFace(i, 1:npts_face) = meshDataa%faces%faces(i)%points
			end if
		end do
		
		! 3. 初始化单元数组
		nCells = maxval(meshDataa%owner)
		allocate(cells%cells(nCells))
		do i = 1, nCells
		    allocate(cells%cells(i)%faceIndices(0))
		    allocate(cells%cells(i)%pointIndices(0))
		end do
		
		! 4. 添加面索引（调用外部的addCellFaceIndices）
		call addCellFaceIndices(meshDataa%owner, cells, nCells)
		call addCellFaceIndices(meshDataa%neighbour, cells, nCells)
		
		! 5. 处理单元点索引（调用内置函数）
		do i = 1, nCells
		    nFacesCell = size(cells%cells(i)%faceIndices)
		    select case(nFacesCell)
		        case(4)
		            call populatePointIndices_Tet(cells%cells(i))
		        case(5)
		            quadFaceCount = 0
		            do j = 1, nFacesCell
		                fi = cells%cells(i)%faceIndices(j)
		                if (faceNpts(fi) == 4) quadFaceCount = quadFaceCount + 1
		            end do
		            if (quadFaceCount == 1) then
		                call populatePointIndices_Pyramid(cells%cells(i))
		            else if (quadFaceCount == 3) then
		                call populatePointIndices_Wedge(cells%cells(i))
		            else
		                print *, "Error: 无法识别的单元类型（单元", i, "）"
		                stop
		            end if
		        case(6)
		        ! 在case(6)中，单元1处理前添加调试
					if (i == 1) then
						print *, "调试：单元1的6个面索引：", cells%cells(i)%faceIndices
						do j = 1, 6
							fi = cells%cells(i)%faceIndices(j)
							print *, "面", fi, "的顶点：", pointIndicesByFace(fi, :)
						end do
						! 检查面1与其他面的公共点（基准面f1应与4个侧面相邻）
						print *, "面1（f1）与其他面的公共点数量："
						do j = 2, 6
							fi = cells%cells(i)%faceIndices(j)
							call get_intersection(pointIndicesByFace(1, :), pointIndicesByFace(fi, :), commonPoints, len_common)
							print *, "面1与面", fi, "：", len_common, "个公共点"
							deallocate(commonPoints)
						end do
					end if
		        ! 新增：检查单元所有面是否为4个顶点（六面体必需）
					do j = 1, size(cells%cells(i)%faceIndices)
						fi = cells%cells(i)%faceIndices(j)
						if (faceNpts(fi) /= 4) then
							print *, "Error: 六面体单元", i, "的面", fi, "顶点数不为4（实际为", faceNpts(fi), "）"
							stop
						end if
					end do
					! 新增：检查面数是否确实为6（防止数据错误）
					if (size(cells%cells(i)%faceIndices) /= 6) then
						print *, "Error: 单元", i, "声明为六面体但实际面数为", size(cells%cells(i)%faceIndices)
						stop
					end if
					call populatePointIndices_Hex(cells%cells(i), ierr)
					if (ierr /= 0) then
						! 增强错误信息：输出单元1的面索引，方便排查
						print *, "Error: 六面体单元", i, "顶点排序失败，错误代码", ierr
						print *, "单元", i, "的面索引列表：", cells%cells(i)%faceIndices
						stop
					end if
		        case default
		            print *, "Error: 不支持的单元面数（单元", i, "，面数", nFacesCell, "）"
		            stop
		    end select
		end do
		
		! 6. 输出点坐标
		if (allocated(pointLocations)) deallocate(pointLocations)
		allocate(pointLocations(size(meshDataa%points,1), 3))
		pointLocations = meshDataa%points


	! 内置所有点索引处理函数（通过contains关联，共享主程序变量）
	contains
		! 四面体单元点索引组织
		subroutine populatePointIndices_Tet(cell)
		    type(Celll), intent(inout) :: cell
		    ! 直接调用内置的addAllNewPoints（无需传递参数，共享pointIndicesByFace）
		    !call addAllNewPoints(cell, cell%faceIndices(1))
		    !call addAllNewPoints(cell, cell%faceIndices(2))
		    print *, "Aa"
		end subroutine populatePointIndices_Tet

		! 金字塔单元点索引组织
		subroutine populatePointIndices_Pyramid(cell)
		    type(Celll), intent(inout) :: cell
		    integer :: i, quadFaceIdx, otherFaceIndex
		    quadFaceIdx = -1
		    do i = 1, size(cell%faceIndices)
		        ! 使用主程序的pointIndicesByFace判断面类型
		        if (faceNpts(cell%faceIndices(i)) == 4) then
		            quadFaceIdx = cell%faceIndices(i)
		            exit
		        end if
		    end do
		    if (quadFaceIdx /= -1) then
		        !call addAllNewPoints(cell, quadFaceIdx)
		        ! 选择另一个面（与Julia逻辑一致）
		        if (quadFaceIdx /= cell%faceIndices(1)) then
		            otherFaceIndex = cell%faceIndices(1)
		        else
		            otherFaceIndex = cell%faceIndices(2)
		        end if
		        !call addAllNewPoints(cell, otherFaceIndex)
		    end if
		    print *, "Bb"
		end subroutine populatePointIndices_Pyramid		
		! 棱柱（楔形）单元：两个三角形底面，其余三个是连接面		
		subroutine populatePointIndices_Wedge(cell)
			type(Celll), intent(inout) :: cell
			integer :: i, j, t1, t2, sideCount, p, score
			integer :: triFaces(2), sideFaces(3)
			integer :: basePts(3), topPts(3), topPerm(3), chosenTop(3)
			logical :: ok
			integer :: perms(6,3)
			triFaces = 0
			sideFaces = 0
			sideCount = 0
			p = 0
			do i = 1, size(cell%faceIndices)
				j = cell%faceIndices(i)
				if (faceNpts(j) == 3) then
					p = p + 1
					if (p <= 2) triFaces(p) = j
				else if (faceNpts(j) == 4) then
					sideCount = sideCount + 1
					if (sideCount <= 3) sideFaces(sideCount) = j
				end if
			end do
			if (p /= 2 .or. sideCount /= 3) then
				print *, "Error: 棱柱单元识别失败（需2个三角面+3个四边面）"
				stop
			end if
			t1 = triFaces(1)
			t2 = triFaces(2)
			basePts = pointIndicesByFace(t1, 1:3)
			topPts = pointIndicesByFace(t2, 1:3)
			perms(1,:) = [1,2,3]
			perms(2,:) = [1,3,2]
			perms(3,:) = [2,1,3]
			perms(4,:) = [2,3,1]
			perms(5,:) = [3,1,2]
			perms(6,:) = [3,2,1]
			ok = .false.
			do i = 1, 6
				topPerm = [topPts(perms(i,1)), topPts(perms(i,2)), topPts(perms(i,3))]
				score = wedge_mapping_score(basePts, topPerm, sideFaces)
				if (score == 3) then
					chosenTop = topPerm
					ok = .true.
					exit
				end if
			end do
			if (.not. ok) then
				print *, "Error: 棱柱顶点配对失败"
				stop
			end if
			if (allocated(cell%pointIndices)) deallocate(cell%pointIndices)
			allocate(cell%pointIndices(6))
			cell%pointIndices(1:3) = basePts
			cell%pointIndices(4:6) = chosenTop
		end subroutine populatePointIndices_Wedge

		integer function wedge_mapping_score(basePts, topPerm, sideFaces) result(score)
			integer, intent(in) :: basePts(3), topPerm(3)
			integer, intent(in) :: sideFaces(3)
			integer :: e, f, a, b, c, d
			logical :: hasA, hasB, hasC, hasD
			score = 0
			do e = 1, 3
				a = basePts(e)
				b = basePts(mod(e,3)+1)
				c = topPerm(mod(e,3)+1)
				d = topPerm(e)
				do f = 1, 3
					hasA = any(pointIndicesByFace(sideFaces(f), 1:faceNpts(sideFaces(f))) == a)
					hasB = any(pointIndicesByFace(sideFaces(f), 1:faceNpts(sideFaces(f))) == b)
					hasC = any(pointIndicesByFace(sideFaces(f), 1:faceNpts(sideFaces(f))) == c)
					hasD = any(pointIndicesByFace(sideFaces(f), 1:faceNpts(sideFaces(f))) == d)
					if (hasA .and. hasB .and. hasC .and. hasD) then
						score = score + 1
						exit
					end if
				end do
			end do
		end function wedge_mapping_score
	    ! 六面体（立方体）单元：按照 VTK 的顺序组织 8 个顶点
        ! 原定义错误：嵌套contains且参数多余，修改为内部子程序（直接访问主程序变量）
		subroutine populatePointIndices_Hex(cell, ierr)
			type(Celll), intent(inout) :: cell
			integer, intent(out) :: ierr
			integer :: f1, f3, lastFace, noResult
			integer, allocatable :: f1Points(:), fiPoints(:)
			integer, allocatable :: unusedFaces(:)
			integer :: i, fi, lastFaceIdx
			logical :: found

			ierr = 0
			found = .false.

			! 1. 复制所有面索引到unusedFaces
			allocate(unusedFaces(size(cell%faceIndices)))
			unusedFaces = cell%faceIndices
			!!print *, "初始面索引: ", unusedFaces

			! 2. 取第一个面作为f1（修复核心错误）
			f1 = unusedFaces(1)
			call remove_element(unusedFaces, 1)  ! 移除第一个元素
			!!print *, "基准面f1 = 面", f1, "，剩余面: ", unusedFaces

			! 3. 获取f1的顶点
			allocate(f1Points(4))
			f1Points = pointIndicesByFace(f1, :)
			!!print *, "f1的顶点: ", f1Points

			! 4. 找对面（与f1无公共点的面，应为面2）
			do i = 1, size(unusedFaces)
				fi = unusedFaces(i)
				allocate(fiPoints(4))
				fiPoints = pointIndicesByFace(fi, :)
				if (disjoint(f1Points, fiPoints)) then
				    !!print *, "找到对面: 面", fi
				    call remove_element(unusedFaces, i)
				    found = .true.
				    deallocate(fiPoints)
				    exit
				end if
				deallocate(fiPoints)
			end do
			if (.not. found) then
				ierr = 1
				!!print *, "错误代码1: 未找到对面"
				return
			end if

			! 检查剩余面是否为4个（关键验证）
			if (size(unusedFaces) /= 4) then
				print *, "错误：剩余面数应为4，实际为", size(unusedFaces)
				ierr = 1
				return
			end if
			!!print *, "删除对面后剩余面（4个侧面）: ", unusedFaces

			! 5. 初始化顶点数组
			if (allocated(cell%pointIndices)) deallocate(cell%pointIndices)
			allocate(cell%pointIndices(8))
			cell%pointIndices(1:4) = f1Points(1:4)
			cell%pointIndices(5:8) = [0,0,0,0]
			!!print *, "初始顶点索引: ", cell%pointIndices

			! 6. 取第一个剩余面作为f3
			f3 = unusedFaces(1)
			call remove_element(unusedFaces, 1)
			!!print *, "处理面f3 = 面", f3, "，剩余面: ", unusedFaces

			! 7. 第一次调用addEdges
			lastFace = addEdges(cell%pointIndices, pointIndicesByFace(f3, :), f1Points, unusedFaces)
			!!print *, "第一次addEdges后顶点: ", cell%pointIndices
			!!print *, "addEdges返回的lastFace: ", lastFace

			! 8. 检查lastFace有效性
			if (lastFace < 1 .or. lastFace > size(unusedFaces)) then
				ierr = 3
				print *, "错误代码3: lastFace无效"
				return
			end if

			! 9. 处理最后一个面
			lastFaceIdx = unusedFaces(lastFace)
			call remove_element(unusedFaces, lastFace)
			!!print *, "处理最后一个面: 面", lastFaceIdx

			! 10. 第二次调用addEdges
			noResult = addEdges(cell%pointIndices, pointIndicesByFace(lastFaceIdx, :), f1Points, unusedFaces)
			!!print *, "第二次addEdges后顶点: ", cell%pointIndices

			! 11. 最终检查
			if (noResult /= -1) then
				ierr = 4
				return
			end if

			deallocate(f1Points, unusedFaces)
		end subroutine populatePointIndices_Hex
		! 模拟Julia的pop!：移除数组最后一个元素，返回被移除的值
		subroutine pop_last(arr, val)
			integer, allocatable, intent(inout) :: arr(:)
			integer, intent(out) :: val
			integer, allocatable :: temp(:)
			integer :: n

			n = size(arr)
			if (n == 0) error stop "pop_last: 数组为空"
			val = arr(n)  ! 取最后一个元素
			if (n == 1) then
				deallocate(arr)  ! 若只有一个元素，直接释放
			else
				allocate(temp(n-1))
				temp(1:n-1) = arr(1:n-1)  ! 保留前n-1个元素
				call move_alloc(temp, arr)
			end if
		end subroutine pop_last

		! 模拟Julia的deleteat!：删除数组第i个元素，后续元素前移
		subroutine deleteat(arr, i)
			integer, allocatable, intent(inout) :: arr(:)
			integer, intent(in) :: i
			integer, allocatable :: temp(:)
			integer :: n, j, k

			n = size(arr)
			if (i < 1 .or. i > n) return  ! 索引无效则不操作
			allocate(temp(n-1))
			k = 1
			do j = 1, n
				if (j /= i) then
				    temp(k) = arr(j)
				    k = k + 1
				end if
			end do
			call move_alloc(temp, arr)
		end subroutine deleteat
		function addEdges(cellPoints, facePoints, endFacePoints, unusedFaces) result(oppositeFaceIndex)
			integer, intent(inout) :: cellPoints(:)
			integer, intent(in) :: facePoints(:), endFacePoints(:)
			integer, allocatable, intent(in) :: unusedFaces(:)
			integer :: oppositeFaceIndex
			integer :: i, fi, p1, p2, pos, len_common
			integer, allocatable :: fiPoints(:), commonPoints(:)
			logical :: p1_in_end
			integer :: processedCount, expectedCount
			logical :: has_opposite  ! 新增：标记是否存在对面

			oppositeFaceIndex = -1
			processedCount = 0
			has_opposite = .false.  ! 初始化为无对面

			do i = 1, size(unusedFaces)
				fi = unusedFaces(i)
				fiPoints = pointIndicesByFace(fi, :)
				call get_intersection(facePoints, fiPoints, commonPoints, len_common)
				!print *, "addEdges: 面", fi, "与当前面公共点数量：", len_common

				if (len_common == 2) then
				    processedCount = processedCount + 1
				    ! 顶点填充逻辑（不变）
				    p1_in_end = any(endFacePoints == commonPoints(1))
				    if (p1_in_end) then
				        p1 = commonPoints(1)
				        p2 = commonPoints(2)
				    else
				        p1 = commonPoints(2)
				        p2 = commonPoints(1)
				    end if
				    pos = findloc(endFacePoints, p1, dim=1)
				    if (pos > 0) cellPoints(pos + 4) = p2

				else if (len_common == 0) then
				    oppositeFaceIndex = i
				    has_opposite = .true.  ! 标记存在对面
				end if
				deallocate(commonPoints)
			end do

			! 关键修改：根据是否有对面计算预期值
			if (has_opposite) then
				expectedCount = size(unusedFaces) - 1  ! 有对面：总面数-1
			else
				expectedCount = size(unusedFaces)      ! 无对面：总面数
			end if

			! 检查是否匹配
			if (processedCount /= expectedCount) then
				!print *, "Error: addEdges处理面数异常，预期", expectedCount, "实际", processedCount
				oppositeFaceIndex = -3
			end if
		end function addEdges
			
		! 辅助函数：计算两个数组的交集
		subroutine get_intersection(a, b, intersection, len)
			integer, intent(in) :: a(:), b(:)
			integer, allocatable, intent(out) :: intersection(:)
			integer, intent(out) :: len
			integer :: i, j, count
			integer, allocatable :: temp(:)  ! 临时数组

			count = 0
			allocate(intersection(min(size(a), size(b))))

			do i = 1, size(a)
				do j = 1, size(b)
				    if (a(i) == b(j)) then
				        count = count + 1
				        intersection(count) = a(i)
				        exit
				    end if
				end do
			end do

			len = count
			if (count < size(intersection)) then
				! 修复：先复制到临时可分配数组，再移动
				allocate(temp(count))
				temp = intersection(1:count)
				call move_alloc(temp, intersection)  ! 现在第一个参数是可分配数组
			end if
		end subroutine get_intersection

		! 辅助函数：检查元素是否在数组中
		logical function is_in(val, arr)
			integer, intent(in) :: val, arr(:)
			integer :: i
			is_in = .false.
			do i = 1, size(arr)
				if (arr(i) == val) then
				    is_in = .true.
				    return
				end if
			end do
		end function is_in

		! 辅助函数：查找元素在数组中的位置
		integer function find_index(val, arr)
			integer, intent(in) :: val, arr(:)
			integer :: i
			find_index = -1
			do i = 1, size(arr)
				if (arr(i) == val) then
				    find_index = i
				    return
				end if
			end do
		end function find_index
		! 辅助函数：检查两个点集是否不相交
		! 辅助函数：检查两个点集是否不相交（同级）
		logical function disjoint(a, b)
			integer, intent(in) :: a(:), b(:)
			integer :: i, j
			disjoint = .true.
			do i = 1, size(a)
			    do j = 1, size(b)
			        if (a(i) == b(j)) then
			            disjoint = .false.
			            return
			        end if
			    end do
			end do
		end function disjoint
		! 辅助子程序：从数组中移除指定索引的元素（同级）
		subroutine remove_element(arr, idx)
		    integer, allocatable, intent(inout) :: arr(:)
		    integer, intent(in) :: idx
		    integer, allocatable :: temp(:)
		    integer :: i, n

		    n = size(arr)
		    if (idx < 1 .or. idx > n) return

		    allocate(temp(n-1))
		    do i = 1, idx-1
		        temp(i) = arr(i)
		    end do
		    do i = idx+1, n
		        temp(i-1) = arr(i)
		    end do
		    call move_alloc(temp, arr)
		end subroutine remove_element		
		subroutine swap(a, b)
			integer, intent(inout) :: a, b
			integer :: temp
			temp = a; a = b; b = temp
		end subroutine swap
		subroutine addCellFaceIndices(adjacentCells, cells, nCells)
			integer, intent(in) :: adjacentCells(:)      ! 面所属单元索引
			integer, intent(in) :: nCells                ! 总单元数
			type(CelllArray), intent(inout) :: cells    ! 单元数组
			
			! 局部变量：统计用
			integer :: f, cellIdx, len, i, num_faces, j
			integer :: count_6_faces                    ! 6个面的单元数量
			integer, allocatable :: other_counts(:), face_nums(:)  ! 其他面数统计
			integer, allocatable :: temp(:)
			logical :: is_duplicate, found
			
			! 初始化统计变量
			count_6_faces = 0
			allocate(other_counts(0), face_nums(0))
			
			! 1. 为单元添加面索引（严格对应Julia逻辑）
			do f = 1, size(adjacentCells)
				cellIdx = adjacentCells(f)
				
				! 检查单元索引有效性（防止越界）
				if (cellIdx < 1 .or. cellIdx > nCells) cycle
				
				! 直接添加面索引（不检查重复，与Julia一致）
				len = size(cells%cells(cellIdx)%faceIndices)
				allocate(temp(len + 1))
				if (len > 0) temp(1:len) = cells%cells(cellIdx)%faceIndices
				temp(len + 1) = f
				call move_alloc(temp, cells%cells(cellIdx)%faceIndices)
			end do
			
			! 2. 统计面数（与Julia逻辑完全一致）
			do cellIdx = 1, nCells
				num_faces = size(cells%cells(cellIdx)%faceIndices)
				
				! 统计6个面的单元
				if (num_faces == 6) then
				    count_6_faces = count_6_faces + 1
				else
				    ! 统计其他面数（动态数组模拟字典）
				    found = .false.
				    do i = 1, size(face_nums)
				        if (face_nums(i) == num_faces) then
				            other_counts(i) = other_counts(i) + 1
				            found = .true.
				            exit
				        end if
				    end do
				    
				    ! 新增面数类型
				    if (.not. found) then
				        allocate(temp(size(face_nums) + 1))
				        if (size(face_nums) > 0) temp(1:size(face_nums)) = face_nums
				        temp(size(face_nums) + 1) = num_faces
				        call move_alloc(temp, face_nums)
				        
				        allocate(temp(size(other_counts) + 1))
				        if (size(other_counts) > 0) temp(1:size(other_counts)) = other_counts
				        temp(size(other_counts) + 1) = 1
				        call move_alloc(temp, other_counts)
				    end if
				end if
			end do
			
			! 3. 输出统计结果（与Julia格式完全一致）
			print *, "包含6个面的单元数量: ", count_6_faces
			if (size(face_nums) > 0) then
				print *, "其他面数的单元统计:"
				! 排序（按面数升序）
				do i = 1, size(face_nums) - 1
				    do j = i + 1, size(face_nums)
				        if (face_nums(j) < face_nums(i)) then
				            call swap(face_nums(i), face_nums(j))
				            call swap(other_counts(i), other_counts(j))
				        end if
				    end do
				end do
				! 打印排序后的结果
				do i = 1, size(face_nums)
				    write(*, '(A, I0, A, I0, A)') "  ", face_nums(i), " 个面的单元: ", other_counts(i), " 个"
				end do
			else
				print *, "所有单元都包含6个面"
			end if
			
			! 释放临时数组
			deallocate(other_counts, face_nums)
		end subroutine addCellFaceIndices		
	end	subroutine OpenFOAMMesh_findCellPts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 读取点文件
    function readOFPointsFile(filePath) result(points)
        character(len=*), intent(in) :: filePath
        real(kind=8), allocatable :: points(:,:)
        character(len=256), allocatable :: lines(:)
        integer :: startLine, pCount, i, nLines, iostat, file_unit
        character(len=256) :: line, bracketsRemoved
        real(kind=8) :: coords(3)
        
        file_unit = get_free_unit()
        open(unit=file_unit, file=filePath, status='old', action='read', iostat=iostat)
        if (iostat /= 0) then
            print *, "Error: Can't open points file: ", trim(filePath)
            allocate(points(0,3))
            return
        end if
        
        nLines = 0
        do
            read(file_unit, '(a)', iostat=iostat)
            if (iostat /= 0) exit
            nLines = nLines + 1
        end do
        rewind(file_unit)
        
        allocate(lines(nLines))
        do i = 1, nLines
            read(file_unit, '(a)') lines(i)
        end do
        close(file_unit)
        
        call OFFile_FindNItems(lines, startLine, pCount)
        if (pCount <= 0) then
            print *, "Error: Invalid point count in ", trim(filePath)
            deallocate(lines)
            allocate(points(0,3))
            return
        end if
        
        allocate(points(pCount, 3))
        do i = 1, pCount
            line = trim(lines(startLine + i - 1))
            bracketsRemoved = line(2:len_trim(line)-1)
            read(bracketsRemoved, *) coords
            points(i,:) = coords
        end do
        deallocate(lines)
    end function readOFPointsFile

    ! 读取面文件
    function readOFFacesFile(filePath) result(faces)
        use iso_fortran_env, only: int64
        implicit none
        character(len=*), intent(in) :: filePath
        type(FaceArray) :: faces
        type(FaceType), allocatable :: tmp_faces(:)
        character(len=1000), allocatable :: lines(:)
        integer :: startLine, fCount, i, nLines, iostat, file_unit, bracketL, bracketR
        character(len=1000) :: line
        integer(int64), allocatable :: pts(:)
        integer :: pt_count
        
        file_unit = get_free_unit()
        open(unit=file_unit, file=filePath, status='old', action='read', iostat=iostat)
        if (iostat /= 0) then
            print *, "Error: Can't open faces file: ", trim(filePath)
            allocate(faces%faces(0))
            return
        end if
        
        nLines = 0
        do
            read(file_unit, '(a)', iostat=iostat)
            if (iostat /= 0) exit
            nLines = nLines + 1
        end do
        rewind(file_unit)
        
        allocate(lines(nLines))
        do i = 1, nLines
            read(file_unit, '(a)') lines(i)
        end do
        close(file_unit)
        
        call OFFile_FindNItems(lines, startLine, fCount)
        if (fCount <= 0) then
            print *, "Error: Invalid face count in ", trim(filePath)
            deallocate(lines)
            allocate(faces%faces(0))
            return
        end if
        
        allocate(tmp_faces(fCount))
        do i = 1, fCount
            line = trim(lines(startLine + i - 1))
            bracketL = index(line, '(')
            bracketR = index(line, ')')
            read(line(1:bracketL-1), *) pt_count
            allocate(pts(pt_count))
            read(line(bracketL+1:bracketR-1), *) pts
            pts = pts + 1  ! 0基转1基
            allocate(tmp_faces(i)%points(pt_count))
            tmp_faces(i)%points = pts
            deallocate(pts)
        end do
        
        allocate(faces%faces(fCount))
        faces%faces = tmp_faces
        deallocate(lines, tmp_faces)
    end function readOFFacesFile

    ! 读取owner文件
    function readOFOwnerFile(filePath) result(owner)
        implicit none
        character(len=*), intent(in) :: filePath
        integer, allocatable :: owner(:)
        character(len=256), allocatable :: lines(:)
        integer :: startLine, oCount, i, nLines, iostat, file_unit
        integer(kind=8) :: val
        
        file_unit = get_free_unit()
        open(unit=file_unit, file=filePath, status='old', action='read', iostat=iostat)
        if (iostat /= 0) then
            print *, "Error: Can't open owner file: ", trim(filePath)
            allocate(owner(0))
            return
        end if
        
        nLines = 0
        do
            read(file_unit, '(a)', iostat=iostat)
            if (iostat /= 0) exit
            nLines = nLines + 1
        end do
        rewind(file_unit)
        
        allocate(lines(nLines))
        do i = 1, nLines
            read(file_unit, '(a)') lines(i)
        end do
        close(file_unit)
        
        call OFFile_FindNItems(lines, startLine, oCount)
        if (oCount <= 0) then
            print *, "Error: Invalid owner count in ", trim(filePath)
            deallocate(lines)
            allocate(owner(0))
            return
        end if
        
        allocate(owner(oCount))
        do i = 1, oCount
            read(lines(startLine + i - 1), *) val
            owner(i) = int(val) + 1
        end do
        deallocate(lines)
    end function readOFOwnerFile

    ! 读取neighbour文件
    function readOFNeighbourFile(filePath) result(neighbour)
		implicit none
		character(len=*), intent(in) :: filePath  ! 输入文件路径
		integer, allocatable :: neighbour(:)      ! 输出邻居单元索引数组
		character(len=256), allocatable :: lines(:)  ! 存储文件所有行
		integer :: startLine, nCount  ! 数据起始行、数据总数
		integer :: i, nLines, iostat, file_unit  ! 循环变量、总行数、I/O状态、文件单元
		integer(kind=8) :: val  ! 临时存储读取的大整数（兼容OpenFOAM格式）

		! 1. 初始化：不预先分配数组（与owner处理一致）
		if (allocated(neighbour)) deallocate(neighbour)

		! 2. 打开文件（与owner逻辑完全一致）
		file_unit = get_free_unit()
		open(unit=file_unit, file=filePath, status='old', action='read', iostat=iostat)
		if (iostat /= 0) then
		    print *, "Error: Can't open neighbour file: ", trim(filePath)  ! 错误信息与owner格式统一
		    allocate(neighbour(0))  ! 错误时返回空数组
		    return
		end if

		! 3. 读取文件总行数（与owner逻辑一致）
		nLines = 0
		do
		    read(file_unit, '(a)', iostat=iostat)  ! 逐行读取计数
		    if (iostat /= 0) exit
		    nLines = nLines + 1
		end do
		rewind(file_unit)  ! 重置文件指针

		! 4. 存储所有行内容（与owner逻辑一致）
		allocate(lines(nLines))  ! 分配存储所有行的数组
		do i = 1, nLines
		    read(file_unit, '(a)') lines(i)  ! 逐行读取内容
		end do
		close(file_unit)  ! 关闭文件（后续使用lines处理）

		! 5. 查找数据起始行和数量（核心：与owner调用同一函数）
		call OFFile_FindNItems(lines, startLine, nCount)

		! 6. 检查数据数量有效性（与owner逻辑一致）
		if (nCount <= 0) then
		    print *, "Error: Invalid neighbour count in ", trim(filePath)  ! 错误信息格式统一
		    deallocate(lines)
		    allocate(neighbour(0))  ! 无效数量时返回空数组
		    return
		end if

		! 7. 分配邻居数组并读取数据（与owner完全一致）
		allocate(neighbour(nCount))  ! 确认数量有效后再分配
		do i = 1, nCount
		    ! 严格按照owner的读取方式：startLine + i - 1索引
		    read(lines(startLine + i - 1), *) val  ! 读取原始值（OpenFOAM格式）
		    neighbour(i) = int(val) + 1  ! 转换为Fortran索引（+1逻辑与owner一致）
		end do

		! 8. 释放临时资源（与owner一致）
		deallocate(lines)

	end function readOFNeighbourFile
    ! 读取边界文件
    subroutine readOFBoundaryFile(filePath, boundaryNames, boundaryNumFaces, boundaryStartFaces)
        implicit none
        character(len=*), intent(in) :: filePath
        character(len=100), allocatable, intent(out) :: boundaryNames(:)
        integer, allocatable, intent(out) :: boundaryNumFaces(:), boundaryStartFaces(:)
        integer :: i, nLines, iostat, startLine, bCount, file_unit
        integer :: bNameLine, bNFacesLine, bStartFaceLine, pos
        character(len=256), allocatable :: bLines(:)
        character(len=256) :: lineStr, dummy
        
        file_unit = get_free_unit()
        open(unit=file_unit, file=filePath, status='old', action='read', iostat=iostat)
        if (iostat /= 0) then
            allocate(boundaryNames(0), boundaryNumFaces(0), boundaryStartFaces(0))
            return
        end if
        
        nLines = 0
        do
            read(file_unit, '(a)', iostat=iostat)
            if (iostat /= 0) exit
            nLines = nLines + 1
        end do
        rewind(file_unit)
        
        allocate(bLines(nLines))
        do i = 1, nLines
            read(file_unit, '(a)') bLines(i)
        end do
        close(file_unit)
        
        call OFFile_FindNItems(bLines, startLine, bCount)
        if (bCount <= 0) then
            deallocate(bLines)
            allocate(boundaryNames(0), boundaryNumFaces(0), boundaryStartFaces(0))
            return
        end if
        
        allocate(boundaryNames(bCount), boundaryNumFaces(bCount), boundaryStartFaces(bCount))
        do i = 1, bCount
            bNameLine = findInLines("{", bLines, startLine) - 1
            boundaryNames(i) = trim(adjustl(bLines(bNameLine)))
            
            bNFacesLine = findInLines("nFaces", bLines, startLine)
            pos = index(bLines(bNFacesLine), "nFaces") + 6
            lineStr = bLines(bNFacesLine)(pos:)
            call removeSemicolon(lineStr)
            read(lineStr, *) boundaryNumFaces(i)
            
            bStartFaceLine = findInLines("startFace", bLines, startLine)
            pos = index(bLines(bStartFaceLine), "startFace") + 9
            lineStr = bLines(bStartFaceLine)
            call removeSemicolon(lineStr)
            read(lineStr, *) dummy, boundaryStartFaces(i)
			boundaryStartFaces(i) = boundaryStartFaces(i) + 1
			!boundaryEndFaces(i) = boundaryStartFaces(i) + boundaryNumFaces(i) - 1
            
            startLine = findInLines("}", bLines, startLine) + 1
        end do
        deallocate(bLines)
    end subroutine readOFBoundaryFile

    ! 辅助函数：去除字符串中的分号（兼容OpenFOAM格式）
    subroutine removeSemicolon(str)
        character(len=*), intent(inout) :: str
        integer :: i
        do i = 1, len_trim(str)
            if (str(i:i) == ';') str(i:i) = ' '
        end do
    end subroutine removeSemicolon

    ! 辅助函数：查找数据数量和起始行
    subroutine OFFile_FindNItems(fileLines, startLine, itemCount)
        implicit none
        character(len=*), intent(in) :: fileLines(:)
        integer, intent(out) :: startLine, itemCount
        integer :: i, pos, iostat
        character(len=256) :: line
        
        itemCount = 0
        startLine = 0
        do i = 1, size(fileLines)
            line = trim(fileLines(i))
            pos = index(line, "nPoints") + index(line, "nFaces") + index(line, "size")
            if (pos == 0) cycle
            
            read(line(pos:), *, iostat=iostat) itemCount
            if (itemCount > 0) then
                startLine = i + 2
                return
            end if
        end do
        
        do i = 1, size(fileLines)
            if (isNumber(trim(fileLines(i)))) then
                read(fileLines(i), *, iostat=iostat) itemCount
                if (itemCount > 0) then
                    startLine = i + 2
                    return
                end if
            end if
        end do
        
        do i = 1, size(fileLines)
            if (index(fileLines(i), "(") > 0) then
                startLine = i + 1
                exit
            end if
        end do
        itemCount = 0
        do i = startLine, size(fileLines)
            if (index(fileLines(i), ")") > 0) exit
            itemCount = itemCount + 1
        end do
    end subroutine OFFile_FindNItems

    ! 辅助函数：判断是否为数字
    logical function isNumber(str)
        character(len=*), intent(in) :: str
        real(kind=8) :: num
        integer :: iostat
        isNumber = .false.
        read(str, *, iostat=iostat) num
        if (iostat == 0) isNumber = .true.
    end function isNumber

    ! 辅助函数：查找包含子串的行
    integer function findInLines(substr, lines, startLine)
        character(len=*), intent(in) :: substr, lines(:)
        integer, intent(in) :: startLine
        integer :: i
        findInLines = 0
        do i = startLine, size(lines)
            if (index(lines(i), substr) > 0) then
                findInLines = i
                return
            end if
        end do
    end function findInLines

    ! 辅助函数：获取空闲单元号
    integer function get_free_unit()
        implicit none
        integer :: i, iostat
        logical :: opened
        do i = 10, 999
            inquire(unit=i, opened=opened, iostat=iostat)
            if (.not. opened .and. iostat == 0) then
                get_free_unit = i
                return
            end if
        end do
        error stop "No free unit available"
    end function get_free_unit

    ! 其他必要函数（简化实现）
    logical function isInArray(pts, pt)
        real(kind=8), intent(in) :: pts(:,:), pt(3)
        integer :: i
        isInArray = .false.
        do i = 1, size(pts,1)
            if (all(abs(pts(i,:)-pt) < 1d-12)) then
                isInArray = .true.
                return
            end if
        end do
    end function isInArray

    subroutine deallocateMeshData(md)
        type(MeshData), intent(inout) :: md
        integer :: i
        if (allocated(md%points)) deallocate(md%points)
        if (allocated(md%owner)) deallocate(md%owner)
        if (allocated(md%neighbour)) deallocate(md%neighbour)
        if (allocated(md%boundaryNames)) deallocate(md%boundaryNames)
        if (allocated(md%boundaryNumFaces)) deallocate(md%boundaryNumFaces)
        if (allocated(md%boundaryStartFaces)) deallocate(md%boundaryStartFaces)
        if (allocated(md%faces%faces)) then
            do i = 1, size(md%faces%faces)
                if (allocated(md%faces%faces(i)%points)) deallocate(md%faces%faces(i)%points)
            end do
            deallocate(md%faces%faces)
        end if
    end subroutine deallocateMeshData


!===============================================================================
! SECTION 7: MODIFIED NUMERICS - RANS versions
!===============================================================================

	!===========================================================================
	! RANS k-omega: Modified Green-Gauss gradient (handles variable number of vars)
	!===========================================================================
	function greenGaussGrad_RANS(mesh, matrix, valuesAtFaces) result(grad)
		implicit none
		type(Meshh), intent(in) :: mesh
		real(kind=8), intent(in) :: matrix(:,:)
		logical, intent(in), optional :: valuesAtFaces
		real(kind=8), allocatable :: grad(:,:,:)
		
		! Same implementation as original, but works with variable dimensions
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4)
		integer :: nVars, f, v, d, ownerCell, neighbourCell
		real(kind=8), allocatable :: faceVals(:,:)
		
		meshInfo = unstructuredMeshInfo(mesh)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		nVars = size(matrix, 2)
		
		allocate(faceVals(nFaces, nVars))
		faceVals = 0.0d0
		
		if (present(valuesAtFaces)) then
			if (valuesAtFaces) then
				faceVals = matrix
			else
				call linInterp_3D_RANS(mesh, matrix, faceVals)
			end if
		else
			call linInterp_3D_RANS(mesh, matrix, faceVals)
		end if
		
		allocate(grad(nCells, nVars, 3))
		grad = 0.0d0
		
		do f = 1, nFaces
			ownerCell = mesh%faces(f,1)
			neighbourCell = mesh%faces(f,2)
			
			do v = 1, nVars
				do d = 1, 3
					grad(ownerCell, v, d) = grad(ownerCell, v, d) + &
					                        mesh%fAVecs(f,d) * faceVals(f,v)
					
					if (neighbourCell > -1) then
						grad(neighbourCell, v, d) = grad(neighbourCell, v, d) - &
						                            mesh%fAVecs(f,d) * faceVals(f,v)
					end if
				end do
			end do
		end do
		
		do d = 1, 3
			do v = 1, nVars
				do f = 1, nCells
					grad(f, v, d) = grad(f, v, d) / mesh%cVols(f)
				end do
			end do
		end do
		
		deallocate(faceVals)
	end function greenGaussGrad_RANS

	!===========================================================================
	! RANS k-omega: Modified linear interpolation
	!===========================================================================
	subroutine linInterp_3D_RANS(mesh, matrix, faceVals)
		implicit none
		type(Meshh), intent(in) :: mesh
		real(kind=8), intent(in) :: matrix(:,:)
		real(kind=8), intent(inout) :: faceVals(:,:)
		
		! Same implementation as original
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4)
		integer :: nVars, f, v, i, c1, c2
		real(kind=8) :: c1Dist, c2Dist, totalDist
		
		meshInfo = unstructuredMeshInfo(mesh)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		nVars = size(matrix, 2)
		
		do f = 1, nFaces - nBdryFaces
			c1 = mesh%faces(f,1)
			c2 = mesh%faces(f,2)
			
			c1Dist = 0.0d0
			c2Dist = 0.0d0
			do i = 1, 3
				c1Dist = c1Dist + (mesh%cCenters(c1,i) - mesh%fCenters(f,i))**2
				c2Dist = c2Dist + (mesh%cCenters(c2,i) - mesh%fCenters(f,i))**2
			end do
			totalDist = c1Dist + c2Dist
			
			do v = 1, nVars
				faceVals(f, v) = matrix(c1, v) * (c2Dist / totalDist) + &
				                 matrix(c2, v) * (c1Dist / totalDist)
			end do
		end do
	end subroutine linInterp_3D_RANS

	!===========================================================================
	! RANS k-omega: Modified face deltas
	!===========================================================================
	function faceDeltas_RANS(mesh, sln) result(deltas)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(in) :: sln
		real(kind=8), allocatable :: deltas(:,:)
		
		! Same implementation as original, works with 7 variables
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4)
		integer :: nVars, f, v, ownerCell, neighbourCell
		
		meshInfo = unstructuredMeshInfo(mesh)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		nVars = size(sln%cellState, 2)
		
		allocate(deltas(nFaces, nVars))
		deltas = 0.0d0
		
		do f = 1, nFaces - nBdryFaces
			ownerCell = mesh%faces(f,1)
			neighbourCell = mesh%faces(f,2)
			
			do v = 1, nVars
				deltas(f, v) = sln%cellState(neighbourCell, v) - sln%cellState(ownerCell, v)
			end do
		end do
	end function faceDeltas_RANS

	!===========================================================================
	! RANS k-omega: Modified decode solution
	!===========================================================================
	subroutine decodeSolution_3D_RANS(sln, fluid)
		implicit none
		type(SolutionState), intent(inout) :: sln
		type(Fluidd), intent(in) :: fluid
		integer :: nCells, i
		
		nCells = size(sln%cellState, 1)
		
		if (.not. allocated(sln%cellPrimitives)) allocate(sln%cellPrimitives(nCells, 7))
		if (.not. allocated(sln%cellFluxes)) allocate(sln%cellFluxes(nCells, 21))
		
		do i = 1, nCells
			call decodePrimitives3D_RANS(sln%cellPrimitives(i,:), sln%cellState(i,:), fluid)
			call calculateFluxes3D_RANS(sln%cellFluxes(i,:), sln%cellPrimitives(i,:), &
			                            sln%cellState(i,:), fluid, sln, i)
		end do
	end subroutine decodeSolution_3D_RANS

	!===========================================================================
	! RANS k-omega: Modified flux integration
	!===========================================================================
	function integrateFluxes_unstructured3D_RANS(mesh, sln) result(fluxResiduals)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(inout) :: sln
		real(kind=8), allocatable :: fluxResiduals(:,:)
		
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4)
		integer :: nVars, f, v, ownerCell, neighbourCell, i1, i2
		real(kind=8) :: flow
		
		meshInfo = unstructuredMeshInfo(mesh)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		nVars = size(sln%cellState, 2)  ! 7 for RANS
		
		if (.not. allocated(sln%fluxResiduals) .or. &
		    size(sln%fluxResiduals,1)/=nCells .or. &
		    size(sln%fluxResiduals,2)/=nVars) then
			if (allocated(sln%fluxResiduals)) deallocate(sln%fluxResiduals)
			allocate(sln%fluxResiduals(nCells, nVars))
		end if
		
		sln%fluxResiduals = 0.0d0
		
		! Integrate face fluxes
		do f = 1, nFaces
			ownerCell = mesh%faces(f,1)
			neighbourCell = mesh%faces(f,2)
			
			do v = 1, nVars
				i1 = (v-1)*3 + 1
				i2 = i1 + 2
				flow = dot_product(sln%faceFluxes(f,i1:i2), mesh%fAVecs(f,:))
				
				sln%fluxResiduals(ownerCell, v) = sln%fluxResiduals(ownerCell, v) - flow
				
				if (neighbourCell > 0 .and. neighbourCell <= nCells) then
					sln%fluxResiduals(neighbourCell, v) = sln%fluxResiduals(neighbourCell, v) + flow
				end if
			end do
		end do
		
		! Add source terms for k and omega equations
		call addTurbulenceSourceTerms(mesh, sln)
		
		! Divide by cell volume
		do v = 1, nVars
			do f = 1, nCells
				sln%fluxResiduals(f, v) = sln%fluxResiduals(f, v) / mesh%cVols(f)
			end do
		end do
		
		allocate(fluxResiduals(nCells, nVars))
		fluxResiduals = sln%fluxResiduals
	end function integrateFluxes_unstructured3D_RANS

	!===========================================================================
	! RANS k-omega: Add turbulence source terms to residuals
	!===========================================================================
	subroutine addTurbulenceSourceTerms(mesh, sln)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(inout) :: sln
		real(kind=8), allocatable :: sources(:,:)
		integer :: nCells, c
		type(Fluidd) :: fluid
		
		nCells = size(sln%cellState, 1)
		allocate(sources(nCells, 2))
		
		fluid = fluiddd()
		call computeTurbulenceSourceTerms(mesh, sln, fluid, sources)
		
		
		do c = 1, nCells
			! integrate sources over the cell before volume normalization
			sln%fluxResiduals(c, 6) = sln%fluxResiduals(c, 6) + sources(c, 1) * mesh%cVols(c)
			sln%fluxResiduals(c, 7) = sln%fluxResiduals(c, 7) + sources(c, 2) * mesh%cVols(c)
		end do

		
		deallocate(sources)
	end subroutine addTurbulenceSourceTerms
	!===========================================================================
	! Original MaxInterp Function and deep_copy_sln
	!===========================================================================

	function maxInterp(mesh, sj, rj) result(faceVals)
		type(Meshh), intent(in) :: mesh
		real(kind=8), intent(in) :: sj(:), rj(:)
		real(kind=8), allocatable :: faceVals(:,:)
		
		integer(kind=8) :: nCells, nFaces, nBoundaries, nBdryFaces, meshInfo(4)
		integer(kind=8) :: f, c1, c2
		integer :: nVars = 2  ! 固定为2个变量（sj和rj）
		
		! 获取网格信息
		meshInfo = unstructuredMeshInfo(mesh)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		
		! 分配内存给输出数组
		allocate(faceVals(nFaces, nVars))
		faceVals = 0.0d0  ! 初始化为0
		
		! 遍历内部面
		do f = 1, nFaces - nBdryFaces
		    c1 = mesh%faces(f, 1)  ! 主单元格
		    c2 = mesh%faces(f, 2)  ! 邻居单元格
		    
		    ! 对每个变量取最大值
		    faceVals(f, 1) = max(sj(c1), sj(c2))  ! 第一列存储sj的插值结果
		    faceVals(f, 2) = max(rj(c1), rj(c2))  ! 第二列存储rj的插值结果
		end do		
	end function maxInterp

	subroutine deep_copy_sln(src, dest)
		type(SolutionState), intent(in) :: src
		type(SolutionState), intent(inout) :: dest

		if (allocated(dest%cellState)) deallocate(dest%cellState)
		allocate(dest%cellState(size(src%cellState,1), size(src%cellState,2)))
		dest%cellState = src%cellState

		if (allocated(dest%cellPrimitives)) deallocate(dest%cellPrimitives)
		allocate(dest%cellPrimitives(size(src%cellPrimitives,1), size(src%cellPrimitives,2)))
		dest%cellPrimitives = src%cellPrimitives

		if (allocated(dest%cellFluxes)) deallocate(dest%cellFluxes)
		allocate(dest%cellFluxes(size(src%cellFluxes,1), size(src%cellFluxes,2)))
		dest%cellFluxes = src%cellFluxes

		if (allocated(dest%fluxResiduals)) deallocate(dest%fluxResiduals)
		allocate(dest%fluxResiduals(size(src%fluxResiduals,1), size(src%fluxResiduals,2)))
		dest%fluxResiduals = src%fluxResiduals

		if (allocated(dest%faceFluxes)) deallocate(dest%faceFluxes)
		allocate(dest%faceFluxes(size(src%faceFluxes,1), size(src%faceFluxes,2)))
		dest%faceFluxes = src%faceFluxes
	end subroutine deep_copy_sln

!===============================================================================
! SECTION 8: TIME INTEGRATION - RANS k-omega modifications
!===============================================================================

	!===========================================================================
	! RANS k-omega: Modified LTS Euler time integration
	!===========================================================================
	subroutine LTSEuler_RANS(mesh, sln, boundaryConditions, fluid, dt)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SolutionState), intent(inout) :: sln
		type(BoundaryCondition), intent(in) :: boundaryConditions(:)
		type(Fluidd), intent(in) :: fluid
		real(kind=8), intent(inout) :: dt(:)
		
		real(kind=8) :: targetCFL
		real(kind=8), allocatable :: fluxResiduals(:,:)
		real(kind=8), allocatable :: cfl_vals(:)
		integer :: nCells, nVars, i, j
		
		targetCFL = dt(1)
		
		! Compute residuals (includes turbulence equations)
		fluxResiduals = unstructured_JSTFlux_RANS(mesh, sln, boundaryConditions, fluid)
		
		nCells = size(sln%cellState, 1)
		nVars = size(fluxResiduals, 2)  ! 7 for RANS
		
		! Compute time step
		!dt = CFL(mesh, sln, fluid, 1.0d0)
		!dt = targetCFL / dt
		allocate(cfl_vals(nCells))
    		cfl_vals = CFL(mesh, sln, fluid, 1.0d0)
    		dt = targetCFL / cfl_vals
    		deallocate(cfl_vals)
		
		! Smooth time step
		dt = smoothTimeStep(dt, mesh, 0.1d0)
		dt = smoothTimeStep(dt, mesh, 0.1d0)
		
		! Update conservative variables
		do j = 1, nVars
			do i = 1, nCells
				sln%cellState(i,j) = sln%cellState(i,j) + fluxResiduals(i,j) * dt(i)
			end do
		end do
		
		!--- RANS k-omega: Apply limiters to turbulence variables ---
		do i = 1, nCells
			! Limit k (must be positive)
			sln%cellState(i, 6) = max(sln%cellState(i, 6), sln%cellState(i, 1) * SMALL_NUM)
			
			! Limit omega (must be positive)
			sln%cellState(i, 7) = max(sln%cellState(i, 7), sln%cellState(i, 1) * SMALL_NUM)
		end do
		
		! Decode solution (compute primitives and fluxes)
		call decodeSolution_3D_RANS(sln, fluid)
		
		deallocate(fluxResiduals)
		do i = 1, nCells
			! Check for invalid values
			if (sln%cellState(i,1) < SMALL_NUM) then
				print *, "ERROR: Negative density in cell", i, " rho=", sln%cellState(i,1)
				print *, "  Primitives:", sln%cellPrimitives(i,:)
				stop
			end if
        
			if (sln%cellState(i,5) < 0.0d0) then
				print *, "WARNING: Negative energy in cell", i, " E=", sln%cellState(i,5)
				! Clip to positive value
				sln%cellState(i,5) = sln%cellState(i,1) * (fluid%Cp * 100.0d0 + &
                                 0.5d0 * sum(sln%cellPrimitives(i,3:5)**2))
			end if
		end do
	end subroutine LTSEuler_RANS

	! Time step smoothing (unchanged from original)
	function smoothTimeStep(dt, mesh, diffusionCoefficient) result(updated_dt)
		real(kind=8), intent(inout) :: dt(:)
		type(Meshh), intent(in) :: mesh
		real(kind=8), intent(in), optional :: diffusionCoefficient
		real(kind=8), allocatable :: updated_dt(:)
		
		integer :: nCells, nFaces, nBoundaries, nBdryFaces, f, ownerCell, neighbourCell, i, meshInfo(4), stat
		real(kind=8) :: coeff, timeFlux
		real(kind=8), allocatable :: timeFluxes(:), surfaceAreas(:)
		
		if (present(diffusionCoefficient)) then
			coeff = diffusionCoefficient
		else
			coeff = 0.2
		end if
		
		meshInfo = unstructuredMeshInfo(mesh)
		nCells = meshInfo(1)
		nFaces = meshInfo(2)
		nBoundaries = meshInfo(3)
		nBdryFaces = meshInfo(4)
		
		allocate(timeFluxes(nCells))
		allocate(surfaceAreas(nCells))
		timeFluxes = 0.0d0
		surfaceAreas = 0.0d0
		
		do f = 1, nFaces-nBdryFaces
			ownerCell = mesh%faces(f, 1)
			neighbourCell = mesh%faces(f, 2)
			
			timeFlux = (dt(ownerCell) - dt(neighbourCell)) * mag(mesh%fAVecs(f,:))
			
			surfaceAreas(ownerCell) = surfaceAreas(ownerCell) + mag(mesh%fAVecs(f,:))
			surfaceAreas(neighbourCell) = surfaceAreas(neighbourCell) + mag(mesh%fAVecs(f,:))
			
			timeFluxes(ownerCell) = timeFluxes(ownerCell) - timeFlux
			timeFluxes(neighbourCell) = timeFluxes(neighbourCell) + timeFlux
		end do
		
		timeFluxes = timeFluxes * (coeff / surfaceAreas)
		
		do i = 1, nCells
			timeFluxes(i) = min(0.0d0, timeFluxes(i))
		end do
		
		dt = dt + timeFluxes
		allocate(updated_dt(size(dt)), stat=stat)
		if (stat /= 0) stop "Error: Allocation updated_dt failed"
		updated_dt = dt
		
		deallocate(timeFluxes, surfaceAreas)
	end function smoothTimeStep

!===============================================================================
! SECTION 9: OUTPUT - RANS k-omega modifications
!===============================================================================

	!===========================================================================
	! RANS k-omega: Modified VTK output (includes turbulence fields)
	!===========================================================================
	subroutine writeVTKFile_RANS(fileName, points, cells, primitives, &
	                             muT, yPlus, nPts, nCells, vtkTypes)
		character(len=*), intent(in) :: fileName
		real(kind=8), intent(in) :: points(:,:)
		type(CelllArray), intent(in) :: cells
		! real(kind=8), intent(in) :: primitives(:,:)  ! Original: 5 vars
		real(kind=8), intent(in) :: primitives(:,:)  ! RANS: 7 vars [P, T, u, v, w, k, omega]
		real(kind=8), intent(in) :: muT(:), yPlus(:)
		integer, intent(in) :: nPts, nCells
		integer, intent(in) :: vtkTypes(8)
		
		integer :: fileUnit, i, j, n, cellType, ptIdx, connOffset
		real(kind=8) :: P, T, Ux, Uy, Uz, k, omega, mu_t, y_plus
		real(kind=8), parameter :: INF_THRESHOLD = 1.0d300
		
		! Check dimensions
		if (size(primitives, 1) /= nCells .or. size(primitives, 2) < 7) then
			print *, "Error: RANS primitives dimension mismatch"
			return
		end if
		
		open(newunit=fileUnit, file=fileName, status='replace', action='write', form='formatted')
		
		! VTK header
		write(fileUnit, '(A)') '# vtk DataFile Version 3.0'
		write(fileUnit, '(A)') 'RANS k-omega Solution from BuFlowModule'
		write(fileUnit, '(A)') 'ASCII'
		write(fileUnit, '(A)') 'DATASET UNSTRUCTURED_GRID'
		
		! Points
		write(fileUnit, '(A,I0,A)') 'POINTS ', nPts, ' double'
		do i = 1, nPts
			write(fileUnit, '(3F16.8)') points(1,i), points(2,i), points(3,i)
		end do
		
		! Cell connectivity
		connOffset = 0
		do i = 1, nCells
			connOffset = connOffset + size(cells%cells(i)%pointIndices) + 1
		end do
		write(fileUnit, '(A,I0,A,I0)') 'CELLS ', nCells, ' ', connOffset
		do i = 1, nCells
			n = size(cells%cells(i)%pointIndices)
			write(fileUnit, '(I4)', advance='no') n
			do j = 1, n
				ptIdx = cells%cells(i)%pointIndices(j) - 1
				write(fileUnit, '(I6)', advance='no') ptIdx
			end do
			write(fileUnit, *)
		end do
		
		! Cell types
		write(fileUnit, '(A,I0)') 'CELL_TYPES ', nCells
		do i = 1, nCells
			n = size(cells%cells(i)%pointIndices)
			cellType = merge(vtkTypes(n), 1, n>=1 .and. n<=8)
			write(fileUnit, '(I4)') cellType
		end do
		
		! Cell data
		write(fileUnit, '(A,I0)') 'CELL_DATA ', nCells
		
		! Pressure
		write(fileUnit, '(A)') 'SCALARS Pressure double 1'
		write(fileUnit, '(A)') 'LOOKUP_TABLE default'
		do i = 1, nCells
			P = primitives(i,1)
			if (P /= P .or. abs(P) > INF_THRESHOLD) P = 0.0d0
			write(fileUnit, '(F16.8)') P
		end do
		
		! Temperature
		write(fileUnit, '(A)') 'SCALARS Temperature double 1'
		write(fileUnit, '(A)') 'LOOKUP_TABLE default'
		do i = 1, nCells
			T = primitives(i,2)
			if (T /= T .or. abs(T) > INF_THRESHOLD) T = 0.0d0
			write(fileUnit, '(F16.8)') T
		end do
		
		! Velocity
		write(fileUnit, '(A)') 'VECTORS Velocity double'
		do i = 1, nCells
			Ux = primitives(i,3)
			Uy = primitives(i,4)
			Uz = primitives(i,5)
			if (Ux /= Ux .or. abs(Ux) > INF_THRESHOLD) Ux = 0.0d0
			if (Uy /= Uy .or. abs(Uy) > INF_THRESHOLD) Uy = 0.0d0
			if (Uz /= Uz .or. abs(Uz) > INF_THRESHOLD) Uz = 0.0d0
			write(fileUnit, '(3F16.8)') Ux, Uy, Uz
		end do
		
		!--- RANS k-omega: Turbulent kinetic energy ---
		write(fileUnit, '(A)') 'SCALARS TKE double 1'
		write(fileUnit, '(A)') 'LOOKUP_TABLE default'
		do i = 1, nCells
			k = primitives(i,6)
			if (k /= k .or. abs(k) > INF_THRESHOLD) k = 0.0d0
			write(fileUnit, '(E16.8)') k
		end do
		
		!--- RANS k-omega: Specific dissipation rate ---
		write(fileUnit, '(A)') 'SCALARS Omega double 1'
		write(fileUnit, '(A)') 'LOOKUP_TABLE default'
		do i = 1, nCells
			omega = primitives(i,7)
			if (omega /= omega .or. abs(omega) > INF_THRESHOLD) omega = 0.0d0
			write(fileUnit, '(E16.8)') omega
		end do
		
		!--- RANS k-omega: Turbulent viscosity ---
		write(fileUnit, '(A)') 'SCALARS EddyViscosity double 1'
		write(fileUnit, '(A)') 'LOOKUP_TABLE default'
		do i = 1, nCells
			mu_t = muT(i)
			if (mu_t /= mu_t .or. abs(mu_t) > INF_THRESHOLD) mu_t = 0.0d0
			write(fileUnit, '(E16.8)') mu_t
		end do
		
		!--- RANS k-omega: Y-plus ---
		write(fileUnit, '(A)') 'SCALARS Yplus double 1'
		write(fileUnit, '(A)') 'LOOKUP_TABLE default'
		do i = 1, nCells
			y_plus = yPlus(i)
			if (y_plus /= y_plus .or. abs(y_plus) > INF_THRESHOLD) y_plus = 0.0d0
			write(fileUnit, '(F16.8)') y_plus
		end do
		
		close(fileUnit)
	end subroutine writeVTKFile_RANS

	!===========================================================================
	! RANS k-omega: Modified output routine
	!===========================================================================
	subroutine outputVTK_RANS(meshPath, cellPrimitives, sln, fileName)
		character(len=*), intent(in) :: meshPath
		real(kind=8), intent(in) :: cellPrimitives(:,:)
		type(SolutionState), intent(in) :: sln
		character(len=*), intent(in), optional :: fileName
		
		real(kind=8), allocatable :: pointLocations(:,:)
		type(CelllArray) :: cells
		character(len=256) :: vtkFileName
		integer :: nCells, nPoints
		
		call OpenFOAMMesh_findCellPts(meshPath, pointLocations, cells)
		nCells = size(cells%cells)
		nPoints = size(pointLocations, 1)
		
		call transposePoints(pointLocations)
		
		if (present(fileName)) then
			vtkFileName = trim(fileName) // ".vtk"
		else
			vtkFileName = "solution_RANS.vtk"
		end if
		
		call writeVTKFile_RANS(vtkFileName, pointLocations, cells, cellPrimitives, &
		                       sln%cellMuT, sln%cellYplus, nPoints, nCells, &
		                       [1,3,5,10,14,13,-1,12])
		
		deallocate(pointLocations)
	end subroutine outputVTK_RANS

	!===========================================================================
	! RANS k-omega: Modified writeOutput wrapper
	!===========================================================================
	subroutine writeOutput_RANS(sln, restartFile, meshPath, createRestartFile, createVTKOutput)
		type(SolutionState), intent(in) :: sln
		character(len=*), intent(in) :: restartFile
		character(len=*), intent(in) :: meshPath
		logical, intent(in) :: createRestartFile
		logical, intent(in) :: createVTKOutput
		
		integer :: maxNum, unit, iostat, num
		character(len=256) :: files(1000)
		integer :: nfiles, i, pos, start, endd
		character(len=256) :: solnName, tmpFile
		
		! Write restart file
		if (createRestartFile) then
			print *, 'Writing RANS Restart File: ', trim(restartFile)
			call writeRestartFile(sln%cellPrimitives, restartFile)
		end if
		
		! Write VTK output
		if (createVTKOutput) then
			! Find maximum solution number
			tmpFile = 'dir_list.tmp'
			call system('ls > '//trim(tmpFile))
			
			nfiles = 0
			open(newunit=unit, file=tmpFile, status='old', action='read', iostat=iostat)
			if (iostat == 0) then
				do while (iostat == 0 .and. nfiles < 1000)
					nfiles = nfiles + 1
					read(unit, '(a)', iostat=iostat) files(nfiles)
				end do
				nfiles = nfiles - 1
				close(unit)
			end if
			call system('rm '//trim(tmpFile))
			
			maxNum = 0
			do i = 1, nfiles
				pos = index(trim(files(i)), 'solution_RANS')
				if (pos == 1) then
					start = len_trim('solution_RANS') + 2
					endd = index(trim(files(i)), '.vtk') - 1
					if (start < endd) then
						read(files(i)(start:endd), *, iostat=iostat) num
						if (iostat == 0 .and. num > maxNum) then
							maxNum = num
						end if
					end if
				end if
			end do
			
			write(solnName, '(a,i0)') 'solution_RANS.', maxNum + 1
			print *, 'Writing ', trim(solnName)
			
			call outputVTK_RANS(meshPath, sln%cellPrimitives, sln, solnName)
		end if
	end subroutine writeOutput_RANS

	! Restart file writer (unchanged)
	subroutine writeRestartFile(cellPrimitives, path)
		real(kind=8), intent(in) :: cellPrimitives(:,:)
		character(len=*), intent(in), optional :: path
		character(len=256) :: filePath
		
		if (present(path)) then
			filePath = path
		else
			filePath = "FvCFDRestart_RANS.txt"
		end if
		
		open(unit=10, file=filePath, status='replace', action='write')
		write(10, *) cellPrimitives
		close(10)
	end subroutine writeRestartFile

	! Transpose points (unchanged)
	subroutine transposePoints(points)
		real(kind=8), allocatable, intent(inout) :: points(:,:)
		real(kind=8), allocatable :: temp(:,:)
		integer :: i, j, nPoints
		
		nPoints = size(points, 1)
		allocate(temp(3, nPoints))
		
		do i = 1, nPoints
			do j = 1, 3
				temp(j, i) = points(i, j)
			end do
		end do
		
		deallocate(points)
		call move_alloc(temp, points)
	end subroutine transposePoints
!===============================================================================
! SECTION 10: Vector Functions
!===============================================================================
!######################### 向量函数 ########################
!#函数：点积，叉积，模，转化为单位向量
	! 情况1：向量 × 向量（均为1维数组）
	function dot_vec_vec(arg1, arg2) result(out)
		real(kind=8), intent(in) :: arg1(:)  ! 1维向量
		real(kind=8), intent(in) :: arg2(:)  ! 1维向量
		real(kind=8), allocatable :: out(:)
		integer :: i, n

		n = size(arg1)
		if (n /= size(arg2)) then
		    print *, "Error: Vectors must have same length"
		    stop
		end if
		allocate(out(1))
		out(1) = 0.0d0
		do i = 1, n
		    out(1) = out(1) + arg1(i) * arg2(i)  ! 合法的1维数组索引
		end do
	end function dot_vec_vec

	! 情况2：向量 × 二维矩阵（含向量矩阵）
	function dot_vec_mat(arg1, arg2) result(out)
		real(kind=8), intent(in) :: arg1(:)    ! 1维向量
		real(kind=8), intent(in) :: arg2(:,:)  ! 二维矩阵（含向量矩阵）
		real(kind=8), allocatable :: out(:), temp(:)
		integer :: i, n_rows, arg1_size, arg2_dim1, arg2_dim2

		arg1_size = size(arg1)
		arg2_dim1 = size(arg2, 1)  ! 矩阵行数
		arg2_dim2 = size(arg2, 2)  ! 矩阵列数

		! 子情况A：向量 × 向量矩阵（N×M矩阵，每行是M维向量）
		if (arg2_dim2 > 1) then
		    if (arg1_size /= arg2_dim2) then
		        print *, "Error: Vector length must match sub-vector length"
		        stop
		    end if
		    allocate(out(arg2_dim1))
		    do i = 1, arg2_dim1
		        temp = dot_vec_vec(arg1, arg2(i,:))   ! 先存储函数返回的数组
				out(i) = temp(1)   ! 调用向量×向量函数
		    end do

		! 子情况B：向量 × 二维矩阵（N×1矩阵，每行是标量）
		else if (arg2_dim2 == 1) then
		    if (arg1_size /= arg2_dim1) then
		        print *, "Error: Vector length must match matrix rows"
		        stop
		    end if
		    allocate(out(1))
		    out(1) = 0.0d0
		    do i = 1, arg2_dim1
		        out(1) = out(1) + arg1(i) * arg2(i,1)  ! 合法的二维数组索引
		    end do
		end if
	end function dot_vec_mat
	! 3D向量叉积（仅适用于3元素向量）
	function cross(v1, v2) result(cross_vec)
		implicit none
		real(kind=8), intent(in) :: v1(3), v2(3)  ! v1=(x1,y1,z1), v2=(x2,y2,z2)
		real(kind=8) :: cross_vec(3)
		
		! 正确的叉乘公式（严格匹配Julia）
		cross_vec(1) = v1(2)*v2(3) - v1(3)*v2(2)        ! y1z2 - z1y2
		cross_vec(2) = v1(3)*v2(1) - v1(1)*v2(3)        ! z1x2 - x1z2（修正符号）
		cross_vec(3) = v1(1)*v2(2) - v1(2)*v2(1)        ! x1y2 - y1x2
	end function cross

	! 向量的模（2-范数）
	real(kind=8) function mag(vec) result(sqrSum)
		real(kind=8), intent(in) :: vec(:)
		integer :: i
		sqrSum = 0.0d0
		do i = 1, size(vec)
		    sqrSum = sqrSum + vec(i)**2  ! 累加元素平方
		end do
		sqrSum = sqrt(sqrSum)  ! 开平方得模长
	end function mag

	! 向量归一化（转化为单位向量）
	function normalize(vec) result(unit_vec)
		real(kind=8), intent(in) :: vec(:)
		real(kind=8), allocatable :: unit_vec(:)
		real(kind=8) :: vec_mag
		vec_mag = mag(vec)  ! 调用mag函数
		! 避免除以零
		if (vec_mag < 1.0d-12) then
		    print *, "Error in normalize: Zero magnitude vector"
		    stop
		end if
		allocate(unit_vec(size(vec)))
		unit_vec = vec / vec_mag  ! 每个元素除以模长
	end function normalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================================
! SECTION 11: EXAMPLE INITIALIZATION - RANS k-omega
!===============================================================================

	!===========================================================================
	! RANS k-omega: Modified compute_CFD for RANS
	!===========================================================================
	subroutine compute_CFD_RANS
		type(Fluidd) :: fluid
		type(Meshh) :: mesh
		real(kind=8), allocatable :: cellPrimitives(:,:)
		real(kind=8) :: P, T, U(3), UunitVec(3), a, machNum, Pt, Tt
		character(256) :: meshPath
		type(BoundaryCondition), allocatable :: boundaryConditions(:)
		logical :: hasNaN
		real(kind=8) :: Pmin, Pmax, Tmin, Tmax, Uxmin, Uxmax
		
		!--- RANS k-omega: Turbulence initialization ---
		real(kind=8) :: k_init, omega_init, turbulence_intensity, length_scale
		
		! Initialize fluid properties
		fluid = Fluidd(Cp=1005.0d0, R=287.05d0, gammaa=1.4d0, &
		               mu_laminar=1.7894d-5, Pr=0.72d0, Pr_turb=0.9d0)
		
		! Flow conditions
		P = 100000.0d0
		T = 300.0d0
		U = [11.0d0, 0.0d0, 0.0d0]  ! 汽车常用速度 ~40 km/h
		UunitVec = normalize(U)
		a = sqrt(fluid%gammaa * fluid%R * T)
		machNum = mag(U)/a
		print *, "Mach number = ", machNum
		
		Pt = P*(1.0d0 + ((fluid%gammaa-1.0d0)/2.0d0)*machNum**2)**(fluid%gammaa/(fluid%gammaa-1.0d0))
		Tt = T*(1.0d0 + ((fluid%gammaa-1.0d0)/2.0d0)*machNum**2)
		
		!--- RANS k-omega: Initialize turbulence quantities ---
		turbulence_intensity = 0.01d0  ! 5% turbulence intensity
		length_scale = 0.01d0          ! Turbulent length scale (m)
		
		k_init = 1.5d0 * (turbulence_intensity * mag(U))**2
		omega_init = sqrt(k_init) / (sqrt(sqrt(BETA_STAR)) * length_scale)
		
		print *, "Initial k = ", k_init
		print *, "Initial omega = ", omega_init
		
		! Boundary conditions
		! 适配 5 个边界: airfoil, empty, inlet, outlet, symmetry
		allocate(boundaryConditions(5))
		
		! Wall
		boundaryConditions(1)%type = wallBoundary
		allocate(boundaryConditions(1)%params(0))
		
		! Empty
		boundaryConditions(2)%type = emptyBoundary
		allocate(boundaryConditions(2)%params(0))
		
		! Inlet (now includes k and omega)
		boundaryConditions(3)%type = InletBoundary
		allocate(boundaryConditions(3)%params(6))
		boundaryConditions(3)%params = [Pt, Tt, UunitVec(1), UunitVec(2), k_init, omega_init]
		
		! Outlet
		boundaryConditions(4)%type = OutletBoundary
		allocate(boundaryConditions(4)%params(1))
		boundaryConditions(4)%params(1) = P

		! Symmetry (Euler/RANS density solver historically treats this as wall)
		boundaryConditions(5)%type = wallBoundary
		allocate(boundaryConditions(5)%params(0))
		
		! Read mesh
		meshPath = "mesh/OFairfoilMesh"
		print *, "=== Loading mesh from: ", trim(meshPath)
		mesh = OpenFOAMMesh(meshPath)
		
		!--- Validate mesh was loaded ---
	if (.not. allocated(mesh%cCenters)) then
		print *, "ERROR: Mesh not loaded properly - cCenters not allocated"
		stop
	end if
	if (.not. allocated(mesh%cVols)) then
		print *, "ERROR: Mesh not loaded properly - cVols not allocated"
		stop
	end if

	print *, "=== Mesh loaded successfully ==="
	print *, "  Number of cells: ", size(mesh%cCenters, 1)
	print *, "  Number of faces: ", size(mesh%faces, 1)
	print *, "  Cell volume range: ", minval(mesh%cVols), maxval(mesh%cVols)
	call flush(6)

	print *, "  P = ", P
	print *, "  T = ", T
	print *, "  U = ", U
	print *, "  k_init = ", k_init
	print *, "  omega_init = ", omega_init

		! Dispatch to selected solver
		if (SOLVER_TYPE == SOLVER_SIMPLE) then
			call configureBoundaryConditionsFromMesh(mesh, boundaryConditions, P, Pt, Tt, &
			                                     UunitVec, k_init, omega_init)
			print *, ''
			print *, '>>> Using SIMPLE pressure-based solver (low Mach) <<<'
			call solve_SIMPLE(mesh, meshPath, boundaryConditions, fluid, &
			                  P, T, U, k_init, omega_init)
		else
			print *, ''
			print *, '>>> Using density-based solver <<<'
			cellPrimitives = initializeUniformSolution3D(mesh, P, T, U(1), U(2), U(3), k_init, omega_init)
			call solve(mesh, meshPath, cellPrimitives, boundaryConditions)
		end if
		
		if (allocated(boundaryConditions)) call freeBoundaryConditions(boundaryConditions)
	end subroutine compute_CFD_RANS


!===============================================================================

!===============================================================================

!-----------------------------------------------------------------------
! 根据 OpenFOAM boundary 名称配置边界条件，避免 SIMPLEC 依赖硬编码顺序
!-----------------------------------------------------------------------
	subroutine configureBoundaryConditionsFromMesh(mesh, boundaryConditions, P, Pt, Tt, &
	                                             UunitVec, k_init, omega_init)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(BoundaryCondition), allocatable, intent(inout) :: boundaryConditions(:)
		real(kind=8), intent(in) :: P, Pt, Tt, UunitVec(3), k_init, omega_init
		integer :: b, nBoundaries
		character(len=100) :: bname

		if (.not. allocated(mesh%boundaryNames)) return
		nBoundaries = size(mesh%boundaryNames)
		if (allocated(boundaryConditions)) call freeBoundaryConditions(boundaryConditions)
		allocate(boundaryConditions(nBoundaries))

		print *, 'Boundary condition mapping from OpenFOAM boundary names:'
		do b = 1, nBoundaries
			bname = trim(adjustl(mesh%boundaryNames(b)))
			do while (len_trim(bname) > 0 .and. iachar(bname(1:1)) == 9)
				bname = trim(adjustl(bname(2:)))
			end do
			select case (trim(bname))
			case ('airfoil', 'car', 'wall', 'walls', 'body')
				boundaryConditions(b)%type = wallBoundary
				allocate(boundaryConditions(b)%params(0))
				print *, '  ', trim(bname), ' -> wall'
			case ('empty', 'frontAndBack')
				boundaryConditions(b)%type = emptyBoundary
				allocate(boundaryConditions(b)%params(0))
				print *, '  ', trim(bname), ' -> empty/slip'
			case ('symmetry', 'symmetryPlane')
				! Symmetry/farfield planes are slip boundaries, not no-slip body walls.
				boundaryConditions(b)%type = emptyBoundary
				allocate(boundaryConditions(b)%params(0))
				print *, '  ', trim(bname), ' -> symmetry/slip'
			case ('inlet', 'Inlet')
				boundaryConditions(b)%type = InletBoundary
				allocate(boundaryConditions(b)%params(6))
				boundaryConditions(b)%params = [Pt, Tt, UunitVec(1), UunitVec(2), k_init, omega_init]
				print *, '  ', trim(bname), ' -> inlet'
			case ('outlet', 'Outlet')
				boundaryConditions(b)%type = OutletBoundary
				allocate(boundaryConditions(b)%params(1))
				boundaryConditions(b)%params(1) = P
				print *, '  ', trim(bname), ' -> outlet'
			case default
				boundaryConditions(b)%type = wallBoundary
				allocate(boundaryConditions(b)%params(0))
				print *, '  ', trim(bname), ' -> wall (default)'
			end select
		end do
		call flush(6)
	end subroutine configureBoundaryConditionsFromMesh

	subroutine freeBoundaryConditions(boundaryConditions)
		implicit none
		type(BoundaryCondition), allocatable, intent(inout) :: boundaryConditions(:)
		integer :: b
		if (.not. allocated(boundaryConditions)) return
		do b = 1, size(boundaryConditions)
			if (allocated(boundaryConditions(b)%params)) deallocate(boundaryConditions(b)%params)
		end do
		deallocate(boundaryConditions)
	end subroutine freeBoundaryConditions

!===============================================================================
! SECTION 12: SIMPLEC 压力基求解器（修正D系数 + 面通量版）
!===============================================================================

	subroutine solve_SIMPLE(mesh, meshPath, boundaryConditions, fluid, &
	                        P_init, T_init, U_init, k_init, omega_init)
		implicit none
		type(Meshh), intent(inout) :: mesh
		character(len=*), intent(in) :: meshPath
		type(BoundaryCondition), intent(in) :: boundaryConditions(:)
		type(Fluidd), intent(in) :: fluid
		real(kind=8), intent(in) :: P_init, T_init, U_init(3), k_init, omega_init
		
		type(SIMPLEState) :: ss
		integer(kind=8) :: meshInfo(4)
		integer :: iter, c, f, nInt
		real(kind=8) :: res_u, res_v, res_p, res_max
		real(kind=8), allocatable :: u_old(:), v_old(:), w_old(:), p_old(:)
		integer :: output_count
		real(kind=8) :: fn(3), rho
		
		meshInfo = unstructuredMeshInfo(mesh)
		ss%nCells = int(meshInfo(1))
		ss%nFaces = int(meshInfo(2))
		ss%nBoundaries = int(meshInfo(3))
		ss%nBdryFaces = int(meshInfo(4))
		nInt = ss%nFaces - ss%nBdryFaces
		
		call computeWallDistances(mesh)
		
		ss%rho_ref = P_init / (fluid%R * T_init)
		ss%p_ref = P_init
		rho = ss%rho_ref
		
		allocate(ss%p(ss%nCells), ss%u(ss%nCells), ss%v(ss%nCells), ss%w(ss%nCells))
		allocate(ss%k(ss%nCells), ss%omega(ss%nCells))
		allocate(ss%mu_t(ss%nCells), ss%mu_l(ss%nCells))
		allocate(ss%p_prime(ss%nCells), ss%aP_u(ss%nCells))
		allocate(ss%rho_field(ss%nCells))
		allocate(ss%phi_f(ss%nFaces))
		allocate(u_old(ss%nCells), v_old(ss%nCells), w_old(ss%nCells), p_old(ss%nCells))
		
		ss%p = 0.0d0
		ss%u = U_init(1)
		ss%v = U_init(2)
		ss%w = U_init(3)
		ss%k = k_init
		ss%omega = omega_init
		ss%rho_field = rho
		ss%mu_l = fluid%mu_laminar
		ss%mu_t = rho * k_init / max(omega_init, SMALL_NUM)
		ss%aP_u = 1.0d0
		ss%p_prime = 0.0d0
		
		! 初始化面通量
		do f = 1, ss%nFaces
			fn = mesh%fAVecs(f,:)
			ss%phi_f(f) = rho * dot_product(U_init, fn)
		end do
		
		print *, ''
		print *, '=========================================='
		print *, ' SIMPLEC 压力基求解器'
		print *, '=========================================='
		write(*,'(A,I8,A,I8)') '  单元:', ss%nCells, '  面:', ss%nFaces
		write(*,'(A,F8.4,A,F6.2)') '  rho=', rho, '  U_ref=', mag(U_init)
		write(*,'(A,F8.2,A)') '  动压=', 0.5d0*rho*mag(U_init)**2, ' Pa'
		print *, ''
		
		output_count = 0
		
		do iter = 1, SIMPLE_MAX_ITER
			u_old = ss%u
			v_old = ss%v
			w_old = ss%w
			p_old = ss%p
			
			! 更新湍流粘性
			do c = 1, ss%nCells
				ss%mu_t(c) = rho * max(ss%k(c),SMALL_NUM) / max(ss%omega(c),SMALL_NUM)
				ss%mu_t(c) = min(ss%mu_t(c), 1.0d4 * ss%mu_l(c))
			end do
			
			! 步骤1：求解动量方程
			call simplec_momentum(mesh, ss, boundaryConditions, fluid, U_init)
			
			! 步骤2：压力修正 + 更新面通量和速度
			call simplec_pressure(mesh, ss, boundaryConditions, U_init)
			
			! 步骤3：全局质量修正
			call simple2_mass_fix(mesh, ss, boundaryConditions, rho, U_init)
			
			! 步骤4：湍流更新（延迟启动）
			if (mod(iter, 100) == 0 .and. iter > 500) then
				call simple_turbulence(mesh, ss, fluid)
			end if
			
			! L2残差
			res_u = sqrt(sum((ss%u-u_old)**2)/ss%nCells) / max(mag(U_init), 1.0d-10)
			res_v = sqrt(sum((ss%v-v_old)**2)/ss%nCells) / max(mag(U_init), 1.0d-10)
			res_p = sqrt(sum((ss%p-p_old)**2)/ss%nCells) / max(0.5d0*rho*mag(U_init)**2, 1.0d-10)
			res_max = max(res_u, res_v, res_p)
			
			if (iter <= 5 .or. mod(iter, 50) == 0) then
				write(*,'(A,I6,A,3ES11.3,A,F8.1,A,F8.1)') &
					' SIMPLEC', iter, '  res(u,v,p)=', res_u, res_v, res_p, &
					'  Umax=', maxval(sqrt(ss%u**2+ss%v**2+ss%w**2)), &
					'  dp=', maxval(ss%p)-minval(ss%p)
				call simplec_report_diagnostics(mesh, ss, boundaryConditions, iter)
				call flush(6)
			end if
			
			if (res_max < SIMPLE_TOLERANCE .and. iter > 50) then
				write(*,'(A,I6,A)') ' === SIMPLEC 收敛于', iter, ' 步 ==='
				exit
			end if
			if (any(ss%u /= ss%u) .or. any(ss%p /= ss%p)) then
				print *, 'ERROR: NaN at iter', iter
				exit
			end if
		end do
		
		print *, ''
		print *, '=== SIMPLEC 结果 ==='
		write(*,'(A,F10.1,A,F10.1)') '  dp(gauge): ', minval(ss%p), ' ~ ', maxval(ss%p)
		write(*,'(A,F10.1,A,F10.1)') '  P(abs)   : ', minval(ss%p)+ss%p_ref, ' ~ ', maxval(ss%p)+ss%p_ref
		write(*,'(A,F10.3,A,F10.3)') '  Ux       : ', minval(ss%u), ' ~ ', maxval(ss%u)
		write(*,'(A,F10.3,A,F10.3)') '  Uy       : ', minval(ss%v), ' ~ ', maxval(ss%v)
		write(*,'(A,F10.3)') '  |V| max  : ', maxval(sqrt(ss%u**2+ss%v**2+ss%w**2))
		
		output_count = 1
		call simple_writeVTK(mesh, ss, meshPath, output_count, fluid)
		
		deallocate(ss%p, ss%u, ss%v, ss%w, ss%k, ss%omega)
		deallocate(ss%mu_t, ss%mu_l, ss%p_prime, ss%aP_u, ss%rho_field)
		deallocate(ss%phi_f, u_old, v_old, w_old, p_old)
	end subroutine solve_SIMPLE

!-----------------------------------------------------------------------
! SIMPLEC 动量方程：返回 aP_orig（原始对角系数）和 H_coeff（邻居系数之和）
!-----------------------------------------------------------------------
	subroutine simplec_momentum(mesh, ss, bcs, fluid, U_in)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SIMPLEState), intent(inout) :: ss
		type(BoundaryCondition), intent(in) :: bcs(:)
		type(Fluidd), intent(in) :: fluid
		real(kind=8), intent(in) :: U_in(3)
		
		integer :: c, f, b, fi, face_idx, oc, nc, gs, idx
		integer :: nInt
		real(kind=8) :: rho, mu_e, d_m, a_d, phi_f
		real(kind=8) :: fn(3), fn_m, nb_s, h_c, a_trans
		integer, parameter :: MNB = 30
		real(kind=8) :: aP(ss%nCells), bU(ss%nCells), bV(ss%nCells), bW(ss%nCells)
		real(kind=8) :: H_sum(ss%nCells)  ! 邻居系数之和，用于SIMPLEC
		real(kind=8) :: aN_c(ss%nCells, MNB)
		integer :: aN_i(ss%nCells, MNB), aN_n(ss%nCells)
		real(kind=8), allocatable :: gP(:,:), pm(:,:), tg(:,:,:)
		
		nInt = ss%nFaces - ss%nBdryFaces
		rho = ss%rho_ref
		
		! 压力梯度
		allocate(pm(ss%nCells, 1)); pm(:,1) = ss%p
		allocate(tg(ss%nCells, 1, 3))
		tg = greenGaussGrad_RANS(mesh, pm, .false.)
		allocate(gP(ss%nCells, 3))
		! Use the full pressure gradient in SIMPLEC momentum.  The previous
		! streamwise-only approximation kept v nearly frozen and produced
		! contour plots dominated by a one-dimensional pressure ramp.
		gP(:,1) = tg(:,1,1)
		gP(:,2) = tg(:,1,2)
		gP(:,3) = tg(:,1,3)
		deallocate(tg, pm)
		
		aP = 0.0d0; bU = 0.0d0; bV = 0.0d0; bW = 0.0d0
		H_sum = 0.0d0
		aN_c = 0.0d0; aN_i = 0; aN_n = 0
		
		! 压力梯度源项
		do c = 1, ss%nCells
			bU(c) = -gP(c,1) * mesh%cVols(c)
			bV(c) = -gP(c,2) * mesh%cVols(c)
			bW(c) = -gP(c,3) * mesh%cVols(c)
		end do
		
		! 内部面：用存储的面通量做upwind
		do f = 1, nInt
			oc = mesh%faces(f,1); nc = mesh%faces(f,2)
			if (oc<1.or.oc>ss%nCells.or.nc<1.or.nc>ss%nCells) cycle
			fn = mesh%fAVecs(f,:); fn_m = mag(fn)
			if (fn_m < 1.0d-30) cycle
			d_m = mag(mesh%cCenters(nc,:)-mesh%cCenters(oc,:))
			if (d_m < 1.0d-30) cycle
			
			mu_e = 0.5d0*(ss%mu_l(oc)+ss%mu_t(oc)+ss%mu_l(nc)+ss%mu_t(nc))
			phi_f = ss%phi_f(f)
			a_d = mu_e * fn_m / d_m
			
			! owner侧
			if (aN_n(oc) < MNB) then
				aN_n(oc) = aN_n(oc) + 1
				aN_i(oc, aN_n(oc)) = nc
				aN_c(oc, aN_n(oc)) = a_d + max(-phi_f, 0.0d0)
				H_sum(oc) = H_sum(oc) + a_d + max(-phi_f, 0.0d0)
			end if
			aP(oc) = aP(oc) + a_d + max(phi_f, 0.0d0)
			
			! neighbour侧
			if (aN_n(nc) < MNB) then
				aN_n(nc) = aN_n(nc) + 1
				aN_i(nc, aN_n(nc)) = oc
				aN_c(nc, aN_n(nc)) = a_d + max(phi_f, 0.0d0)
				H_sum(nc) = H_sum(nc) + a_d + max(phi_f, 0.0d0)
			end if
			aP(nc) = aP(nc) + a_d + max(-phi_f, 0.0d0)
		end do
		
		! 边界面
		do b = 1, min(ss%nBoundaries, size(bcs))
			do fi = 1, size(mesh%boundaryFaces, 2)
				face_idx = mesh%boundaryFaces(b, fi)
				if (face_idx == 0) exit
				if (face_idx<1.or.face_idx>ss%nFaces) cycle
				oc = mesh%faces(face_idx, 1)
				if (oc<1.or.oc>ss%nCells) cycle
				fn = mesh%fAVecs(face_idx,:); fn_m = mag(fn)
				if (fn_m < 1.0d-30) cycle
				d_m = max(mag(mesh%fCenters(face_idx,:)-mesh%cCenters(oc,:)), 1.0d-10)
				mu_e = ss%mu_l(oc) + ss%mu_t(oc)
				
				select case(bcs(b)%type)
				case(wallBoundary)
					a_d = mu_e * fn_m / d_m
					aP(oc) = aP(oc) + a_d
				case(InletBoundary)
					phi_f = ss%phi_f(face_idx)
					! Fixed-value inlet: include incoming convective flux in the
					! known-value source instead of treating the inlet as weak diffusion only.
					a_d = mu_e * fn_m / d_m + max(-phi_f, 0.0d0)
					aP(oc) = aP(oc) + a_d
					bU(oc) = bU(oc) + a_d * U_in(1)
					bV(oc) = bV(oc) + a_d * U_in(2)
					bW(oc) = bW(oc) + a_d * U_in(3)
				case(OutletBoundary)
					phi_f = ss%phi_f(face_idx)
					if (phi_f > 0.0d0) then
						aP(oc) = aP(oc) + phi_f
					else
						! Stabilise outlet backflow as a bounded known-value inflow.
						aP(oc) = aP(oc) + max(-phi_f, 0.0d0)
						bU(oc) = bU(oc) + max(-phi_f, 0.0d0) * U_in(1)
						bV(oc) = bV(oc) + max(-phi_f, 0.0d0) * U_in(2)
						bW(oc) = bW(oc) + max(-phi_f, 0.0d0) * U_in(3)
					end if
				case(emptyBoundary)
					continue
				end select
			end do
		end do
		
		! 保存原始aP（用于SIMPLEC的D系数）
		! 伪瞬态对角项：限制小/畸变单元中的压力修正D=V/aP，
		! 同时不直接截断求解后的速度。
		do c = 1, ss%nCells
			h_c = max(mesh%cVols(c)**(1.0d0/3.0d0), 1.0d-8)
			if (allocated(mesh%wallDistance)) h_c = min(h_c, max(mesh%wallDistance(c), 1.0d-8))
			a_trans = rho * mesh%cVols(c) * max(mag(U_in), 1.0d0) / (max(SIMPLEC_PSEUDO_CFL, 1.0d-6) * h_c)
			aP(c) = max(aP(c) + a_trans, 1.0d-20)
			bU(c) = bU(c) + a_trans * ss%u(c)
			bV(c) = bV(c) + a_trans * ss%v(c)
			bW(c) = bW(c) + a_trans * ss%w(c)
			! SIMPLEC pressure-correction coefficient uses the neighbour-stripped
			! momentum diagonal (aP - H).  The floor prevents excessive D on nearly
			! balanced diffusion cells without clipping the solved fields.
			ss%aP_u(c) = max(aP(c) - H_sum(c), SIMPLEC_D_DIAG_FLOOR * aP(c), 1.0d-20)
		end do
		
		! under-relaxation：aP_tilde = aP / alpha_u
		! 注意：这里不再截断历史速度，避免用人为限幅污染物理解。
		do c = 1, ss%nCells
			bU(c) = bU(c) + (1.0d0-SIMPLE_ALPHA_U) * (aP(c)/SIMPLE_ALPHA_U) * ss%u(c)
			bV(c) = bV(c) + (1.0d0-SIMPLE_ALPHA_U) * (aP(c)/SIMPLE_ALPHA_U) * ss%v(c)
			bW(c) = bW(c) + (1.0d0-SIMPLE_ALPHA_U) * (aP(c)/SIMPLE_ALPHA_U) * ss%w(c)
		end do
		
		! GS求解（用 aP_tilde = aP/alpha_u）
		do gs = 1, SIMPLE_MAX_INNER_U
			do c = 1, ss%nCells
				nb_s = 0.0d0
				do idx = 1, aN_n(c)
					nb_s = nb_s + aN_c(c,idx) * ss%u(aN_i(c,idx))
				end do
				ss%u(c) = (bU(c) + nb_s) / (aP(c)/SIMPLE_ALPHA_U)
			end do
			do c = 1, ss%nCells
				nb_s = 0.0d0
				do idx = 1, aN_n(c)
					nb_s = nb_s + aN_c(c,idx) * ss%v(aN_i(c,idx))
				end do
				ss%v(c) = (bV(c) + nb_s) / (aP(c)/SIMPLE_ALPHA_U)
			end do
			do c = 1, ss%nCells
				nb_s = 0.0d0
				do idx = 1, aN_n(c)
					nb_s = nb_s + aN_c(c,idx) * ss%w(aN_i(c,idx))
				end do
				ss%w(c) = (bW(c) + nb_s) / (aP(c)/SIMPLE_ALPHA_U)
			end do
		end do
		
		deallocate(gP)
	end subroutine simplec_momentum

!-----------------------------------------------------------------------
! SIMPLEC 压力修正：PCG求解压力方程并保持Rhie-Chow面通量一致
!-----------------------------------------------------------------------
	subroutine simplec_pressure(mesh, ss, bcs, U_in)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SIMPLEState), intent(inout) :: ss
		type(BoundaryCondition), intent(in) :: bcs(:)
		real(kind=8), intent(in) :: U_in(3)
		
		integer :: c, f, b, fi, face_idx, oc, nc, nInt
		real(kind=8) :: rho, fn(3), fn_m, d_m, D_f, flux_s, fv(3), D_c
		real(kind=8) :: alpha_p_eff, p_step_limit, p_prime_max
		integer, parameter :: MNB = 30
		real(kind=8) :: aP_pp(ss%nCells), src(ss%nCells)
		real(kind=8) :: aN_pp(ss%nCells, MNB)
		integer :: nb_i(ss%nCells, MNB), nb_n(ss%nCells)
		real(kind=8), allocatable :: gP(:,:), gPp(:,:), pm(:,:), tg(:,:,:)
		
		nInt = ss%nFaces - ss%nBdryFaces
		rho = ss%rho_ref
		
		! 当前压力梯度
		allocate(pm(ss%nCells, 1)); pm(:,1) = ss%p
		allocate(tg(ss%nCells, 1, 3))
		tg = greenGaussGrad_RANS(mesh, pm, .false.)
		allocate(gP(ss%nCells, 3))
		! Rhie-Chow interpolation must use the same full pressure gradient as
		! the momentum equation to suppress checkerboard modes in every direction.
		gP(:,1) = tg(:,1,1)
		gP(:,2) = tg(:,1,2)
		gP(:,3) = tg(:,1,3)
		deallocate(tg, pm)
		
		aP_pp = 0.0d0; src = 0.0d0
		aN_pp = 0.0d0; nb_i = 0; nb_n = 0
		ss%p_prime = 0.0d0
		
		! 内部面
		do f = 1, nInt
			oc = mesh%faces(f,1); nc = mesh%faces(f,2)
			if (oc<1.or.oc>ss%nCells.or.nc<1.or.nc>ss%nCells) cycle
			fn = mesh%fAVecs(f,:); fn_m = mag(fn)
			if (fn_m < 1.0d-30) cycle
			d_m = mag(mesh%cCenters(nc,:)-mesh%cCenters(oc,:))
			if (d_m < 1.0d-30) cycle
			
			! D_f 用未松弛的aP_u，不使用aP/alpha_u
			D_f = SIMPLEC_RC_DAMPING * 0.5d0*(mesh%cVols(oc)/ss%aP_u(oc) + mesh%cVols(nc)/ss%aP_u(nc))
			
			! Rhie-Chow 预测面通量
			fv(1) = 0.5d0*(ss%u(oc)+ss%u(nc))
			fv(2) = 0.5d0*(ss%v(oc)+ss%v(nc))
			fv(3) = 0.5d0*(ss%w(oc)+ss%w(nc))
			flux_s = rho * dot_product(fv, fn) &
				- rho * D_f * ((ss%p(nc)-ss%p(oc))/d_m * fn_m &
				  - dot_product(0.5d0*(gP(oc,:)+gP(nc,:)), fn))
			
			! 存储预测面通量
			ss%phi_f(f) = flux_s
			
			! 连续性源项
			src(oc) = src(oc) - flux_s
			src(nc) = src(nc) + flux_s
			
			! 压力方程系数
			D_f = rho * D_f * fn_m / d_m
			aP_pp(oc) = aP_pp(oc) + D_f
			aP_pp(nc) = aP_pp(nc) + D_f
			
			if (nb_n(oc)<MNB) then
				nb_n(oc)=nb_n(oc)+1; nb_i(oc,nb_n(oc))=nc; aN_pp(oc,nb_n(oc))=D_f
			end if
			if (nb_n(nc)<MNB) then
				nb_n(nc)=nb_n(nc)+1; nb_i(nc,nb_n(nc))=oc; aN_pp(nc,nb_n(nc))=D_f
			end if
		end do
		
		! 边界面
		do b = 1, min(ss%nBoundaries, size(bcs))
			do fi = 1, size(mesh%boundaryFaces, 2)
				face_idx = mesh%boundaryFaces(b, fi)
				if (face_idx == 0) exit
				if (face_idx<1.or.face_idx>ss%nFaces) cycle
				oc = mesh%faces(face_idx, 1)
				if (oc<1.or.oc>ss%nCells) cycle
				fn = mesh%fAVecs(face_idx,:); fn_m = mag(fn)
				if (fn_m < 1.0d-30) cycle
				
				select case(bcs(b)%type)
				case(wallBoundary, emptyBoundary)
					ss%phi_f(face_idx) = 0.0d0
				case(InletBoundary)
					ss%phi_f(face_idx) = rho * dot_product(U_in, fn)
					src(oc) = src(oc) - ss%phi_f(face_idx)
				case(OutletBoundary)
					fv = (/ss%u(oc), ss%v(oc), ss%w(oc)/)
					flux_s = rho * dot_product(fv, fn)
					ss%phi_f(face_idx) = flux_s
					src(oc) = src(oc) - flux_s
					! 出口 p'=0 贡献
					d_m = max(mag(mesh%fCenters(face_idx,:)-mesh%cCenters(oc,:)), 1.0d-10)
					D_f = SIMPLEC_RC_DAMPING * mesh%cVols(oc) / ss%aP_u(oc)
					D_f = rho * D_f * fn_m / d_m
					aP_pp(oc) = aP_pp(oc) + D_f
				end select
			end do
		end do
		
		do c = 1, ss%nCells
			if (aP_pp(c) < 1.0d-30) aP_pp(c) = 1.0d0
		end do
		
		! PCG 求解 p'。压力矩阵是 aP*p' - sum(aN*p'_nb) = src，
		! 对角预处理共轭梯度比逐点GS更快消除全局压力误差。
		call simplec_pressure_pcg(ss%nCells, MNB, aP_pp, aN_pp, nb_i, nb_n, src, &
		                          ss%p_prime, SIMPLE_MAX_INNER_P, SIMPLE_PCG_TOLERANCE)
		
		! SIMPLEC pressure update uses adaptive equation under-relaxation.
		! This limits a single pressure-correction step relative to dynamic pressure
		! without clipping the accumulated pressure field.
		p_step_limit = SIMPLEC_MAX_PRESSURE_STEP_FACTOR * max(0.5d0*rho*mag(U_in)**2, 1.0d0)
		p_prime_max = maxval(abs(ss%p_prime))
		alpha_p_eff = SIMPLE_ALPHA_P
		if (SIMPLE_ALPHA_P * p_prime_max > p_step_limit) then
			alpha_p_eff = p_step_limit / max(p_prime_max, 1.0d-30)
		end if
		ss%p = ss%p + alpha_p_eff * ss%p_prime
		
		! 减去均值防漂移
		ss%p = ss%p - sum(ss%p) / dble(ss%nCells)
		
		! p'梯度用于速度修正
		allocate(pm(ss%nCells, 1)); pm(:,1) = ss%p_prime
		allocate(tg(ss%nCells, 1, 3))
		tg = greenGaussGrad_RANS(mesh, pm, .false.)
		allocate(gPp(ss%nCells, 3))
		! Correct every velocity component with the full pressure-correction gradient.
		gPp(:,1) = tg(:,1,1)
		gPp(:,2) = tg(:,1,2)
		gPp(:,3) = tg(:,1,3)
		deallocate(tg, pm)
		
		! 速度修正：用未松弛的 aP_u
		do c = 1, ss%nCells
			D_c = SIMPLEC_RC_DAMPING * mesh%cVols(c) / ss%aP_u(c)
			ss%u(c) = ss%u(c) - alpha_p_eff * D_c * gPp(c,1)
			ss%v(c) = ss%v(c) - alpha_p_eff * D_c * gPp(c,2)
			ss%w(c) = ss%w(c) - alpha_p_eff * D_c * gPp(c,3)
		end do
		
		! 面通量增量修正
		do f = 1, nInt
			oc = mesh%faces(f,1); nc = mesh%faces(f,2)
			if (oc<1.or.oc>ss%nCells.or.nc<1.or.nc>ss%nCells) cycle
			fn = mesh%fAVecs(f,:); fn_m = mag(fn)
			if (fn_m < 1.0d-30) cycle
			d_m = mag(mesh%cCenters(nc,:)-mesh%cCenters(oc,:))
			if (d_m < 1.0d-30) cycle
			D_f = SIMPLEC_RC_DAMPING * 0.5d0*(mesh%cVols(oc)/ss%aP_u(oc) + mesh%cVols(nc)/ss%aP_u(nc))
			D_f = rho * D_f * fn_m / d_m
			ss%phi_f(f) = ss%phi_f(f) - alpha_p_eff * D_f * (ss%p_prime(nc) - ss%p_prime(oc))
		end do
		
		! 出口面通量更新
		do b = 1, min(ss%nBoundaries, size(bcs))
			if (bcs(b)%type /= OutletBoundary) cycle
			do fi = 1, size(mesh%boundaryFaces, 2)
				face_idx = mesh%boundaryFaces(b, fi)
				if (face_idx == 0) exit
				if (face_idx<1.or.face_idx>ss%nFaces) cycle
				oc = mesh%faces(face_idx, 1)
				if (oc<1.or.oc>ss%nCells) cycle
				fv = (/ss%u(oc), ss%v(oc), ss%w(oc)/)
				ss%phi_f(face_idx) = rho * dot_product(fv, mesh%fAVecs(face_idx,:))
			end do
		end do
		
		deallocate(gP, gPp)
	end subroutine simplec_pressure

!-----------------------------------------------------------------------
! 对角预处理共轭梯度求解压力修正方程
!-----------------------------------------------------------------------
	subroutine simplec_pressure_pcg(nCells, maxNb, aP, aN, nb_i, nb_n, rhs, x, maxIter, relTol)
		implicit none
		integer, intent(in) :: nCells, maxNb, maxIter
		real(kind=8), intent(in) :: aP(nCells), aN(nCells, maxNb), rhs(nCells), relTol
		integer, intent(in) :: nb_i(nCells, maxNb), nb_n(nCells)
		real(kind=8), intent(inout) :: x(nCells)
		integer :: iter, c
		real(kind=8), allocatable :: r(:), z(:), pvec(:), Avec(:)
		real(kind=8) :: alpha_cg, beta_cg, rz_old, rz_new, pAp, bnorm, rnorm
		real(kind=8) :: diag

		allocate(r(nCells), z(nCells), pvec(nCells), Avec(nCells))

		call simplec_apply_pressure_matrix(nCells, maxNb, aP, aN, nb_i, nb_n, x, Avec)
		r = rhs - Avec
		bnorm = sqrt(sum(rhs*rhs))
		if (bnorm < 1.0d-30) then
			x = 0.0d0
			deallocate(r, z, pvec, Avec)
			return
		end if

		do c = 1, nCells
			diag = max(aP(c), 1.0d-30)
			z(c) = r(c) / diag
		end do
		pvec = z
		rz_old = sum(r*z)

		do iter = 1, maxIter
			call simplec_apply_pressure_matrix(nCells, maxNb, aP, aN, nb_i, nb_n, pvec, Avec)
			pAp = sum(pvec*Avec)
			if (abs(pAp) < 1.0d-300) exit
			alpha_cg = rz_old / pAp
			x = x + alpha_cg * pvec
			r = r - alpha_cg * Avec
			rnorm = sqrt(sum(r*r))
			if (rnorm <= relTol * bnorm) exit
			do c = 1, nCells
				diag = max(aP(c), 1.0d-30)
				z(c) = r(c) / diag
			end do
			rz_new = sum(r*z)
			if (abs(rz_old) < 1.0d-300) exit
			beta_cg = rz_new / rz_old
			pvec = z + beta_cg * pvec
			rz_old = rz_new
		end do

		deallocate(r, z, pvec, Avec)
	end subroutine simplec_pressure_pcg

	subroutine simplec_apply_pressure_matrix(nCells, maxNb, aP, aN, nb_i, nb_n, x, Ax)
		implicit none
		integer, intent(in) :: nCells, maxNb
		real(kind=8), intent(in) :: aP(nCells), aN(nCells, maxNb), x(nCells)
		integer, intent(in) :: nb_i(nCells, maxNb), nb_n(nCells)
		real(kind=8), intent(out) :: Ax(nCells)
		integer :: c, idx

		do c = 1, nCells
			Ax(c) = aP(c) * x(c)
			do idx = 1, nb_n(c)
				Ax(c) = Ax(c) - aN(c,idx) * x(nb_i(c,idx))
			end do
		end do
	end subroutine simplec_apply_pressure_matrix

!-----------------------------------------------------------------------
! SIMPLEC 诊断：定位压力/速度极值和各边界状态，避免只看全局极值
!-----------------------------------------------------------------------
	subroutine simplec_report_diagnostics(mesh, ss, bcs, iter)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SIMPLEState), intent(in) :: ss
		type(BoundaryCondition), intent(in) :: bcs(:)
		integer, intent(in) :: iter
		integer :: c, b, fi, face_idx, oc, pmin_cell, pmax_cell, umax_cell, nFacesB
		integer :: n_wall, n_front, n_roof, n_rear, n_wake, n_upstream
		real(kind=8) :: speed, pmin_val, pmax_val, umax_val, p_range
		real(kind=8) :: b_pmin, b_pmax, b_psum, b_umax, b_flux
		real(kind=8) :: v_abs_max, v_rms, x_min, x_max, y_min, y_max, xw_min, xw_max, yw_min, yw_max
		real(kind=8) :: wall_dx, wall_dy, front_p, roof_p, rear_p, wake_u, upstream_u, wake_deficit
		logical, allocatable :: wall_cell(:)
		character(len=100) :: bname

		allocate(wall_cell(ss%nCells)); wall_cell = .false.
		pmin_cell = 1; pmax_cell = 1; umax_cell = 1
		pmin_val = ss%p(1); pmax_val = ss%p(1)
		umax_val = sqrt(ss%u(1)**2 + ss%v(1)**2 + ss%w(1)**2)
		v_abs_max = abs(ss%v(1)); v_rms = ss%v(1)**2
		x_min = mesh%cCenters(1,1); x_max = mesh%cCenters(1,1)
		y_min = mesh%cCenters(1,2); y_max = mesh%cCenters(1,2)
		do c = 2, ss%nCells
			x_min = min(x_min, mesh%cCenters(c,1)); x_max = max(x_max, mesh%cCenters(c,1))
			y_min = min(y_min, mesh%cCenters(c,2)); y_max = max(y_max, mesh%cCenters(c,2))
			v_abs_max = max(v_abs_max, abs(ss%v(c)))
			v_rms = v_rms + ss%v(c)**2
			if (ss%p(c) < pmin_val) then
				pmin_val = ss%p(c); pmin_cell = c
			end if
			if (ss%p(c) > pmax_val) then
				pmax_val = ss%p(c); pmax_cell = c
			end if
			speed = sqrt(ss%u(c)**2 + ss%v(c)**2 + ss%w(c)**2)
			if (speed > umax_val) then
				umax_val = speed; umax_cell = c
			end if
		end do
		v_rms = sqrt(v_rms / dble(ss%nCells))

		p_range = pmax_val - pmin_val
		write(*,'(A,I6,A,I8,A,3F10.4,A,I8,A,3F10.4,A,I8,A,F10.4,A,ES11.3)') &
			'  diag', iter, ' pMinCell=', pmin_cell, ' xyz=', mesh%cCenters(pmin_cell,:), &
			' pMaxCell=', pmax_cell, ' xyz=', mesh%cCenters(pmax_cell,:), &
			' uMaxCell=', umax_cell, ' |U|=', umax_val, ' rawDp=', p_range
		if (p_range > SIMPLE_MAX_PRESSURE_RANGE .or. umax_val > SIMPLE_MAX_SPEED) then
			write(*,'(A,ES11.3,A,F10.3,A)') '  diag warning: rawDp=', p_range, &
				' Umax=', umax_val, ' exceeds diagnostic target; not clipped.'
		end if

		do b = 1, min(ss%nBoundaries, size(bcs))
			nFacesB = 0; b_psum = 0.0d0; b_umax = 0.0d0; b_flux = 0.0d0
			b_pmin = huge(1.0d0); b_pmax = -huge(1.0d0)
			bname = 'boundary'
			if (allocated(mesh%boundaryNames) .and. b <= size(mesh%boundaryNames)) then
				bname = trim(adjustl(mesh%boundaryNames(b)))
			end if
			do fi = 1, size(mesh%boundaryFaces, 2)
				face_idx = mesh%boundaryFaces(b, fi)
				if (face_idx == 0) exit
				if (face_idx<1 .or. face_idx>ss%nFaces) cycle
				oc = mesh%faces(face_idx, 1)
				if (oc<1 .or. oc>ss%nCells) cycle
				nFacesB = nFacesB + 1
				b_psum = b_psum + ss%p(oc)
				b_pmin = min(b_pmin, ss%p(oc)); b_pmax = max(b_pmax, ss%p(oc))
				speed = sqrt(ss%u(oc)**2 + ss%v(oc)**2 + ss%w(oc)**2)
				b_umax = max(b_umax, speed)
				b_flux = b_flux + ss%phi_f(face_idx)
			end do
			if (nFacesB > 0) then
				write(*,'(A,I2,1X,A,A,ES11.3,A,ES11.3,A,ES11.3,A,F8.3,A,ES11.3)') &
					'    bc', b, trim(bname), ' pAvg=', b_psum/dble(nFacesB), &
					' pMin=', b_pmin, ' pMax=', b_pmax, ' uMax=', b_umax, ' flux=', b_flux
			end if
		end do


		! Spatial sanity diagnostics for contour plots.  A physically meaningful
		! car/airfoil low-Mach solution should not look like a pure one-dimensional
		! inlet-to-outlet ramp: wall sectors and wake/upstream samples should expose
		! stagnation, roof/rear suction and wake velocity deficit trends.
		n_wall = 0
		xw_min = huge(1.0d0); xw_max = -huge(1.0d0)
		yw_min = huge(1.0d0); yw_max = -huge(1.0d0)
		do b = 1, min(ss%nBoundaries, size(bcs))
			if (bcs(b)%type /= wallBoundary) cycle
			bname = 'boundary'
			if (allocated(mesh%boundaryNames) .and. b <= size(mesh%boundaryNames)) then
				bname = trim(adjustl(mesh%boundaryNames(b)))
			end if
			! Symmetry/slip patches are not the car/airfoil body for contour-physics
			! sector diagnostics.
			if (index(bname, 'symmetry') > 0) cycle
			do fi = 1, size(mesh%boundaryFaces, 2)
				face_idx = mesh%boundaryFaces(b, fi)
				if (face_idx == 0) exit
				if (face_idx<1 .or. face_idx>ss%nFaces) cycle
				oc = mesh%faces(face_idx, 1)
				if (oc<1 .or. oc>ss%nCells) cycle
				if (.not. wall_cell(oc)) then
					wall_cell(oc) = .true.; n_wall = n_wall + 1
					xw_min = min(xw_min, mesh%cCenters(oc,1)); xw_max = max(xw_max, mesh%cCenters(oc,1))
					yw_min = min(yw_min, mesh%cCenters(oc,2)); yw_max = max(yw_max, mesh%cCenters(oc,2))
				end if
			end do
		end do

		front_p = 0.0d0; roof_p = 0.0d0; rear_p = 0.0d0; wake_u = 0.0d0; upstream_u = 0.0d0
		n_front = 0; n_roof = 0; n_rear = 0; n_wake = 0; n_upstream = 0
		if (n_wall > 0) then
			wall_dx = max(xw_max - xw_min, 1.0d-12)
			wall_dy = max(yw_max - yw_min, 1.0d-12)
			do c = 1, ss%nCells
				if (wall_cell(c)) then
					if (mesh%cCenters(c,1) <= xw_min + 0.15d0*wall_dx) then
						front_p = front_p + ss%p(c); n_front = n_front + 1
					end if
					if (mesh%cCenters(c,2) >= yw_min + 0.70d0*wall_dy) then
						roof_p = roof_p + ss%p(c); n_roof = n_roof + 1
					end if
					if (mesh%cCenters(c,1) >= xw_max - 0.15d0*wall_dx) then
						rear_p = rear_p + ss%p(c); n_rear = n_rear + 1
					end if
				end if
				if (mesh%cCenters(c,1) > xw_max + 0.05d0*wall_dx .and. &
				    mesh%cCenters(c,2) >= yw_min .and. mesh%cCenters(c,2) <= yw_max) then
					wake_u = wake_u + ss%u(c); n_wake = n_wake + 1
				end if
				if (mesh%cCenters(c,1) < xw_min - 0.05d0*wall_dx .and. &
				    mesh%cCenters(c,2) >= yw_min .and. mesh%cCenters(c,2) <= yw_max) then
					upstream_u = upstream_u + ss%u(c); n_upstream = n_upstream + 1
				end if
			end do
			if (n_front > 0) front_p = front_p / dble(n_front)
			if (n_roof > 0) roof_p = roof_p / dble(n_roof)
			if (n_rear > 0) rear_p = rear_p / dble(n_rear)
			if (n_wake > 0) wake_u = wake_u / dble(n_wake)
			if (n_upstream > 0) upstream_u = upstream_u / dble(n_upstream)
			wake_deficit = upstream_u - wake_u
			write(*,'(A,I6,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3)') &
				'  phys', iter, ' vMax=', v_abs_max, ' vRms=', v_rms, ' frontP=', front_p, &
				' roofP=', roof_p, ' rearP=', rear_p, ' wakeUx=', wake_u, ' wakeDef=', wake_deficit
		end if

		deallocate(wall_cell)
	end subroutine simplec_report_diagnostics

!-----------------------------------------------------------------------
! 全局质量通量修正
!-----------------------------------------------------------------------
	subroutine simple2_mass_fix(mesh, ss, bcs, rho, U_in)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SIMPLEState), intent(inout) :: ss
		type(BoundaryCondition), intent(in) :: bcs(:)
		real(kind=8), intent(in) :: rho, U_in(3)
		integer :: b, fi, face_idx, oc
		real(kind=8) :: fn(3), flux_in, flux_out, corr
		
		flux_in = 0.0d0; flux_out = 0.0d0
		do b = 1, min(ss%nBoundaries, size(bcs))
			do fi = 1, size(mesh%boundaryFaces, 2)
				face_idx = mesh%boundaryFaces(b, fi)
				if (face_idx == 0) exit
				if (face_idx<1.or.face_idx>ss%nFaces) cycle
				select case(bcs(b)%type)
				case(InletBoundary)
					flux_in = flux_in + ss%phi_f(face_idx)
				case(OutletBoundary)
					flux_out = flux_out + ss%phi_f(face_idx)
				end select
			end do
		end do
		
		if (abs(flux_out) > 1.0d-20) then
			corr = -flux_in / flux_out
			if (corr > 0.5d0 .and. corr < 2.0d0) then
				do b = 1, min(ss%nBoundaries, size(bcs))
					if (bcs(b)%type /= OutletBoundary) cycle
					do fi = 1, size(mesh%boundaryFaces, 2)
						face_idx = mesh%boundaryFaces(b, fi)
						if (face_idx == 0) exit
						if (face_idx<1.or.face_idx>ss%nFaces) cycle
						oc = mesh%faces(face_idx, 1)
						if (oc<1.or.oc>ss%nCells) cycle
						ss%phi_f(face_idx) = ss%phi_f(face_idx) * corr
					end do
				end do
			end if
		end if
	end subroutine simple2_mass_fix

!-----------------------------------------------------------------------
! 湍流更新
!-----------------------------------------------------------------------
	subroutine simple_turbulence(mesh, ss, fluid)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SIMPLEState), intent(inout) :: ss
		type(Fluidd), intent(in) :: fluid
		integer :: c, i, j
		real(kind=8) :: rho, S_mag, P_k, k_new, omega_new
		real(kind=8), allocatable :: gradU(:,:,:), vel(:,:)
		real(kind=8) :: gU(3,3), S(3,3)
		
		rho = ss%rho_ref
		allocate(vel(ss%nCells, 3))
		vel(:,1)=ss%u; vel(:,2)=ss%v; vel(:,3)=ss%w
		allocate(gradU(ss%nCells, 3, 3))
		gradU = greenGaussGrad_RANS(mesh, vel, .false.)
		deallocate(vel)
		
		do c = 1, ss%nCells
			gU = gradU(c,:,:)
			S_mag = 0.0d0
			do i = 1,3; do j = 1,3
				S(i,j) = 0.5d0*(gU(i,j)+gU(j,i))
				S_mag = S_mag + S(i,j)**2
			end do; end do
			S_mag = sqrt(2.0d0 * S_mag)
			
			P_k = ss%mu_t(c)*S_mag**2
			P_k = min(P_k, 10.0d0*BETA_STAR*rho*ss%omega(c)*max(ss%k(c),SMALL_NUM))
			k_new = max(P_k / max(BETA_STAR*rho*ss%omega(c),SMALL_NUM), SMALL_NUM)
			ss%k(c) = max(ss%k(c) + SIMPLE_ALPHA_K*(k_new-ss%k(c)), SMALL_NUM)
			omega_new = max(ALPHA*S_mag**2 / max(BETA*ss%omega(c),SMALL_NUM), SMALL_NUM)
			ss%omega(c) = max(ss%omega(c) + SIMPLE_ALPHA_OMEGA*(omega_new-ss%omega(c)), SMALL_NUM)
		end do
		deallocate(gradU)
	end subroutine simple_turbulence

!-----------------------------------------------------------------------
! VTK 输出
!-----------------------------------------------------------------------
	subroutine simple_writeVTK(mesh, ss, meshPath, output_count, fluid)
		implicit none
		type(Meshh), intent(in) :: mesh
		type(SIMPLEState), intent(in) :: ss
		character(len=*), intent(in) :: meshPath
		integer, intent(in) :: output_count
		type(Fluidd), intent(in) :: fluid
		real(kind=8), allocatable :: prims(:,:)
		type(SolutionState) :: st
		character(len=256) :: vn
		allocate(prims(ss%nCells, 7))
		prims(:,1)=ss%p+ss%p_ref
		prims(:,2)=(ss%p+ss%p_ref)/(fluid%R*ss%rho_ref)
		prims(:,3)=ss%u; prims(:,4)=ss%v; prims(:,5)=ss%w
		prims(:,6)=ss%k; prims(:,7)=ss%omega
		allocate(st%cellMuT(ss%nCells), st%cellYplus(ss%nCells))
		st%cellMuT=ss%mu_t; st%cellYplus=0.0d0
		write(vn,'(A,I0)') 'solution_SIMPLE.', output_count
		print *, 'Writing ', trim(vn)
		call outputVTK_RANS(meshPath, prims, st, vn)
		deallocate(prims, st%cellMuT, st%cellYplus)
	end subroutine simple_writeVTK

! SECTION 13: TURKEL LOW-MACH PRECONDITIONING (for density-based solver)
!===============================================================================

	function preconditioned_spectral_radius(V_mag, a, gamma) result(rho_p)
		implicit none
		real(kind=8), intent(in) :: V_mag, a, gamma
		real(kind=8) :: rho_p
		real(kind=8) :: M_local, M_ref, alpha2, cprime
		
		if (.not. LOW_MACH_PRECOND) then
			rho_p = V_mag + a
			return
		end if
		
		M_local = V_mag / max(a, SMALL_NUM)
		M_ref = max(M_local, M_CUTOFF)
		M_ref = min(M_ref, 1.0d0)
		alpha2 = 1.0d0 - M_ref**2
		cprime = sqrt(((1.0d0 - alpha2) * V_mag)**2 + 4.0d0 * alpha2 * a**2)
		rho_p = 0.5d0 * ((1.0d0 + alpha2) * V_mag + cprime)
	end function preconditioned_spectral_radius

end module BuFlowModule

