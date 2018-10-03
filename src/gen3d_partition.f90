!**********************************************
!*gen3d_partition.F                           *
!*                                            *
!*Parition the mesh                           *
!*                                            * 
!*   Partit v 1.0.0                           * 
!*                                            *
!**********************************************
!*Description:                                *
!* File controls the datastructure            *
!* partioning for the pre-processor.          *
!**********************************************

    module partition

! Include the definitions of other modules. 
! These could be changed if other types are 
! wanted
    use BoundaryRegister ! handles boundary faces
    use TetrahedralElements ! handles elements
    use PrismElements
    use PyramidElements
    use HexahedralElements
    use CoordinateRegister ! provides coordinates for nodes in grid
    use Grid ! handles a grid   

    IMPLICIT NONE ! to disallow implicit variable declarations

include "mpif.h"
!    use mpi

    integer :: nparts        ! number of partitions 

    character(len=32)  :: arg
    character(len=100) :: infileName
    character(len=100) :: outfileName
    logical :: FILEEXISTS, ISOPEN


    integer,parameter:: NSD = 3 ! number of space dimensions

    type(CoordinateRegisterData),pointer :: crp
    type(CoordinateRegisterData),pointer :: crpPrev
    type(CoordinateRegisterData),pointer :: crpPrev2
    type(TetrahedralElementData),pointer :: tep
    type(PrismElementData),pointer :: prp 
    type(PyramidElementData),pointer :: pyp 
    type(HexahedralElementData),pointer :: hep
    type(GridData),pointer :: grp(:)
    type(BoundaryRegisterData), pointer :: brp

    character :: problemName*80

    integer :: numberOfNodes,numberOfElements,numberOfSurfaceSegments,numberOfLineSegments
    integer :: numberOfBoundaryNodes,numberOfBoundarySides,numberOfBoundaryFaces
    integer :: numberOfTriangularBFs,numberOfRectangularBFs,numberOfSurfaceNormals
    integer :: numberOfBoundarySurfaceBlocks,numberOfBoundarySegmentBlocks                         
    integer :: numberOfHexEl,numberOfPriEl,numberOfPyrEl,numberOfTetEl
    integer :: numberOfGrids ! number of grids to produce (1: just fine grid)
    integer :: numberOfSides

    integer :: numberOfDomains ! for parallelization

    integer :: numberOfCycles,timestepNumber,stepInCycle,numberOfStepsPerCycle
    integer :: inputNumber,inputNumber2

    integer :: problemNameLength 
    integer :: InvVis
    real :: directionalityParameter,minimumAspectRatio,groundAngle
    real :: normalSmoothingFactor
    real :: wallDistanceLimit
    logical :: doVisualization,readPart,mergeLonesomeNodes,isHybrid,isRollGround
    integer :: visualizationMode,numberOfDirectionalAggl,normalSmoothingIterations

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: coords
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: elemNodeConn

    ! for METIS partitioning
    INTEGER, pointer     :: vwgt=>null(), vsize=>null()
    INTEGER :: options_metis(40)

    INTEGER :: ee, ii, jj, kk, ind, objval, errpar
    INTEGER :: n1, n2, n3, n4, n5, n6, n7, n8

      INTEGER, DIMENSION(:), ALLOCATABLE :: elem_proc_id
      INTEGER, DIMENSION(:), ALLOCATABLE :: elem_proc_id_local
      INTEGER, DIMENSION(:), ALLOCATABLE :: node_proc_id
      INTEGER, DIMENSION(:), ALLOCATABLE :: elemdist
      INTEGER, DIMENSION(:), ALLOCATABLE :: eptr
      INTEGER, DIMENSION(:), ALLOCATABLE :: eind
      INTEGER, DIMENSION(:), ALLOCATABLE :: elmwgt
      INTEGER, DIMENSION(:), ALLOCATABLE :: intVecTemp
      INTEGER, DIMENSION(:), ALLOCATABLE :: recvcounts
      INTEGER, DIMENSION(:), ALLOCATABLE :: displs

      !DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: tpwgts
      REAL*4, DIMENSION(:), ALLOCATABLE :: tpwgts
      REAL*4, DIMENSION(:), ALLOCATABLE :: ubvec
      DOUBLE PRECISION :: fact

      INTEGER :: wgtflag, numflag, ncon, ncommonnodes, edgecut

      INTEGER ::  n_mpi_procs     ! total number of processors in the group
      INTEGER ::  this_mpi_proc   ! rank of the current processor

      INTEGER ::  nNode_global   ! number of nodes in the whole mesh
      INTEGER ::  nNode_local    ! number of nodes owned by the local processor
      INTEGER ::  nElem_global   ! number of global elements (in the whole model)
      INTEGER ::  nElem_local    ! number of local  elements (owned by the local processor)
      INTEGER ::  npElem         ! number of nodes per element
      INTEGER ::  elem_num_start ! starting element in the current processor
      INTEGER ::  elem_num_end   ! end element in the current processor
    INTEGER :: metisType     ! metis algorithm type - nodal/dual
    INTEGER :: io            ! for checking input/output (from files)


    contains

!--------------------------------------------------------------------------
    subroutine main()
    ! runs the show, calls other functions

    call MPI_Comm_size(MPI_COMM_WORLD, n_mpi_procs, errpar);
    call MPI_Comm_rank(MPI_COMM_WORLD, this_mpi_proc, errpar);

    call communicate() ! communicates with user - calls procedures to get and
                     ! set up data from file
    call readFineGrid() ! reads fine grid from file and prepares data structure

    call partitionFineGrid() ! prepares fine grid data structure for agglomeration

    call writevtk() ! prepares fine grid data structure for agglomeration

    call cleanup()

    end subroutine main

!--------------------------------------------------------------------------
    subroutine communicate()

      IMPLICIT NONE
 
      integer :: j,allocateStatus,RUNFILE

      write(*,*) "" 
      write(*,'(A)') "************************************************"
      write(*,'(A)') "    WELCOME TO PARTIT - Mesh Partioner    "  
      write(*,'(A)') "************************************************"
      write(*,*) "" 

      !Set file names
      !The file names are specified as inputs from the command line
      IF( iargc() < 2 ) THEN
        WRITE(*,*) "Number of input entries is not sufficient "
        STOP "Aborting..."
      END IF

      CALL getarg(1, arg)
      READ(arg,*) problemName
      problemNameLength = nameLen(problemName)

      CALL getarg(2, arg)
      READ(arg,*) nparts

      !CALL getarg(3, arg)
      !READ(arg,*) isHybrid

    end subroutine communicate  

!-------------------------------------------------------------------------- 
!-------------------------------------------------------------------------- 
    subroutine readFineGrid()
! communicates with user and reads fine grid

    IMPLICIT NONE

    integer, parameter :: SURFINFILE = 25,VOLINFILE = 26,SEGINFILE = 27,WALLFILE = 28
    integer :: i,j,ip,allocateStatus,dummy,numberOfBoundaryFacesSurf
    real :: dummyr
    character :: text*10

    integer :: DCOORFILE,DCOORFILE2,itLen,itLen2
    logical :: coorMov,readWALLFile
    character :: itExt*10,itForm*4,itExt2*10,itForm2*4
  
    infileName = trim(problemName)//'.plt'

    write(*,*) "Reading input file ", infileName
  
    INQUIRE(file=infileName, EXIST=FILEEXISTS)
    WRITE(*,*) " FILEEXISTS = ", FILEEXISTS
    IF(FILEEXISTS .NEQV. .TRUE.) THEN
      WRITE(*,*) "File ... ", infileName, "does not exist"
      call EXIT(1)
    END IF

    open(VOLINFILE,file=problemName(1:problemNameLength)//'.plt',form='unformatted',status='old')
    write(*,*) 'File '//problemName(1:problemNameLength)//'.plt'//' opened...' 

    DCOORFILE = -1
    DCOORFILE2= -1

    ! Start reading files

    isHybrid = .false.

    if(isHybrid) then 
      read(VOLINFILE) numberOfElements, numberOfNodes, numberOfBoundaryFaces,&
                      numberOfHexEl,numberOfPriEl,numberOfPyrEl,numberOfTetEl,&
                      numberOfRectangularBFs,numberOfTriangularBFs,numberOfSides 
    else
      read(VOLINFILE) numberOfElements, numberOfNodes, numberOfBoundaryFaces
      numberOfTetEl = numberOfElements
      numberOfTriangularBFs = numberOfBoundaryFaces
      numberOfRectangularBFs = 0
      numberOfHexEl = 0
      numberOfPriEl = 0
      numberOfPyrEl = 0
    end if


    write(*,*) "" 
    write(*,'(''Elements         :'',i8)') numberOfElements  
    write(*,'(''Nodes            :'',i8)') numberOfNodes 
    write(*,'(''Boundary faces   :'',i8)') numberOfBoundaryFaces 
    write(*,'(''Hexahedra        :'',i8)') numberOfHexEl
    write(*,'(''Prisms           :'',i8)') numberOfPriEl
    write(*,'(''Pyramids         :'',i8)') numberOfPyrEl
    write(*,'(''Tetrahedra       :'',i8)') numberOfTetEl
    write(*,'(''Triangular faces :'',i8)') numberOfTriangularBFs 
    write(*,'(''Rectangular faces:'',i8)') numberOfRectangularBFs 
    write(*,'(''Number of edges  :'',i8)') numberOfSides
    write(*,*) ""
  
    ! initiates data structure for input

    ALLOCATE(coords(numberOfNodes,3))

    npElem=4
    ALLOCATE(elemNodeConn(numberOfTetEl,npElem))

    ! read element data 
    write(*,*) "Reading element data..." 

    !DO ii=1,numberOfTetEl
      !READ(VOLINFILE,*, iostat=io) elemNodeConn(ii,1), elemNodeConn(ii,2), elemNodeConn(ii,3), elemNodeConn(ii,4)
    !END DO
    !read(VOLINFILE) ((tep%elements(ip)%pointIndexes(i),ip=1,nel),i=1,4)

    read(VOLINFILE) ((elemNodeConn(i,j),i=1,numberOfTetEl),j=1,4)

    ! read coordinate data 
    write(*,*) "Reading coordinate data..."
    !call readCoordinateData(VOLINFILE,crp) ! calls point register module to read data 

    !do ii=1,numberOfNodes
      !read(VOLINFILE,*, iostat=io) coords(ii,1), coords(ii,2), coords(ii,3)
    !end do
    !read(INFILE) ((crp%points(ip,i),ip=1,NP),i=1,NSD)

    read(VOLINFILE) ((coords(i,j),i=1,numberOfNodes),j=1,3)

    write(*,*) " "
    write(*,*) "Finished reading the input file ..."
    write(*,*) " "
  end subroutine readFineGrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine partitionFineGrid

    IMPLICIT NONE



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !! Partition the mesh. Here METIS is used.
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    WRITE(*,*) " Partitioning the mesh "

      nElem_global = numberOfElements
      nNode_global = numberOfNodes

      write(*,*) "nElem_global    = ", nElem_global, this_mpi_proc, n_mpi_procs

      ALLOCATE(elemdist(n_mpi_procs+1))
      elemdist = 0
      ALLOCATE(elem_proc_id(nElem_global))
      elem_proc_id = 0
      ALLOCATE(node_proc_id(nNode_global))
      node_proc_id = 0

      fact = nElem_global/n_mpi_procs
      ind = ceiling(fact)

      elem_num_start = ind*this_mpi_proc+1

      if(this_mpi_proc .EQ. (n_mpi_procs-1) ) then
        elem_num_end   = nElem_global
      else
        elem_num_end   = ind*(this_mpi_proc+1)
      end if

      nElem_local = elem_num_end - elem_num_start + 1
      ALLOCATE(elem_proc_id_local(nElem_local))
      elem_proc_id_local = 0

      write(*,*) "elem_num_start = ", elem_num_start, this_mpi_proc, n_mpi_procs
      write(*,*) "elem_num_end   = ", elem_num_end, this_mpi_proc, n_mpi_procs
      write(*,*) "nElem_local    = ", nElem_local, this_mpi_proc, n_mpi_procs

      do ii=1,n_mpi_procs
        n1 = ind*(ii-1)+1
        n2 = ind*ii

        if(ii .EQ. n_mpi_procs ) then
          n2  = nElem_global
        end if

        elemdist(ii+1) = elemdist(ii) + n2-n1+1
      end do

      write(*,*) elemdist

      call MPI_Barrier(MPI_COMM_WORLD, errpar)
      WRITE(*,*) " Building eptr & eind arrays ", this_mpi_proc
      ! array of size (nElem_local+1) which stores
      ! number of nodes per each element,
      ! with 0 as the first entry and
      ! nElem_local*npElem as the last entry
      ind = nElem_local+1
      ALLOCATE( eptr(ind) )

      ! array of size 'nElem_local*npElem'
      ! which stores eleme<->node connectivity information
      npElem=4
      ind = nElem_local*npElem
      ALLOCATE( eind(ind) )

      ALLOCATE( intVecTemp(npElem) )

      ee=elem_num_start
      DO ii=1, nElem_local
        kk = (ii-1)*npElem

        eptr(ii) = kk

        intVecTemp = elemNodeConn(ee,:)

        DO jj=1, npElem
          eind(kk+jj) = intVecTemp(jj)-1
        END DO
        ee = ee+1
      END DO

      eptr(nElem_local+1) = nElem_local*npElem


      ! Call the METIS's options and change them if necessary
      !call ParMETIS_SetDefaultOptions(options_metis)
      call MPI_Barrier(MPI_COMM_WORLD, errpar)
      WRITE(*,*) " Before Metis "

      ALLOCATE(elmwgt(nElem_local))
      elmwgt = 1
      !
      wgtflag = 0
      !
      numflag = 0
      !
      ncon    = 1
      !
      !nparts  = 4 !n_mpi_procs
      !
      ALLOCATE(tpwgts(ncon*nparts))
      fact = 1.0/float(nparts)
      tpwgts = fact
      write(*,*) tpwgts
      !
      ALLOCATE(ubvec(ncon))
      ubvec = 1.05


      options_metis(1) = 0

      ! METIS partition routine
      call ParMETIS_V3_PartMeshKway(elemdist, eptr, eind, elmwgt, wgtflag, numflag, &
      & ncon, ncommonnodes, nparts, tpwgts, ubvec, &
      & options_metis, edgecut, elem_proc_id_local, MPI_COMM_WORLD, errpar)

      call MPI_Barrier(MPI_COMM_WORLD, errpar)
      WRITE(*,*) " After Metis "

      ! IF(ind == METIS_OK) THEN
      !   WRITE(*,*) " METIS partition routine successful "
      ! ELSE
      !   STOP " METIS partition routine FAILED "
      ! END IF

      ALLOCATE(recvcounts(n_mpi_procs))
      ALLOCATE(displs(n_mpi_procs))

      call MPI_Allgather(nElem_local, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD, errpar)

      displs(1) = 0;
      do ii=2, n_mpi_procs
        displs(ii) = elemdist(ii)
      end do

      call MPI_Allgatherv(elem_proc_id_local, nElem_local, MPI_INT, elem_proc_id, &
      & recvcounts, displs, MPI_INT, MPI_COMM_WORLD, errpar);

    WRITE(*,*) " After Metis "

    ! IF(ind == METIS_OK) THEN
    !   WRITE(*,*) " METIS partition routine successful "
    ! ELSE
    !   STOP " METIS partition routine FAILED "
    ! END IF
    !DO ii=1,numberOfTetEl
      !WRITE(*,*) ii, elem_proc_id(ii)
    !ENDDO

    end subroutine partitionFineGrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine writevtk

      IMPLICIT NONE


      WRITE(*,*) "Writing VTK file"

      outfileName="partition-mesh.vtk"

      OPEN(11, file=outfileName, STATUS="UNKNOWN", ACTION="WRITE")

      ! Directives
      WRITE(11,'(A)') "# vtk DataFile Version 4.0"
      WRITE(11,'(A)') "Mesh parition"

      ! ASCII or Binary
      WRITE(11,*) "ASCII"

      ! Type of dataset : Structured/Unstructured/Rectilinear/Polydata...
      WRITE(11,'(A)') "DATASET UNSTRUCTURED_GRID"
      
      ! Coordinates of the points (nodes)
      WRITE(11,'(A,I8,A)') "POINTS ", numberOfNodes, " float"
      DO ii=1,numberOfNodes
        WRITE(11,'(F12.6,F12.6,F12.6)') coords(ii,1), coords(ii,2), coords(ii,3)
      END DO

      ! Element<->Nodes connectivity
      ! In VTK terminology, Cell<->Points
      ! <number of nodes> <n1> <n2> <n3> ...
      ! Starting index is 0
      ind = numberOfTetEl*(npElem+1)
      WRITE(11,'(A,I8,I8)') "CELLS ", numberOfTetEl, ind

      ! Tetrahedral element
      IF(npElem == 4) THEN
          n1 = 10
          DO ee=1,numberOfTetEl
            WRITE(11,'(I6,I12,I12,I12,I12)') npElem, elemNodeConn(ee,1)-1, elemNodeConn(ee,2)-1, elemNodeConn(ee,3)-1, elemNodeConn(ee,4)-1
          END DO
      END IF

      ! Cell type, as per VTK
      WRITE(11,'(A,I8)') "CELL_TYPES", numberOfTetEl
      DO ee=1,numberOfTetEl
        WRITE(11,'(I3)') n1
      END DO

      ! Cell data. Processor ID
      WRITE(11,'(A,I8)') "CELL_DATA", numberOfTetEl
      WRITE(11,'(A)') "SCALARS procid int 1"
      WRITE(11,'(A)') "LOOKUP_TABLE default"
      DO ee=1,numberOfTetEl
        WRITE(11,'(I3)') elem_proc_id(ee)
      END DO

      ! close the file
      CLOSE(11)
  
    end subroutine writevtk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cleanup

      IMPLICIT NONE
  
        WRITE(*,*) "Deallocating the memory"
  
        IF( ALLOCATED(coords) )        DEALLOCATE(coords)
        IF( ALLOCATED(elemNodeConn) )  DEALLOCATE(elemNodeConn)
        IF( ALLOCATED(intVecTemp) )    DEALLOCATE(intVecTemp)
        IF( ALLOCATED(elem_proc_id) )  DEALLOCATE(elem_proc_id)
        IF( ALLOCATED(node_proc_id) )  DEALLOCATE(node_proc_id)
        IF( ALLOCATED(eptr) )          DEALLOCATE(eptr)
        IF( ALLOCATED(eind) )          DEALLOCATE(eind)
    
    end subroutine cleanup
!------------------------------------------------------------------
!------------------------------------------------------------------
integer function nameLen(fn)
    IMPLICIT NONE
 
    character*80 :: fn
    integer :: i

    do i = 80,1,-1
      nameLen = i
      if(fn(i:i)/=' ') GOTO 77 ! EXIT
    end do 
    nameLen = 0
77  end function nameLen
!--------------------------------------------------------------------------
 subroutine membreak()
 IMPLICIT NONE
  
 character*20 :: tex
 
 write(*,*) "Type enter to continue"
 read(*,*) tex
 
 end subroutine membreak
!--------------------------------------------------------------------------

    end module partition

