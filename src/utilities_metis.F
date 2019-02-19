!
! Mesh partitioning using METIS
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 13-Nov-2017
! Place : Swansea, UK
!
!
!

      MODULE Utilities_metis

      CONTAINS

!***********************************************************************************************************************
!> \brief Partition the Finite Element mesh using Metis
!> \param MetisAlgoType    ! IN,  Metis algorithm type --> Nodal (==1) or Dual (else)
!> \param nElem,           ! IN,  total number of elements in the mesh
!> \param nNode,           ! IN,  total number of nodes in the mesh
!> \param nodes_per_elem   ! IN,  nodes per element. (Needs to change this to array if the mesh contains mixed elements)
!> \param elemNodeConn,    ! IN,  element<->node connectivity
!> \param nodes_per_side,  ! IN,  nodes per side, for partitioning with Dual (element) graph
!> \param nParts,          ! IN,  number of partitions
!> \param elem_proc_id,    ! OUT, processor ID of the elements
!> \param node_proc_id     ! OUT, processor ID of the nodes
!> \author C. Kadapa
!> \date 13 Nov 2017
!>
!> extra, for future coding purposes
!> \param MixedElemsFlag   ! IN,  Flag to specify if the mesh contains mixed elements
!***********************************************************************************************************************

      SUBROUTINE PartitionMeshMetis(
     1     MetisAlgoType,   ! IN,  Metis algorithm type --> Nodal or Dual
     2     nElem,           ! IN,  total number of elements in the mesh
     3     nNode,           ! IN,  total number of nodes in the mesh
     4     nodes_per_elem   ! IN,  nodes per element
     5     elemNodeConn,    ! IN,  element<->node connectivity
     6     nodes_per_side,  ! IN,  nodes per side, for partitioning with Dual (element) graph
     7     nParts,          ! IN,  number of partitions
     8     elem_proc_id,    ! OUT, processor ID of the elements
     9     node_proc_id     ! OUT, processor ID of the nodes
      )

        IMPLICIT NONE

        ! function arguments
        INTEGER, INTENT(IN) :: MetisAlgoType, nElem, nNode, nodes_per_elem, nodes_per_side, nParts
        INTEGER, DIMENSION(:), INTENT(IN)   :: elemNodeConn
        INTEGER, DIMENSION(:), INTENT(OUT)  :: elem_proc_id, node_proc_id

        INTEGER, DIMENSION(:), ALLOCATABLE :: eptr, eind

        ! local data
        INTEGER :: ind, ee
        INTEGER, pointer     :: vwgt=>null(), vsize=>null()
        DOUBLE PRECISION, pointer     :: tpwgts=>null()
        INTEGER :: options_metis(100)


        WRITE(*,*) " Partitioning the mesh ......."

        ! array of size (nElem_global+1) which stores
        ! number of nodes per each element,
        ! with 0 as the first entry and
        ! nElem_global*nodes_per_elem as the last entry
        ind = nElem+1
        ALLOCATE( eptr(ind) )

        ! array of size 'nElem_global*nodes_per_elem'
        ! which stores eleme<->node connectivity information
        ind = nElem_global*nodes_per_elem
        ALLOCATE( eind(ind) )

        WRITE(*,*) " Before Metis ...... "

          DO ee=1,nElem
            eptr(ee) = (ee-1)*nodes_per_elem

            kk = (ee-1)*nodes_per_elem;

            DO ii=1, nodes_per_elem
              eind(kk+ii) = elemNodeConn(ee,ii)-1
            END DO
          END DO

          eptr(nElem+1) = nElem*nodes_per_elem

          ! METIS options
          call METIS_SetDefaultOptions(options_metis)

        WRITE(*,*) " Before Metis ...... "

          ! METIS partition routine
        IF(MetisAlgoType == 1) THEN
          call METIS_PartMeshNodal(
     1      nElem_global, nNode_global, eptr, eind, vwgt, vsize, 
     2      nParts, tpwgts, options_metis,
     3      objval, elem_proc_id, node_proc_id);
        ELSE
          call METIS_PartMeshDual(
     1      nElem_global, nNode_global, eptr, eind, vwgt, vsize, 
     2      nodes_per_side, nParts, tpwgts, options_metis,
     3      objval, elem_proc_id, node_proc_id)
       END IF

      END SUBROUTINE Utilities_metis

      END MODULE Utilities_metis
