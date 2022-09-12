module ADTProjections
    !
    !      Module, which defines the derived data types and the arrays to
    !      store multiple ADT's. An array is chosen to store multiple
    !      ATD's rather than a linked list, because this is more
    !      convenient during the search. When the ADT's are built there
    !      is some additional work due to reallocation. However this is
    !      negligible due to the usage of pointers.
    !
    use precision, only: realType
    implicit none
    save
    !
    !      Define the functions needed for the sorting of the derived
    !      data types to be private, i.e. they can only be accessed
    !      within this module.
    !
    public
    private :: adtBBoxTargetTypeLessEqual
    private :: adtBBoxTargetTypeLess
    private :: adtTypeAssign
    !
    !      Definition of the derived data type store a leaf of an ADT.
    !      The ADT itself is an array of these leaves.
    !
    type adtLeafType

        ! children(2): Children of the parent. If negative it means that
        !              it is a terminal leaf and the absolute values
        !              indicate the bounding box id's. Note that it is
        !              allowed that 1 child is negative and the other
        !              positive.
        ! xMin(6):     The minimum coordinates of the leaf.
        ! xMax(6):     The maximum coordinates of the leaf.

        integer, dimension(2) :: children
        real(kind=realType), dimension(6) :: xMin, xMax

    end type adtLeafType
    !
    !      The definition of adtBBoxTargetType, which stores the data of
    !      a possible bounding box which minimizes the distances to the
    !      given coordinate.
    !
    type adtBBoxTargetType

        ! ID:       The id of the bounding box in the list.
        ! posDist2: the possible minimum distance squared to the active
        !           coordinate.

        integer :: ID
        real(kind=realType) :: posDist2

    end type adtBBoxTargetType

    ! Interfaces for the extension of the operators <= and <.
    ! These are needed for the sorting of BBoxTargetType. Note
    ! that the = operator does not need to be defined, because
    ! BBoxTargetType only contains primitive types.

    interface operator(<=)
        module procedure adtBBoxTargetTypeLessEqual
    end interface operator(<=)

    interface operator(<)
        module procedure adtBBoxTargetTypeLess
    end interface operator(<)

    !
    !      Definition of the derived data type to store an ADT.
    !
    type adtType

        ! adtType:  Type of ADT. Possibilities are adtSurfaceADT and
        !           adtVolumeADT.
        ! adtID:    The given ID of the ADT.
        ! isActive: Whether or not the ADT is active. If not, this
        !           entry could be used during a reallocation.

        integer :: adtType
        character(len=64) :: adtID
        logical :: isActive

        ! nNodes:  Number of local nodes in the given grid.
        ! nTria:   Number of local triangles in the given grid.
        ! nQuads:  Idem for the quadrilaterals.
        ! nTetra:  Idem for the tetrahedra.
        ! nPyra:   Idem for the pyramids.
        ! nPrisms: Idem for the prisms.
        ! nHexa:   Idem for the hexahedra.

        integer :: nNodes, nTria, nQuads
        integer :: nTetra, nPyra, nPrisms, nHexa

        ! coor(3,nNodes): Nodal coordinates of the local grid.
        !                 To save memory this pointer is not
        !                 allocated, but set to the data given.
        real(kind=realType), dimension(:, :), pointer :: coor

        ! triaConn(3,nTria):     Local connectivity of the triangles.
        !                        To save memory this pointer is not
        !                        allocated, but set to the data given.
        ! quadsConn(4,nQuads):   Idem for the quadrilaterals.
        ! tetraConn(4,nTetra):   Idem for the tetrahedra.
        ! pyraConn(5,nPyra):     Idem for the pyramids.
        ! prismsConn(6,nPrisms): Idem for the prisms.
        ! hexaConn(8,nHexa):     Idem for the hexahedra.

        integer, dimension(:, :), pointer :: triaConn
        integer, dimension(:, :), pointer :: quadsConn
        integer, dimension(:, :), pointer :: tetraConn
        integer, dimension(:, :), pointer :: pyraConn
        integer, dimension(:, :), pointer :: prismsConn
        integer, dimension(:, :), pointer :: hexaConn

        ! nRootLeaves:        Number of non-empty root leaves.
        !                     This number is of course less than or
        !                     equal to nProcs.
        ! myEntryInRootProcs: If the local tree is not empty, this
        !                     contains the entry in rootLeavesProcs.
        ! rootLeavesProcs(:): The corresponding processor ID's.
        ! rootBBoxes(3,2,:):  The 3D bounding boxes of the non-empty
        !                     root leaves.

        integer :: nRootLeaves, myEntryInRootProcs
        integer, dimension(:), pointer :: rootLeavesProcs
        real(kind=realType), dimension(:, :, :), pointer :: rootBBoxes

        ! nBBoxes:              Number of bounding boxes stored in
        !                       the ADT.
        ! elementType(nBBoxes): The corresponding element type of the
        !                       bounding box.
        ! elementID(nBBoxes):   The corresponding entry in the element
        !                       connectivity of the bounding box.
        ! xBBox(6,nBBoxes):     The coordinates of the bounding boxes
        !                       of the elements stored in the ADT.

        integer :: nBBoxes

        integer, dimension(:), pointer :: elementType
        integer, dimension(:), pointer :: elementID
        real(kind=realType), dimension(:, :), pointer :: xBBox

        ! nLeaves:         Number of present in the ADT. Due to the
        !                  variable splitting the tree is optimally
        !                  balanced and therefore nLeaves = nBBoxes -1.
        ! ADTree(nLeaves): The alternating digital tree.

        integer :: nLeaves
        type(adtLeafType), dimension(:), pointer :: ADTree

    end type adtType

    !
    !  Variables stored in this module.
    !
    ! nStack:   Number of elements allocated in the stack array;
    !           needed for a more efficient implementation of the
    !           local qsort routines.
    ! stack(:): The corresponding array to store the stack.
    !           This is a pointer, such that the reallocation
    !           is easier.

    integer :: nStack
    integer, dimension(:), pointer :: stack

    integer, parameter :: adtSurfaceADT = 1
    integer, parameter :: adtVolumeADT = 2

    integer, parameter :: adtTriangle = 1
    integer, parameter :: adtQuadrilateral = 2

    ! Interface for the extension of the operator =.

    interface assignment(=)
        module procedure adtTypeAssign
    end interface assignment(=)

contains

    !===============================================================

    logical function adtBBoxTargetTypeLessEqual(g1, g2)
        !
        !        This function returns .true. if g1 <= g2. The comparison is
        !        firstly based on the possible minimum distance such that the
        !        most likely candidates are treated first. In case of ties
        !        the boundary box ID is considered.
        !        Function intent(in) arguments.
        !        ------------------------------
        !        g1, g2: The two instances of the derived datatype that most
        !                be compared.
        !
        implicit none
        !
        !       Function arguments.
        !
        type(adtBBoxTargetType), intent(in) :: g1, g2

        ! Compare the possible minimum distances.

        if (g1%posDist2 < g2%posDist2) then
            adtBBoxTargetTypeLessEqual = .true.
            return
        else if (g1%posDist2 > g2%posDist2) then
            adtBBoxTargetTypeLessEqual = .false.
            return
        end if

        ! Compare the bounding box ID's.

        if (g1%ID < g2%ID) then
            adtBBoxTargetTypeLessEqual = .true.
            return
        else if (g1%ID > g2%ID) then
            adtBBoxTargetTypeLessEqual = .false.
            return
        end if

        ! g1 and g2 are identical. Return .true.

        adtBBoxTargetTypeLessEqual = .true.

    end function adtBBoxTargetTypeLessEqual

    !===============================================================

    logical function adtBBoxTargetTypeLess(g1, g2)
        !
        !        This function returns .true. if g1 < g2. The comparison is
        !        firstly based on the possible minimum distance such that the
        !        most likely candidates are treated first. In case of ties
        !        the boundary box ID is considered.
        !        Function intent(in) arguments.
        !        ------------------------------
        !        g1, g2: The two instances of the derived datatype that most
        !                be compared.
        !
        implicit none
        !
        !       Function arguments.
        !
        type(adtBBoxTargetType), intent(in) :: g1, g2

        ! Compare the possible minimum distances.

        if (g1%posDist2 < g2%posDist2) then
            adtBBoxTargetTypeLess = .true.
            return
        else if (g1%posDist2 > g2%posDist2) then
            adtBBoxTargetTypeLess = .false.
            return
        end if

        ! Compare the bounding box ID's.

        if (g1%ID < g2%ID) then
            adtBBoxTargetTypeLess = .true.
            return
        else if (g1%ID > g2%ID) then
            adtBBoxTargetTypeLess = .false.
            return
        end if

        ! g1 and g2 are identical. Return .false.

        adtBBoxTargetTypeLess = .false.

    end function adtBBoxTargetTypeLess

    !===============================================================

    subroutine adtTypeAssign(g1, g2)
        !
        !        This subroutine defines the generic assignment operator for
        !        the derived datatype adtType. The contents of g1 is copied
        !        into g2, where pointers just point to the other pointers,
        !        i.e. no additional allocation takes place.
        !        Subroutine intent(in) arguments.
        !        --------------------------------
        !        g2: Entity whose data should be copied.
        !        Subroutine intent(out) arguments.
        !        ---------------------------------
        !        g1: Entity whose data should be assigned.
        !
        implicit none
        !
        !       Subroutine arguments.
        !
        type(adtType), intent(in) :: g2
        type(adtType), intent(out) :: g1

        g1%adtType = g2%adtType
        g1%adtID = g2%adtID
        g1%isActive = g2%isActive

        g1%nNodes = g2%nNodes
        g1%nTria = g2%nTria
        g1%nQuads = g2%nQuads
        g1%nTetra = g2%nTetra
        g1%nPyra = g2%nPyra
        g1%nPrisms = g2%nPrisms
        g1%nHexa = g2%nHexa

        g1%coor => g2%coor
        g1%triaConn => g2%triaConn
        g1%quadsConn => g2%quadsConn
        g1%tetraConn => g2%tetraConn
        g1%pyraConn => g2%pyraConn
        g1%prismsConn => g2%prismsConn
        g1%hexaConn => g2%hexaConn

        g1%nRootLeaves = g2%nRootLeaves
        g1%myEntryInRootProcs = g2%myEntryInRootProcs
        g1%rootLeavesProcs => g2%rootLeavesProcs
        g1%rootBBoxes => g2%rootBBoxes

        g1%nBBoxes = g2%nBBoxes
        g1%elementType => g2%elementType
        g1%elementID => g2%elementID
        g1%xBBox => g2%xBBox

        g1%nLeaves = g2%nLeaves
        g1%ADTree => g2%ADTree

    end subroutine adtTypeAssign

    subroutine adtTerminate(ADT, routineName, errorMessage)
        !
        !        This routine writes the given error message to standard
        !        output and terminates the executation of the program.
        !        Subroutine intent(in) arguments.
        !        --------------------------------
        !        routineName: Name of the routine where the error occured.
        !        ADT:         Currently active ADT.
        !        Subroutine intent(inout) arguments.
        !        -----------------------------------
        !        errorMessage: On input it contains the error message to be
        !                      written. It is modified in this routine, such
        !                      that it fits on one line. On output its
        !                      contents is undefined, which does not matter
        !                      a whole lot.
        !
        implicit none
        !
        !       Subroutine arguments
        !
        type(adtType), intent(in) :: ADT

        character(len=*), intent(in) :: routineName
        character(len=*), intent(in) :: errorMessage
        !
        !       Local parameter
        !
        integer, parameter :: maxCharLine = 55
        !
        !       Local variables
        !
        integer :: ierr, len, i2
        logical :: firstTime

        character(len=len_trim(errorMessage)) :: message
        character(len=8) :: integerString

        ! Copy the errorMessage into message. It is not possible to work
        ! with errorMessage directly, because it is modified in this
        ! routine. Sometimes a constant string is passed to this routine
        ! and some compilers simply fail then.

        message = errorMessage

        ! Print a nice error message. In case of a parallel executable
        ! also the processor ID is printed.

        print "(a)", "#"
        print "(a)", "#=========================== !!! Error !!! &
             &============================"

        print "(2a)", "#* Run-time error in procedure ", &
            trim(routineName)

        ! Loop to write the error message. If the message is too long it
        ! is split over several lines.

        firstTime = .true.
        do
            ! Determine the remaining error message to be written.
            ! If longer than the maximum number of characters allowed
            ! on a line, it is attempted to split the message.

            message = adjustl(message)
            len = len_trim(message)
            i2 = min(maxCharLine, len)

            if (i2 < len) i2 = index(message(:i2), " ", .true.) - 1
            if (i2 < 0) i2 = index(message, " ") - 1
            if (i2 < 0) i2 = len

            ! Write this part of the error message. If it is the first
            ! line of the message some additional stuff is printed.

            if (firstTime) then
                print "(2a)", "#* Error message: ", trim(message(:i2))
                firstTime = .false.
            else
                print "(2a)", "#*                ", trim(message(:i2))
            end if

            ! Exit the loop if the entire message has been written.

            if (i2 == len) exit

            ! Adapt the string for the next part to be written.

            message = message(i2 + 1:)

        end do

        ! Write the trailing message.

        print "(a)", "#*"
        print "(a)", "#* Now exiting"
        print "(a)", "#==========================================&
             &============================"
        print "(a)", "#"

        ! Call abort and stop the program. This stop should be done in
        ! abort, but just to be sure.

        stop

    end subroutine adtTerminate

    subroutine qsortBBoxes(arr, nn, ADT, dir)
        !
        !        This routine sorts the integer array arr, such that the
        !        coordinate of the corresponding bounding box in the
        !        direction dir is in increasing order. Note that the array to
        !        store the stack is stored in this module. The reason is that
        !        this routine is called quite often and in this way constant
        !        allocation, deallocation and reallocation of stack is
        !        avoided.
        !        Subroutine intent(in) arguments.
        !        --------------------------------
        !        nn:  Size of the array to be sorted.
        !        ADT: The ADTfrom which the coordinate of
        !             the bounding box must be taken.
        !        dir: Index of the coordinate, which must be sorted.
        !        Subroutine intent(inout) arguments.
        !        -----------------------------------
        !        arr(:): On input it contains the bounding box ID's which
        !                must be sorted. On output these ID's are sorted,
        !                such that the given coordinate is in increasing
        !                order.
        !
        implicit none
        !
        !       Subroutine arguments.
        !
        type(adtType), intent(in) :: ADT
        integer, intent(in) :: nn, dir

        integer, dimension(:), intent(inout) :: arr
        !
        !       Local parameters.
        !
        integer, parameter :: m = 7
        !
        !       Local variables.
        !
        integer :: i, j, k, r, l, jStack
        integer :: a, tmp

        real(kind=realType) :: ra
        real(kind=realType), dimension(:, :), pointer :: xBBox

        ! Set the pointer for the coordinates of the bounding boxes.

        xBBox => ADT%xBBox

        ! Initialize the variables that control the sorting.

        jStack = 0
        l = 1
        r = nn

        ! Start of the algorithm.

        sortLoop: do

            ! Check for the size of the subarray.

            testInsertion: if ((r - l) < m) then

                ! Perform the insertion sort.

                do j = (l + 1), r
                    a = arr(j)
                    ra = xBBox(dir, a)
                    do i = (j - 1), l, -1
                        if (xBBox(dir, arr(i)) <= ra) exit
                        arr(i + 1) = arr(i)
                    end do
                    arr(i + 1) = a
                end do

                ! In case there are no more elements on the stack, exit from
                ! the outermost do-loop. Algorithm has finished.

                if (jStack == 0) exit sortLoop

                ! Pop stack and begin a new round of partitioning.

                r = stack(jStack)
                l = stack(jStack - 1)
                jStack = jStack - 2

                else testInsertion

                ! Subarray is larger than the threshold for a linear sort.
                ! Choose median of left, center and right elements as
                ! partitioning element a. Also rearrange so that
                ! (l) <= (l+1) <= (r).

                k = (l + r) / 2
                tmp = arr(k)      ! Swap the elements
                arr(k) = arr(l + 1)    ! k and l+1.
                arr(l + 1) = tmp

                if (xBBox(dir, arr(r)) < xBBox(dir, arr(l))) then
                    tmp = arr(l)             ! Swap the elements
                    arr(l) = arr(r)             ! r and l.
                    arr(r) = tmp
                end if

                if (xBBox(dir, arr(r)) < xBBox(dir, arr(l + 1))) then
                    tmp = arr(l + 1)         ! Swap the elements
                    arr(l + 1) = arr(r)           ! r and l+1.
                    arr(r) = tmp
                end if

                if (xBBox(dir, arr(l + 1)) < xBBox(dir, arr(l))) then
                    tmp = arr(l + 1)         ! Swap the elements
                    arr(l + 1) = arr(l)           ! l and l+1.
                    arr(l) = tmp
                end if

                ! Initialize the pointers for partitioning.

                i = l + 1
                j = r
                a = arr(l + 1)
                ra = xBBox(dir, a)

                ! The innermost loop.

                innerLoop: do

                    ! Scan up to find element >= a.
                    do
                        i = i + 1
                        if (ra <= xBBox(dir, arr(i))) exit
                    end do

                    ! Scan down to find element <= a.
                    do
                        j = j - 1
                        if (xBBox(dir, arr(j)) <= ra) exit
                    end do

                    ! Exit the loop in case the pointers i and j crossed.

                    if (j < i) exit innerLoop

                    ! Swap the element i and j.

                    tmp = arr(i)
                    arr(i) = arr(j)
                    arr(j) = tmp

                end do innerLoop

                ! Swap the entries j and l+1. Remember that a equals
                ! arr(l+1).

                arr(l + 1) = arr(j)
                arr(j) = a

                ! Push pointers to larger subarray on stack; process smaller
                ! subarray immediately. Check if enough memory is available.
                ! If not reallocate it.

                jStack = jStack + 2

                if (jStack > nStack) call reallocPlus(stack, nStack, 100, ADT)

                if ((r - i + 1) >= (j - l)) then
                    stack(jStack) = r
                    r = j - 1
                    stack(jStack - 1) = j
                else
                    stack(jStack) = j - 1
                    stack(jStack - 1) = l
                    l = j
                end if

            end if testInsertion
        end do sortLoop

    end subroutine qsortBBoxes

    subroutine qsortBBoxTargets(arr, nn, ADT)
        !
        !        This routine sorts the given number of bounding box targets
        !        in increasing order, based on the generalized < operator.
        !
        implicit none
        !
        !       Subroutine arguments
        !
        type(adtType), intent(in) :: ADT
        integer, intent(in) :: nn

        type(adtBBoxTargetType), dimension(:), pointer :: arr
        !
        !       Local variables
        !
        integer, parameter :: m = 7

        integer :: i, j, k, r, l, jStack

        type(adtBBoxTargetType) :: a, tmp

        ! Initialize the variables that control the sorting.

        jStack = 0
        l = 1
        r = nn

        ! Start of the algorithm

        sortLoop: do

            ! Check for the size of the subarray.

            testInsertion: if ((r - l) < m) then

                ! Perform insertion sort

                do j = l + 1, r
                    a = arr(j)
                    do i = (j - 1), l, -1
                        if (arr(i) <= a) exit
                        arr(i + 1) = arr(i)
                    end do
                    arr(i + 1) = a
                end do

                ! In case there are no more elements on the stack, exit from
                ! the outermost do-loop. Algorithm has finished.

                if (jStack == 0) exit sortLoop

                ! Pop stack and begin a new round of partitioning.

                r = stack(jStack)
                l = stack(jStack - 1)
                jStack = jStack - 2

                else testInsertion

                ! Subarray is larger than the threshold for a linear sort.
                ! Choose median of left, center and right elements as
                ! partitioning element a. Also rearrange so that
                ! (l) <= (l+1) <= (r).

                k = (l + r) / 2
                tmp = arr(k)      ! Swap the elements
                arr(k) = arr(l + 1)    ! k and l+1.
                arr(l + 1) = tmp

                if (arr(r) < arr(l)) then
                    tmp = arr(l)             ! Swap the elements
                    arr(l) = arr(r)             ! r and l.
                    arr(r) = tmp
                end if

                if (arr(r) < arr(l + 1)) then
                    tmp = arr(l + 1)         ! Swap the elements
                    arr(l + 1) = arr(r)           ! r and l+1.
                    arr(r) = tmp
                end if

                if (arr(l + 1) < arr(l)) then
                    tmp = arr(l + 1)         ! Swap the elements
                    arr(l + 1) = arr(l)           ! l and l+1.
                    arr(l) = tmp
                end if

                ! Initialize the pointers for partitioning.

                i = l + 1
                j = r
                a = arr(l + 1)

                ! The innermost loop

                innerLoop: do

                    ! Scan up to find element >= a.
                    do
                        i = i + 1
                        if (a <= arr(i)) exit
                    end do

                    ! Scan down to find element <= a.
                    do
                        j = j - 1
                        if (arr(j) <= a) exit
                    end do

                    ! Exit the loop in case the pointers i and j crossed.

                    if (j < i) exit innerLoop

                    ! Swap the element i and j.

                    tmp = arr(i)
                    arr(i) = arr(j)
                    arr(j) = tmp

                end do innerLoop

                ! Swap the entries j and l+1. Remember that a equals
                ! arr(l+1).

                arr(l + 1) = arr(j)
                arr(j) = a

                ! Push pointers to larger subarray on stack; process smaller
                ! subarray immediately. Check if enough memory is available.
                ! If not reallocate it.

                jStack = jStack + 2

                if (jStack > nStack) call reallocPlus(stack, nStack, 100, ADT)

                if ((r - i + 1) >= (j - l)) then
                    stack(jStack) = r
                    r = j - 1
                    stack(jStack - 1) = j
                else
                    stack(jStack) = j - 1
                    stack(jStack - 1) = l
                    l = j
                end if

            end if testInsertion
        end do sortLoop

    end subroutine qsortBBoxTargets

    subroutine reallocBBoxTargetTypePlus(arr, nSize, nInc, ADT)
        !
        !        This routine reallocates the memory of the given
        !        adtBBoxTargetType pointer array.
        !        Subroutine intent(in) arguments.
        !        --------------------------------
        !        ADT:         Currently active ADT.
        !        nInc: Increment of the size of the array.
        !        Subroutine intent(inout) arguments.
        !        -----------------------------------
        !        nSize: On input it contains the size of the given array.
        !               On output this value is incremented by nInc.
        !        Subroutine pointer arguments.
        !        -----------------------------
        !        arr: Array to be reallocated.
        !
        implicit none
        !
        !       Subroutine arguments.
        !
        type(adtType), intent(in) :: ADT
        integer, intent(in) :: nInc
        integer, intent(inout) :: nSize

        type(adtBBoxTargetType), dimension(:), pointer :: arr
        !
        !       Local variables.
        !
        integer :: ierr
        integer :: i, nOld

        type(adtBBoxTargetType), dimension(:), pointer :: tmp

        ! Store the input value of nSize and set the pointer tmp to the
        ! original array.

        nOld = nSize
        tmp => arr

        ! Allocate the new memory for the array.

        nSize = nSize + nInc
        allocate (arr(nSize), stat=ierr)
        if (ierr /= 0) &
            call adtTerminate(ADT, "reallocBBoxTargetTypePlus", &
                              "Memory allocation failure for arr.")

        ! Copy the data from the original array into arr.

        nOld = min(nOld, nSize)
        do i = 1, nOld
            arr(i) = tmp(i)
        end do

        ! Release the memory of tmp, which points to the original
        ! memory of the given array.

        deallocate (tmp, stat=ierr)
        if (ierr /= 0) &
            call adtTerminate(ADT, "reallocBBoxTargetTypePlus", &
                              "Deallocation failure for tmp.")

    end subroutine reallocBBoxTargetTypePlus

    subroutine reallocPlus(arr, nSize, nInc, ADT)
        !
        !        This internal routine reallocates the memory of the given
        !        pointer array.
        !        Subroutine intent(in) arguments.
        !        --------------------------------
        !        ADT:         Currently active ADT.
        !        nInc: Increment of the size of the array.
        !        Subroutine intent(inout) arguments.
        !        -----------------------------------
        !        nSize: On input it contains the size of the given array.
        !               On output this value is incremented by nInc.
        !        Subroutine pointer arguments.
        !        -----------------------------
        !        arr: Array to be reallocated.
        !
        implicit none
        !
        !       Subroutine arguments.
        !
        type(adtType), intent(in) :: ADT
        integer, intent(in) :: nInc
        integer, intent(inout) :: nSize

        integer, dimension(:), pointer :: arr
        !
        !       Local variables.
        !
        integer :: ierr
        integer :: i, nOld

        integer, dimension(:), pointer :: tmp

        ! Store the input value of nSize and set the pointer tmp to the
        ! original array.

        nOld = nSize
        tmp => arr

        ! Allocate the new memory for the array.

        nSize = nSize + nInc
        allocate (arr(nSize), stat=ierr)
        if (ierr /= 0) &
            call adtTerminate(ADT, "reallocPlus", &
                              "Memory allocation failure for arr.")

        ! Copy the data from the original array into arr.

        nOld = min(nOld, nSize)
        do i = 1, nOld
            arr(i) = tmp(i)
        end do

        ! Release the memory of tmp, which points to the original
        ! memory of the given array.

        deallocate (tmp, stat=ierr)
        if (ierr /= 0) &
            call adtTerminate(ADT, "reallocPlus", &
                              "Deallocation failure for tmp.")

    end subroutine reallocPlus

    subroutine buildADT(ADT)
        !
        !        This routine builds the 6 dimensional ADT for the given
        !        ADT. When this routine is called it is assumed that the
        !        bounding boxes of the grid have already been computed; the
        !        ADT for these bounding boxes is built here.
        !        Subroutine intent(inout) arguments.
        !        --------------------------------
        !        ADT: adt derived type to build
        !
        implicit none
        !
        !       Subroutine arguments.
        !
        type(adtType), intent(inout) :: ADT
        !
        !       Local variables.
        !
        integer :: ierr

        integer :: i, j, k, ii, kk, ll, mm, nn, nfl, nfr
        integer :: nLeaves, nBBoxes, splitDir
        integer :: nLeavesToDivide, nLeavesToDivideNew
        integer :: nLeavesTot

        integer, dimension(:), pointer :: BB_IDs
        integer, dimension(:), pointer :: BB_IDsNew
        integer, dimension(:), pointer :: nBB_IDs
        integer, dimension(:), pointer :: nBB_IDsNew
        integer, dimension(:), pointer :: curLeaf
        integer, dimension(:), pointer :: curLeafNew
        integer, dimension(:), pointer :: tmpIntPointer
        real(kind=realType), dimension(:, :), pointer :: xBBox
        type(adtLeafType), dimension(:), pointer :: ADTree

        ! Initialize nStack and allocate the corresponding array stack.
        ! These are used in the qsort routine for the bounding boxes.
        ! As this routine is called quite often it is more efficient to
        ! have the stack array available rather than allocate it over
        ! and over again.

        nStack = 100
        allocate (stack(nStack), stat=ierr)
        if (ierr /= 0) &
            call adtTerminate(ADT, "buildADT", &
                              "Memory allocation failure for stack.")

        ! Determine the number of leaves of the adt. It can be proved
        ! that nLeaves equals nBBoxes - 1 for an optimally balanced
        ! tree. Take the exceptional case of nBBoxes == 0 and
        ! nBBoxes == 1 into account.

        nBBoxes = ADT%nBBoxes
        nLeaves = nBBoxes - 1
        if (nBBoxes <= 1) nLeaves = nLeaves + 1

        ADT%nLeaves = nLeaves

        ! Allocate the memory for the adt.

        allocate (ADT%ADTree(nLeaves), stat=ierr)
        if (ierr /= 0) &
            call adtTerminate(ADT, "buildADT", &
                              "Memory allocation failure for ADTree.")

        ! Set some pointers to make the code more readable.

        xBBox => ADT%xBBox
        ADTree => ADT%ADTree

        ! Allocate the memory for the arrays which control the
        ! subdivision of the leaves.

        nn = (nBBoxes + 1) / 2
        nn = max(nn, 1)

        allocate (BB_IDs(nBBoxes), BB_IDsNew(nBBoxes), &
                  nBB_IDs(0:nn), nBB_IDsNew(0:nn), &
                  curLeaf(nn), curLeafNew(nn), stat=ierr)
        if (ierr /= 0) &
             call adtTerminate(ADT, "buildADT", &
             "Memory allocation failure for the arrays &
             &used in the subdivision.")

        ! Initialize the arrays BB_IDs, nBB_IDs and curLeaf, such that
        ! all bounding boxes belong to the root leaf. Also set the
        ! counters nLeavesToDivide and nLeavesTot, depending on the
        ! situation

        nBB_IDs(0) = 0; nBB_IDs(1) = nBBoxes
        curLeaf(1) = 1

        do i = 1, nBBoxes
            BB_IDs(i) = i
        end do

        nLeavesToDivide = min(nLeaves, 1)
        nLeavesTot = nLeavesToDivide

        ! Initialize splitDir to 0, such that the first time it will
        ! split in direction 1.

        splitDir = 0

        ! Loop to subdivide the leaves. The division is such that the
        ! adt is optimally balanced.

        leafDivision: do

            ! Criterion to exit the loop.

            if (nLeavesToDivide == 0) exit

            ! Initializations for the next round of subdivisions and
            ! increment splitDir.

            nLeavesToDivideNew = 0
            nBB_IDsNew(0) = 0

            splitdir = splitDir + 1
            if (splitDir > 6) splitDir = 1

            ! Loop over the current number of leaves to be divided.

            currentLeavesLoop: do i = 1, nLeavesToDivide

                ! Store the number of bounding boxes present in the leaf
                ! in nn, the current leaf number in mm and i-1 in ii.

                ii = i - 1
                nn = nBB_IDs(i) - nBB_IDs(ii)
                mm = curLeaf(i)

                ! Determine the bounding box coordinates of this leaf.

                ll = BB_IDs(nBB_IDs(ii) + 1)
                ADTree(mm)%xMin(1) = xBBox(1, ll)
                ADTree(mm)%xMin(2) = xBBox(2, ll)
                ADTree(mm)%xMin(3) = xBBox(3, ll)
                ADTree(mm)%xMin(4) = xBBox(4, ll)
                ADTree(mm)%xMin(5) = xBBox(5, ll)
                ADTree(mm)%xMin(6) = xBBox(6, ll)

                ADTree(mm)%xMax(1) = xBBox(1, ll)
                ADTree(mm)%xMax(2) = xBBox(2, ll)
                ADTree(mm)%xMax(3) = xBBox(3, ll)
                ADTree(mm)%xMax(4) = xBBox(4, ll)
                ADTree(mm)%xMax(5) = xBBox(5, ll)
                ADTree(mm)%xMax(6) = xBBox(6, ll)

                do j = (nBB_IDs(ii) + 2), nBB_IDs(i)
                    ll = BB_IDs(j)

                    ADTree(mm)%xMin(1) = min(ADTree(mm)%xMin(1), xBBox(1, ll))
                    ADTree(mm)%xMin(2) = min(ADTree(mm)%xMin(2), xBBox(2, ll))
                    ADTree(mm)%xMin(3) = min(ADTree(mm)%xMin(3), xBBox(3, ll))
                    ADTree(mm)%xMin(4) = min(ADTree(mm)%xMin(4), xBBox(4, ll))
                    ADTree(mm)%xMin(5) = min(ADTree(mm)%xMin(5), xBBox(5, ll))
                    ADTree(mm)%xMin(6) = min(ADTree(mm)%xMin(6), xBBox(6, ll))

                    ADTree(mm)%xMax(1) = max(ADTree(mm)%xMax(1), xBBox(1, ll))
                    ADTree(mm)%xMax(2) = max(ADTree(mm)%xMax(2), xBBox(2, ll))
                    ADTree(mm)%xMax(3) = max(ADTree(mm)%xMax(3), xBBox(3, ll))
                    ADTree(mm)%xMax(4) = max(ADTree(mm)%xMax(4), xBBox(4, ll))
                    ADTree(mm)%xMax(5) = max(ADTree(mm)%xMax(5), xBBox(5, ll))
                    ADTree(mm)%xMax(6) = max(ADTree(mm)%xMax(6), xBBox(6, ll))
                end do

                ! Determine the situation of the leaf. It is either a
                ! terminal leaf or a leaf that must be subdivided.

                terminalTest: if (nn <= 2) then

                    ! Terminal leaf. Store the ID's of the bounding boxes with
                    ! negative numbers in children.

                    ADTree(mm)%children(1) = -BB_IDs(nBB_IDs(ii) + 1)
                    ADTree(mm)%children(2) = -BB_IDs(nBB_IDs(i))

                    else terminalTest

                    ! Leaf must be divided. Sort the bounding boxes of the
                    ! current leaf in increasing order; the sorting is based
                    ! on the coordinate in the split direction.

                    call qsortBBoxes(BB_IDs(nBB_IDs(ii) + 1:), nn, ADT, splitDir)

                    ! Determine the number of bounding boxes in the left leaf.
                    ! This number is at least 2. The actual number stored in
                    ! kk is this number plus an offset. Also initialize the
                    ! counter nfl, which is used to store the bounding boxes
                    ! in the arrays for the new round.

                    kk = (nn + 1) / 2 + nBB_IDs(ii)
                    nfl = nBB_IDsNew(nLeavesToDivideNew)

                    ! Copy the ID's of the left bounding boxes into BB_IDsNew.
                    ! Also update nLeavesToDivideNew and the corresponding
                    ! entry in nBB_IDsNew.

                    do k = (nBB_IDs(ii) + 1), kk
                        nfl = nfl + 1
                        BB_IDsNew(nfl) = BB_IDs(k)
                    end do

                    nLeavesToDivideNew = nLeavesToDivideNew + 1
                    nBB_IDsNew(nLeavesToDivideNew) = nfl

                    ! Update the total number of leaves and store this number
                    ! in child 1 of the current leaf and in the current leaves
                    ! for the next round.

                    nLeavesTot = nLeavesTot + 1
                    ADTree(mm)%children(1) = nLeavesTot
                    curLeafNew(nLeavesToDivideNew) = nLeavesTot

                    ! The right leaf will only be created if it has more than
                    ! one bounding box in it, i.e. if the original leaf has
                    ! more than three bounding boxes. If the new leaf only has
                    ! one bounding box in it, it is not created; instead the
                    ! bounding box is stored in the current leaf.

                    if (nn == 3) then

                        ! Only three bounding boxes present in the current leaf.
                        ! The right leaf is not created and the last bounding
                        ! box is stored as the second child of the current leaf.

                        ADTree(mm)%children(2) = -BB_IDs(nBB_IDs(i))

                    else

                        ! More than 3 bounding boxes are present and thus the
                        ! right leaf is created. Copy the ID's from BB_IDs into
                        ! BB_IDsNew and update the counters for the new round.

                        nfr = nBB_IDsNew(nLeavesToDivideNew)
                        do k = (kk + 1), nBB_IDs(i)
                            nfr = nfr + 1
                            BB_IDsNew(nfr) = BB_IDs(k)
                        end do

                        nLeavesToDivideNew = nLeavesToDivideNew + 1
                        nBB_IDsNew(nLeavesToDivideNew) = nfr

                        ! Update the total number of leaves and store this number
                        ! in child 2 of the current leaf and in the current
                        ! leaves for the next round.

                        nLeavesTot = nLeavesTot + 1
                        ADTree(mm)%children(2) = nLeavesTot
                        curLeafNew(nLeavesToDivideNew) = nLeavesTot

                    end if

                end if terminalTest

            end do currentLeavesLoop

            ! Swap the pointers for the next round.

            nLeavesToDivide = nLeavesToDivideNew

            tmpIntPointer => BB_IDs
            BB_IDs => BB_IDsNew
            BB_IDsNew => tmpIntPointer

            tmpIntPointer => nBB_IDs
            nBB_IDs => nBB_IDsNew
            nBB_IDsNew => tmpIntPointer

            tmpIntPointer => curLeaf
            curLeaf => curLeafNew
            curLeafNew => tmpIntPointer

        end do leafDivision

        ! Deallocate the arrays used to build the local tree.

        deallocate (stack, BB_IDs, BB_IDsNew, nBB_IDs, nBB_IDsNew, &
                    curLeaf, curLeafNew, stat=ierr)
        if (ierr /= 0) &
            call adtTerminate(ADT, "buildADT", &
                              "Deallocation failure for the local arrays.")

    end subroutine buildADT

    subroutine searchQuads(pts, conn, searchPts, nPts, nConn, nSearchPts, faceID, uv)
        !
        !        This routine searches for the closest point on a set of
        !        quads for each searchPt. An ADT tree is built and used
        !        for the search and subsequently destroyed.
        !
        implicit none

        ! Input
        integer, intent(in) :: nPts, nConn, nSearchPts
        real(kind=realType), intent(in), target :: pts(3, nPts), searchPts(3, nSearchPts)
        integer, intent(in), target :: conn(4, nConn)

        ! Ouput
        real(kind=realType), intent(out) :: uv(2, nSearchPts)
        integer, intent(out) :: faceID(nSearchPts)

        ! Working Variables
        integer :: ierr, ll, nNPE, intInfo(3)
        integer :: i, j, mm
        real(kind=realType), dimension(3) :: xMin, xMax
        real(kind=realType) :: uvw(5), coor(4)
        type(adtType) :: ADT
        type(adtBBoxTargetType), dimension(:), pointer :: BB
        integer, dimension(:), pointer :: frontLeaves
        integer, dimension(:), pointer :: frontLeavesNew
        real(kind=realType), dimension(3, 2) :: dummy
        integer nInterpol

        nInterpol = 0
        ! Set the ADT type, which is a surface ADT.
        ADT%adtType = adtSurfaceADT

        ! Copy the number of nodes and volume elements and set the number
        ! of surface elements to 0; only a volume grid has been given.

        ADT%nNodes = nPts
        ADT%nHexa = 0
        ADT%nTetra = 0
        ADT%nPyra = 0
        ADT%nPrisms = 0
        ADT%nTria = 0
        ADT%nQuads = nConn

        ! Set the pointers for the coordinates and the
        ! volume connectivities.

        ADT%coor => pts
        ADT%quadsConn => conn
        nullify (ADT%triaConn)
        ADT%nBBoxes = nConn

        ! Allocate the memory for the bounding box coordinates, the
        ! corresponding element type and the index in the connectivity.

        allocate (ADT%xBBox(6, nConn))
        allocate (ADT%elementType(nConn))
        allocate (ADT%elementID(nConn))

        ! All quads
        ADT%elementType = adtQuadrilateral

        ! Loop over the number of elements and store the bounding
        ! box info.
        nNPE = 4

        do i = 1, nConn

            mm = conn(1, i)
            xMin(1) = pts(1, mm); xMax(1) = pts(1, mm)
            xMin(2) = pts(2, mm); xMax(2) = pts(2, mm)
            xMin(3) = pts(3, mm); xMax(3) = pts(3, mm)

            do j = 2, nNPE
                mm = conn(j, i)

                xMin(1) = min(xMin(1), pts(1, mm))
                xMin(2) = min(xMin(2), pts(2, mm))
                xMin(3) = min(xMin(3), pts(3, mm))

                xMax(1) = max(xMax(1), pts(1, mm))
                xMax(2) = max(xMax(2), pts(2, mm))
                xMax(3) = max(xMax(3), pts(3, mm))
            end do

            ADT%xBBox(1, i) = xMin(1)
            ADT%xBBox(2, i) = xMin(2)
            ADT%xBBox(3, i) = xMin(3)

            ADT%xBBox(4, i) = xMax(1)
            ADT%xBBox(5, i) = xMax(2)
            ADT%xBBox(6, i) = xMax(3)

            ! elementID is just sequential since we only have 1 element type
            ADT%elementID(i) = i

        end do

        ! Build the ADT from the now known boundary boxes.

        call buildADT(ADT)

        ! Allocate the (pointer) memory that may be resized as necessary for
        ! the singlePoint search routine.
        allocate (BB(10), frontLeaves(25), frontLeavesNew(25), stack(100))

        ! Now do the searches
        do i = 1, nSearchPts
            coor(1:3) = searchPts(:, i)
            coor(4) = 1e20
            call minDistanceTreeSearchSinglePoint(ADT, coor, intInfo, uvw, dummy, &
                                                  nInterpol, BB, frontLeaves, frontLeavesNew)
            faceID(i) = intInfo(3)
            uv(:, i) = uvw(1:2)
        end do

        ! Release the memory
        deallocate (BB, frontLeaves, frontLeavesNew, stack)
        deallocate (ADT%xBBox)
        deallocate (ADT%elementType)
        deallocate (ADT%elementID)
        deallocate (ADT%ADTree)

    end subroutine searchQuads

    subroutine minDistanceTreeSearchSinglePoint(ADT, coor, intInfo, &
                                                uvw, arrDonor, nInterpol, BB, frontLeaves, frontLeavesNew)
        !
        !        This routine performs the actual minimum distance search for
        !        a single point on the local tree. It is local in the sens
        !        that no communication is involved. This routine does the
        !        actual search. The minDistanceTreeSearch is just a wrapper
        !        around this routine. The reason for the split is that the
        !        overset mesh connectivity requires efficient calling with
        !        a single coordinate. Therefore, this rouine does not
        !        allocate/deallocate any variables.
        !        Subroutine intent(in) arguments.
        !        --------------------------------
        !        ADT:       ADT type whose ADT must be searched
        !        coor:      The coordinates and the currently stored minimum
        !                   distance squared of these points:
        !                   coor(1): Coordinate 1.
        !                   coor(2): Coordinate 2.
        !                   coor(3): Coordinate 3.
        !                   coor(4): The currently stored minimum distance
        !                   squared.
        !        nInterpol: Number of variables to be interpolated.
        !        arrDonor:  Array with the donor data; needed to obtain the
        !                   interpolated data.
        !        Subroutine intent(out) arguments.
        !        ---------------------------------
        !        intInfo: 1D integer array, in which the following output
        !                 will be stored:
        !                 intInfo(1): processor ID of the processor where
        !                               the element is stored. This of course
        !                               is myID. If no element is found this
        !                               value is set to -1.
        !                 intInfo(2): The element type of the element.
        !                 intInfo(3): The element ID of the element in the
        !                               the connectivity.
        !        uvw:     2D floating point array to store the parametric
        !                 coordinates of the point in the transformed element
        !                 as well as the new distance squared and the
        !                 interpolated solution:
        !                 uvw(1): Parametric u-weight.
        !                 uvw(2): Parametric v-weight.
        !                 uvw(3): Parametric w-weight.
        !                 uvw(4): The new distance squared.
        !                 uvw(5): Interpolated solution, if desired. It is
        !                            possible to call this routine with
        !                            nInterpol == 0.
        !
        implicit none
        !
        !       Subroutine arguments.
        !
        type(adtType), intent(inout) :: ADT
        integer, intent(in) :: nInterpol

        real(kind=realType), dimension(4), intent(in) :: coor
        real(kind=realType), dimension(:, :), intent(in) :: arrDonor

        integer, dimension(3), intent(out) :: intInfo
        real(kind=realType), dimension(5), intent(out) :: uvw
        integer, dimension(:), pointer :: frontLeaves
        integer, dimension(:), pointer :: frontLeavesNew
        type(adtBBoxTargetType), dimension(:), pointer :: BB
        !
        !       Local parameters used in the Newton algorithm.
        !
        integer, parameter :: iterMax = 15
        real(kind=realType), parameter :: adtEps = 1.e-25_realType
        real(kind=realType), parameter :: thresConv = 1.e-10_realType
        !
        !       Local variables.
        !
        integer :: ierr

        integer :: ii, kk, ll, mm, nn, activeLeaf
        integer :: nBB, nFrontLeaves, nFrontLeavesNew
        integer :: nAllocBB, nAllocFront, nNodeElement
        integer :: i, kkk

        integer, dimension(8) :: n, m

        real(kind=realType) :: dx, dy, dz, d1, d2, invLen, val
        real(kind=realType) :: u, v, w, uv, uold, vold, vn, du, dv
        real(kind=realType) :: uu, vv, ww

        real(kind=realType), dimension(2) :: dd
        real(kind=realType), dimension(3) :: x1, x21, x41, x3142, xf
        real(kind=realType), dimension(3) :: vf, vt, a, b, norm, an, bn
        real(kind=realType), dimension(3) :: chi
        real(kind=realType), dimension(8) :: weight

        real(kind=realType), dimension(:, :), pointer :: xBBox

        logical :: elementFound
        type(adtLeafType), dimension(:), pointer :: ADTree

        ! Set some pointers to make the code more readable.

        xBBox => ADT%xBBox
        ADTree => ADT%ADTree

        ! Initial allocation of the arrays for the tree traversal as well
        ! as the stack array used in the qsort routine. The latter is
        ! done, because the qsort routine is called for every coordinate
        ! and therefore it is more efficient to allocate the stack once
        ! rather than over and over again. The disadvantage of course is
        ! that an essentially local variable, stack, is now stored in
        ! adtData.

        nAllocBB = size(BB)
        nAllocFront = size(frontLeaves)
        nStack = size(stack)

        ! Initialize the processor ID to -1 to indicate that no
        ! corresponding volume element is found and the new minimum
        ! distance squared to the old value.

        intInfo(1) = -1
        uvw(4) = coor(4)
        !
        !          Part 1. Determine the possible minimum distance squared to
        !                  the root leaf. If larger than the current distance
        !                  there is no need to search this tree.
        !
        if (coor(1) < ADTree(1)%xMin(1)) then
            dx = coor(1) - ADTree(1)%xMin(1)
        else if (coor(1) > ADTree(1)%xMax(4)) then
            dx = coor(1) - ADTree(1)%xMax(4)
        else
            dx = 0.0_realType
        end if

        if (coor(2) < ADTree(1)%xMin(2)) then
            dy = coor(2) - ADTree(1)%xMin(2)
        else if (coor(2) > ADTree(1)%xMax(5)) then
            dy = coor(2) - ADTree(1)%xMax(5)
        else
            dy = 0.0_realType
        end if

        if (coor(3) < ADTree(1)%xMin(3)) then
            dz = coor(3) - ADTree(1)%xMin(3)
        else if (coor(3) > ADTree(1)%xMax(6)) then
            dz = coor(3) - ADTree(1)%xMax(6)
        else
            dz = 0.0_realType
        end if

        ! Continue with the next coordinate if the possible distance
        ! squared to the root leaf is larger than the currently stored
        ! value.

        if ((dx * dx + dy * dy + dz * dz) >= uvw(4)) return
        !
        !          Part 2. Find a likely bounding box, which minimizes the
        !                  guaranteed distance.
        !
        activeLeaf = 1

        ! Traverse the tree until a terminal leaf is found.

        treeTraversal1: do

            ! Exit the loop when a terminal leaf has been found.
            ! This is indicated by a negative value of activeLeaf.

            if (activeLeaf < 0) exit treeTraversal1

            ! Determine the guaranteed distance squared for both children.

            do mm = 1, 2

                ! Determine whether the child contains a bounding box or
                ! another leaf of the tree.

                ll = ADTree(activeLeaf)%children(mm)
                if (ll > 0) then

                    ! Child contains a leaf. Determine the guaranteed distance
                    ! vector to the leaf.

                    d1 = abs(coor(1) - ADTree(ll)%xMin(1))
                    d2 = abs(coor(1) - ADTree(ll)%xMax(4))
                    dx = max(d1, d2)

                    d1 = abs(coor(2) - ADTree(ll)%xMin(2))
                    d2 = abs(coor(2) - ADTree(ll)%xMax(5))
                    dy = max(d1, d2)

                    d1 = abs(coor(3) - ADTree(ll)%xMin(3))
                    d2 = abs(coor(3) - ADTree(ll)%xMax(6))
                    dz = max(d1, d2)

                else

                    ! Child contains a bounding box. Determine the guaranteed
                    ! distance vector to it.

                    ll = -ll

                    d1 = abs(coor(1) - xBBox(1, ll))
                    d2 = abs(coor(1) - xBBox(4, ll))
                    dx = max(d1, d2)

                    d1 = abs(coor(2) - xBBox(2, ll))
                    d2 = abs(coor(2) - xBBox(5, ll))
                    dy = max(d1, d2)

                    d1 = abs(coor(3) - xBBox(3, ll))
                    d2 = abs(coor(3) - xBBox(6, ll))
                    dz = max(d1, d2)

                end if

                ! Compute the guaranteed distance squared for this child.

                dd(mm) = dx * dx + dy * dy + dz * dz

            end do

            ! Determine which will be the next leaf in the tree traversal.
            ! This will be the leaf which has the minimum guaranteed
            ! distance. In case of ties take the left leaf, because this
            ! leaf may have more children.

            if (dd(1) <= dd(2)) then
                activeLeaf = ADTree(activeLeaf)%children(1)
            else
                activeLeaf = ADTree(activeLeaf)%children(2)
            end if

        end do treeTraversal1

        ! Store the minimum of the just computed guaranteed distance
        ! squared and the currently stored value in uvw.

        uvw(4) = min(uvw(4), dd(1), dd(2))
        !
        !          Part 3. Find the bounding boxes whose possible minimum
        !                  distance is less than the currently stored value.
        !
        ! In part 1 it was already tested that the possible distance
        ! squared of the root leaf was less than the current value.
        ! Therefore initialize the current front to the root leaf and
        ! set the number of bounding boxes to 0.

        nBB = 0

        nFrontLeaves = 1
        frontLeaves(1) = 1

        ! Second tree traversal. Now to find all possible bounding
        ! box candidates.

        treeTraversal2: do

            ! Initialize the number of leaves for the new front, i.e.
            ! the front of the next round, to 0.

            nFrontLeavesNew = 0

            ! Loop over the leaves of the current front.

            currentFrontLoop: do ii = 1, nFrontLeaves

                ! Store the ID of the leaf a bit easier and loop over
                ! its two children.

                ll = frontLeaves(ii)

                childrenLoop: do mm = 1, 2

                    ! Determine whether this child contains a bounding box
                    ! or a leaf of the next level.

                    kk = ADTree(ll)%children(mm)
                    terminalTest: if (kk < 0) then

                        ! Child contains a bounding box. Determine the possible
                        ! minimum distance squared to this bounding box.

                        kk = -kk

                        if (coor(1) < xBBox(1, kk)) then
                            dx = coor(1) - xBBox(1, kk)
                        else if (coor(1) > xBBox(4, kk)) then
                            dx = coor(1) - xBBox(4, kk)
                        else
                            dx = 0.0_realType
                        end if

                        if (coor(2) < xBBox(2, kk)) then
                            dy = coor(2) - xBBox(2, kk)
                        else if (coor(2) > xBBox(5, kk)) then
                            dy = coor(2) - xBBox(5, kk)
                        else
                            dy = 0.0_realType
                        end if

                        if (coor(3) < xBBox(3, kk)) then
                            dz = coor(3) - xBBox(3, kk)
                        else if (coor(3) > xBBox(6, kk)) then
                            dz = coor(3) - xBBox(6, kk)
                        else
                            dz = 0.0_realType
                        end if

                        d2 = dx * dx + dy * dy + dz * dz

                        ! If this distance squared is less than the current
                        ! value, store this bounding box as a target.

                        testStoreBBox: if (d2 < uvw(4)) then

                            ! Check if the memory must be reallocated.

                            if (nBB == nAllocBB) &
                                call reallocBBoxTargetTypePlus(BB, nAllocBB, &
                                                               100, ADT)

                            ! Update the counter and store the data.

                            nBB = nBB + 1
                            BB(nBB)%ID = kk
                            BB(nBB)%posDist2 = d2

                            ! Although in step 2, i.e. the first tree traversal,
                            ! the guaranteed distance squared to a bounding box
                            ! has already been computed, this has been done only
                            ! for a likely candidate and not for all the possible
                            ! candidates. As this test is relatively cheap, do it
                            ! now for this bounding box.

                            d1 = abs(coor(1) - xBBox(1, kk))
                            d2 = abs(coor(1) - xBBox(4, kk))
                            dx = max(d1, d2)

                            d1 = abs(coor(2) - xBBox(2, kk))
                            d2 = abs(coor(2) - xBBox(5, kk))
                            dy = max(d1, d2)

                            d1 = abs(coor(3) - xBBox(3, kk))
                            d2 = abs(coor(3) - xBBox(6, kk))
                            dz = max(d1, d2)

                            d2 = dx * dx + dy * dy + dz * dz
                            uvw(4) = min(uvw(4), d2)

                        end if testStoreBBox

                        else terminalTest

                        ! Child contains a leaf. Compute the possible minimum
                        ! distance squared to the current coordinate.

                        if (coor(1) < ADTree(kk)%xMin(1)) then
                            dx = coor(1) - ADTree(kk)%xMin(1)
                        else if (coor(1) > ADTree(kk)%xMax(4)) then
                            dx = coor(1) - ADTree(kk)%xMax(4)
                        else
                            dx = 0.0_realType
                        end if

                        if (coor(2) < ADTree(kk)%xMin(2)) then
                            dy = coor(2) - ADTree(kk)%xMin(2)
                        else if (coor(2) > ADTree(kk)%xMax(5)) then
                            dy = coor(2) - ADTree(kk)%xMax(5)
                        else
                            dy = 0.0_realType
                        end if

                        if (coor(3) < ADTree(kk)%xMin(3)) then
                            dz = coor(3) - ADTree(kk)%xMin(3)
                        else if (coor(3) > ADTree(kk)%xMax(6)) then
                            dz = coor(3) - ADTree(kk)%xMax(6)
                        else
                            dz = 0.0_realType
                        end if

                        d2 = dx * dx + dy * dy + dz * dz

                        ! If this distance squared is less than the current
                        ! value, store this leaf in the new front.

                        testStoreLeave: if (d2 < uvw(4)) then

                            ! Check if enough memory has been allocated and
                            ! store the leaf.

                            if (nFrontLeavesNew == nAllocFront) then
                                i = nAllocFront
                                call reallocPlus(frontLeavesNew, i, 250, ADT)
                                call reallocPlus(frontLeaves, nAllocFront, 250, ADT)
                            end if

                            nFrontLeavesNew = nFrontLeavesNew + 1
                            frontLeavesNew(nFrontLeavesNew) = kk

                            ! Compute the guaranteed distance squared to this leaf.
                            ! It may be less than the currently stored value.

                            d1 = abs(coor(1) - ADTree(kk)%xMin(1))
                            d2 = abs(coor(1) - ADTree(kk)%xMax(4))
                            dx = max(d1, d2)

                            d1 = abs(coor(2) - ADTree(kk)%xMin(2))
                            d2 = abs(coor(2) - ADTree(kk)%xMax(5))
                            dy = max(d1, d2)

                            d1 = abs(coor(3) - ADTree(kk)%xMin(3))
                            d2 = abs(coor(3) - ADTree(kk)%xMax(6))
                            dz = max(d1, d2)

                            d2 = dx * dx + dy * dy + dz * dz
                            uvw(4) = min(uvw(4), d2)

                        end if testStoreLeave

                    end if terminalTest
                end do childrenLoop
            end do currentFrontLoop

            ! End of the loop over the current front. If the new front
            ! is empty the entire tree has been traversed and an exit is
            ! made from the corresponding loop.

            if (nFrontLeavesNew == 0) exit treeTraversal2

            ! Copy the data of the new front leaves into the current
            ! front for the next round.

            nFrontLeaves = nFrontLeavesNew
            do ll = 1, nFrontLeaves
                frontLeaves(ll) = frontLeavesNew(ll)
            end do

        end do treeTraversal2

        ! Sort the target bounding boxes in increasing order such that
        ! the one with the smallest possible distance is first.

        call qsortBBoxTargets(BB, nBB, ADT)
        !
        !          Part 4: Loop over the selected bounding boxes and check if
        !                  the corresponding element minimizes the distance.
        !
        elementFound = .false.

        BBoxLoop: do mm = 1, nBB

            ! Exit the loop if the possible minimum distance of this
            ! bounding box is not smaller than the current value.
            ! Remember that BB has been sorted in increasing order.

            if (uvw(4) <= BB(mm)%posDist2) exit BBoxLoop

            ! Determine the element type stored in this bounding box.

            kk = BB(mm)%ID
            select case (ADT%elementType(kk))

            case (adtTriangle)
                call adtTerminate(ADT, "minDistanceTreeSearch", &
                     "Minimum distance search for &
                     &triangles not implemented yet")

                !=========================================================

            case (adtQuadrilateral)

                ! Temporary implementation. I'm waiting for Juan to come
                ! up with his more sophisticated algorithm.

                ! This is a surface element, so set the parametric weight
                ! w to zero.

                w = 0.0_realType

                ! Determine the 4 vectors which completely describe
                ! the quadrilateral face

                ll = ADT%elementID(kk)
                n(1) = ADT%quadsConn(1, ll)
                n(2) = ADT%quadsConn(2, ll)
                n(3) = ADT%quadsConn(3, ll)
                n(4) = ADT%quadsConn(4, ll)

                x1(1) = ADT%coor(1, n(1))
                x1(2) = ADT%coor(2, n(1))
                x1(3) = ADT%coor(3, n(1))

                x21(1) = ADT%coor(1, n(2)) - x1(1)
                x21(2) = ADT%coor(2, n(2)) - x1(2)
                x21(3) = ADT%coor(3, n(2)) - x1(3)

                x41(1) = ADT%coor(1, n(4)) - x1(1)
                x41(2) = ADT%coor(2, n(4)) - x1(2)
                x41(3) = ADT%coor(3, n(4)) - x1(3)

                x3142(1) = ADT%coor(1, n(3)) - x1(1) - x21(1) - x41(1)
                x3142(2) = ADT%coor(2, n(3)) - x1(2) - x21(2) - x41(2)
                x3142(3) = ADT%coor(3, n(3)) - x1(3) - x21(3) - x41(3)

                ! Initialize u and v to 0.5 and determine the
                ! corresponding coordinates on the face, which is the
                ! centroid.

                u = 0.5_realType
                v = 0.5_realType
                uv = u * v

                xf(1) = x1(1) + u * x21(1) + v * x41(1) + uv * x3142(1)
                xf(2) = x1(2) + u * x21(2) + v * x41(2) + uv * x3142(2)
                xf(3) = x1(3) + u * x21(3) + v * x41(3) + uv * x3142(3)

                ! Newton loop to determine the point on the surface,
                ! which minimizes the distance to the given coordinate.

                NewtonQuads: do ll = 1, iterMax

                    ! Store the current values of u and v for a stop
                    ! criterion later on.

                    uold = u
                    vold = v

                    ! Determine the vector vf from xf to given coordinate.

                    vf(1) = coor(1) - xf(1)
                    vf(2) = coor(2) - xf(2)
                    vf(3) = coor(3) - xf(3)

                    ! Determine the tangent vectors in u- and v-direction.
                    ! Store these in a and b respectively.

                    a(1) = x21(1) + v * x3142(1)
                    a(2) = x21(2) + v * x3142(2)
                    a(3) = x21(3) + v * x3142(3)

                    b(1) = x41(1) + u * x3142(1)
                    b(2) = x41(2) + u * x3142(2)
                    b(3) = x41(3) + u * x3142(3)

                    ! Determine the normal vector of the face by taking the
                    ! cross product of a and b. Afterwards this vector will
                    ! be scaled to a unit vector.

                    norm(1) = a(2) * b(3) - a(3) * b(2)
                    norm(2) = a(3) * b(1) - a(1) * b(3)
                    norm(3) = a(1) * b(2) - a(2) * b(1)

                    invLen = 1.0_realType / max(adtEps, sqrt(norm(1) * norm(1) &
                                                             + norm(2) * norm(2) &
                                                             + norm(3) * norm(3)))

                    norm(1) = norm(1) * invLen
                    norm(2) = norm(2) * invLen
                    norm(3) = norm(3) * invLen

                    ! Determine the projection of the vector vf onto
                    ! the face.

                    vn = vf(1) * norm(1) + vf(2) * norm(2) + vf(3) * norm(3)
                    vt(1) = vf(1) - vn * norm(1)
                    vt(2) = vf(2) - vn * norm(2)
                    vt(3) = vf(3) - vn * norm(3)

                    ! The vector vt points from the current point on the
                    ! face to the new point. However this new point lies on
                    ! the plane determined by the vectors a and b, but not
                    ! necessarily on the face itself. The new point on the
                    ! face is obtained by projecting the point in the a-b
                    ! plane onto the face. this can be done by determining
                    ! the coefficients du and dv, such that vt = du*a + dv*b.
                    ! To solve du and dv the vectors normal to a and b
                    ! inside the plane ab are needed.

                    an(1) = a(2) * norm(3) - a(3) * norm(2)
                    an(2) = a(3) * norm(1) - a(1) * norm(3)
                    an(3) = a(1) * norm(2) - a(2) * norm(1)

                    bn(1) = b(2) * norm(3) - b(3) * norm(2)
                    bn(2) = b(3) * norm(1) - b(1) * norm(3)
                    bn(3) = b(1) * norm(2) - b(2) * norm(1)

                    ! Solve du and dv. the clipping of vn should not be
                    ! active, as this would mean that the vectors a and b
                    ! are parallel. This corresponds to a quad degenerated
                    ! to a line, which should not occur in the surface mesh.

                    vn = a(1) * bn(1) + a(2) * bn(2) + a(3) * bn(3)
                    vn = sign(max(adtEps, abs(vn)), vn)
                    du = (vt(1) * bn(1) + vt(2) * bn(2) + vt(3) * bn(3)) / vn

                    vn = b(1) * an(1) + b(2) * an(2) + b(3) * an(3)
                    vn = sign(max(adtEps, abs(vn)), vn)
                    dv = (vt(1) * an(1) + vt(2) * an(2) + vt(3) * an(3)) / vn

                    ! Determine the new parameter values uu and vv. These
                    ! are limited to 0 <= (uu,vv) <= 1.

                    u = u + du; u = min(1.0_realType, max(0.0_realType, u))
                    v = v + dv; v = min(1.0_realType, max(0.0_realType, v))

                    ! Determine the final values of the corrections.

                    du = abs(u - uold)
                    dv = abs(v - vold)

                    ! Determine the new coordinates of the point xf.

                    uv = u * v
                    xf(1) = x1(1) + u * x21(1) + v * x41(1) + uv * x3142(1)
                    xf(2) = x1(2) + u * x21(2) + v * x41(2) + uv * x3142(2)
                    xf(3) = x1(3) + u * x21(3) + v * x41(3) + uv * x3142(3)

                    ! Exit the loop if the update of the parametric
                    ! weights is below the threshold

                    val = sqrt(du * du + dv * dv)
                    if (val <= thresConv) exit NewtonQuads

                end do NewtonQuads

                ! Compute the distance squared between the given
                ! coordinate and the point xf.

                dx = coor(1) - xf(1)
                dy = coor(2) - xf(2)
                dz = coor(3) - xf(3)

                val = dx * dx + dy * dy + dz * dz

                ! If the distance squared is less than the current value
                ! store the wall distance and interpolation info and
                ! indicate that an element was found.

                if (val < uvw(4)) then
                    uvw(4) = val
                    nNodeElement = 4
                    elementFound = .true.

                    kkk = kk; uu = u; vv = v; ww = w
                    m(1) = n(1); m(2) = n(2); m(3) = n(3); m(4) = n(4)

                    weight(1) = (1.0_realType - u) * (1.0_realType - v)
                    weight(2) = u * (1.0_realType - v)
                    weight(3) = u * v
                    weight(4) = (1.0_realType - u) * v
                end if
            end select
        end do BBoxLoop

        ! Check if an element was found. As all the minimum distance
        ! searches are initialized by the calling routine (to support
        ! periodic searches) this is not always the case.

        if (elementFound) then

            ! Store the interpolation information for this point.  First
            ! the integer info, i.e. the processor ID (always 0), element
            ! type and local element ID.

            intInfo(1) = 0
            intInfo(2) = ADT%elementType(kkk)
            intInfo(3) = ADT%elementID(kkk)

            ! The parametric weights. Note that the wall distance
            ! squared, stored in the 4th position of uvw, already
            ! contains the correct value.

            uvw(1) = uu
            uvw(2) = vv
            uvw(3) = ww

            ! The interpolated solution, if needed.

            do ll = 1, nInterpol
                ii = 4 + ll
                uvw(ii) = weight(1) * arrDonor(ll, m(1))
                do i = 2, nNodeElement
                    uvw(ii) = uvw(ii) + weight(i) * arrDonor(ll, m(i))
                end do
            end do

        end if
    end subroutine minDistanceTreeSearchSinglePoint

end module ADTProjections
