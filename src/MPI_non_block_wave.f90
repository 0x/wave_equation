    
   program wavedd
    
    use mpi
    
    implicit none
    
    integer :: error, nprocs, rank
    integer :: istatus(MPI_STATUS_SIZE)
	  
    integer, parameter :: iprocs = 3, jprocs = 3
    integer, parameter :: isrc=14, jsrc=14
    real(8), parameter :: Lx=1000, Ly=1000, t=2
    real(8), parameter :: CFL=0.5, c=1
    real(8), parameter :: dx=0.1, dy=0.1
    
    integer :: i, j
    integer :: itable(-1:iprocs, -1:jprocs)
    integer :: irank, myranki, myrankj
    integer :: inext, jnext, iprev, jprev
    real(8), dimension(:,:), allocatable :: wn, wnm1, wnp1
    real(8), dimension(:), allocatable :: row1s, row2s, row1r, row2r
    integer :: isend1, isend2, jsend1, jsend2
    integer :: irecv1, irecv2, jrecv1, jrecv2
    integer :: ista, iend, jsta, jend
    real(8) :: tt
    integer, parameter :: nxlocal=Lx/dx, nylocal=Ly/dy
    integer, parameter :: nx=nxlocal+2, ny=nylocal+2
    real(8), parameter :: dt=CFL*dx/c
    ! MPI I/O
    ! integer(kind=MPI_OFFSET_KIND), parameter :: zero_off = 0
    ! integer, dimension(mpi_status_size) :: wstatus
    
    
    call mpi_init(error)
    call mpi_comm_size(MPI_COMM_WORLD,nprocs,error)
    call mpi_comm_rank(MPI_COMM_WORLD,rank,error)
    ! write(*,'("np = ",i2,2x,"id = ",i2)') nprocs,rank
    if (nprocs /= iprocs*jprocs) then
        write(*, *) 'ERROR, mpi_comm_size: ', nprocs
        stop
    endif
    
    ! Itable interface 
    do j = -1, jprocs
        do i = -1, iprocs
            itable(i,j) = MPI_PROC_NULL
        enddo
    enddo
    irank = 0
    do i = 0, iprocs-1
        do j = 0, jprocs-1
            itable(i,j) = irank
            if (rank == irank) then
                myranki = i;
                myrankj = j
            endif
            irank = irank + 1
        enddo
    enddo
    jnext = itable(myranki, myrankj + 1)
    jprev = itable(myranki, myrankj - 1)
    inext = itable(myranki+1, myrankj)
    iprev = itable(myranki-1, myrankj)
    
    allocate(wn(nx,ny))
    allocate(wnp1(nx,ny))
    allocate(wnm1(nx,ny))
    allocate(row1s(ny))
    allocate(row2s(ny))
    allocate(row1r(ny))
    allocate(row2r(ny))
    
    ! Calculation of the wave equation
    wn=0
    wnp1=wn
    tt=0
    
    ! MPI I/O
    ! istart = myrank*nxlocal+1
    ! iend = (myrank + 1)*nxlocal
    ! create subarrays
    ! call mpi_type_create_subarray( 2, [iprocs*nxlocal,jprocs*nylocal], [iend-istart+1,jend-jstart+1], &
    !                           [istart,jstart], mpi_order_fortran, mpi_real, localarray, error )
    ! call mpi_type_commit( localarray, error )

do while (tt<t)
        if (myranki /= 0) then
            row1s=wnp1(2,:)
        endif
        if (myranki /= iprocs-1) then
            row2s=wnp1(nx-1,:)
        endif

        call mpi_isend(wnp1(:,ny-1), nx, MPI_REAL8, jnext, 1, MPI_COMM_WORLD, isend1, error)
        call mpi_isend(wnp1(:,2), nx, MPI_REAL8, jprev, 1, MPI_COMM_WORLD, isend2, error)
        call mpi_isend(row2s, ny, MPI_REAL8, inext, 1, MPI_COMM_WORLD, jsend1, error)
        call mpi_isend(row1s, ny, MPI_REAL8, iprev, 1, MPI_COMM_WORLD, jsend2, error)
        call mpi_irecv(wnp1(:,1), nx, MPI_REAL8, jprev, 1, MPI_COMM_WORLD, irecv1, error)
        call mpi_irecv(wnp1(:,ny), nx, MPI_REAL8, jnext, 1, MPI_COMM_WORLD, irecv2, error)
        call mpi_irecv(row1r, ny, MPI_REAL8, iprev, 1, MPI_COMM_WORLD, jrecv1, error)
        call mpi_irecv(row2r, ny, MPI_REAL8, inext, 1, MPI_COMM_WORLD, jrecv2, error)

        ista=3
        jsta=3
        iend=nx-2
        jend=ny-2
        if (myranki == 0) then
            ista=ista+1
        endif
        if (myranki == iprocs - 1) then
            iend=iend-1
        endif
        if (myrankj == 0) then
            jsta=jsta+1
        endif
        if (myrankj == jprocs - 1) then
            jend=jend-1
        endif
        wnm1=wn
        wn=wnp1
        do i=ista,iend
            do j=jsta,jend
                wnp1(i,j)=2*wn(i,j)-wnm1(i,j)+CFL*CFL*(wn(i+1,j)+wn(i,j+1)-4*wn(i,j)+wn(i-1,j)+wn(i,j-1))
                if (isrc == myranki*nxlocal+i-1 .and. jsrc == myrankj*nylocal+j-1) then
                    wnp1(i,j)=wnp1(i,j)+dt*dt*20*sin(30*3.14*tt/20)
                endif
            enddo
        enddo

        call mpi_wait(isend1, istatus, error)
        call mpi_wait(isend2, istatus, error)
        call mpi_wait(jsend1, istatus, error)
        call mpi_wait(jsend2, istatus, error)
        call mpi_wait(irecv1, istatus, error)
        call mpi_wait(irecv2, istatus, error)
        call mpi_wait(jrecv1, istatus, error)
        call mpi_wait(jrecv2, istatus, error)

        if (myranki /= 0) then
            wnp1(1,:) = row1r
        endif
        if (myranki /= iprocs-1) then
            wnp1(nx,:) = row2r
        endif

        ista=2
        jsta=2
        iend=nx-1
        jend=ny-1
        if (myranki == 0) then
            wn(ista,:)=0
            ista=ista+1
        endif
        if (myranki == iprocs - 1) then
            wn(iend,:)=0
            iend=iend-1
        endif
        if (myrankj == 0) then
            wn(:,jsta)=0
            jsta=jsta+1
        endif
        if (myrankj == jprocs - 1) then
            wn(:,jend)=0
            jend=jend-1
        endif

        wn(1,:)=wnp1(1,:)
        wn(nx,:)=wnp1(nx,:)
        wn(:,1)=wnp1(:,1)
        wn(:,ny)=wnp1(:,ny)

        i=ista
        do j=jsta,jend
            wnp1(i,j)=2*wn(i,j)-wnm1(i,j)+CFL*CFL*(wn(i+1,j)+wn(i,j+1)-4*wn(i,j)+wn(i-1,j)+wn(i,j-1))
            if (isrc == myranki*nxlocal+i-1 .and. jsrc == myrankj*nylocal+j-1) then
                wnp1(i,j)=wnp1(i,j)+dt*dt*20*sin(30*3.14*tt/20)
            endif
        enddo
        i=iend
        do j=jsta,jend
            wnp1(i,j)=2*wn(i,j)-wnm1(i,j)+CFL*CFL*(wn(i+1,j)+wn(i,j+1)-4*wn(i,j)+wn(i-1,j)+wn(i,j-1))
            if (isrc == myranki*nxlocal+i-1 .and. jsrc == myrankj*nylocal+j-1) then
                wnp1(i,j)=wnp1(i,j)+dt*dt*20*sin(30*3.14*tt/20)
            endif
        enddo
        j=jsta
        do i=ista,iend
            wnp1(i,j)=2*wn(i,j)-wnm1(i,j)+CFL*CFL*(wn(i+1,j)+wn(i,j+1)-4*wn(i,j)+wn(i-1,j)+wn(i,j-1))
            if (isrc == myranki*nxlocal+i-1 .and. jsrc == myrankj*nylocal+j-1) then
                wnp1(i,j)=wnp1(i,j)+dt*dt*20*sin(30*3.14*tt/20)
            endif
        enddo
        j=jend
        do i=ista,iend
            wnp1(i,j)=2*wn(i,j)-wnm1(i,j)+CFL*CFL*(wn(i+1,j)+wn(i,j+1)-4*wn(i,j)+wn(i-1,j)+wn(i,j-1))
            if (isrc == myranki*nxlocal+i-1 .and. jsrc == myrankj*nylocal+j-1) then
                wnp1(i,j)=wnp1(i,j)+dt*dt*20*sin(30*3.14*tt/20)
            endif
        enddo

        tt=tt+dt
    enddo
	
    ! MPI I/O
    ! write to file
    ! call mpi_file_open( mpi_comm_world, 'test.txt', IOR(MPI_mode_create,MPI_mode_wronly), &
    !              mpi_info_null, fileno, error )
    ! call mpi_file_set_view( fileno, zero_off, mpi_real, localarray, "native", mpi_info_null, error )
    ! call mpi_file_write_all( fileno, wn, (jend-jstart+1)*(iend-istart+1), MPI_real, wstatus, error )
    ! call mpi_file_close( fileno, error )
    
    ! Display and write the result 
    !do irank=0, nprocs-1
    !    if (rank==irank) then
    !        write(*,*) "rank: ", rank
    !        write(*,*) wn(2:nx-1,2:ny-1)
    !        
    !        open(unit = 1 , file = "wavedd.txt", position="append")
    !        write(1,*) wn(2:nx-1,2:ny-1)
    !        close(1)        
    !    endif
    !    call mpi_barrier(MPI_COMM_WORLD, error)
    !enddo
    
    call mpi_finalize(error)
    end program wavedd
