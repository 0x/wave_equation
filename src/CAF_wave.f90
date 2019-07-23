    program wavedd
   
    implicit none
    
    integer, parameter :: jprocs = 3, iprocs = 3
    integer, parameter :: isrc=14, jsrc=14
    real(8), parameter :: Lx=1000, Ly=1000, t=2
    real(8), parameter :: CFL=0.5, c=1
    real(8), parameter :: dx=0.1, dy=0.1
    
    integer :: i, j
    integer :: irank
    real(8), dimension(:,:), allocatable :: wn, wnm1
    integer :: ista, iend, jsta, jend
    real(8) :: tt
    integer, parameter :: nxlocal=Lx/dx, nylocal=Ly/dy
    integer, parameter :: nx=nxlocal+2, ny=nylocal+2
    real(8), parameter :: dt=CFL*dx/c
    integer :: curimage(2)
    ! alloc coarray
    real(8), allocatable :: wnp1(:,:)[:,:]
    
    
    allocate(wn(nx,ny))
    allocate(wnp1(nx,ny)[iprocs,*])
    allocate(wnm1(nx,ny))
    
    curimage=this_image(wnp1)
    ! write(*,*) this_image(wnp1)
    sync all
  
    !! Calculation of the wave equation
    wn=0
    wnp1=wn
    tt=0
    
    do while (tt<t)  
    
        ! sync all
        
        if (curimage(1) > 1) &
            wnp1(1,:)=wnp1(nx-1,:)[curimage(1)-1,curimage(2)]
        if (curimage(1) < iprocs) &
            wnp1(nx,:)=wnp1(2,:)[curimage(1)+1,curimage(2)]
        if (curimage(2)>1) &
            wnp1(:,1)=wnp1(:,ny-1)[curimage(1),curimage(2)-1]
        if (curimage(2) < jprocs) &
            wnp1(:,ny)=wnp1(:,2)[curimage(1),curimage(2)+1] 
        
        sync all
        
        ista=2
        jsta=2
        iend=nx-1
        jend=ny-1
        if (curimage(1) == 1) then            
            wn(ista,:)=0           
            ista=ista+1
        endif
        if (curimage(1) == iprocs) then           
            wn(iend,:)=0           
            iend=iend-1
        endif
        if (curimage(2) == 1) then           
            wn(:,jsta)=0           
            jsta=jsta+1
        endif
        if (curimage(2) == jprocs) then           
            wn(:,jend)=0
            jend=jend-1
        endif       
       
        wnm1=wn
        wn=wnp1
        do i=ista,iend
            do j=jsta,jend
                wnp1(i,j)=2*wn(i,j)-wnm1(i,j)+CFL*CFL*(wn(i+1,j)+wn(i,j+1)-4*wn(i,j)+wn(i-1,j)+wn(i,j-1))
                if (isrc == (curimage(1)-1)*nxlocal+i-1 .and. jsrc == (curimage(2)-1)*nylocal+j-1) then
                    wnp1(i,j)=wnp1(i,j)+dt*dt*20*sin(30*3.14*tt/20)
                endif
            enddo
        enddo
        tt=tt+dt
    enddo
     
    ! Display and write the result 
    !do irank=1, num_images()
    !    if (irank==this_image()) then
    !        write(*,*) "image: ", this_image()
    !        write(*,*) wn(2:nx-1,2:ny-1)
    !        open(unit = 1 , file = "wavedd.txt", position="append")
    !        write(1,*) wn(2:nx-1,2:ny-1)
    !        close(1)
    !   endif
    !   sync all
    !enddo
    
    end program wavedd
