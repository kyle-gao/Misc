        program kyle_makegamma
        implicit none
        logical more
        integer ( kind = 4 ) dim_num,level,j
	integer ( kind = 4 ), allocatable:: compositions(:,:),gamma(:,:)
  
      level=4
      dim_num=3
!calls and writes make_gamma
      allocate(gamma(1,dim_num))
      call make_gamma(level,dim_num,gamma)
      do j=1,(size(gamma)/dim_num)
  	write (*,*) gamma (j,:) 
      	end do
      contains
!subrotine make_gamma(p,m,gam)
!p=level,m=dimension number,gam(:,:)=input/output
!makes coefficient table gamma by appending compositions(i=1:p,n)
      subroutine make_gamma(p,m,gam)
	implicit none
	integer (kind=4) m,p,l1,l2,i
	integer ( kind = 4 ), allocatable::comp_list(:,:),temp(:,:)
	integer (kind=4),allocatable::gam(:,:)
	do i=1,p
	  allocate(comp_list(0,0))
	  call make_compositions(i,m,comp_list)
	  l1=size(gam)/m
	  l2=size(comp_list)/m
!	  create a temporary list of size(gam)+size(comp_list)
	  allocate (temp(l1+l2,m))
	  temp(:l1,:)=gam
	  temp(l1+1:l2,:)=comp_list
	  deallocate(comp_list)
!	  switch between gam and temp
	  call move_alloc(temp,gam) 
	  end do 
      end subroutine
!subroutine make_composition (n,k,comp_list)
!n=integer to be composed,k=how many parts,comp_list(:,:)=input/output
!composition of integer n into k parts like in mathemtica 
!returns comp_list, list of compositions of lenght n+k-1 choose k-1
!eg compositions(6,3) returns comp_list(28,3) the 28 compositions of 6 into 3 parts
      subroutine make_compositions (n,k,comp_list)
        implicit none
        logical more
        integer ( kind = 4 ) n,i,k,t,h
        integer ( kind = 4 ), dimension (k) :: level_1d
	integer ( kind = 4 ), allocatable::comp_list(:,:), temp(:,:)
      more = .false.
      h=0
      t=0
      i=0
      level_1d=(/0,0,0/)
      do while (level_1d(k) /= n)
	i=i+1
	allocate (temp(i,k))
	temp(:i-1,:)=comp_list(:i-1,:)
        call comp_next ( n, k, level_1d, more, h, t )        
	temp(i,:)=level_1d
	call move_alloc(temp,comp_list)
        end do
      return

      end subroutine
      
!subroutine comp_next as defined in 
      subroutine comp_next ( n, k, a, more, h, t )
      implicit none
      
      integer ( kind = 4 ) k,h,n,t,a(k)
      logical more

!
!  The first computation.
!
      if ( .not. more ) then
         
         t = n
         h = 0
         a(1) = n
         a(2:k) = 0
          !
!  The next computation.
!
      else
         
         if ( 1 < t ) then
            h = 0
         end if
         
         h = h + 1
         t = a(h)
         a(h) = 0
         a(1) = t - 1
         a(h+1) = a(h+1) + 1
         
      end if
!     
!     This is the last element of the sequence if all the
!     items are in the last slot.
!     
      more = ( a(k) /= n )
      
      return
      end subroutine
      end program
      
