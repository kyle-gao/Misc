        program kyle_makegamma
        implicit none
        logical more
        integer ( kind = 4 ) dim_num,level,j
	integer(kind=4),allocatable::gamma(:,:),box(:,:)

	integer (kind=4),allocatable::test_less(:)
	integer ( kind = 4 ):: k(3),jtest(4,3)
!test
!      k=[4,5,1]
!      jtest=transpose(reshape([3,4,1,1,2,3,5,1,1,1,1,1],[3,4]))

      level=4
      dim_num=3
      allocate(gamma(1,dim_num))
      gamma=0
      call make_gamma(level,dim_num,gamma)
      allocate(test_less(0))

      write(*,*) 'GAMMA'
      do j=1,size(gamma,1)
	write (*,*) gamma(j,:)
      end do

      allocate(box(0,0))
      call make_box(gamma,box)
      write(*,*) 'box(8,:)'
      write(*,*)  box(8,:)
      write(*,*) 'gamma(8,:)'
      write(*,*) gamma(8,:)
      contains
!subroutine make_box builds a coefficient box out of a gamma 
      subroutine make_box(gamma,out_box)
	implicit none
	integer (kind=4)::gamma(:,:),i,l,l_out,l_less
	integer (kind=4),allocatable::temp(:,:),less_list(:),out_box(:,:)
	l=size(gamma,1)
	do i=1,l
	  l_out=size(out_box,1)
	  allocate(less_list(0))
          less_list=0
	  call less_than_k(gamma,gamma(i,:),less_list)
	  allocate(temp(i,l))
	  temp=0
	  l_less=size(less_list)
	  temp(:i-1,:)=out_box(:i-1,:)
	  temp(i,:)=less_list 
	  call move_alloc(temp,out_box)
          deallocate(less_list)	  
	end do
      end subroutine
!subroutine less_than_K(j_list,k,out_less) takes input list of vectors j_list like k, and vector 
!k, returns a vector with 1 in the i'th component if the i'th vector of j_list 
!has vector component all less or equal to the corresponding components of k.
!out_less is passed as output
      subroutine less_than_k(j_list,k,out_less)
	implicit none
	integer (kind=4) i,j,len_j,len_k,len_out,j_list(:,:),k(:)
	integer (kind=4),allocatable::temp1(:),temp2(:),out_less(:)
	len_k=size(k)
	len_j=size(j_list,1)
	do i=1,len_j
	  allocate(temp1(len_j))
	  temp1=0
	  do j=1,len_k
	    if (j_list(i,j)<=k(j)) then
		temp1(j)=0
	    else 
		temp1(j)=1
	    end if
	    end do
          
	  allocate(temp2(len_j))
          temp2=0
          temp2(:i-1)=out_less(:i-1)
	  if (sum(temp1)==0) then	    
            temp2(i)=1
	  else 
	    temp2(i)=0
	  end if	  
	  call move_alloc(temp2,out_less)
	  deallocate(temp1)
	  end do 
      end subroutine
!calls and writes make_gamma
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
          comp_list=0
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
        temp=0
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
!  The first computation.
      if ( .not. more ) then
         
         t = n
         h = 0
         a(1) = n
         a(2:k) = 0
!  The next computation.
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
!     This is the last element of the sequence if all the
!     items are in the last slot.
!     
      more = ( a(k) /= n )
      
      return
      end subroutine
      end program
      
