module maths
  use common
  implicit none

  interface median
     module procedure median_s, median_d
  end interface median

  interface mean
     module procedure mean_s, mean_d,mean_sw, mean_dw
  end interface mean
  
  interface chi2_norm
     module procedure chi2_norm_s, chi2_norm_d
  end interface chi2_norm
  
  interface percentile
     module procedure percentile_s, percentile_d
  end interface percentile

  private
  public :: median,mean,chi2_norm,percentile,qselect,dqselect
  
contains
  SUBROUTINE median_s (vec_in, med, mederr)
    !     median the vector vec, return the median &med and the 68°.27 
    !     percentile median error mederr. Now powered by QuickSelect! 
    !     (see F77 Numerical Recipes

    INTEGER :: veclen, i,j
    REAL, dimension(:), intent(in) :: vec_in
    REAL, dimension(size(vec_in,1)) :: vec, err
    REAL  :: med, mederr, temp
    
    vec=vec_in    
    veclen = size(vec)
    
    if (veclen .lt. 1) then
       write(*,*) 'ERROR: attempting to median a <0 long vector'
       STOP
    endif

    if (veclen .eq. 1) then
       med = vec(1) ; mederr = -1. ; return
    endif

    if (mod(veclen,2) .eq. 1)  med = qselect(veclen/2 + 1, veclen, vec)    
    if (mod(veclen,2) .eq. 0)  med = (qselect(veclen/2,veclen,vec) + qselect(veclen/2+1,veclen,vec) )/ 2.
    err = abs(vec - med)
    mederr = qselect(ifix(0.6827*float(veclen)+0.500),veclen,err)
    
    return
  END SUBROUTINE median_s
  
  SUBROUTINE median_d (vec_in, med, mederr)
    !     median the vector vec, return the median &med and the 68�.27 
    !     percentile median error mederr. Now powered by QuickSelect! 
    !     (see F77 Numerical Recipes
    
    INTEGER :: veclen, i,j
    REAL (kind=PR), dimension(:),intent(in) :: vec_in  
    REAL (kind=PR), dimension(size(vec_in,1)) :: vec, err
    REAL (kind=PR) :: med, mederr, temp

    veclen = size(vec_in,1)
    if (veclen .lt. 1) then
       write(*,*) 'ERROR: attempting to median a <0 long vector'
       STOP
    endif

    if (veclen .eq. 1) then
       med = vec_in(1) ;  mederr = -1. ; return
    endif
    
    veclen = size(vec)
    vec=vec_in

!!$    ! I'm not sure that qselect is still working with double precision numbers....
!!$    !     OLD (EXPENSIVE) METHOD
!!$
!!$    do i=1, veclen - 1
!!$       do j=i+1, veclen 
!!$          if (vec(j) .lt. vec(i)) then
!!$             temp   = vec(i) ;  vec(i) = vec(j) ;  vec(j) = temp
!!$          endif
!!$       enddo
!!$    enddo
!!$    
!!$    !     veclen odd: median is the central element
!!$    if (mod(veclen,2) .eq. 1) med = vec(veclen/2 + 1)
!!$
!!$    !     veclen even: median is the average of the 2 middle elements
!!$    if (mod(veclen,2) .eq. 0) med = (vec(veclen/2) + vec(veclen/2+1)) / 2.
!!$
!!$    do i=1, veclen
!!$       err(i) = abs(vec(i) - med)
!!$    enddo
!!$
!!$    do i=1, veclen - 1
!!$       do j=i+1, veclen
!!$          if (err(j) .lt. err(i)) then
!!$             temp   = err(i) ;  err(i) = err(j) ;  err(j) = temp
!!$          endif
!!$       enddo
!!$    enddo
!!$
!!$    !     the median error contain the 68.27% of the deviations nint(number) rounds $number to the nearest integer
!!$    mederr = err(nint(0.6827*float(veclen)))
    
    if (mod(veclen,2) .eq. 1)  med = dqselect(veclen/2 + 1, veclen, vec)    
    if (mod(veclen,2) .eq. 0)  med = (dqselect(veclen/2,veclen,vec) + dqselect(veclen/2+1,veclen,vec) )/ 2.
    err = abs(vec - med)
    mederr = dqselect(ifix(0.6827*float(veclen)+0.500),veclen,err)

    return
  END SUBROUTINE median_d

  subroutine mean_s (vec, avg, p68, std)
    REAL, dimension(:), intent(in) :: vec
    REAL, intent(out) :: avg, p68, std
    REAL (kind=PR) :: avg_out, p68_out, std_out

    call mean_d(dble(vec),avg_out,p68_out,std_out)
    avg = real(avg_out)
    p68 = real(p68_out)
    std = real(std_out)
    return
  end subroutine mean_s

  SUBROUTINE mean_d (vec, avg, p68, std)
    !     mean the vector vec, return the mean &avg and the 68°.27 
    !     percentile mean error avgerr. Now powered by QuickSelect! 
    !     (see F77 Numerical Recipes, page 334) 

    INTEGER :: veclen, i
    REAL (kind=PR), dimension(:), intent(in) :: vec
    REAL (kind=PR), intent(out) :: avg, p68, std
    real, dimension(size(vec,1)) :: terr
    REAL (kind=PR) :: A,Q,A_prv,Q_prv

    veclen = size(vec,1)
    if (veclen .lt. 1) then
       write(*,*) 'ERROR: attempting to mean a <0 long vector'
       STOP
    endif
    
    avg = 0.0D0 ;  std = 0.0D0 ;  p68 = 0.0D0

    if (veclen .eq. 1) then
       avg = vec(1) ;  std = -1.0D0 ;  p68 = -1.0D0
       return
    endif
    
    A_prv = 0.0D0 ;  A = 0.0D0
    Q_prv = 0.0D0 ;  Q = 0.0D0  
    do i=1,veclen
       A = A_prv + (vec(i)-A_prv)/real(i)
       Q = Q_prv + (vec(i)-A_prv)*(vec(i)-A)
       A_prv = A
       Q_prv = Q
    end do
    avg = A
    std = sqrt(Q/real(veclen))
    
    terr = real(abs(vec - avg))
    p68 = dble(qselect(nint(0.6827*real(veclen)),veclen,terr))
    
    return
  END SUBROUTINE mean_d

  SUBROUTINE mean_sw(vec_in, weight_in, avg, p68, std)
    !     mean the vector vec, return the mean &avg and the 68.27 
    !     percentile mean error avgerr. Now powered by QuickSelect! 
    !     (see F77 Numerical Recipes, page 334) 

    REAL, dimension(:), intent(in) :: vec_in, weight_in
    REAL, intent(out) :: avg, p68, std
    REAL (kind=PR) :: avg_out, p68_out, std_out
    
    call mean_dw(dble(vec_in), dble(weight_in), avg_out, p68_out, std_out)
    avg = real(avg_out)
    p68 = real(p68_out)
    std = real(std_out)
    return
  end SUBROUTINE mean_sw
  
  SUBROUTINE mean_dw(vec_in, weight_in, avg, p68, std)
    !     mean the vector vec, return the mean &avg and the 68°.27 
    !     percentile mean error avgerr. Now powered by QuickSelect! 
    !     (see F77 Numerical Recipes, page 334) 


    INTEGER :: veclen, i, n
    REAL (kind=PR), dimension(:), intent(in) :: vec_in, weight_in
    REAL (kind=PR), dimension(size(vec_in,1)) :: vec, weight
    REAL , dimension(size(vec_in,1)) :: terr
    REAL (kind=PR), intent(out)  :: avg, p68, std
    REAL (kind=PR) :: A,Q,W,A_prv,Q_prv,W_prv
    
    veclen = size(vec_in,1)
    if (size(weight_in).ne.veclen) stop
    vec=vec_in
    weight = weight_in
    
    avg = 0.0D0 ;  std = 0.0D0 ;  p68 = 0.0D0
    
    if (veclen .lt. 1) then
       write(*,*) 'ERROR: attempting to mean a <0 long vector'
       STOP
    endif

    if (veclen .eq. 1) then
       avg = vec(1) ;  std = -1.0D0 ;  p68 = -1.0D0
       return
    end if
    
    A_prv = 0.0D0 ;  A = 0.0D0
    Q_prv = 0.0D0 ;  Q = 0.0D0  
    W_prv = 0.0D0 ;  W = 0.0D0
    do i = 1,veclen
       if (weight(i).le. 0.00d0) cycle
       W = W_prv + weight(i)
       A = A_prv + weight(i)/W*(vec(i)-A_prv)
       Q = Q_prv + weight(i)*(vec(i)-A_prv)*(vec(i)-A)
       W_prv = W
       A_prv = A
       Q_prv = Q
    end do
    avg = A
    std = sqrt(Q/W)
    
    n=0
    do i=1,veclen
       if (weight(i).gt.0.0) then
          n = n+1
          terr(n) = abs(vec_in(i)-avg)
       end if
    end do
    p68 = dble(qselect(nint(0.6827*real(n)),n,terr(1:n)))
    
    return
  END SUBROUTINE mean_dw
  
  SUBROUTINE chi2_norm_s (vec_in, err_in, mod_in, chi2_out, chi2_norm_out)

    REAL, dimension(:), intent(in) :: vec_in, mod_in, err_in
    REAL, dimension(size(vec_in,1)) :: vec, mod, err
    REAL, intent(out)  :: chi2_out, chi2_norm_out
    REAL :: n_tot
    
    chi2_out = -9.9
    chi2_norm_out = -9.9

    mod = 1._4
    n_tot = sum(mod, dim=1, mask=(err_in.gt.0._4))

    if (n_tot.lt.3) return
    
    vec = vec_in
    err = err_in
    mod = mod_in
    
    where (err_in.le.0._4)  
       vec = mod
       err = 999999._4
    end where
    
    chi2_out = sum( ((vec-mod)/err)**2,dim=1)
    chi2_norm_out =  sum( ((vec-mod)/err)**2,dim=1)/(n_tot-1._4)
    return
  
  END SUBROUTINE chi2_norm_s

  SUBROUTINE chi2_norm_d (vec_in, err_in, mod_in, chi2_out, chi2_norm_out)

    REAL (kind=PR), dimension(:), intent(in) :: vec_in, mod_in, err_in
    REAL (kind=PR), dimension(size(vec_in,1)) :: vec, mod, err
    REAL (kind=PR), intent(out)  :: chi2_out, chi2_norm_out
    REAL (kind=PR) :: n_tot
    
    chi2_out = -9.9
    chi2_norm_out = -9.9

    mod = 1._PR
    n_tot = sum(mod, dim=1, mask=(err_in.gt.0._PR))

    if (n_tot.lt.3) return
    
    vec = vec_in
    err = err_in
    mod = mod_in
    
    where (err_in.le.0._PR)  
       vec = mod
       err = 999999._PR
    end where
    
    chi2_out = sum( ((vec-mod)/err)**2,dim=1)
    chi2_norm_out =  sum( ((vec-mod)/err)**2,dim=1)/(n_tot-1._PR)
    
    return
  
  END SUBROUTINE chi2_norm_d
  
  real FUNCTION qselect(k,n,arr)
    !     Returns the kth smallest value in the array arr(1:n). The input 
    !     array will be rearranged to have this value in location arr(k), 
    !     with all smaller elements moved to arr(1:k-1) (in arbitrary order)
    !     and all larger elements in arr[k+1..n] (also in arbitrary order).

    INTEGER k,n
    REAL arr(*)
    INTEGER i,ir,j,l,mid
    REAL a,temp

    l=1
    ir=n

1   if(ir-l.le.1)then 
       if(ir-l.eq.1)then 
          if(arr(ir).lt.arr(l))then
             temp=arr(l)
             arr(l)=arr(ir)
             arr(ir)=temp
          endif
       endif
       qselect=arr(k)
       return
    else
       mid=(l+ir)/2 
       temp=arr(mid)
       arr(mid)=arr(l+1)
       arr(l+1)=temp
       if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
       endif
       if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
       endif
       if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
       endif
       i=l+1 
       j=ir
       a=arr(l+1) 
3      continue 
       i=i+1 
       if(arr(i).lt.a)goto 3
4      continue
       j=j-1 
       if(arr(j).gt.a)goto 4
       if(j.lt.i)goto 5 
       temp=arr(i) 
       arr(i)=arr(j)
       arr(j)=temp
       goto 3 
5      arr(l+1)=arr(j) 
       arr(j)=a
       if(j.ge.k)ir=j-1 
       if(j.le.k)l=i
    endif
    goto 1


  end FUNCTION qselect

  real function percentile_s(vec_in,n)
    real, dimension(:), intent(in) :: vec_in
    REAL, dimension(size(vec_in,1)) :: vec
    integer, intent(in) :: n
    integer :: size_tot,pos_perc
    vec = vec_in
    size_tot = size(vec,dim=1)
    pos_perc = nint(real(n)/100._4*real(size_tot))
    if (pos_perc.le.0) pos_perc=1
    if (pos_perc.gt.size_tot) pos_perc=size_tot
    if (n.le.1) pos_perc = 1
    if (n.eq.100) pos_perc = size_tot
    percentile_s = qselect(pos_perc,size_tot,vec)
    return
  end function percentile_s
  
  real (kind=PR) function percentile_d(vec_in,n)
    real (kind=PR), dimension(:), intent(in) :: vec_in
    REAL (kind=PR), dimension(size(vec_in,1)) :: vec
    integer, intent(in) :: n
    integer :: size_tot,pos_perc
    vec = vec_in
    size_tot = size(vec,dim=1)
    pos_perc = nint(real(n)/100._4*real(size_tot))
    if (pos_perc.le.0) pos_perc=1
    if (pos_perc.gt.size_tot) pos_perc=size_tot
    if (n.le.1) pos_perc = 1
    if (n.eq.100) pos_perc = size_tot
    percentile_d = dqselect(pos_perc,size_tot,vec)
    return
  end function percentile_d
  
  real (kind=PR) FUNCTION dqselect(k,n,arr)
    !     Returns the kth smallest value in the array arr(1:n). The input 
    !     array will be rearranged to have this value in location arr(k), 
    !     with all smaller elements moved to arr(1:k-1) (in arbitrary order)
    !     and all larger elements in arr[k+1..n] (also in arbitrary order).

    INTEGER k,n
    REAL (kind=PR) :: arr(*)
    INTEGER i,ir,j,l,mid
    REAL (kind=PR) :: a,temp

    l=1
    ir=n

1   if(ir-l.le.1)then 
       if(ir-l.eq.1)then 
          if(arr(ir).lt.arr(l))then
             temp=arr(l)
             arr(l)=arr(ir)
             arr(ir)=temp
          endif
       endif
       dqselect=arr(k)
       return
    else
       mid=(l+ir)/2 
       temp=arr(mid)
       arr(mid)=arr(l+1)
       arr(l+1)=temp
       if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
       endif
       if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
       endif
       if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
       endif
       i=l+1 
       j=ir
       a=arr(l+1) 
3      continue 
       i=i+1 
       if(arr(i).lt.a)goto 3
4      continue
       j=j-1 
       if(arr(j).gt.a)goto 4
       if(j.lt.i)goto 5 
       temp=arr(i) 
       arr(i)=arr(j)
       arr(j)=temp
       goto 3 
5      arr(l+1)=arr(j) 
       arr(j)=a
       if(j.ge.k)ir=j-1 
       if(j.le.k)l=i
    endif
    goto 1


  end FUNCTION dqselect




!!$   ROUTINE NEVER TESTED
!!$
!!$      SUBROUTINE mode(pdata,ndata,xhistmax)
!!$      IMPLICIT NONE
!!$      INTEGER ndata,nbin,i,k
!!$      REAL pdata(_VRES_*2)
!!$      REAL med, mederr
!!$      REAL pmin,pmax
!!$
!!$      REAL bin
!!$      REAL xhist(100),yhist(100), yhistmax, xhistmax, kmax
!!$
!!$      pmin=pdata(1)
!!$      pmax=pdata(1)
!!$      do i=2,ndata
!!$         if (pdata(i).gt.pmax) pmax=pdata(i)
!!$         if (pdata(i).lt.pmin) pmin=pdata(i)
!!$      enddo
!!$
!!$      call median(pdata(1:ndata),ndata,med,mederr)
!!$
!!$      nbin=int(sqrt(ndata*1.0))
!!$      nbin=20
!!$      if (nbin.gt.100) nbin=100
!!$      bin=6*mederr/nbin
!!$      
!!$!      write(*,*) 'MED:',med,mederr, nbin,bin
!!$
!!$      xhist(1)=med-3*mederr-bin/2
!!$      do k=2,nbin
!!$         xhist(k)=xhist(1)+(k-1)*bin
!!$      enddo
!!$      yhist=0
!!$      ! very stupoid way to create histogram
!!$      do k=1,nbin
!!$         do i=1,ndata       
!!$            if ( pdata(i).ge. xhist(k)-bin/2 .and.
!!$     *           pdata(i).lt. xhist(k)+bin/2 )
!!$     *           yhist(k)=yhist(k)+1
!!$         enddo
!!$      enddo
!!$      
!!$      kmax=1
!!$      yhistmax=yhist(1)
!!$      do k=1,nbin
!!$!         write(*,*) 'x-y hist', xhist(k),yhist(k)
!!$
!!$         if (yhist(k).gt.yhistmax) then
!!$            yhistmax=yhist(k)
!!$            kmax=k
!!$         endif
!!$      enddo
!!$      xhistmax=xhist(kmax)
!!$!      write(*,*) xhistmax
!!$
!!$      RETURN
!!$      END SUBROUTINE
!!$      


end module maths
