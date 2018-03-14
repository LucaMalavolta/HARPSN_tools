
module spectrum_rebin
  use common
  implicit none
contains
  subroutine rebin_keepnorm(WL,fl,wl_rb,fl_rb,verbose)
    ! assuming assumed-shape variables
    integer, optional, intent(in) :: verbose
    integer ::  n1,n2,k,kk,n,k1,k2,i1,i2,kt
    real (PD), dimension(:)  ::  fl, fl_rb
    real (PD), dimension(:)  ::  wl, wl_rb
    real (PD) :: x1,x2,x3,y1,y2,y3

    
    n1 = size(wl,1);  n2 = size(wl_rb,1)
    fl_rb = 0.0
    
    do n=1,n2
       ! write(*,*) n,n2
       if (wl_rb(n).lt.wl(1) .or. wl_rb(n).gt.wl(n1)) cycle 
       k=0
       do kk=1,n1-1
          if (present(verbose)) write(verbose,*) n,n2,kk,n1

          if (WL_RB(n).ge.WL(kk).and.WL_RB(n).lt.wl(kk+1)) then
             
             k1=kk ; k2=kk+1
             i1=k1 ; i2=k2
             kt=0
             
             do while (k1.gt.i1-5 .and. k2.le.i2+5) 
                kt=kt+1
                if ( fl(k1).ge.0.0 .and. fl(k2).ge.0.0) exit
                if ( fl(k1).le.0.0 .or. .not.(fl(k1).gt.-10e18 .and. fl(k1).lt.10e18)) k1=k1-1
                if ( fl(k2).le.0.0 .or. .not.(fl(k1).gt.-10e18 .and. fl(k1).lt.10e18)) k2=k2+1
                if (k1.lt.1 .or. k2.gt.n1) exit
                if (kt.gt.100) exit ! something wrong happened while iterating, we force the exit from the cycle
             enddo
             if  (k1.le.i1-5 .or. k2.gt.i2+5 .or. kt.gt.100) cycle
             
             x1=wl(k1) ; x2=wl(k2)
             y1=dble(fl(k1)) ; y2=dble(fl(k2))
             fl_rb(n)=real((y2-y1)/sngl(x2-x1)*sngl(wl_rb(n)-x1)+y1)
             if ( fl_rb(n).gt.1.0D8 .or. fl_rb(n).lt.-1.0D8) fl_rb(n) = 0.0
             exit
          endif
       enddo

    enddo
    
    return
  end subroutine rebin_keepnorm

  subroutine rebin_keepflux(wl,dw,fl,wl_rb,dw_rb,fl_rb,no_norm)
    !normalization is on by default
    real (PD), dimension(:), intent(in ) :: wl, dw, fl, wl_rb, dw_rb
    real (PD), dimension(:), intent(out) :: fl_rb
    real (PD) :: wlb, wle, pix_e, pix_s, fl_nm
    integer ::  ns, ns_prv,  nl, n1, n2,  i, j, ir, ib, norm
    integer, intent(in), optional :: no_norm
    
    norm = 1
    if (present(no_norm)) norm = 0
    
    fl_rb = 0.000d0
    n1 = size(wl,1);  n2 = size(wl_rb,1)
    ns_prv=1
    do i=1,n2
       !Starting and ending point of the bin
       wlb =  wl_rb(i) - dw_rb(i)/2.000d0
       wle =  wl_rb(i) + dw_rb(i)/2.000d0

       !Normalized flux value within the bin
       fl_nm = 0.000d0

       !! b->blue and r->red side of the original spectrum which include the bin
       !! ib and ie are initialized with values close to the ones of the last iteration to save time
       ib = ns_prv ;  ir = ns_prv
       
       do ns=ns_prv,n1-1
          if (wl(ib)+dw(ib)*0.500d0.lt.wlb) ib = ib+1
          if (wl(ir)+dw(ir)*0.500d0.lt.wle) ir = ir+1
          
!!$          write(*,*)  wlb, wle, ib, ir, ns, n1, n2, wl(ib)+dw(ib)*0.500d0, wl(ir)+dw(ir)*0.500d0
          if (ir.lt.ns-3) exit
       end do

       !! Fail-safe checks
       if (ib.gt.ns_prv) ns_prv = ib - 3
       if (ib.le.0 .or. ir.gt.n1) cycle
       if (ib.gt.ir) cycle
       if (ns_prv.le.0) ns_prv=1
       
       if (ib.eq.ir)  then
          pix_s    = (wle-wlb ) / dw(ib)     ! fraction 
          fl_rb(i)   = fl_rb(i)   + pix_s*fl(ib)
          fl_nm      = fl_nm      + pix_s
       else if (ib+1.eq.ir) then
          pix_s     = (wl(ib)+dw(ib)*0.5-wlb) / dw(ib)
          pix_e     = (wle-(wl(ir)-dw(ir)*0.5))/ dw(ir)
!!$          write(*,*) dw(ib), dw(ir), pix_s, pix_e, fl_rb(i), fl_nm, (pix_s*fl(ib) + pix_e*fl(ir))
          fl_rb(i)   = fl_rb(i)   + (pix_s*fl(ib) + pix_e*fl(ir))
          fl_nm      = fl_nm   + (pix_s + pix_e)
!!$          stop 
       else 
          pix_s     = (wl(ib)+dw(ib)*0.5-wlb)/dw(ib)
          pix_e     = (wle-(wl(ir)-dw(ir)*0.5))/ dw(ir)
          fl_rb(i)  = fl_rb(i)   + (pix_s*fl(ib)+pix_e*fl(ir))
          fl_nm      = fl_nm   + (pix_s + pix_e)
          do j=ib+1,ir-1
             fl_rb(i)  = fl_rb(i)   + fl(j)
             fl_nm      = fl_nm   +  1.00d0
          end do
       end if
       if (norm.eq.1 .and. fl_nm.gt.0.000d0) fl_rb(i) = fl_rb(i) / fl_nm
       !pay attention to the ir-i index, 
    enddo
    return
  end subroutine rebin_keepflux
  
end module spectrum_rebin
