module Spectrum_linefit
  use common
  use voigt_fit
  use levenberg_marquardt_mod5V
  implicit none

  private
  public :: spectral_lines_fit, get_spectral_lines, get_spectral_line_EW, get_single_spectral_lines, spectral_lines_singlefit

contains

  real (PR) function d1mach(I)
    integer,intent(IN):: I
    !
    !  DOUBLE-PRECISION MACHINE CONSTANTS
    !  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
    !  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
    !  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
    !  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
    !  D1MACH( 5) = LOG10(B)
    !
    real(PR),dimension(5),save::X
    logical,save:: ran = .false.

    if(.not. ran)then
       X(1) = tiny(X(1))
       X(2) = huge(X(2))
       X(3) = epsilon(X(3))/radix(X(3))
       X(4) = epsilon(X(4))
       X(5) = log10(real(radix(X(5)),PR))
       ran = .true.
    end if

    if (I < 1 .or. I > 5) then
       write(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
       stop
    end if
    D1MACH = X(I)
    return
  end function d1mach

  real (PR) function get_spectral_lines(x,pams_val,back,cond_lines)
    real (PR), dimension(:,:), intent(in) :: pams_val
    real (PR), intent(in) :: x, back
    logical, dimension(:), intent(in), optional :: cond_lines
    logical, dimension(:), allocatable :: int_cond
    REAL (PR) :: x1,x2,x3,x4, line
    
    integer :: k, nline
    
    !! PAMS legend
    !! PAMS_xxx(1:,) = Center of the line
    !! PAMS_xxx(2:,) = FWHM of the Gaussian profile
    !! PAMS_xxx(3:,) = FWHM of the Lorentian profile (when the Hjerting function is used)
    !! PAMS_xxx(4:,) = AMPL of the line 
    nline = size(pams_val,dim=2)

    allocate(int_cond(nline))
    int_cond = .true.
    if (present(cond_lines)) int_cond = cond_lines
    
    get_spectral_lines = back
    do k=1,nline
       if (int_cond(k) .and. pams_val(3,k).le.0._PR) &
            get_spectral_lines=get_spectral_lines + &
            pams_val(4,k)*exp((-(x-pams_val(1,k))**2)/(2._PR*(pams_val(2,k)/sgm2fwhm)**2)) 
       if (int_cond(k) .and. pams_val(3,k).gt.0._PR) &
            get_spectral_lines=get_spectral_lines + & 
            pams_val(4,k)*get_ext_pseudo_voigt(x-pams_val(1,k),pams_val(2,k),pams_val(3,k))
    end do
    deallocate(int_cond)
    
    return

  end function get_spectral_lines


  real (PR) function get_single_spectral_lines(x,pams_val_1l,back)
    real (PR), dimension(:), intent(in) :: pams_val_1l
    real (PR), intent(in) :: x, back
    REAL (PR) :: x1,x2,x3,x4, line
    
    integer :: k, nline
    
    !! PAMS legend
    !! PAMS_xxx(1:,) = Center of the line
    !! PAMS_xxx(2:,) = FWHM of the Gaussian profile
    !! PAMS_xxx(3:,) = FWHM of the Lorentian profile (when the Hjerting function is used)
    !! PAMS_xxx(4:,) = AMPL of the line 

    
    if ( pams_val_1l(3).le.0._PR) &
         get_single_spectral_lines= back + &
         pams_val_1l(4)*exp((-(x-pams_val_1l(1))**2)/(2._PR*(pams_val_1l(2)/sgm2fwhm)**2)) 
    if (pams_val_1l(3).gt.0._PR) &
         get_single_spectral_lines= back + &
         pams_val_1l(4)*get_ext_pseudo_voigt(x-pams_val_1l(1),pams_val_1l(2),pams_val_1l(3))
    
    return

  end function get_single_spectral_lines


  real (PR) function get_spectral_line_EW(pams_val_1l)
    real (PR), dimension(4) :: pams_val_1l
    real (PR) :: x_step,xx, y_add
    integer, parameter :: half_size=10000
    integer :: ii

    x_step = 20*pams_val_1l(2)/real(half_size)
    get_spectral_line_EW = 0._PR
    if (pams_val_1l(3).le.0._PR) then
       get_spectral_line_EW  =  sqrt(2._PR*pid)*abs(pams_val_1l(4))*pams_val_1l(2)/sgm2fwhm*1000._PR
    else
       do ii=1,half_size
          xx = ii*x_step
          y_add = abs(pams_val_1l(4))*(get_ext_pseudo_voigt(xx,pams_val_1l(2),pams_val_1l(3))*1000._PR)
          if (xx.gt.10*pams_val_1l(2) .or.  y_add .lt. 0.001) exit
          get_spectral_line_EW=get_spectral_line_EW + x_step*y_add
       end do
       xx = 0._PR
       get_spectral_line_EW = 2._PR * get_spectral_line_EW + &
            x_step *abs(pams_val_1l(4))*(get_ext_pseudo_voigt(xx,pams_val_1l(2),pams_val_1l(3))*1000._PR)
    end if
    return
    
  end function get_spectral_line_EW


  subroutine spectral_lines_fit(xdata,ydata,wdata,nline,pams_val,back_val,pams_fix,back_fix,cond_lines)
    integer :: m, k, nline, i, j
    integer, dimension(:,:), intent(in), optional :: pams_fix
    integer, optional, intent(in) :: back_fix
    logical, dimension(:)  , intent(in), optional :: cond_lines
    
    real (PR), dimension(:), intent(in) ::xdata,ydata,wdata

    real (PR), dimension(:,:), intent(inout) :: pams_val
    real (PR), intent(inout) :: back_val 

    real (PR), dimension(size(xdata,dim=1)) :: xdata_rc
    real (PR), dimension(:), allocatable :: PAMSfix_cond
    real (PR), dimension(2) :: BACKfix_cond

    real (PR) :: cntr_tmp

    integer :: N,INFO
    integer, dimension(4*nline+1) ::  IWA
    real (PR) :: TOL
    REAL (PR),dimension(4*nline+1) ::x,diag,sig
    REAL (PR),dimension(size(xdata)) ::fvec


    if (size(ydata).ne.size(xdata)) stop
    if (size(wdata).ne.size(xdata)) stop

    !! PAMS legend
    !! PAMS_xxx(1,:) = Center of the line
    !! PAMS_xxx(2,:) = FWHM of the Gaussian profile
    !! PAMS_xxx(3,:) = FWHM of the Lorentian profile (when the Hjerting function is used)
    !! PAMS_xxx(4,:) = AMPL of the line 
    
    !! Legend for cond_PR
    !! cond_PR .eq. 1 .or. 2 -> line is included in the fit, gaussian profile
    !! cond_PR .eq. 3 .or. 4 -> line is included in the fit, Voig profile
    !! cond_PR .eq. 2 .or. cond_PR .eq. 4 .or. cond_PR .eq. -2 -> background is kept at a fixed value
    !!                                      (regardless if the line is considered in the gaussian fit)

    cntr_tmp = xdata(size(xdata,dim=1)/2)

    allocate(PAMSfix_cond(9*nline))
    PAMSfix_cond = 0._PR
    
    do i=1,nline       
       k = (i-1)*4 + 1 ; PAMSfix_cond(k) = pams_val(1,i)-cntr_tmp
       k = (i-1)*4 + 2 ; PAMSfix_cond(k) = pams_val(2,i)
       k = (i-1)*4 + 3 ; PAMSfix_cond(k) = pams_val(3,i)
       k = (i-1)*4 + 4 ; PAMSfix_cond(k) = pams_val(4,i)/(1._PR - PAMS_val(4,i))
    end do
    if (present(pams_fix)) then
       do i=1,nline       
          k = (i-1)*4 + 1 + nline*4 ; PAMSfix_cond(k) = real(pams_fix(1,i))
          k = (i-1)*4 + 2 + nline*4 ; PAMSfix_cond(k) = real(pams_fix(2,i))
          k = (i-1)*4 + 3 + nline*4 ; PAMSfix_cond(k) = real(pams_fix(3,i))
          k = (i-1)*4 + 4 + nline*4 ; PAMSfix_cond(k) = real(pams_fix(4,i))
       end do
    end if
       
    if (present(cond_lines)) then
       k = nline*8
       where (cond_lines) PAMSfix_cond(k+1:k+nline) = 1._PR
    end if
    
    BACKfix_cond = 0._PR
    BACKfix_cond(1) = back_val
    if (present(back_fix)) then
       BACKfix_cond(2) = real(back_fix)
    end if
    
    
    do i=1,nline
       k = (i-1)*4 + 1 ;  X(k)=  PAMS_val(1,i)-cntr_tmp
       k = (i-1)*4 + 2 ;  X(k)=  PAMS_val(2,i)
       k = (i-1)*4 + 3 ;  X(k)=  PAMS_val(3,i)
       k = (i-1)*4 + 4 ;  X(k)=  PAMS_val(4,i)/(1._PR - PAMS_val(4,i))
    end do
    x(nline*4+1) = back_val
    
    xdata_rc = xdata - cntr_tmp
    m = size(xdata)

    fvec = 0._PR; diag = 0._PR; sig = 0._PR
    TOL = D1MACH(4) 
    INFO=0 
    N = nline*4 + 1
    sig = 0._PR
    call lmdif2_mod5V(spectral_lines_f,m,n,x,xdata_rc,ydata,wdata,PAMSfix_cond,BACKfix_cond,fvec,tol,diag,sig,info,iwa)

    do i=1,nline
       k = (i-1)*4 + 1 ;   PAMS_val(1,i) = X(k)+cntr_tmp 
       k = (i-1)*4 + 2 ;   PAMS_val(2,i) = abs(X(k))           
       k = (i-1)*4 + 3 ;   PAMS_val(3,i) = abs(X(k))     
       k = (i-1)*4 + 4 ;   PAMS_val(4,i) = X(k)/(1._PR+X(k))  
    end do
    back_val = x(nline*4+1) 
    
    do i=1,nline       
       do j = 1, 4
          k = (i-1)*4 + j + nline*4 ;
          if (PAMSfix_cond(k).lt.0._PR) pams_val(j,i) = 0._PR
          if (PAMSfix_cond(k).gt.0._PR) then
             if (j.eq.1) pams_val(j,i) = PAMSfix_cond(k-nline*4)+cntr_tmp
             if (j.eq.2) pams_val(j,i) = PAMSfix_cond(k-nline*4)
             if (j.eq.3) pams_val(j,i) = PAMSfix_cond(k-nline*4)
             if (j.eq.4) pams_val(j,i) = PAMSfix_cond(k-nline*4)/(1._PR + PAMSfix_cond(k-nline*4))
          end if
       end do
    end do
    
    if (BACKfix_cond(2).gt.0) back_val = BACKfix_cond(1)
    deallocate(PAMSfix_cond)
    
    return
  end subroutine spectral_lines_fit
  

  SUBROUTINE spectral_lines_f(m,n,x,xd,yd,wd,PAMSfix_cond,BACKfix_cond,fvec,iflag)
    use common
    IMPLICIT NONE
    !INTEGER,PARAMETER::dp = SELECTED_REAL_KIND(12,60)
    INTEGER,INTENT(IN)::m,n
    REAL (PR),INTENT(IN)::x(:)
    REAL (PR),INTENT(OUT)::fvec(:)
    INTEGER,INTENT(inout)::iflag
    REAL (PR), intent(in), dimension(:) :: xd,yd,wd,BACKfix_cond
    REAL (PR), intent(in), dimension(:) :: PAMSfix_cond
    INTEGER :: k,p,q,i, nline
    REAL (PR) :: x1,x2,x3,x4, line
    
    nline = (size(X,dim=1)-1)/4
    do k=1,M
       line = X(n)
       if (BACKfix_cond(2) .gt. 0._PR)  line = BACKfix_cond(1)
       
       do i=1,nline
          q = i + 8*nline
          if (PAMSfix_cond(q) .lt. 1._PR) cycle

          q = (i-1)*4 + 3 
          if (PAMSfix_cond(q) .gt. 0._PR) then

             p = (i-1)*4 + 1 ;   x1 = X(p)  
             if (PAMSfix_cond(p+4*nline) .gt.0._PR)  x1 = PAMSfix_cond(p)

             p = (i-1)*4 + 2 ;   x2 = X(p)
             if (PAMSfix_cond(p+4*nline) .gt.0._PR)  x2 = PAMSfix_cond(p)

             p = (i-1)*4 + 3 ;   x3 = abs(X(p))
             if (PAMSfix_cond(p+4*nline) .gt.0._PR)  x3 = PAMSfix_cond(p)

             p = (i-1)*4 + 4 ;   x4 = X(p)/(1._PR + X(p))
             if (PAMSfix_cond(p+4*nline) .gt.0._PR)  x4 = PAMSfix_cond(p)
             
             line=line + x4*(get_ext_pseudo_voigt(xd(k)-X1,x2,x3))
          else
             p = (i-1)*4 + 1 ;   x1 = X(p)
             if (PAMSfix_cond(p+4*nline) .gt.0._PR)  x1 = PAMSfix_cond(p)

             p = (i-1)*4 + 2 ;   x2 = X(p)/sgm2fwhm
             if (PAMSfix_cond(p+4*nline) .gt.0._PR)  x2 = PAMSfix_cond(p)

             p = (i-1)*4 + 4 ;   x4 = X(p)/(1._PR + X(p))
             if (PAMSfix_cond(p+4*nline) .gt.0._PR)  x4 = PAMSfix_cond(p)
             
             line=line + X4*exp((-(xd(k)-X1)**2)/(2._PR*X2**2)) 
             
          end if
             
          
          
       end do
       FVEC(k)= (yd(k) - line)*wd(k)!*(abs((X(3)-yd(k))/(X(3)-X(4))))**2
    enddo
    
    return
  end subroutine spectral_lines_f


  subroutine spectral_lines_singlefit(xdata,ydata,wdata,pams_val,back_val,pams_fix,back_fix)
    integer :: m, k, i, j
    integer, dimension(:), intent(in), optional :: pams_fix
    integer, optional, intent(in) :: back_fix
    
    real (PR), dimension(:), intent(in) ::xdata,ydata,wdata

    real (PR), dimension(:), intent(inout) :: pams_val
    real (PR), intent(inout) :: back_val 

    real (PR), dimension(size(xdata,dim=1)) :: xdata_rc
    real (PR), dimension(8) :: PAMSfix_cond
    real (PR), dimension(2) :: BACKfix_cond

    real (PR) :: cntr_tmp

    integer :: N,INFO
    integer, dimension(5) ::  IWA
    real (PR) :: TOL
    REAL (PR),dimension(5) ::x,diag,sig
    REAL (PR),dimension(size(xdata)) ::fvec


    if (size(ydata).ne.size(xdata)) stop
    if (size(wdata).ne.size(xdata)) stop

    !! PAMS legend
    !! PAMS_xxx(1) = Center of the line
    !! PAMS_xxx(2) = FWHM of the Gaussian profile
    !! PAMS_xxx(3) = FWHM of the Lorentian profile (when the Hjerting function is used)
    !! PAMS_xxx(4) = AMPL of the line 
    
    !! Legend for cond_PR
    !! cond_PR .eq. 1 .or. 2 -> line is included in the fit, gaussian profile
    !! cond_PR .eq. 3 .or. 4 -> line is included in the fit, Voig profile
    !! cond_PR .eq. 2 .or. cond_PR .eq. 4 .or. cond_PR .eq. -2 -> background is kept at a fixed value
    !!                                      (regardless if the line is considered in the gaussian fit)

    cntr_tmp = xdata(size(xdata,dim=1)/2)

    PAMSfix_cond = 0._PR
    
    PAMSfix_cond(1) = pams_val(1)-cntr_tmp
    PAMSfix_cond(2) = pams_val(2)
    PAMSfix_cond(3) = pams_val(3)
    PAMSfix_cond(4) = pams_val(4)/(1._PR - PAMS_val(4))

    if (present(pams_fix)) then
       PAMSfix_cond(5) = real(pams_fix(1))
       PAMSfix_cond(6) = real(pams_fix(2))
       PAMSfix_cond(7) = real(pams_fix(3))
       PAMSfix_cond(8) = real(pams_fix(4))
    end if
      
    
    BACKfix_cond = 0._PR
    BACKfix_cond(1) = back_val
    if (present(back_fix)) then
       BACKfix_cond(2) = real(back_fix)
    end if
    
    
    
    X(1)=  PAMS_val(1)-cntr_tmp
    X(2)=  PAMS_val(2)
    X(3)=  PAMS_val(3)
    X(4)=  PAMS_val(4)/(1._PR - PAMS_val(4))
    X(5) = back_val
    
    xdata_rc = xdata - cntr_tmp
    m = size(xdata)
    
    fvec = 0._PR; diag = 0._PR; sig = 0._PR
    TOL = D1MACH(4) 
    INFO=0 
    N =5
    sig = 0._PR
    call lmdif2_mod5V(spectral_lines_singlefit_f,m,n,x,xdata_rc,ydata,wdata,PAMSfix_cond,BACKfix_cond,fvec,tol,diag,sig,info,iwa)
    
    
    PAMS_val(1) = X(1)+cntr_tmp 
    PAMS_val(2) = abs(X(2))           
    PAMS_val(3) = abs(X(3))     
    PAMS_val(4) = X(4)/(1._PR+X(4))  
    back_val = x(5) 
    
    do j = 1, 4
       k =  j + 4 ;
       if (PAMSfix_cond(k).lt.0._PR) pams_val(j) = 0._PR
       if (PAMSfix_cond(k).gt.0._PR) then
          if (j.eq.1) pams_val(j) = PAMSfix_cond(j)+cntr_tmp
          if (j.eq.2) pams_val(j) = PAMSfix_cond(j)
          if (j.eq.3) pams_val(j) = PAMSfix_cond(j)
          if (j.eq.4) pams_val(j) = PAMSfix_cond(j)/(1._PR + PAMSfix_cond(j))
       end if
    end do
    
    if (BACKfix_cond(2).gt.0) back_val = BACKfix_cond(1)
    
    return
  end subroutine spectral_lines_singlefit
  

  SUBROUTINE spectral_lines_singlefit_f(m,n,x,xd,yd,wd,PAMSfix_cond,BACKfix_cond,fvec,iflag)
    use common
    IMPLICIT NONE
    INTEGER,INTENT(IN)::m,n
    REAL (PR),INTENT(IN)::x(:)
    REAL (PR),INTENT(OUT)::fvec(:)
    INTEGER,INTENT(inout)::iflag
    REAL (PR), intent(in), dimension(:) :: xd,yd,wd,BACKfix_cond
    REAL (PR), intent(in), dimension(:) :: PAMSfix_cond
    INTEGER :: k,i
    REAL (PR) :: x1,x2,x3,x4, line

    do k=1,M
       line = X(5)
       if (BACKfix_cond(2) .gt. 0._PR)  line = BACKfix_cond(1)

       if (PAMSfix_cond(3) .gt. 0._PR) then

          x1 = X(1)  
          if (PAMSfix_cond(5) .gt.0._PR)  x1 = PAMSfix_cond(1)

          x2 = abs(X(2))
          if (PAMSfix_cond(6) .gt.0._PR)  x2 = PAMSfix_cond(2)

          x3 = abs(X(3))
          if (PAMSfix_cond(7) .gt.0._PR)  x3 = PAMSfix_cond(3)

          x4 = X(4)/(1._PR + X(4))
          if (PAMSfix_cond(8) .gt.0._PR)  x4 = PAMSfix_cond(4)

          line=line + x4*(get_ext_pseudo_voigt(xd(k)-X1,x2,x3))
       else
          x1 = X(1)
          if (PAMSfix_cond(5) .gt.0._PR)  x1 = PAMSfix_cond(1)

          x2 = X(2)/sgm2fwhm
          if (PAMSfix_cond(6) .gt.0._PR)  x2 = PAMSfix_cond(2)

          x4 = X(4)/(1._PR + X(4))
          if (PAMSfix_cond(8) .gt.0._PR)  x4 = PAMSfix_cond(4)

          line=line + X4*exp((-(xd(k)-X1)**2)/(2._PR*X2**2)) 

       end if


       FVEC(k)= (yd(k) - line)*wd(k)!*(abs((X(3)-yd(k))/(X(3)-X(4))))**2
    enddo

    return
  end subroutine spectral_lines_singlefit_f

  

end module Spectrum_Linefit

