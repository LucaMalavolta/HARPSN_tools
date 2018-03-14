program fits_coaddspectra_e2ds
  use Common
  use Fits
  use Instrument_HARPS
  use Filepaths
  use Spectrum_Rebin
  use Maths
  implicit none


  integer :: i, j, n,ij, ii, j_is, j_ie
  integer :: lun_input, lun_fits, lun_out1, lun_out2
  integer :: opt_output

  character (len = nch_file) :: file_input, file_date, file_rad, file_s1d_name, file_s2d_name, file_blaze, &
       file_s1d_output, star_name, file_output
  character (len = 2) :: mask_sel
  character (len = 1) :: fiber_sel

  integer :: read_iostat
  real (kind=8) :: rv_fin, berv

  integer :: x_axis, y_axis
  integer, allocatable, dimension (:) :: is, ie
  real (kind=8), allocatable, dimension (:,:)  :: r_pix, b_pix, s_pix, dw_pix

  real :: crval, crpix, cdelt
  real :: nr_joint, ccd_gain

  real (kind=8) :: rb_bin

  real (kind=8), allocatable, dimension (:)  :: r_rb_pix, j_rb_pix, dw_rb_pix,d1_rb_pix
  real (kind=8), allocatable, dimension(:,:) :: s_rb_pix, b_rb_pix

  real (kind=8), allocatable, dimension(:) ::  sum_rb_pix, snr_rb_pix
  integer :: n_rebin


  integer :: sts, sts_out

  character (len=nch_key) :: ctype,bunit


  real (kind=8) :: wl_min, wl_max , mjd_mean


  integer :: n_com, date_opt
  character (len = nch_file) :: arch_red, opt_in

  call get_lun(lun_input)
  call get_lun(lun_out1)
  call get_lun(lun_out2)
  call get_lun(lun_fits )

  n_com = iargc()

  if (n_com.lt.4) then
     write(*,*) 'Input 1: input file list'
     write(*,*) 'Input 2: ouput rad name'
     write(*,*) 'Input 3: rebinning step'
     write(*,*) 'Input 4: 1: extensive output files'
     write(*,*) 'Input 5 (optional): Archive directory '
     write(*,*) 'Input 6 (optional): 0 = do not use date as prefix '
     write(*,*) 'Input 6: to be used if all the files are in the same directory (as Yabi reprocessed files are) '
     !!write(*,*) 'Input 8 (optional): ESO flag (1 for HARPS) '
     stop
  end if

  call getarg(1,file_input)
  call getarg(2,file_output)
  call getarg(3,opt_in)
  read(opt_in,*) rb_bin
  call getarg(4,opt_in)
  read(opt_in,*) opt_output

  arch_red = archive_red
  date_opt = 1
  if (n_com.eq.6) then
     call getarg(5,arch_red)
     call getarg(6,opt_in); read (opt_in,*) date_opt
  end if



  open(lun_input,file=trim(file_input), status='old')

  n_rebin = nint(wl_range/rb_bin)+1

  wl_min =    0._PR
  wl_max = 9900._PR


  !! Allocating array for e2ds
  allocate(r_rb_pix(n_rebin), j_rb_pix(n_rebin), snr_rb_pix(n_rebin), dw_rb_pix(n_rebin), sum_rb_pix(n_rebin))
  allocate(d1_rb_pix(n_rebin))

  do i=1,n_rebin
     r_rb_pix(i) = wl_start + (i-1)*rb_bin
  end do
  dw_rb_pix = rb_bin

  read_iostat = 0 ; n = 0

  mjd_mean = 0._PR

  sum_rb_pix = 0._PR
  snr_rb_pix = 0._PR

  if (opt_output.gt.0) then
     file_s1d_output =  trim(file_output) //'_s1d_all.dat'
     open(lun_out1,file=trim(file_s1d_output),status='new')
     file_s1d_output =  trim(file_output) //'_rebin_all.dat'
     open(lun_out2,file=trim(file_s1d_output),status='new')
  end if

  do while (read_iostat.eq.0)
     read(lun_input,*,iostat=read_iostat) file_date, file_rad, mask_sel, fiber_sel, star_name
     if (read_iostat.ne.0) exit
     n = n + 1
     write(*,*)

     if (date_opt.eq.1) then
        file_s2d_name = trim(arch_red)//trim(file_date)//'/'//trim(file_rad)//'_e2ds_'//fiber_sel//'.fits'
     else
        file_s2d_name = trim(arch_red)//'/'//trim(file_rad)//'_e2ds_'//fiber_sel//'.fits'
     end if

     call fits_check(file_s2d_name,sts_out)
     write(*,*) n, trim(file_date),'   ', trim(file_rad),'   ', trim(file_s2d_name), '   FOUND= ',sts_out,'(1=true)'

     if (sts_out.eq.0) cycle
     call fits_size_2d(file_s2d_name,x_axis,y_axis)
     allocate(s_pix(x_axis,y_axis),r_pix(x_axis,y_axis), b_pix(x_axis,y_axis), dw_pix(x_axis,y_axis))
     allocate(is(y_axis), ie(y_axis))
     allocate(s_rb_pix(n_rebin,y_axis))
     allocate(b_rb_pix(n_rebin,y_axis))

     s_pix = 0._PR ; r_pix = 0._PR ; b_pix = 0._PR ; dw_pix = 0._PR
     is = 0; ie = 0
     s_rb_pix = 0._PR ; b_rb_pix = 0._PR

     call fits_open_s2d (lun_fits,file_s2d_name,s_pix,verbose=6)

     call key2wave(lun_fits,r_pix)
     call singleheader_CCD_GAIN(lun_fits,ccd_gain,sts)
     call singleheader_BLAZE(lun_fits,file_blaze,sts)
     call fits_close(lun_fits)
     !! WARNING: e2ds are already corrected for gain!!
     !! s_pix = s_pix*ccd_gain

     !! CCF file is udes to read additional informations of spectra
     if (date_opt.eq.1) then
        file_s1d_name = trim(arch_red)//trim(file_date)//'/'//trim(file_rad)//'_ccf_'//mask_sel//'_'//fiber_sel//'.fits'
     else
        file_s1d_name = trim(arch_red)//'/'//trim(file_rad)//'_ccf_'//mask_sel//'_'//fiber_sel//'.fits'
     end if

     call fits_check(file_s1d_name,sts_out)
     write(*,*) n, trim(file_date),'   ', trim(file_rad),'   ', trim(file_s1d_name), '   FOUND= ',sts_out,'(1=true)'

     call fits_actv(lun_fits,file_s1d_name)
     call singleheader_RVC(lun_fits,rv_fin,sts)
     call singleheader_BERV(lun_fits,berv,sts)
     call fits_header_getinfo(lun_fits)
     call fits_close(lun_fits)

     if (date_opt.eq.1) then
        file_s2d_name = trim(arch_red)//trim(file_date)//'/'//trim(file_blaze)
     else
        file_s2d_name = trim(arch_red)//'/'//trim(file_blaze)
     end if

     call fits_check(file_s2d_name,sts_out)
     write(*,*) n, trim(file_date),'   ', trim(file_blaze),'   ', trim(file_s2d_name), '   FOUND= ',sts_out,'(1=true)'
     call fits_open_s2d (lun_fits,file_s2d_name,b_pix,verbose=6)
     call fits_close(lun_fits)



     !!!!!!!!!r_pix = (1._PR)*((berv-rv_drift-rv_fin)/cc+1.00000D0)*r_pix
     write(*,*) ' BERV = ', berv , 'RV = ', rv_fin
     r_pix = (1._PR)*((berv-rv_fin)/cc+1.00000D0)*r_pix

     if ( r_pix(1,1) .gt. wl_min ) wl_min =  r_pix(1,1)
     if ( r_pix(x_axis,y_axis) .lt. wl_max ) wl_max =  r_pix(x_axis,y_axis)

     !! First derivative, needed for a correct rebinning of the spectra
     dw_pix(2:x_axis-1,:) = (r_pix(3:x_axis,:) - r_pix(1:x_axis-2,:))/2
     dw_pix(1,:) =  dw_pix(2,:)
     dw_pix(x_axis,:) = dw_pix(x_axis-1,:)


     !! Rebinning of the individual orders. We keep track of the starting end ending point of each order in the rebinned solution
     b_rb_pix = 0._PR
     do ii=1,y_axis
        is(ii) = minloc(abs(r_pix(     1,ii)-r_rb_pix(:)),dim=1) + 1
        ie(ii) = minloc(abs(r_pix(x_axis,ii)-r_rb_pix(:)),dim=1) - 1
        i = is(ii)
        j = ie(ii)
        call rebin_keepflux(r_pix(:,ii),dw_pix(:,ii),s_pix(:,ii),r_rb_pix(i:j),dw_rb_pix(i:j),b_rb_pix(i:j,ii))
        s_pix(:,ii) = s_pix(:,ii)/b_pix(:,ii)!/(dw_pix(:,ii)/0.01_PR) -> substitued with no_norm=1 option

        call rebin_keepflux(r_pix(:,ii),dw_pix(:,ii),s_pix(:,ii),r_rb_pix(i:j),dw_rb_pix(i:j),s_rb_pix(i:j,ii),no_norm=1)

     end do


     !! Now: joining the orders

     !! fix for the intial part of the spectra
     j_rb_pix(is(1):ie(1)) =  s_rb_pix(is(1):ie(1),1)
     snr_rb_pix = snr_rb_pix + b_rb_pix(:,1)

     do ii=2,y_axis
        j_ie = ie(ii-1)
        j_is = is(ii)
!!$        write(*,*) ii,  j_ie, j_is, r_rb_pix(j_ie), r_rb_pix(j_is)
        if (j_ie .eq. 0 .or. j_is.eq.0) cycle

        if (j_ie .gt. j_is) then
           nr_joint = real(j_ie - j_is -1)
           do ij = j_is, j_ie
              j_rb_pix(ij) = s_rb_pix(ij,ii-1) * real( j_ie - ij,PR)/nr_joint + s_rb_pix(ij,ii) * real(ij-j_is,PR)/nr_joint
           end do
           j_rb_pix(j_ie:ie(ii)) =  s_rb_pix(j_ie:ie(ii),ii)
        else
           j_rb_pix(j_is:ie(ii)) =  s_rb_pix(j_is:ie(ii),ii)
        end if
        snr_rb_pix = snr_rb_pix + b_rb_pix(:,ii)


     end do

     sum_rb_pix = sum_rb_pix + j_rb_pix


     if (opt_output.gt.0) then
        do j=1,y_axis
           do i = 1,x_axis
!!$              if (r_pix(i,j).lt.3700 .or. r_pix(i,j).gt.3910) cycle
              if (r_pix(i,j).lt.5700 .or. r_pix(i,j).gt.6400) cycle
              write(lun_out1,'(I3,2I7,F10.4,3E14.6)') n, i, j, r_pix(i,j), b_pix(i,j), s_pix(i,j), dw_pix(i,j)
           end do
        end do
        do i = 1,n_rebin
!!$           if (r_rb_pix(i).lt.3700 .or. r_rb_pix(i).gt.3910) cycle
           if (r_rb_pix(i).lt.5700 .or. r_rb_pix(i).gt.6400) cycle
           write(lun_out2,'(I3,I7,F10.4,3E14.6)') n, i, r_rb_pix(i), j_rb_pix(i), snr_rb_pix(i), dw_rb_pix(i)
        end do

     end if

     mjd_mean = mjd_mean + head_mjdobs

     deallocate(s_pix,r_pix, b_pix, dw_pix)
     deallocate(is, ie, s_rb_pix, b_rb_pix)

  end do
  close(lun_input)

  if (opt_output.gt.0) then
     close(lun_out1)
     close(lun_out2)
  end if

  !!! first derivative
  d1_rb_pix = 0.00000000d0
  d1_rb_pix(2:n_rebin-1) = (sum_rb_pix(3:n_rebin)-sum_rb_pix(1:n_rebin-2))/(r_rb_pix(3:n_rebin)-r_rb_pix(1:n_rebin-2))
  d1_rb_pix(1) = d1_rb_pix(2)
  d1_rb_pix(n_rebin) = d1_rb_pix(n_rebin-1)

  file_s1d_output = trim(file_output) //'_coadd.dat'
  open(lun_out1,file=file_s1d_output,status='new')
  do i = 1,n_rebin
     write(lun_out1,'(I9,F11.4,5E15.7)') i, r_rb_pix(i), sum_rb_pix(i), dw_rb_pix(i), snr_rb_pix(i), d1_rb_pix(i), j_rb_pix(i)
  end do
  close(lun_out1)

  mjd_mean = mjd_mean / real(n)

  ctype = 'Angstrom'
  bunit = 'Relative Flux'
  crpix = 1.0
  crval = real(wl_start)
  cdelt = real(rb_bin)
  head_radvel = 0._PR
  head_mjdobs = mjd_mean

  file_s1d_output = trim(file_output) //'_coadd.fits'
  call fits_write_1d(lun_out1,file_s1d_output,sum_rb_pix)
  call scale_to_wcs(lun_out1,1,ctype=ctype,crval=crval,crpix=crpix,cdelt=cdelt,bunit=bunit)
  call fits_header_putinfo(lun_out1)
  call fits_close(lun_out1)


  file_s1d_output = trim(file_output) //'_snr.fits'
  call fits_write_1d(lun_out1,file_s1d_output,snr_rb_pix)
  call scale_to_wcs(lun_out1,1,ctype=ctype,crval=crval,crpix=crpix,cdelt=cdelt,bunit=bunit)
  call fits_header_putinfo(lun_out1)
  call fits_close(lun_out1)



  !! dellocating array for e2ds
  deallocate(r_rb_pix, j_rb_pix, snr_rb_pix, dw_rb_pix, sum_rb_pix)
  deallocate(d1_rb_pix)

  write(*,*) 'PROGRAM HAS ENDED'


  call free_lun(lun_out1)
  call free_lun(lun_out2)
  call free_lun(lun_input)
  call free_lun(lun_fits )
end program fits_coaddspectra_e2ds
