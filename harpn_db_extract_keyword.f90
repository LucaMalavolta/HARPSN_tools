program db_extract_keyword
  use Common
  use Fits
  use Filepaths
  use Instrument_HARPN
  implicit none


  integer :: i, j, n
  integer :: lun_input, lun_fits, lun_output
  character (len = nch_file) :: file_input, file_date, file_rad, fits_name, file_output, target_name
  character (len = nch_file) :: mask_sel
  character (len = 1) :: fiber_sel

  integer :: read_iostat

!!$ KEYWORDS
  real (kind=8) :: bis_span, bis_rv, ccf_contrast, ccf_rv, ccf_fwhm, ccf_rvc, ccf_noise
  integer (kind=8)  ::  drift_nbcos, ccf_maxccp, ccf_lines
  character (len=nch_key) :: ccf_mask, drift_qc
  character (len=nch_file) :: drift_ref_ccf
  real (kind=8) ::  berv, bjd, bervmax, airmass
  real :: drift_spe_rv, drift_rflux, drift_noise, drift_ref_rv, &
       drift_ccf_rv, drift_rv_used, drift_algo,ddvrms, &
       sn1, sn46, sn68, tel_mag
  real (kind=8):: ra_deg, dec_deg

  integer :: sts, sts_out
  character (len=64) :: comment
  character (len=1) :: nc

  integer :: n_com, date_opt
  character (len = nch_file) :: arch_red, opt_in

  call get_lun(lun_input)
  call get_lun(lun_output)
  call get_lun(lun_fits ) 
  
  n_com = iargc()
  
  if (n_com.lt.2 .or. n_com.gt.4) then
     write(*,*) 'Input 1:  file_input'
     write(*,*) 'Input 2:  file_output'
     write(*,*) 'Input 3 (opt):  arch_red'
     write(*,*) 'Input 4 (opt):  date flag (0 to ignore date in filename)'
     stop
  end if

  call getarg(1,file_input)
  call getarg(2,file_output)
  open(lun_input ,file=trim(file_input) , status='old')
  open(lun_output,file=trim(file_output), status='new')

  arch_red = archive_red
  date_opt = 1
  if (n_com.eq.4) then
     call getarg(3,arch_red)
     call getarg(4,opt_in); read (opt_in,*) date_opt
  end if
  
  read_iostat = 0 ; n = 0

  write(lun_output,'(A)') '#name 1.s  date 2.s  file 3.s  mask 4.s  fiber 5.  BJD 6  ' // &
       'BJD_short 7  RV_ccf 8  RV_ccf_drift 9  CCF_noise 10  DVR_noise 11 CCF_contrast 12  '//&
       'CCF_fwhm 13  BIS_rv 14 BIS_span15  berv 16  drift_ccf_rv 17  drift_rv_used 18  '//&
       'sn1 19  sn46 20  sn68 21  airmass 22 ra_deg 23 dec_deg 24 tel_mag 25'

  write(*,*) file_input
  do while (read_iostat.eq.0)
     read(lun_input,*,iostat=read_iostat) file_date, file_rad, mask_sel, fiber_sel, target_name
     if (read_iostat.ne.0) exit
     
     if (date_opt.eq.1) then
        fits_name = trim(arch_red)//trim(file_date)//'/'//trim(file_rad)//'_bis_'//trim(mask_sel)//'_'//fiber_sel//'.fits'
     else
        fits_name = trim(arch_red)//'/'//trim(file_rad)//'_bis_'//trim(mask_sel)//'_'//fiber_sel//'.fits'
     end if
     
     write(*,*) trim(fits_name)
     call fits_check(fits_name,sts_out)
     if (sts_out.eq.0) cycle
     
     call fits_actv(lun_fits,fits_name)
!!$     write(nc,'(I1)') 1 ; sts = 0
     ra_deg = 0._PR
     dec_deg = 0._PR
     
     call FTGKYD(lun_fits,'RA'   ,ra_deg     ,comment,sts);  sts = 0 !Right ascension (degrees as float)
     call FTGKYD(lun_fits,'DEC'  ,dec_deg    ,comment,sts);  sts = 0 !Declination (degrees as float)
     
     call FTGKYD(lun_fits,'HIERARCH TNG DRS BIS SPAN'     ,bis_span     ,comment,sts);  sts = 0 !Bisector velocity span [km/s]
     call FTGKYD(lun_fits,'HIERARCH TNG DRS BIS RV'       ,bis_rv       ,comment,sts);  sts = 0 !Bisector mean velocity [km/s]
     call FTGKYJ(lun_fits,'HIERARCH TNG DRS CCF MAXCPP'   ,ccf_maxccp   ,comment,sts);  sts = 0 !max count/pixel of CCF (e-)
     call FTGKYD(lun_fits,'HIERARCH TNG DRS CCF CONTRAST' ,ccf_contrast ,comment,sts);  sts = 0 !Contrast of  CCF (%)
     call FTGKYD(lun_fits,'HIERARCH TNG DRS CCF RV'       ,ccf_rv       ,comment,sts);  sts = 0 !Baryc RV (no drift correction) (km
     call FTGKYD(lun_fits,'HIERARCH TNG DRS CCF FWHM'     ,ccf_fwhm     ,comment,sts);  sts = 0 !FWHM of CCF (km/s)
     call FTGKYD(lun_fits,'HIERARCH TNG DRS CCF RVC'      ,ccf_rvc      ,comment,sts);  sts = 0 !Baryc RV (drift corrected) (km/s)
     call FTGKYJ(lun_fits,'HIERARCH TNG DRS CCF LINES'    ,ccf_lines    ,comment,sts);  sts = 0 !nr lines used
     call FTGKYS(lun_fits,'HIERARCH TNG DRS CCF MASK'     ,ccf_mask     ,comment,sts);  sts = 0 !Mask filename
     call FTGKYD(lun_fits,'HIERARCH TNG DRS CCF NOISE'    ,ccf_noise    ,comment,sts);  sts = 0 !Photon noise on CCF RV (km/s)

!!$     call FTGKYL(lun_fits,'TIME'                          ,time        ,comment,sts);  sts = 0 ! Time FITS file was written (hh:mm:ss)
!!$     call FTGKYL(lun_fits,'EXPSTART'                      ,expstart    ,comment,sts);  sts = 0 !Exposure start time (UT) (hh:mm:ss)

     call FTGKYD(lun_fits,'HIERARCH TNG DRS BERV'         ,berv         ,comment,sts);  sts = 0 !Barycentric Earth Radial Velocity
     call FTGKYD(lun_fits,'HIERARCH TNG DRS BJD'          ,bjd          ,comment,sts);  sts = 0 !Barycentric Julian Day
     call FTGKYD(lun_fits,'HIERARCH TNG DRS BERVMX'       ,bervmax      ,comment,sts);  sts = 0 !Maximum BERV

!! these are REAL single precision variables (if not differently specified)
     call FTGKYE(lun_fits,'HIERARCH TNG DRS DRIFT SPE RV' ,drift_spe_rv ,comment,sts);  sts = 0 !Th RV Drift [m/s]
     call FTGKYJ(lun_fits,'HIERARCH TNG DRS DRIFT NBCOS'  ,drift_nbcos  ,comment,sts);  sts = 0 !Th Drift Nbr cosmic detected (integer)
     call FTGKYE(lun_fits,'HIERARCH TNG DRS DRIFT RFLUX'  ,drift_rflux  ,comment,sts);  sts = 0 !Th Drift Flux ratio
     call FTGKYE(lun_fits,'HIERARCH TNG DRS DRIFT NOISE'  ,drift_noise  ,comment,sts);  sts = 0 !Th Drift photon noise [m/s]
     !call FTGKYS(lun_fits,'HIERARCH TNG DRS DRIFT REF CCF',drift_ref_ccf,comment,sts);  sts = 0 !fits 
     call FTGKYE(lun_fits,'HIERARCH TNG DRS DRIFT REF RV' ,drift_ref_rv ,comment,sts);  sts = 0 !RV of Th ref CCF [m/s]
     call FTGKYE(lun_fits,'HIERARCH TNG DRS DRIFT CCF RV' ,drift_ccf_rv ,comment,sts);  sts = 0 !CCF RV Drift [m/s]
     call FTGKYE(lun_fits,'HIERARCH TNG DRS DRIFT RV USED',drift_rv_used,comment,sts);  sts = 0 !Used RV Drift [m/s]
     !call FTGKYS(lun_fits,'HIERARCH TNG DRS DRIFT ALGO'   ,drift_algo   ,comment,sts);  sts = 0 !Algo used for Drift meas. (tex)
     call FTGKYE(lun_fits,'HIERARCH TNG DRS DRIFT QC'     ,drift_qc     ,comment,sts);  sts = 0 !QC on drift (text)
     call FTGKYE(lun_fits,'HIERARCH TNG DRS DVRMS'        ,ddvrms       ,comment,sts);  sts = 0 !Estimated RV uncertainty [m/s]

     call FTGKYE(lun_fits,'HIERARCH TNG DRS SPE EXT SN1'   ,sn1         ,comment,sts);  sts = 0 !S_N order center1
     call FTGKYE(lun_fits,'HIERARCH TNG DRS SPE EXT SN46' ,sn46         ,comment,sts);  sts = 0 !S_N order center46
     call FTGKYE(lun_fits,'HIERARCH TNG DRS SPE EXT SN68' ,sn68         ,comment,sts);  sts = 0 !S_N order center68
     call FTGKYE(lun_fits,'HIERARCH TNG TEL TARG MAG'     ,tel_mag      ,comment,sts);  sts = 0
     call FTGKYD(lun_fits,'AIRMASS'                       ,airmass      ,comment,sts);  sts = 0 !Air mass
     call fits_close(lun_fits)
     
  !!   write(lun_output,'(5(A,5x),2F18.8,11F13.6,3F9.1)') & 
     write(lun_output,'(5(A,5x),F16.6,F13.6,11F10.4,3F9.1,F8.4,2F12.6,F8.2)') & 
          trim(target_name),trim(file_date), trim(file_rad),  trim(mask_sel), trim(fiber_sel), &
          bjd,bjd-2450000._PR,ccf_rv,ccf_rvc,ccf_noise,ddvrms/1000._PR, &
          ccf_contrast, ccf_fwhm, bis_rv, bis_span, berv, &
          drift_ccf_rv ,drift_rv_used, &
          sn1,sn46,sn68, airmass, ra_deg, dec_deg, tel_mag
     
     !write(*,*) bis_span,bis_rv,ccf_maxccp,ccf_contrast,ccf_rv
     !write(*,*) ccf_fwhm,ccf_rvc,ccf_lines ,ccf_mask,ccf_noise ,berv 
     !write(*,*) bjd,bervmax,drift_spe_rv ,drift_nbcos,drift_rflux
     !write(*,*) drift_noise,drift_ref_ccf,drift_ref_rv ,drift_ccf_rv ,drift_rv_used
     !write(*,*) drift_algo ,drift_qc ,ddvrms
     
     
  end do
  close(lun_input)
  close(lun_output)
  


  
  call free_lun(lun_output)
  call free_lun(lun_input)
  call free_lun(lun_fits ) 
end program db_extract_keyword
