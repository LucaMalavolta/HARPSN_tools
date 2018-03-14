module Spectrum_Libraries
  use common
  use filepaths
  implicit none
  
  !> Solar spectra from Kurucz
  
  integer, parameter :: solkur_size = 762078  
  character (len=*), parameter :: solkur_file = tables_path // 'sol_kurucz.dat'
  real (kind=PR), dimension(solkur_size) :: solkur_wl = 0._PR, solkur_fl =0._Pr, solkur_dw =0._PR
  
  integer, parameter :: coelho2005_teff_nit = 15, coelho2005_logg_nit = 11, coelho2005_gfeh_nit = 8
  

  character(len=*), parameter, private :: cohelo2005_dir = stellib_path // 'Coelho2005/'
  
  integer, dimension(coelho2005_teff_nit), parameter, private :: &
       coelho2005_Teff_list = (/3500,3750,4000,4250,4500,4750,5000,5250,5500,5750,6000,6250,6500,6750,7000/)
  real , dimension(coelho2005_logg_nit), parameter, private :: &
       coelho2005_logg_list = (/0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0/)
  real , dimension(coelho2005_gfeh_nit), parameter, private :: &
       coelho2005_gfeh_list = (/-2.5, -2.0, -1.5, -1.0, -0.5, 0.00, 0.2, .5/)
  
  character (len=4), dimension(coelho2005_teff_nit), parameter, private  ::  coelho2005_Teff_char = & 
       (/'3500','3750','4000','4250','4500','4750','5000','5250','5500','5750','6000','6250','6500','6750','7000'/)
  character (len=2), dimension(coelho2005_logg_nit), parameter, private  :: &
       coelho2005_logg_char = (/'00','05','10','15','20','25','30','35','40','45','50'/)
  character (len=3), dimension(coelho2005_gfeh_nit), parameter, private  :: &
       coelho2005_gfeh_char = (/'m25','m20','m15','m10','m05','p00','p02','p05'/)
  
  


contains
  subroutine read_solar_kurucz
    integer :: io_unit, ii
    call get_lun(io_unit)
    write(*,*) tables_path //solkur_file
    open(io_unit, file=solkur_file, status='old')
    do ii=1,solkur_size
       read(io_unit,*) solkur_wl(ii), solkur_fl(ii)
    end do
    close(io_unit)
    call free_lun(io_unit)
    
    solkur_dw(2:solkur_size-1) = (solkur_wl(3:solkur_size) - solkur_wl(1:solkur_size-2))/2.000d0
    solkur_dw(1) = solkur_dw(2) ; solkur_dw(solkur_size) = solkur_dw(solkur_size-1)
    
    return
  end subroutine read_solar_kurucz
  
  subroutine coelho2005_select_template(teff_inp,logg_inp,gfeh_inp,type,outfile)
    
    !! Synthetic stellar library from Coelho et al. (2005) (A&A, 443, 735). 
    !! http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:astro-ph/0505511
    !! http://www.mpa-garching.mpg.de/PUBLICATIONS/DATA/SYNTHSTELLIB/synthetic_stellar_spectra.html
    
    !! type = 1 -> Solar scaled  ([alpha/Fe] = 0.0)
    !| type = 2 -> Alpha enhanced ([alpha/Fe] = 0.4) 
    !! type = 3 ->  Alpha enhanced with solar Calcium ([alpha/Fe] = 0.4 && [Ca/Fe] = 0.0)
    !! For given temperature, gravity and overall metallicty, the template file name is given
    
    integer :: type, i_t, i_g, i_m
    character (len=nch_file) :: outfile, dir, ext
    real, intent(in) :: teff_inp,logg_inp,gfeh_inp
    
     real   , dimension(coelho2005_teff_nit) :: adiff_t
     real   , dimension(coelho2005_logg_nit) :: adiff_g
     real   , dimension(coelho2005_gfeh_nit) :: adiff_m
     
    if (type.eq.1) then
       dir = 'solar_scaled/'
       ext = 'p00.ms.fits.gz'
    end if
    if (type.eq.2) then
       dir = 'alpha_enhanced/'
       ext = 'p04.ms.fits.gz'
    end if
    if (type.eq.3) then
       dir = 'alpha_enhanced_casolar/'
       ext = 'p04_casolar.ms.fits.gz'
    end if
    
    !! there is a bug with minloc, for xx_inp = xx_list(1) I get minloc=2 O__o 
    !! if I don't use the adiff escamotage...
    
    adiff_t = abs(coelho2005_Teff_list-Teff_inp);  i_t = minloc(adiff_t,dim=1)
    adiff_g = abs(coelho2005_logg_list-logg_inp);  i_g = minloc(adiff_g,dim=1)
    adiff_m = abs(coelho2005_gfeh_list-gfeh_inp);  i_m = minloc(adiff_m,dim=1)

    outfile = cohelo2005_dir//trim(dir)//coelho2005_teff_char(i_t)//'_'// &
         coelho2005_logg_char(i_g)//'_'//coelho2005_gfeh_char(i_m)//trim(ext)
    
    return
  end subroutine coelho2005_select_template
  
  real function  coelho2005_templates_vmic(logg_inp)
    ! from Coelho et al. 2005
    ! former get_vmic_of_templates(logg_inp)
    real :: logg_inp, logg
    real   , dimension(coelho2005_logg_nit) :: adiff_g
    integer :: i_g
    
    adiff_g = abs(coelho2005_logg_list-logg_inp)
    i_g = minloc(adiff_g,dim=1)
    logg =  coelho2005_logg_list(i_g)
    coelho2005_templates_vmic = 1.8
    if (logg.ge.3.0) coelho2005_templates_vmic = 1.0
    if (logg.le.2.5 .and. logg.ge. 1.5)  coelho2005_templates_vmic = 1.8
    if (logg.le. 1.0)  coelho2005_templates_vmic = 2.5
    return
    
  end function coelho2005_templates_vmic
  

end module Spectrum_Libraries
