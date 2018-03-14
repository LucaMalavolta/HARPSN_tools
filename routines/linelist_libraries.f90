module Linelist_Libraries
  use Common
  use Filepaths

  implicit none

  integer :: llist_n 
  real (PR), dimension(:), allocatable :: llist_wl, llist_ep, llist_lgf, llist_elm
  character (len=nch_file) :: llist_name, llist_file
  
  !Nave et al. 1994: table dimensions
  integer, parameter :: N941_size =  846, N942_size = 9759, N944_size = 9759, N945_size =  125

  character (len=*), parameter :: tmp_file = './temp.dat'

  !Nave et al. 1994 table positions:
  character (len=*), parameter :: &
       N941_file = code_path // 'tables/nave1994_table1', &
       N942_file = code_path // 'tables/nave1994_table2', &
       N944_file = code_path // 'tables/nave1994_table4', &
       N945_file = code_path // 'tables/nave1994_table5'

  !Nave et al. 1994 table 1
  character (len=1),  dimension(N941_size) :: N941_n_J
  character (len=10), dimension(N941_size) :: N941_Term
  character (len=15), dimension(N941_size) :: N941_Config
  integer, dimension(N941_size) :: N941_J
  real (PR), dimension(N941_size) :: N941_Level,N941_e_Level

  !Nave et al. 1994 table 2
  character (len=1), dimension(N942_size) :: N942_n_Term1,N942_n_J1,N942_q_Sigma, N942_r_log_gf
  character (len=5), dimension(N942_size) :: N942_Blen
  character (len=6), dimension(N942_size) :: N942_I
  character (len=10), dimension(N942_size) :: N942_Term1,N942_Term2
  integer, dimension(N942_size) :: N942_Num,N942_J1,N942_J2,N942_So_R
  real (PR), dimension(N942_size) :: N942_Lvac,N942_Lair,N942_Lo_R,N942_Sigma,N942_El,N942_Eu,N942_Elc, N942_Euc, N942_log_gf

  !Nave et al. 1994 table 4
  character (len=4), dimension(N944_size) :: N944_Num
  real (PR), dimension(N944_size) :: N944_Lambda,N944_Sigma

  !Nave et al. 1994 table 5
  real (PR), dimension(N945_size) :: N945_I,N945_Lambda,N945_Sigma


  ! Schuler 2011b line list
  character (len=*), parameter :: schuler2011b_name = 'schuler2011b_lines', &
    schuler2011b_file = tables_path // 'schuler2011b_lines'
  integer, parameter :: schuler2011_n = 190
  real (PR), dimension(schuler2011_n):: schuler2011_wl, schuler2011_ep,schuler2011_lgf,schuler2011_sew,schuler2011_elm


  ! Sousa line list
  character (len=*), parameter :: sousa2010_allelm_name = 'sousa2010_allelm_table1', &
  sousa2010_allelm_file = tables_path // 'sousa2010_allelm_table1'
  integer, parameter :: sousa_n = 498
  real (PR), dimension(sousa_n):: sousa_wl, sousa_ep,sousa_lgf,sousa_sew,sousa_elm
    
  ! Sousa log gf line list
  character (len=*) , parameter:: sousa2010_THsel_name = 'sousa2010_original_thsel', &
  sousa2010_THsel_file = tables_path // 'sousa2010_original_thsel'
  integer, parameter :: south_n = 230
  real (PR), dimension(south_n):: south_wl, south_ep,south_lgf,south_sew,south_elm

  ! Ivans+ 2003 line list
  character (len=*) , parameter:: Ivans2003_name = 'Ivans2003_all' , &
   Ivans2003_file = tables_path // 'Ivans2003_all' 
  integer, parameter :: Ivans2003_n = 437
  real (PR), dimension(Ivans2003_n):: Ivans2003_wl, Ivans2003_ep,Ivans2003_lgf,Ivans2003_sew,Ivans2003_elm

  ! Ramirez & Cohen 2003 line list
  character (len=*) , parameter:: RamCo2003_name = 'RamirezCohen2003_all', &
    RamCo2003_file = tables_path // 'RamirezCohen2003_all'
  integer, parameter :: RamCo2003_n = 334
  real (PR), dimension(RamCo2003_n):: RamCo2003_wl, RamCo2003_ep,RamCo2003_lgf,RamCo2003_sew,RamCo2003_elm

  ! Lai+ 2011 line list
  character (len=*) , parameter:: Lai2011_name = 'Lai2011_all', &
    Lai2011_file = tables_path // 'Lai2011_all'
  integer, parameter :: Lai2011_n = 334
  real (PR), dimension(Lai2011_n):: Lai2011_wl, Lai2011_ep,Lai2011_lgf,Lai2011_sew,Lai2011_elm

  ! Ivans 2006 line list
  character (len=*), parameter ::ivans_name =  'ivans2006_table1' , &
    ivans_file = tables_path // 'ivans2006_table1' 
  integer, parameter :: ivans_n = 645
  real (PR), dimension(ivans_n) :: ivans_wl,ivans_ep,ivans_lgf,ivans_elm


  ! TiI from Lawler2013+ and TiII from Wood2013+
  character (len=*), parameter ::Ti_LW2013_name =  'Ti_LawlerWood2013_join.dat' ,&
    Ti_LW2013_file = tables_path // 'Ti_LawlerWood2013_join.dat'  
  integer, parameter :: Ti_LW2013_n = 69
  real (PR), dimension(Ti_LW2013_n) :: Ti_LW2013_wl,Ti_LW2013_ep,Ti_LW2013_lgf,Ti_LW2013_elm,Ti_LW2013_sew  
  
  ! Sousa 2011 Stable website original 
  character (len=*), parameter ::sousa_org_name =  'sousa_webstable_orig' , &
    sousa_org_file = tables_path // 'sousa_webstable_orig' 
  integer, parameter :: sousa_org_n = 298
  real (PR), dimension(sousa_org_n) :: sousa_org_wl,sousa_org_ep,sousa_org_lgf,sousa_org_elm,sousa_org_sew  
  
    ! Sousa 2011 Stable website FAKE 
  character (len=*), parameter ::sousa_fake_name =  'sousa_webstable_fake' , &
    sousa_fake_file = tables_path // 'sousa_webstable_fake' 
  integer, parameter :: sousa_fake_n = 298
  real (PR), dimension(sousa_fake_n) :: sousa_fake_wl,sousa_fake_ep,sousa_fake_lgf,sousa_fake_elm,sousa_fake_sew  
  
  ! Sousa 2011 Stable website refined 
  ! Removed line: 123   26.1 5427.81 6.72    -1.481   4.4
  ! From the original linelist
  character (len=*), parameter ::sousa_ref_name =  'sousa_webstable_ref' , &
    sousa_ref_file = tables_path // 'sousa_webstable_ref' 
  integer, parameter :: sousa_ref_n = 298
  real (PR), dimension(sousa_ref_n) :: sousa_ref_wl,sousa_ref_ep,sousa_ref_lgf,sousa_ref_elm,sousa_ref_sew  

  ! Sousa 2011 Stable website refined obtained using MOOG2010-2013
  ! Removed line: 123   26.1 5427.81 6.72    -1.481   4.4
  ! From the original linelist
  character (len=*), parameter :: sousa_ref_2010_name =  'sousa_webstable_ref_2010' , &
    sousa_ref_2010_file = tables_path // 'sousa_webstable_ref_2010'
  integer, parameter :: sousa_ref_2010_n = 298
  real (PR), dimension(sousa_ref_2010_n) ::&
  sousa_ref_2010_wl,sousa_ref_2010_ep,sousa_ref_2010_lgf,sousa_ref_2010_elm,sousa_ref_2010_sew  

  !Tsantaki 2013 cool-stars original
  character (len=*), parameter ::tsan2013_org_name =  'tsantaki2013_orig' , &
    tsan2013_org_file = tables_path // 'tsantaki2013_orig' 
  integer, parameter :: tsan2013_org_n = 137
  real (PR), dimension(tsan2013_org_n) :: tsan2013_org_wl,tsan2013_org_ep,tsan2013_org_lgf,tsan2013_org_elm,tsan2013_org_sew  
  
  !Tsantaki 2013 cool-stars refined
  character (len=*), parameter ::tsan2013_ref_name =  'tsantaki2013_ref' , &
    tsan2013_ref_file = tables_path // 'tsantaki2013_ref' 
  integer, parameter :: tsan2013_ref_n = 137
  real (PR), dimension(tsan2013_ref_n) :: tsan2013_ref_wl,tsan2013_ref_ep,tsan2013_ref_lgf,tsan2013_ref_elm,tsan2013_ref_sew  

  !Neves 2009 refractory elements refined
  character (len=*), parameter ::neves2009_ref_name =  'neves2009_ref' , &
    neves2009_ref_file = tables_path // 'neves2009_ref' 
  integer, parameter :: neves2009_ref_n = 180
  real (PR), dimension(neves2009_ref_n) :: neves2009_ref_wl,neves2009_ref_ep,neves2009_ref_lgf,neves2009_ref_elm,neves2009_ref_sew  
  
  ! Bocek Topcu Iron list
  character (len=*), parameter :: BocekTopcu_Fe_name =  'BocekTopcu_Fe_orig' , &
     BocekTopcu_Fe_file = tables_path // 'BocekTopcu_Fe_orig' 
  integer, parameter :: BocekTopcu_Fe_n = 73
  real (PR), dimension(BocekTopcu_Fe_n) :: BocekTopcu_Fe_wl,BocekTopcu_Fe_ep,BocekTopcu_Fe_lgf,BocekTopcu_Fe_elm,BocekTopcu_Fe_sew  
  
  
contains
  subroutine nave1994_read

    character (len=*), parameter :: fmt_table1 = '(A15,A10,I1,A1,F10.3,F4.3)', &
         fmt_table2 = '(I4,A10,A1,A10,I1,A1,I1,A6,F10.4,F10.4,F6.1,F10.3,A1,I4,F5.3,F5.3,F9.3,F9.3,F6.2,A1,A5)', &
         fmt_table4 = '(F10.4,F9.3,A4)', &
         fmt_table5 = '(F5.2,F10.4,F10.4)'
    integer :: unit, i
    !I'm following the same name convention used in the Vizier descritions of the tables
    ! http://cdsarc.u-strasbg.fr/viz-bin/Cat?cat=J%2FApJS%2F94%2F221&target=readme&
    !N94n -> N for Nace, 94 for the year (1994), n for the table number
    !reading from table 1

    call get_lun(unit)

    open(unit, file=N941_file, status='old')
    do i = 1,N941_size
       read(unit,fmt_table1) N941_Config(i), N941_Term(i), N941_J(i), N941_n_J(i), N941_Level(i), N941_e_Level(i)
    end do
    close(unit)

    !reading from table 2
    open(unit, file=N942_file, status='old')
    do i = 1,N942_size   
       read(unit,fmt_table2) N942_Num(i), N942_Term1(i), N942_n_Term1(i), N942_Term2(i), &
            N942_J1(i), N942_n_J1(i), N942_J2(i), N942_I(i), N942_Lvac(i), N942_Lair(i), N942_Lo_R(i), &
            N942_Sigma(i), N942_q_Sigma(i), N942_So_R(i), N942_El(i), N942_Eu(i), N942_Elc(i), &
            N942_Euc(i), N942_log_gf(i), N942_r_log_gf(i), N942_Blen(i)
    end do
    close(unit)

    !reading from table 4
    open(unit, file=N944_file, status='old')
    do i = 1,N944_size   
       read(unit,fmt_table4) N944_Lambda(i), N944_Sigma(i), N944_Num(i)
    end do
    close(unit)

    !reading from Table 5
    open(unit, file=N945_file, status='old')
    do i = 1,N945_size
       read(unit,fmt_table5) N945_I(i), N945_Lambda(i), N945_Sigma(i)
    end do
    close(unit)
    call free_lun(unit)

    return

  end subroutine nave1994_read


  subroutine schuler2011_read
    integer :: unit,nl

    call get_lun(unit)
    open(unit,file=schuler2011b_file,status='old')

    schuler2011_sew = 0.00
    do nl = 1,schuler2011_n
       read(unit,*) schuler2011_elm(nl),schuler2011_wl(nl),schuler2011_ep(nl), schuler2011_lgf(nl), schuler2011_sew(nl)
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine schuler2011_read

  subroutine sousa2010_read
    integer :: unit,nl
    character (len=6) :: sousa_ch6
    real :: elm

    call get_lun(unit)
    open(unit,file=sousa2010_allelm_file,status='old')

    sousa_sew = 0.00
    do nl = 1,sousa_n
       read(unit,*) sousa_wl(nl),sousa_ep(nl), sousa_lgf(nl), sousa_ch6
       call elm_char2real(sousa_ch6,elm)
       sousa_elm(nl) = elm
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine sousa2010_read


  subroutine sousa2010_THsel_read
    integer :: unit,nl, nl_int

    call get_lun(unit)
    open(unit,file=sousa2010_THsel_file,status='old')
    south_sew = 0._PR
    do nl = 1,south_n
       read(unit,*) nl_int,south_elm(nl),south_wl(nl),south_ep(nl), &
            south_lgf(nl),south_sew(nl)
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine sousa2010_THsel_read

  subroutine RamirezCohen2003_read
    integer :: unit,nl, nl_int

    call get_lun(unit)
    open(unit,file=RamCo2003_file,status='old')
    RamCo2003_sew = 0._PR
    do nl = 1,RamCo2003_n
       read(unit,*) nl_int,RamCo2003_elm(nl),RamCo2003_wl(nl),RamCo2003_ep(nl), &
            RamCo2003_lgf(nl),RamCo2003_sew(nl)
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine RamirezCohen2003_read

  subroutine Ivans2003_read
    integer :: unit,nl, nl_int

    call get_lun(unit)
    open(unit,file=Ivans2003_file,status='old')
    Ivans2003_sew = 0._PR
    do nl = 1,Ivans2003_n
       read(unit,*) nl_int,Ivans2003_elm(nl),Ivans2003_wl(nl),Ivans2003_ep(nl), &
            Ivans2003_lgf(nl),Ivans2003_sew(nl)
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine Ivans2003_read

  subroutine Lai2011_read
    integer :: unit,nl, nl_int

    call get_lun(unit)
    open(unit,file=Lai2011_file,status='old')
    Lai2011_sew = 0._PR
    do nl = 1,Lai2011_n
       read(unit,*) nl_int,Lai2011_elm(nl),Lai2011_wl(nl),Lai2011_ep(nl), &
            Lai2011_lgf(nl),Lai2011_sew(nl)
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine Lai2011_read


  subroutine elm_real2char(elm,ch,space_on)
    character(len=6) :: ch
    character(len=3) :: ch_ion
    integer, intent(in), optional :: space_on
    real :: elm, elm_ion
    integer :: elm_int
    ch = '      '

    elm_int = int(elm)
    elm_ion = elm-elm_int


    if (elm_int.eq. 1) ch = 'H'
    if (elm_int.eq. 2) ch = 'He' 
    if (elm_int.eq. 3) ch = 'Li' 
    if (elm_int.eq. 4) ch = 'Be' 
    if (elm_int.eq. 5) ch = 'B' 
    if (elm_int.eq. 6) ch = 'C' 
    if (elm_int.eq. 7) ch = 'N' 
    if (elm_int.eq. 8) ch = 'O' 
    if (elm_int.eq. 9) ch = 'F' 
    if (elm_int.eq. 10) ch = 'Ne' 
    if (elm_int.eq. 11) ch = 'Na' 
    if (elm_int.eq. 12) ch = 'Mg' 
    if (elm_int.eq. 13) ch = 'Al' 
    if (elm_int.eq. 14) ch = 'Si' 
    if (elm_int.eq. 15) ch = 'P' 
    if (elm_int.eq. 16) ch = 'S' 
    if (elm_int.eq. 17) ch = 'Cl' 
    if (elm_int.eq. 18) ch = 'Ar' 
    if (elm_int.eq. 19) ch = 'K' 
    if (elm_int.eq. 20) ch = 'Ca' 
    if (elm_int.eq. 21) ch = 'Sc' 
    if (elm_int.eq. 22) ch = 'Ti' 
    if (elm_int.eq. 23) ch = 'V' 
    if (elm_int.eq. 24) ch = 'Cr' 
    if (elm_int.eq. 25) ch = 'Mn' 
    if (elm_int.eq. 26) ch = 'Fe' 
    if (elm_int.eq. 27) ch = 'Co' 
    if (elm_int.eq. 28) ch = 'Ni' 
    if (elm_int.eq. 29) ch = 'Cu'     

    if (elm_ion.eq.0.0) ch_ion = 'I  ' 
    if (elm_ion.eq.0.1) ch_ion = 'II ' 
    if (elm_ion.eq.0.2) ch_ion = 'III' 
    if (elm_ion.eq.0.3) ch_ion = 'IV ' 
    if (elm_ion.eq.0.4) ch_ion = 'V  ' 

    if (present(space_on))  then
       ch = trim(ch) // ' ' // ch_ion
    else
       ch = trim(ch) // ch_ion
    end if
    return
  end subroutine elm_real2char

  subroutine elm_char2real(ch,elm)
    character(len=6) :: ch
    character(len=4) :: ch_ion,ch_elm
    real :: elm, elm_ion
    integer :: elm_int, i, nn

    elm_int = 0
    elm_ion = 0.0
    nn = 1
    do i=2,6
       !' Iodium starts with I, so we must skip the first letter...'
       if (ch(i:i).eq.'I') exit
       nn=nn+1
    end do
    ch_elm = trim(ch(1:nn))
    ch_ion = trim(ch(nn+1:6))

    if (ch_elm.eq. 'H'  ) elm_int =1 
    if (ch_elm.eq. 'He' ) elm_int =2 
    if (ch_elm.eq. 'Li' ) elm_int =3 
    if (ch_elm.eq. 'Be' ) elm_int =4 
    if (ch_elm.eq. 'B'  ) elm_int =5 
    if (ch_elm.eq. 'C'  ) elm_int =6 
    if (ch_elm.eq. 'N'  ) elm_int =7 
    if (ch_elm.eq. 'O'  ) elm_int =8 
    if (ch_elm.eq. 'F'  ) elm_int =9 
    if (ch_elm.eq. 'Ne' ) elm_int =10
    if (ch_elm.eq. 'Na' ) elm_int =11
    if (ch_elm.eq. 'Mg' ) elm_int =12
    if (ch_elm.eq. 'Al' ) elm_int =13
    if (ch_elm.eq. 'Si' ) elm_int =14
    if (ch_elm.eq. 'P'  ) elm_int =15
    if (ch_elm.eq. 'S'  ) elm_int =16
    if (ch_elm.eq. 'Cl' ) elm_int =17
    if (ch_elm.eq. 'Ar' ) elm_int =18
    if (ch_elm.eq. 'K'  ) elm_int =19
    if (ch_elm.eq. 'Ca' ) elm_int =20
    if (ch_elm.eq. 'Sc' ) elm_int =21
    if (ch_elm.eq. 'Ti' ) elm_int =22
    if (ch_elm.eq. 'V'  ) elm_int =23
    if (ch_elm.eq. 'Cr' ) elm_int =24
    if (ch_elm.eq. 'Mn' ) elm_int =25
    if (ch_elm.eq. 'Fe' ) elm_int =26
    if (ch_elm.eq. 'Co' ) elm_int =27
    if (ch_elm.eq. 'Ni' ) elm_int =28
    if (ch_elm.eq. 'Cu' ) elm_int =29   

    if (ch_ion .eq. 'I'  ) elm_ion=0.0
    if (ch_ion .eq. 'II' ) elm_ion=0.1
    if (ch_ion .eq. 'III') elm_ion=0.2
    if (ch_ion .eq. 'IV' ) elm_ion=0.3
    if (ch_ion .eq. 'V'  ) elm_ion=0.4

    elm = real(elm_int)+elm_ion

    return
  end subroutine elm_char2real

  subroutine ivans_read
    integer :: unit, sts, n

    character (len=*), parameter :: ivans_fmt = '(A2,F5.2,I3,1x,A2,F8.2,4x,F4.2,F8.3)'
    character (len=4) :: ivans_ch2_1, ivans_ch_2
    real :: asa_el, log_gf
    integer :: elm_i

!!$ From Ivans et al. 2006 for HD221170
!!$   1-  2  A2    ---     El        Element
!!$   4-  7  F4.2  ---     Abd       Adopted Solar Abundance for El (H=12)
!!$   9- 10  I2    ---     Z         Atomic number of El
!!$  12- 13  A2    ---     Ion       Ionization State
!!$  15- 21  F7.2  0.1nm   Wave      ? Wavelength in Angstroms
!!$  23- 24  A2    ---   f_Wave     Flag on Wavelength (1) - SUPPRESSED
!!$  26- 29  F4.2  eV      ExtPot    ? Excitation Potential
!!$  31- 37  F7.3  [-]     log(gf)   ? Log of the Oscillator Strength

    call get_lun(unit)
    sts = 0
    open(unit,file=ivans_file,status='old')

    !! we skip the first two lines since they contain anomalous values...
    read(unit, *, iostat=sts) ivans_ch2_1
    read(unit, *, iostat=sts) ivans_ch2_1
    do n=1,ivans_n
       read(unit,ivans_fmt) ivans_ch2_1, asa_el, elm_i, ivans_ch_2, ivans_wl(n), ivans_ep(n), log_gf
       ivans_elm(n) = real(elm_i)
       if (ivans_ch_2.eq.'II') ivans_elm(n) = ivans_elm(n) + 0.1
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine ivans_read

  !!----------------

  subroutine Ti_LawlerWood_read
    integer :: unit,nl, nl_int

    call get_lun(unit)
    open(unit,file=Ti_LW2013_file,status='old')
    Ti_LW2013_sew = 0._PR
    do nl = 1,Ti_LW2013_n
       read(unit,*) nl_int,Ti_LW2013_elm(nl),Ti_LW2013_wl(nl),Ti_LW2013_ep(nl), &
            Ti_LW2013_lgf(nl),Ti_LW2013_sew(nl)
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine Ti_LawlerWood_read


  !!----------------
   subroutine Sousa_WebStable_orig_read
    integer :: unit,nl, nl_int

    call get_lun(unit)
    open(unit,file=sousa_org_file,status='old')
    sousa_org_sew = 0._PR
    do nl = 1,sousa_org_n
       read(unit,*) nl_int,sousa_org_elm(nl),sousa_org_wl(nl),sousa_org_ep(nl), &
            sousa_org_lgf(nl),sousa_org_sew(nl)
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine Sousa_WebStable_orig_read
  
  subroutine Sousa_WebStable_fake_read
    integer :: unit,nl, nl_int

    call get_lun(unit)
    open(unit,file=sousa_fake_file,status='old')
    sousa_fake_sew = 0._PR
    do nl = 1,sousa_fake_n
       read(unit,*) nl_int,sousa_fake_elm(nl),sousa_fake_wl(nl),sousa_fake_ep(nl), &
            sousa_fake_lgf(nl),sousa_fake_sew(nl)
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine Sousa_WebStable_fake_read
  
  
  subroutine Sousa_WebStable_ref_read
    integer :: unit,nl, nl_int

    call get_lun(unit)
    open(unit,file=sousa_ref_file,status='old')
    sousa_ref_sew = 0._PR
    do nl = 1,sousa_ref_n
       read(unit,*) nl_int,sousa_ref_elm(nl),sousa_ref_wl(nl),sousa_ref_ep(nl), &
            sousa_ref_lgf(nl),sousa_ref_sew(nl)
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine Sousa_WebStable_ref_read


  subroutine Sousa_WebStable_ref_2010_read
    integer :: unit,nl, nl_int

    call get_lun(unit)
    open(unit,file=sousa_ref_2010_file,status='old')
    sousa_ref_2010_sew = 0._PR
    do nl = 1,sousa_ref_2010_n
       read(unit,*) nl_int,sousa_ref_2010_elm(nl),sousa_ref_2010_wl(nl),sousa_ref_2010_ep(nl), &
            sousa_ref_2010_lgf(nl),sousa_ref_2010_sew(nl)
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine Sousa_WebStable_ref_2010_read

    subroutine Tsantaki2013_orig_read
    integer :: unit,nl, nl_int

    call get_lun(unit)
    open(unit,file=tsan2013_org_file,status='old')
    tsan2013_org_sew = 0._PR
    do nl = 1,tsan2013_org_n
       read(unit,*) nl_int,tsan2013_org_elm(nl),tsan2013_org_wl(nl),tsan2013_org_ep(nl), &
            tsan2013_org_lgf(nl),tsan2013_org_sew(nl)
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine Tsantaki2013_orig_read

  
  subroutine Tsantaki2013_ref_read
    integer :: unit,nl, nl_int

    call get_lun(unit)
    open(unit,file=tsan2013_ref_file,status='old')
    tsan2013_ref_sew = 0._PR
    do nl = 1,tsan2013_ref_n
       read(unit,*) nl_int,tsan2013_ref_elm(nl),tsan2013_ref_wl(nl),tsan2013_ref_ep(nl), &
            tsan2013_ref_lgf(nl),tsan2013_ref_sew(nl)
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine Tsantaki2013_ref_read

  subroutine Neves2009_ref_read
    integer :: unit,nl, nl_int

    call get_lun(unit)
    open(unit,file=neves2009_ref_file,status='old')
    neves2009_ref_sew = 0._PR
    do nl = 1,neves2009_ref_n
       read(unit,*) nl_int,neves2009_ref_elm(nl),neves2009_ref_wl(nl),neves2009_ref_ep(nl), &
            neves2009_ref_lgf(nl),neves2009_ref_sew(nl)
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine Neves2009_ref_read
  
   
  subroutine BocekTopcu_Fe_read
    integer :: unit,nl, nl_int

    call get_lun(unit)
    open(unit,file=BocekTopcu_Fe_file,status='old')
    BocekTopcu_Fe_sew = 0._PR
    do nl = 1,BocekTopcu_Fe_n
       read(unit,*) nl_int,BocekTopcu_Fe_elm(nl),BocekTopcu_Fe_wl(nl),BocekTopcu_Fe_ep(nl), &
            BocekTopcu_Fe_lgf(nl),BocekTopcu_Fe_sew(nl)
    end do
    close(unit)
    call free_lun(unit)

    return
  end subroutine BocekTopcu_Fe_read
  
  


  subroutine linelist_select(list_input, no_read)
    character (len=nch_file), intent(in) :: list_input
    integer :: lun_tmp, nl, io_stat, ii, nl_int
    logical, intent(in), optional :: no_read
    logical :: nr_f
    llist_n = 0

    if (allocated(llist_wl) .or. allocated(llist_ep) .or. allocated(llist_lgf) .or. allocated(llist_elm)) llist_n = -1
    if (allocated(llist_wl)) deallocate(llist_wl)
    if (allocated(llist_ep)) deallocate(llist_ep)
    if (allocated(llist_lgf)) deallocate(llist_lgf)
    if (allocated(llist_elm)) deallocate(llist_elm)
    if (llist_n.le.-1) then
       llist_n = 0
       return
    end if

    nr_f = .false.
    if (present(no_read)) nr_f = no_read
    
    
    llist_n = 0
    name_case: select case (trim(list_input))
    case ('sousa2010')
       write(*,*) 'Using Sousa et al. 2010 log(gf)'
       llist_name = sousa2010_allelm_name
       llist_file = sousa2010_allelm_file
       llist_n = sousa_n
       if (nr_f) exit name_case
       
       call sousa2010_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = sousa_wl  
       llist_ep  = sousa_ep  
       llist_lgf = sousa_lgf 
       llist_elm = sousa_elm 

    case ('sousa2010_THsel')
       write(*,*) 'Using of Sousa et al. 2010 log(gf), selected against Th-Ar cont'
       llist_name = sousa2010_THsel_name
       llist_file = sousa2010_THsel_file
       llist_n = south_n
       if (nr_f) exit name_case

       call sousa2010_THsel_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = south_wl  
       llist_ep  = south_ep  
       llist_lgf = south_lgf 
       llist_elm = south_elm 

    case ('schuler2011')
       llist_name = schuler2011b_name
       llist_file = schuler2011b_file
       llist_n = schuler2011_n
       if (nr_f) exit name_case
       
       call schuler2011_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = schuler2011_wl  
       llist_ep  = schuler2011_ep  
       llist_lgf = schuler2011_lgf 
       llist_elm = schuler2011_elm 

    case ('ramirez_cohen2003')
       llist_name = RamCo2003_name
       llist_file = RamCO2003_file
       llist_n = RamCo2003_n
       if (nr_f) exit name_case
  
       call RamirezCohen2003_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = RamCo2003_wl  
       llist_ep  = RamCo2003_ep  
       llist_lgf = RamCo2003_lgf 
       llist_elm = RamCo2003_elm 

    case ('ivans2003')
       llist_name = Ivans2003_name
       llist_file = Ivans2003_file
       llist_n = Ivans2003_n
       if (nr_f) exit name_case

       call Ivans2003_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = Ivans2003_wl  
       llist_ep  = Ivans2003_ep  
       llist_lgf = Ivans2003_lgf 
       llist_elm = Ivans2003_elm 

    case ('lai2011')
       llist_name = Lai2011_name
       llist_file = Lai2011_file
       llist_n = Lai2011_n
       if (nr_f) exit name_case

       call Lai2011_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = Lai2011_wl  
       llist_ep  = Lai2011_ep  
       llist_lgf = Lai2011_lgf 
       llist_elm = Lai2011_elm 


    case ('ti_lawlerwood')
       llist_name = Ti_LW2013_name
       llist_file = Ti_LW2013_file
       llist_n = Ti_LW2013_n
       if (nr_f) exit name_case

       write(*,*) 'Using Ti experimental log(gf) from Lawler2013+ and Wood2013+'
       call Ti_LawlerWood_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = Ti_LW2013_wl  
       llist_ep  = Ti_LW2013_ep  
       llist_lgf = Ti_LW2013_lgf 
       llist_elm = Ti_LW2013_elm 

    case ('sousa_webstable_orig')
       llist_name = sousa_org_name
       llist_file = sousa_org_file
       llist_n = sousa_org_n
       if (nr_f) exit name_case

       write(*,*) 'Using original values of Sousa et al. 2010 log(gf)'
       call Sousa_WebStable_orig_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = sousa_org_wl  
       llist_ep  = sousa_org_ep  
       llist_lgf = sousa_org_lgf 
       llist_elm = sousa_org_elm 
       
    case ('sousa_webstable_fake')
       llist_name = sousa_fake_name
       llist_file = sousa_fake_file
       llist_n = sousa_fake_n
       if (nr_f) exit name_case

       write(*,*) 'Using original values of Sousa et al. 2010 log(gf)'
       call Sousa_WebStable_fake_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = sousa_fake_wl  
       llist_ep  = sousa_fake_ep  
       llist_lgf = sousa_fake_lgf 
       llist_elm = sousa_fake_elm 
    
    case ('sousa_webref')
       llist_name = sousa_ref_name
       llist_file = sousa_ref_file
       llist_n = sousa_ref_n
       if (nr_f) exit name_case

       write(*,*) 'Using Revisited values of Sousa et al. 2010 log(gf)'
       call Sousa_WebStable_ref_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = sousa_ref_wl  
       llist_ep  = sousa_ref_ep  
       llist_lgf = sousa_ref_lgf 
       llist_elm = sousa_ref_elm 

    case ('sousa_webref_2010')
       llist_name = sousa_ref_2010_name
       llist_file = sousa_ref_2010_file
       llist_n = sousa_ref_2010_n
       if (nr_f) exit name_case

       write(*,*) 'Using NEW Revisited values of Sousa et al. 2010 log(gf) with MOOG2010'
       call Sousa_WebStable_ref_2010_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = sousa_ref_2010_wl  
       llist_ep  = sousa_ref_2010_ep  
       llist_lgf = sousa_ref_2010_lgf 
       llist_elm = sousa_ref_2010_elm 

    case ('tsantaki2013_orig')
       llist_name = tsan2013_org_name
       llist_file = tsan2013_org_file
       llist_n = tsan2013_org_n
       if (nr_f) exit name_case

       write(*,*) 'Using Revisited values of Sousa et al. 2010 log(gf)'
       call Tsantaki2013_orig_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = tsan2013_org_wl  
       llist_ep  = tsan2013_org_ep  
       llist_lgf = tsan2013_org_lgf 
       llist_elm = tsan2013_org_elm      
       
       
    case ('tsantaki2013_ref')
       llist_name = tsan2013_ref_name
       llist_file = tsan2013_ref_file
       llist_n = tsan2013_ref_n
       if (nr_f) exit name_case

       write(*,*) 'Using Revisited values of Sousa et al. 2010 log(gf)'
       call Tsantaki2013_ref_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = tsan2013_ref_wl  
       llist_ep  = tsan2013_ref_ep  
       llist_lgf = tsan2013_ref_lgf 
       llist_elm = tsan2013_ref_elm      

    case ('neves2009_ref')
       llist_name = neves2009_ref_name
       llist_file = neves2009_ref_file
       llist_n = neves2009_ref_n
       if (nr_f) exit name_case

       write(*,*) 'Using Revisited values of Neves et al. 2009 log(gf)'
       call Neves2009_ref_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = neves2009_ref_wl  
       llist_ep  = neves2009_ref_ep  
       llist_lgf = neves2009_ref_lgf 
       llist_elm = neves2009_ref_elm     
       
    case ('bocektopcu_fe')
       llist_name = BocekTopcu_Fe_name
       llist_file = BocekTopcu_Fe_file
       llist_n = BocekTopcu_Fe_n
       if (nr_f) exit name_case

       write(*,*) 'Using NEW Revisited values of Sousa et al. 2010 log(gf)'
       call BocekTopcu_Fe_read
       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))

       llist_wl  = BocekTopcu_Fe_wl  
       llist_ep  = BocekTopcu_Fe_ep  
       llist_lgf = BocekTopcu_Fe_lgf 
       llist_elm = BocekTopcu_Fe_elm 
       
    case default
       call get_lun(lun_tmp)
       open(lun_tmp,file=trim(list_input),status='old')
       llist_n = 0
       io_stat = 0
       do while (io_stat.eq.0) 
          read(lun_tmp,*,iostat=io_stat) ii
          if (io_stat.eq.0) llist_n = llist_n+1
       end do
       rewind(lun_tmp)

       allocate(llist_wl(llist_n),llist_ep(llist_n),llist_lgf(llist_n),llist_elm(llist_n))
       llist_wl = 1._PR
       do nl = 1,llist_n
          read(lun_tmp,*) nl_int, llist_elm(nl), llist_wl(nl),llist_ep(nl), llist_lgf(nl)
       end do
       close (lun_tmp)
       call free_lun(lun_tmp)
    end select name_case

    if (llist_n .eq.0 ) then 
       write(*,*) ' ERROR: no lines have been read'
       stop
    end if
    return
  end subroutine linelist_select


end module Linelist_Libraries

