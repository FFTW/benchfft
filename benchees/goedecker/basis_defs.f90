!{\src2tex{textfont=tt}}
!!****f* ABINIT/basis_defs
!! NAME
!! basis_defs
!!
!! FUNCTION
!! This module contains definitions for a number of named constants and
!! physical constants
!!
!! COPYRIGHT
!! Copyright (C) 2000-2001 ABINIT group (HM, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! Of the named constants,
!! by far the most important are those that define the 'kind' types of
!! virtually all the variables used in a (well-written) FORTRAN 90 code
!! the content of this file is derived from 'Numerical Recipes in Fortran 90'
!! W.H. Press et al., volume 2 of 'Fortran Numerical Recipes', Cambridge
!! University Press, Second Edition (1996), p. 937 and 1361
!!
!! TODO
!! Unit numbers from abinit.f should be transferred to here.
!!
!! SOURCE

 module basis_defs

 implicit none

!Keyword 'integer' stands for default integer type
!and may be used whenever integer are presumed to be small

!nb of bytes related to an integer subtype n such as -10^9 < n < 10^9
 integer, parameter :: i4b=selected_int_kind(9)

!Idem for smaller integer subtypes
 integer, parameter :: i2b=selected_int_kind(4)
 integer, parameter :: i1b=selected_int_kind(2)

!nb of bytes related to default simple-precision real/complex subtypes
!(= 4 for many machine architectures, = 8 for Cray T3E for instance)
!integer, parameter :: sp=kind(1.0)          ! Single precision should not be used
!integer, parameter :: spc=kind((1.0,1.0))

!nb of bytes related to default double-precision real/complex subtypes
!(= 8 for many machine architectures)
 integer, parameter :: dp=kind(1.0d0)
 integer, parameter :: dpc=kind((1.0d0,1.0d0))  ! Complex should not be used presently
                                                ! except for use of libraries

!Example:
! integer, parameter :: urp=selected_real_kind((p=)12,(r=)50)
! real((kind=)urp) :: d
! d=5.04876_urp   ! for a real d with 12 significative digits
! and such as 10^-50 < |d| < 10^50

!To modify sp/spc and / or dp/dpc, insert instructions such as 'dp='
! but do not modify the other declarations in this module

!Default logical type
 integer, parameter :: lgt=kind(.true.)

!Some constants:
 integer, parameter :: integer_not_used=0
 logical, parameter :: logical_not_used=.true.

!UNIX standard input, standard output and error output units
 integer, parameter :: std_in_un=5
 integer, parameter :: std_out_un=6
 integer, parameter :: err_out_un=0

!Real constants
 real(dp), parameter :: zero=0._dp
 real(dp), parameter :: one=1._dp
 real(dp), parameter :: two=2._dp
 real(dp), parameter :: three=3._dp
 real(dp), parameter :: four=4._dp
 real(dp), parameter :: five=5._dp
 real(dp), parameter :: six=6._dp
 real(dp), parameter :: seven=7._dp
 real(dp), parameter :: eight=8._dp
 real(dp), parameter :: nine=9._dp
 real(dp), parameter :: ten=10._dp

!Fractionary real constants
 real(dp), parameter :: quarter=0.25_dp
 real(dp), parameter :: third=0.3333333333333333333333333333333333333333333_dp
 real(dp), parameter :: half=0.50_dp
 real(dp), parameter :: two_thirds=0.66666666666666666666666666666666666666_dp
 real(dp), parameter :: three_quarters=0.75_dp
 real(dp), parameter :: four_thirds=1.3333333333333333333333333333333333333_dp
 real(dp), parameter :: five_thirds=1.6666666666666666666666666666666666666_dp

!Real constants derived from pi
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: two_pi=6.283185307179586476925286766559012_dp
 real(dp), parameter :: four_pi=12.566370614359172953850573533118024_dp
 real(dp), parameter :: rad_to_deg=57.29577951308232087679815481410494_dp ! (= 180./pi)
 real(dp), parameter :: deg_to_rad=0.01745329251994329576923690768488609_dp ! (= pi/180.)
!The following are not used
!real(dp), parameter :: half_pi=1.570796326794896619231321691639753_dp
!real(dp), parameter :: third_pi=1.047197551196597746154214461093169_dp
!real(dp), parameter :: quarter_pi=7.853981633974483096156608458198765_dp
!real(dp), parameter :: two_thirds_pi=2.094395102393195492308428922186337_dp


!Real precision
 real(dp), parameter :: smallest_positive_real = epsilon(one)
 real(dp), parameter :: greatest_real = huge(one)
 real(dp), parameter :: smallest_real = -greatest_real

!Real physical constants 
!Revised fundamental constants from Physics Today August 1989 p.8.
!(from 1986 least squares adjustment)
 real(dp), parameter :: Bohr_Ang=0.529177249_dp    ! 1 Bohr, in Angstrom
 real(dp), parameter :: Ha_cmm1=219474.6306726_dp  ! 1 Hartree, in cm^-1 
 real(dp), parameter :: Ha_eV=27.2113961_dp ! 1 Hartree, in eV
 real(dp), parameter :: Ha_Hz=6579.683899978_dp ! 1 Hartree, in Hz
 real(dp), parameter :: e_Cb=1.60217733d-19 ! minus the electron charge, in Coulomb
 real(dp), parameter :: kb_HaK=3.166829d-06 ! Boltzmann constant in Ha/K
 real(dp), parameter :: amu_emass=1822.88851_dp ! 1 atomic mass unit, in electronic mass
!This value is 1Ha/bohr^3 in J/m^3 : used 1Ha=27.2113961eV ;
!1eV=1.60217733d-19J ; 1Bohr=0.529177249d-10m.
 real(dp), parameter :: HaBohr3_GPa=29421.033_dp ! 1 Ha/Bohr^3, in GPa
 real(dp), parameter :: Avogadro=6.0221367d23 ! per mole

!Character constants
 character*1, parameter :: ch10 = '\n'

 end module basis_defs
!!***
