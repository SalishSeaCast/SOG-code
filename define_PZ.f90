SUBROUTINE define_PZ(M, PZ_bins, D_bins, M2, Pmicro, Pnano, NO, NH, Detritus, &
                      PZ)

! This subroutine takes all the separate variables (microplankton,
! nanoplankton, nitrate, ammonium and detritus and loads them
! sequentially into the PZ vector for the ODE solver to use.

  USE mean_param, only : bins, snow 


  IMPLICIT NONE

  INTEGER, INTENT (IN) :: M
  TYPE (bins), INTENT (IN) :: PZ_bins
  INTEGER, INTENT (IN) :: D_bins
  INTEGER, INTENT (IN) :: M2
  DOUBLE PRECISION, DIMENSION(0:M+1), INTENT(IN)::Pmicro, Pnano, NO, NH
  TYPE(snow), DIMENSION(D_bins), INTENT(IN)::Detritus
  DOUBLE PRECISION, DIMENSION (M2), INTENT (OUT) :: PZ

! Local Variables

  INTEGER :: bPZ, ePZ ! start position and end position in PZ array
  INTEGER :: j

  PZ = 0.

! Microplankton

  bPz = (PZ_bins%micro-1) * M + 1
  ePZ = PZ_bins%micro * M
  PZ(bPZ:ePZ) = Pmicro(1:M)

! Nanoplankton

  bPz = (PZ_bins%nano-1) * M + 1
  ePZ = PZ_bins%nano * M
  PZ(bPZ:ePZ) = Pnano(1:M)

! Nitrate

  bPz = (PZ_bins%NO-1) * M + 1
  ePZ = PZ_bins%NO * M
  PZ(bPZ:ePZ) = NO(1:M)

! Ammonium

  bPz = (PZ_bins%NH-1) * M + 1
  ePZ = PZ_bins%NH * M
  PZ(bPZ:ePZ) = NH(1:M)

! Detritus
  do j=1,D_bins
     bPz = (PZ_bins%Quant + (j-1) ) * M + 1
     ePz = (PZ_bins%Quant + j) * M
     PZ(bPz:ePz) = Detritus(j)%D%new(1:M)
  enddo


END SUBROUTINE define_PZ

