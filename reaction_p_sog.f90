SUBROUTINE reaction_p_sog(M, PZ_bins, D_bins, PZ_quant, PZ_det, Pmicro, Pnano, NO, &
     NH, Detritus, Gvector)

      USE mean_param

      IMPLICIT NONE

      INTEGER, INTENT(IN)::M
      TYPE (bins), INTENT(IN) :: PZ_bins
      INTEGER, INTENT(IN)::D_bins  !
      DOUBLE PRECISION, DIMENSION(PZ_bins%Quant*M), INTENT(IN)::PZ_quant
      DOUBLE PRECISION, DIMENSION(D_bins*M),INTENT(IN)::PZ_det
      DOUBLE PRECISION, DIMENSION(0:M+1), INTENT(IN)::Pmicro,Pnano, NO, NH
      TYPE(snow), DIMENSION(D_bins), INTENT(IN)::Detritus
      TYPE(UVST), INTENT(OUT)::Gvector
! Local variables
      INTEGER :: bPZ, ePZ ! start position and end position in PZ array
      INTEGER::j

      ! Micro plankton
      bPZ = (PZ_bins%micro-1) * M + 1
      ePZ = PZ_bins%micro * M
      Gvector%p%micro = PZ_quant(bPZ:ePZ) - Pmicro(1:M)

      ! Nano plankton
      bPZ = (PZ_bins%nano-1) * M + 1
      ePZ = PZ_bins%nano * M
      Gvector%p%nano = PZ_quant(bPZ:ePZ) - Pnano(1:M)

      ! Nitrate
      bPZ = (PZ_bins%NO-1) * M + 1
      ePZ = PZ_bins%NO * M
      Gvector%N%O = PZ_quant(bPZ:ePZ) - NO(1:M)

      ! Ammonimum
      bPZ = (PZ_bins%NH-1) * M + 1
      ePZ = PZ_bins%NH * M
      Gvector%N%H = PZ_quant(bPZ:ePZ) - NH(1:M)

      ! Detritus
      DO j = 1,D_bins
         Gvector%d(j)%bin = PZ_det((j-1)*M+1:j*M) - Detritus(j)%D%old(1:M)
      END DO

END SUBROUTINE reaction_p_sog
