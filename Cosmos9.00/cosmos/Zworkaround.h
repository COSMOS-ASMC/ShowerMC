!        This is used when next 3 items are not recognized
!      by ld when they are in module.
!      (happens for cheavyInt.f and chAcol.f)
!      values are transferred here via cworkaround in cGetXsec.f
!
      integer::TargetNucleonNo, TargetProtonNo, colElemNo
      real(8)::TargetXs
      common /workaround/TargetXs, TargetNucleonNo, TargetProtonNo,
     *          colElemNo

