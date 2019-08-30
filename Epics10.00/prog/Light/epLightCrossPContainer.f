      subroutine epLightCrossPContainer(cnx, icon)
!       crossing with parent container whic partially contain the current comp.

!       we don't know what is outside of cnx; so canoot know the refraction index
!       there.
!       This case will not happen since  partial conatinment has not been used so
!       far. so propaget light until crossed point


      integer,intent(inout):: cnx ! outside is this comp. #. 
                         ! if reflection happens, will become Cn
      integer,intent(out)::icon
   
      write(0,*) ' partially contained component to the paretent' 
      write(0,*) ' not well defined for light transpprt '
      icon =0
      end
