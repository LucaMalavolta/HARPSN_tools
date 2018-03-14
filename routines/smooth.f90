module Smooth
  use common
  implicit none
  
contains
  subroutine smooth_boxcar(input, width, output)
    real (PR), dimension(:), intent(in) :: input
    integer, intent(in) :: width
    real (PR), dimension(:), intent(out) :: output
    
    integer :: i,j,i1, i2, hwidth
    
    output = input
    i1 =  int((width-1)/2) + 1 
    i2 =  size(input,dim=1) - int((width-1)/2) - 1 
    do i = i1, i2 
       hwidth = width/2
       if (hwidth.le.0) hwidth = 1
       output(i) = sum(input( i-hwidth:i+hwidth-1),dim=1) /(2._PR*hwidth)
    end do
    return
  end subroutine smooth_boxcar
end module Smooth
    
