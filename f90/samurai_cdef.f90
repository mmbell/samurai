  ! C functions declaration
  ! Include this file in your f90 code so that the compiler knows about the function signatures
  
  interface

     function create_vardriver3d_c(config) bind(C, name="create_vardriver3D")
       use iso_c_binding
       use samModule
       implicit none
       type(samurai_config) :: config
       type(c_ptr) :: create_vardriver3d_c
     end function create_vardriver3d_c

     function create_vardriver3d_from_xml_c(xmlfile) bind(C, name="create_vardriver3D_From_File")
       use iso_c_binding
       implicit none
       character(kind=c_char) :: xmlfile(*)
       type(c_ptr) :: create_vardriver3d_from_xml_c
     end function create_vardriver3d_from_xml_c

     subroutine delete_vardriver3d_c(driver) bind(C, name="delete_vardriver3D")
       use iso_c_binding
       implicit none
       type(c_ptr), value :: driver
     end subroutine delete_vardriver3d_c

     subroutine run_vardriver3d_c( &
          driver, &
          nx, ny, nsigma, &
          dx, dy, &
          sigmas, latitude, longitude, &
          u1,v1, w1, th1, p1, &
          usam, vsam, wsam, thsam, psam ) bind(C, name="run_vardriver3D")
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: driver
       integer(c_int), value, intent(in)  :: nx, ny, nsigma
       real(c_float),  value, intent(in)  :: dx, dy
       
       real(c_float),  intent(in)  :: sigmas(nsigma)
       real(c_float),  intent(in)  :: latitude(nx, ny), longitude(nx, ny)
       
       real(c_float),  intent(in)  :: u1(nx, ny, nsigma), v1(nx, ny, nsigma)
       real(c_float),  intent(in)  :: w1(nx, ny, nsigma)
       real(c_float),  intent(in)  :: th1(nx, ny, nsigma), p1(nx, ny, nsigma)
       
       real(c_float),  intent(out) :: usam(nx, ny, nsigma), vsam(nx, ny, nsigma)
       real(c_float),  intent(out) :: wsam(nx, ny, nsigma)	
       real(c_float),  intent(out) :: thsam(nx, ny, nsigma), psam(nx, ny, nsigma)

     end subroutine run_vardriver3d_c

  end interface
