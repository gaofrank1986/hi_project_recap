program main
        use green_funcs,only : mirror,DGFunc
        implicit none
        real(8) :: a(3),b(3)
        real(8) :: h
        a = (/0.d0,0.d0,-1.d0/)
      h = 10.d0
      print *,h
      print *,a
      print *,DGFunc(b,a)
      print *,DGFunc(b,a)
      print *,DGFunc(b,a)
      b=DGFunc(b,a)
      print *,b
      end program
