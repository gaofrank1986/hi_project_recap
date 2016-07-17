program tes
    use io
    type(Ostream) :: test
    !real(4) a
    integer :: a
    !call test%initialize("test",6)

    a=2
    !print *,test%title
    test=Ostream("test",6)
    call test%fout(test%toString(a,'(i9)'))
    !call test<<'2'
    !call test%fout(test%toString(a,'(f6.2)'),'2')
    !print *,test%toString(a)
    
    !call tos('string')
    !character(len=:),allocatable :: string
    !string ="A"
    !print *,string
contains
    subroutine tos(this)
        character(len=*) :: this
        print *,this
    end subroutine
end program
