module io
    implicit none
    public

    type Ostream
        private
        integer port;
        character(len=:),allocatable :: header
    contains
        procedure :: initialize
        procedure ::  fout
        procedure,private:: toString_int
        !procedure,private :: toString_real8
        !procedure,private :: toString_real4
        procedure,private :: toString_real
        generic :: toString =>toString_int,toString_real
        !generic :: toString =>toString_int,toString_real8,toString_real4,toString_realk
        !generic :: operator(<) => fout

    end type

    interface Ostream 
        procedure :: constructor
    end interface

contains
    function constructor(string,port)
        type(Ostream) ::constructor
        character(len=*),intent(in)::string
        integer,intent(in) :: port
        call constructor%initialize(string,port)
    end function


    subroutine initialize(this,string,port)
        class (Ostream) ::  this
        character(len=*),intent(in)::string
        integer,intent(in) :: port
        this%header='['//string//']'
        this%port=port
    end subroutine initialize

    subroutine fout(this,string)
    !subroutine fout(this,string,ctl)
        class (Ostream),intent(in) ::  this
        character(len=*),intent(in)::string
        !character(len=*),optional::ctl
        write(this%port,fmt='(a)') this%header//' '//string
    end subroutine
   
    ! this need to be polymorphism for integer/logical/double/complex etc.
    function toString_int(this,arg) result(string)
        class (Ostream),intent(in) ::  this
        integer,intent(in)::arg
        character(:),allocatable :: string
        character(len=1024) :: tmp
        write(tmp,fmt='(i4)') arg 
        string=trim(tmp)
    end function

    !function toString_real8(this,arg,ctl) result(string)
        !class (Ostream) ::  this
        !character(:),allocatable :: string
        !character(*),optional :: ctl
        !character(len=1024) :: tmp
        !real(8)::arg
        
        !if(present(ctl)) then
            !write(tmp,fmt=ctl) arg
        !else
            !write(tmp,fmt='(f10.5)') arg 
        !end if
        !string=trim(tmp)
    !end function

    !function toString_real4(this,arg,ctl) result(string)
        !class (Ostream) ::  this
        !character(:),allocatable :: string
        !character(*),optional :: ctl
        !character(len=1024) :: tmp
        !real(4)::arg
        
        !if(present(ctl)) then
            !write(tmp,fmt=ctl) arg
        !else
            !write(tmp,fmt='(f8.4)') arg 
        !end if
        !string=trim(tmp)
    !end function


    function toString_real(this,arg,ctl) result(string)
        class (Ostream) ::  this
        character(:),allocatable :: string
        character(*),optional :: ctl
        character(len=1024) :: tmp
        real::arg
        
        if(present(ctl)) then
            write(tmp,fmt=ctl) arg
        else
            write(tmp,fmt='(f10.5)') arg 
        end if
        string=trim(tmp)
    end function
end module io

