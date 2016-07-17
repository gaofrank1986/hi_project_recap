module io
    implicit none
    public

    type Ostream
        private
        integer port;
        character(len=:),allocatable :: header
    contains
        procedure,private :: initialize
        procedure ::  fout
        procedure,private:: toString_int
        !procedure,private :: toString_real8
        !procedure,private :: toString_real4
        procedure,private :: toString_real
        generic :: toString =>toString_int,toString_real
        !generic :: toString =>toString_int,toString_real8,toString_real4,toString_realk
        !generic :: operator(<) => fout

    end type
    ! @ interface : this interface bound constructor with type Ostream to mimic a constructor behavior
    interface Ostream 
        procedure :: constructor
    end interface

contains
    ! @func : same as C++ constructor
    function constructor(string,port)
        type(Ostream) ::constructor
        character(len=*),intent(in)::string
        integer,intent(in) :: port
        call constructor%initialize(string,port)
    end function

    ! @func : used in constructor to intialized object 
    ! @param: [string] header appear at beginning
    ! @param: [port] which termianl to write output, 6 for stdout, other for file
    subroutine initialize(this,string,port)
        class (Ostream) ::  this
        character(len=*),intent(in)::string
        integer,intent(in) :: port
        this%header='['//string//']'
        this%port=port
    end subroutine initialize

    ! @func : output given string to designated port
    subroutine fout(this,string)
        class (Ostream),intent(in) ::  this
        character(len=*),intent(in)::string
        !character(len=*),optional::ctl
        write(this%port,fmt='(a)') this%header//' '//string
    end subroutine
   
    ! @func : [toString_int] convert a integer into string 
    ! @ param: ctl optional, format of integer output
    ! todo this need to be polymorphism for integer/logical/double/complex etc.
    function toString_int(this,arg,ctl) result(string)
        class (Ostream),intent(in) ::  this
        character(*),optional :: ctl
        integer,intent(in)::arg
        character(:),allocatable :: string
        character(len=1024) :: tmp
        if(present(ctl)) then
            write(tmp,fmt=ctl) arg
        else
            write(tmp,fmt='(i4)') arg ! output to string
        end if
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

    ! @func : [toString_real] convert a given real number to a string
    ! @param: [ctl] optional, control the format
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

