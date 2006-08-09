	subroutine getpar_init(flag)
	implicit	none
	integer		flag

	character	fdate*24

	include		"gpcom.h"

	do_rep = flag

	if (do_rep==1) &
     		write(*,'(a50," = ",a)') "Date", fdate()

	return
      end subroutine getpar_init

      logical function getparl(name,flag)
	character*(*)	name
	integer		flag

	include		"gpcom.h"

	logical		enable
	character	label*25, desc*50, str*1

	read(5,*) label, enable, desc
	if (name .ne. label) then
		print *, "GETPARI: Expecting variable ", name, &
     			" but got ",trim(label), " instead."
		stop
	end if

	if ((do_rep .ne. 0) .and. (flag .ne. 0)) then
		write(str,'(L1)') enable
		str = AdjustL(str)
		write(*,'(a50," = ",a1)') trim(desc), trim(str)
	endif

	getparl = enable

	return
      end function getparl


	double precision function getpard(name,flag)
	character*(*)	name
	integer		flag

	include		"gpcom.h"

	double precision		x
	character	label*25, desc*50, str*16

	read(5,*) label, x, desc
	if (name .ne. label) then
		print *, "GETPARD: Expecting variable ", name, &
     			" but got ", trim(label), " instead."
		stop
	end if

	if ((do_rep .ne. 0) .and. (flag .ne. 0)) then
		write(str,'(1PG16.8)') x
		str = AdjustL(str)
		write(*,'(a50," = ",a16)') trim(desc), str
	endif

	getpard = x

	return
	end


	subroutine getpariv(name,v,n,flag)
	implicit	none
	character*(*)	name
	integer		v(*), n, flag


	include		"gpcom.h"

	integer		i
	character	label*25, desc*50, str*16

	read(5,*) label, (v(i), i=1,n), desc
	if (name .ne. label) then
		print *, "GETPARIV: Expecting variable ", name, &
     			" but got ",trim(label), " instead."
		stop
	end if

	if ((do_rep .ne. 0) .and. (flag .ne. 0)) then
		write(str,'(I16)') v(1)
		str = AdjustL(str)
		write(*,'(a46," [1] = ",a16)') trim(desc), str

		do i = 2, n
			write(str,'(I16)') v(i)
			str = AdjustL(str)
			write(*,'(46X," [",I1,"] = ",a16)') i, str
		enddo
	endif

	return
      end subroutine getpariv


	integer function getpari(name,flag)
	character*(*)	name
	integer		flag

	include		"gpcom.h"

	integer		i
	character	label*25, desc*50, str*16

	read(5,*) label, i, desc
	if (name .ne. label) then
		print *, "GETPARI: Expecting variable ", name, &
     			" but got ",trim(label), " instead."
		stop
	end if

	if ((do_rep .ne. 0) .and. (flag .ne. 0)) then
		write(str,'(I16)') i
		str = AdjustL(str)
		write(*,'(a50," = ",a16)') trim(desc), str
	endif

	getpari = i

	return
	end


	subroutine getpardv(name,v,n,flag)
	implicit	none
	double precision		v(*)
	character*(*)	name
	integer		n, flag

	include		"gpcom.h"

	integer		i
	character	label*25, desc*50, str*16

	read(5,*) label, (v(i), i=1,n), desc
	if (name .ne. label) then
		print *, "GETPARDV: Expecting variable ", name, &
     			" but got ",trim(label), " instead."
		stop
	end if

	if ((do_rep .ne. 0) .and. (flag .ne. 0)) then
		write(str,'(1PG16.8)') v(1)
		str = AdjustL(str)
		write(*,'(a46," [1] = ",a16)') trim(desc), str

		do i = 2, n
			write(str,'(1PG16.8)') v(i)
			str = AdjustL(str)
			write(*,'(46X," [",I1,"] = ",a16)') i, str
		enddo
	endif

	return
	end


	character*80 function getpars(name,flag)
	character*(*)	name
	integer		flag

	include		"gpcom.h"

	character	label*25, s*80, desc*50


	read(5,*) label, s, desc

	if (name .ne. label) then
		print *, "GETPARS: Expecting ", name, " but got ", &
     			trim(label), " instead."
		stop
	end if

	if ((do_rep .ne. 0) .and. (flag .ne. 0)) &
     		write(*,'(a50," = ",a)') trim(desc), trim(s)

	getpars = s

	return
      end function getpars
