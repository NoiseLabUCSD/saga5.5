      logical function blanks(msg)
      integer l1
      character*(*) msg

c**   Returns true if the input string is composed of all blanks

      blanks = .false.
      l1     = lenstr(msg)
      ifs    = index(msg,char(32))

c**   First handle leading spaces from option string.

      i = 1
      if (ifs.eq.1) then
 5       i=i+1
         if (ichar(msg(i:i)).eq.32) goto 5
         ifs = index(msg(i:l1),char(32))+i-1
      end if

c**   Next check for empty message string
      
      if (ifs.eq.l1) blanks = .true.
      return
      end
