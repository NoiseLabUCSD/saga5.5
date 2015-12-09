      subroutine mname(chstr,leng)
c
c: Subroutine to determine the length of a character
c: string without blanks
c: chstr -- the string (length 64)
c: leng -- the length without blanks
c
      implicit integer*4(i-n)
      character*64 chstr
      character*1 char1
c
      do 40 i=64,1,-1
         char1 = chstr(i:i)
       	 if(char1 .ne. ' ') then
            leng = i
            goto 45
       	 endif
40    continue
45    continue
c
      return
      end
