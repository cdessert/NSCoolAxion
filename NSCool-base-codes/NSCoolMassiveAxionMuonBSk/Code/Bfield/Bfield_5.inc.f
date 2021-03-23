      if (ifield.eq.2) then
       do i=1,imax,2
        dels=abs(nstoke(i)-onstoke(i))
        rats=dels/(abs(nstoke(i))+1.d-100)
        if(rats.gt.ratios)then
         ratios=rats
         irats=i
        end if
       end do
      end if

      do i=1,imax,2
       onstoke(i)=nstoke(i)
      end do
