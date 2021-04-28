c THIS USES STOKE TO CONTROL THE TIME STEP:
c STOKE variation:

      if (m_dot.eq.0.d0) then
       ichs=0
       mdstoke=0.d0
       do i=1,imax-2,2
        mds=abs(stoke(i)-nstoke(i))/(abs(nstoke(i))+1.d-5)
        if (mds.gt.mdstoke)then
         mdstoke=mds
         ichs=i
        end if
       end do
       chstoke=0.d0
       if((svar-1.d0).lt.mdstoke)then
        dt=dt*(svar-1.d0)/mdstoke
        if (time.lt.1.d5) then
         dt=min(dt,dt1)
        else
         dt=min(dt,dt0)
        end if
        chstoke=1.d0
       end if
      else
       chstoke=0.d0
      end if
