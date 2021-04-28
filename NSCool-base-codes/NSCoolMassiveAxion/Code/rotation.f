      subroutine spin_down(p1,dtime,beta_rot,bfield,p2)
c CHECKED on Nov.5, 1998
       implicit real*8(a-h,k-z)
        p22=p1**2+beta_rot*bfield**2*dtime
        p2=sqrt(p22)
       return
      end
c *******************************************************
c *******************************************************
      subroutine get_beta(mass,radius,beta_rot,m_i)
c CHECKED on Nov.5, 1998
       implicit real*8(a-h,k-z)
       parameter(pi=3.14159)
       parameter(g=6.67e-8,c=2.99e10,msun=2.e33)
        lambda=1./(1.-2.*g*mass*msun/radius/c**2)
        m_i= 0.21 *mass*msun*radius**2 * lambda
        beta_rot= 16.*pi**2/3./c**3*radius**6/m_i
       return
      end
c *******************************************************
c *******************************************************
