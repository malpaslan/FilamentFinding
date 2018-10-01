c
       subroutine cross(x1,n1,x2,n2,binmin,binmax,binwid,binout,
     + binn,weights1,weights2,nd)

       integer nd,n1,n2,binloc
       integer i,j, ic, binn
       
       real*8 x1(n1,nd), x2(n2,nd), dtemp, binout(binn)
       real*8 binmin,binmax,binwid,weights1(n1),weights2(n2),binmax2
       binmax2=binmax**2
            do  15 i= 1, n1
               do 10 j =1,n2

c** accumulate squared differences
               dtemp= 0.0
               do 5 ic= 1, nd 
                    dtemp= dtemp + (x1(i,ic) - x2(j,ic))**2
                    if( dtemp.gt.binmax2) goto 10
 5             continue
c squareroot to get distance 
            dtemp=sqrt(dtemp)
c check this will fit in the histogram
            if((dtemp.ge.binmin).and.(dtemp.le.binmax)) then
c if it fits, find where to put it
                binloc=int((dtemp-binmin)/binwid)+1
c add to the correct hist element the product of both weights
                binout(binloc)=binout(binloc)+(weights1(i)*weights2(j))
            endif

10             continue
15          continue

       return
       end

