      subroutine aqwait(aqp, istatu)

c 30/3/2004:
cJAC	use dflib
	      
      integer aqp, istatu, iuni

c 30/3/2004: Compaq routine to complete buffer output:
cJAC	if(.not. commitqq(aqp))istatu = 1
	!!!if(.not. int(FLUSH(aqp)))istatu = 1      
        !!!if(int(FLUSH(aqp)).ne.1)istatu=1
        call FLUSH(aqp)
        istatu=1 
      return
      end
