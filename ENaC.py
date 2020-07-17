import numpy as np
from values import *
from defs import *

def ENaC(cell,i,memb_id,hNaMP,area,jvol):
	
	facNaMP=(30.0/(30.0+cell.conc[0,0]))*(50.0/(50.0+cell.conc[0,1]))
	if cell.segment =='DCT':
		LzD=i-1
		x1=LzD+1
		xn=N
		xdct2=2.0/3.0

		if x1/xn<xdct2:
			#NCCexp=1.0
			ENaCexp=0.0
		else:
			#NCCexp=2*(1-(x1/xn-xdct2)/(1-xdct2))
			ENaCexp=(x1/xn-xdct2)/(1-xdct2)
		facphMP=1.0
		if cell.sex == 'female':
			hENaC=ENaCexp*hNaMP*facNaMP*facphMP*2.0
		else:
			hENaC=ENaCexp*hNaMP*facNaMP*facphMP
	elif cell.segment == 'CNT':
		facCaMP=1.0
		facphMP=1.0
		NaMPq0=cell.vol_init[0]-(2.0e-6)/60/Vref
		facFvMP=max(0.01,1+3*((cell.vol[0]/NaMPq0)-1))
		if cell.sex == 'female':
			hENaC=hNaMP*facNaMP*facCaMP*facphMP*facFvMP*1.3
		else:
			hENaC=hNaMP*facNaMP*facCaMP*facphMP*facFvMP
	elif cell.segment == 'CCD':
		facphMP=1.0
		NaMPq0=cell.vol_init[0]-(0.1e-6)/60/Vref
		facFvMP=max(0.01,1+3*((cell.vol[0]/NaMPq0)-1))
		if cell.sex == 'female':
			hENaC=hNaMP*facNaMP*facphMP*facFvMP*1.5
		else:
			hENaC=hNaMP*facNaMP*facphMP*facFvMP
	elif cell.segment == 'OMCD':
		facphMP=1.0
		if cell.sex == 'female':
			hENaC=hNaMP*facNaMP*facphMP*1.2
		else:
			hENaC=hNaMP*facNaMP*facphMP

	XI=zval[0]*F*EPref/RT*(cell.ep[0]-cell.ep[1])
	dint=np.exp(-XI)
	if (abs(1-dint)<1e-6):
		dJENaC=cell.area[0,1]*hENaC*(cell.conc[0,0]-cell.conc[0,1])
	else:
		dJENaC=cell.area[0,1]*hENaC*XI*(cell.conc[0,0]-cell.conc[0,1]*dint)/(1-dint)

	if cell.segment=='CNT' or cell.segment == 'CCD':
		concdiff=cell.conc[0,0]-cell.conc[0,1]
		if abs(concdiff)>1e-6:
			concmean=(cell.conc[0,0]-cell.conc[0,1])/np.log(abs(cell.conc[0,0]/cell.conc[0,1]))
			dimless=(Pfref*Vwbar*Cref)/href
			convect=(1-cell.sig[0,0,1])*concmean*jvol[0,1]*dimless
			dJENaC=dJENaC+convect

	return [0],[dJENaC]
