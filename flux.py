from defs import *
from values import *
import electrochemical 
import water
import glucose
import cotransport
import NHE3
import ATPase
import NKCC
import KCC
import NCC
import ENaC
import Pendrin
import AE1
import NHE1
import NaKCl2
import math

def compute_fluxes (cell,j):

	# update LIS-bath surface area, based on LIS volume
	cell.area[4][5] = 0.02*max(cell.vol[4]/cell.volref[4],1.0)
	cell.area[5][4] = cell.area[4][5]  
	
	# compute water fluxes
	jvol = water.compute_water_fluxes(cell)
	# if cell.segment =='PT':
	# 	print('Water Fluxes:')
	# 	print(jvol)
	# 	pause = input()
	
	# compute electrochemical convective and diffusive fluxes
	jsol,delmu = electrochemical.compute_ecd_fluxes(cell,jvol)
	
	for i in range(len(cell.trans)):
		
		transporter_type = cell.trans[i].type
		memb_id = cell.trans[i].membrane_id

		if transporter_type == 'SGLT1':
			solute_id,fluxsglt1=glucose.sglt1(cell,cell.ep,memb_id,cell.trans[i].act,cell.area)
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] += fluxsglt1[i]
#            print('SGLT1 returns',solute_id,fluxsglt1)

		elif transporter_type == 'SGLT2':
			solute_id,fluxsglt2=glucose.sglt2(cell,cell.ep,memb_id,cell.trans[i].act,cell.area)        
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] += fluxsglt2[i]
#            print('SGLT2 returns',solute_id,fluxsglt2)

		elif transporter_type == 'GLUT1':
			solute_id,fluxglut1=glucose.glut1(cell,cell.ep,memb_id,cell.trans[i].act,cell.area)           
			jsol[solute_id][memb_id[0]][memb_id[1]] += fluxglut1
#            print('GLUT1 returns',fluxglut1)

		elif transporter_type == 'GLUT2':
			solute_id,fluxglut2=glucose.glut2(cell,cell.ep,memb_id,cell.trans[i].act,cell.area)           
			jsol[solute_id][memb_id[0]][memb_id[1]] += fluxglut2
#            print('GLUT2 returns',fluxglut2,solute_id)
			
		elif transporter_type == 'NHE3':
			solute_id,fluxnhe3=NHE3.nhe3(cell,cell.ep,memb_id,cell.trans[i].act,cell.area)           
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] += fluxnhe3[i]             
#            print('NHE3 returns',fluxnhe3)

		elif transporter_type == 'NaKATPase':
			solute_id,fluxnakatpase=ATPase.nakatpase(cell,cell.ep,memb_id,cell.trans[i].act,cell.area)           
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] += fluxnakatpase[i]           
#            print('Na-K-ATPase returns',fluxnakatpase)

		elif transporter_type == 'HATPase':
			solute_id,fluxhatpase=ATPase.hatpase(cell,cell.ep,memb_id,cell.trans[i].act,cell.area)           
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] += fluxhatpase[i]     
#            print('H-ATPase returns',fluxhatpase,solute_id)
		elif transporter_type == 'NKCC2A':
			solute_id,fluxnkcc2a=NKCC.nkcc2(cell,memb_id,cell.trans[i].act,cell.area,'A')
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] += fluxnkcc2a[i]
		elif transporter_type == 'NKCC2B':
			solute_id,fluxnkcc2b=NKCC.nkcc2(cell,memb_id,cell.trans[i].act,cell.area,'B')
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] += fluxnkcc2b[i]
		elif transporter_type == 'NKCC2F':
			solute_id,fluxnkcc2f=NKCC.nkcc2(cell,memb_id,cell.trans[i].act,cell.area,'F')
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] += fluxnkcc2f[i]        
		elif transporter_type == 'KCC4':
			solute_id,fluxkcc4=KCC.kcc4(cell.conc,memb_id,cell.trans[i].act,cell.area)
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] += fluxkcc4[i]
		elif transporter_type == 'ENaC':
			solute_id,fluxENaC=ENaC.ENaC(cell,j,memb_id,cell.trans[i].act,cell.area,jvol)
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] +=fluxENaC[i]
		elif transporter_type == 'NCC':
			solute_id,fluxncc=NCC.NCC(cell,j,memb_id,cell.trans[i].act,cell.area)
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] +=fluxncc[i]
		elif transporter_type == 'Pendrin':
			solute_id,fluxPendrin=Pendrin.Pendrin(cell,memb_id,cell.trans[i].act,cell.area)
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] +=fluxPendrin[i]
		elif transporter_type =='AE1':
			solute_id,fluxAE1=AE1.AE1(cell,memb_id,cell.trans[i].act,cell.area)
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] +=fluxAE1[i]
		elif transporter_type == 'HKATPase':
			solute_id,fluxhkatpase=ATPase.hkatpase(cell,memb_id,cell.trans[i].act,cell.area)
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] +=fluxhkatpase[i]
		elif transporter_type == 'NHE1':
			solute_id,fluxnhe1=NHE1.NHE1(cell,memb_id,cell.trans[i].act,cell.area)
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] +=fluxnhe1[i]
		elif transporter_type == 'NKCC1':
			solute_id,fluxnakcl2=NKCC.nkcc1(cell,memb_id,cell.trans[i].act,delmu)
			for i in range(len(solute_id)):
				jsol[solute_id[i]][memb_id[0]][memb_id[1]] +=fluxnakcl2[i]
		else:
			print('What is this?',transporter_type)
			
	jsol = cotransport.compute_cotransport(cell,delmu,jsol)

	
	if cell.segment=='PT' or cell.segment == 'S3':

		if cell.segment == 'PT':
			TS = 1.3
			scaleT = 1.0
		elif cell.segment == 'S3':
			TS = 1.3
			scaleT = 0.5

		#torque-modulated effects
		
		PM=cell.pres[0]
	
		Radref = 0.0025/2.0e0
		fac1 = 8.0*visc*(cell.vol_init[0]*Vref)*torqL/(Radref**2)
		fac2 = 1.0 + (torqL+torqd)/Radref + 0.50*((torqL/Radref)**2)
		TM0= fac1*fac2
	
		RMtorq = torqR*(1.0e0+torqvm*(PM - PbloodPT))
		factor1 = 8.0*visc*(cell.vol[0]*Vref)*torqL/(RMtorq**2)
		factor2 = 1.0 + (torqL+torqd)/RMtorq + 0.50*((torqL/RMtorq)**2)
		Torque = factor1*factor2
	
		Scaletorq = 1.0 + TS*scaleT*(Torque/TM0-1.0)
		#print(TS,scaleT,Torque,TM0)
		for i in range(NS):
			jsol[i][0][1]=Scaletorq*jsol[i][0][1]
			jsol[i][1][4]=Scaletorq*jsol[i][1][4]
			jsol[i][1][5]=Scaletorq*jsol[i][1][5]
   
		jvol[0][1]=Scaletorq*jvol[0][1]
		jvol[1][4]=Scaletorq*jvol[1][4]
		jvol[1][5]=Scaletorq*jvol[1][5]
	# if cell.segment == 'LDL':
	# 	print(jsol)
	# 	input('pause')
		
		# print(TS,scaleT,Torque,TM0)
		# print(visc,cell.vol[0],Vref,torqL,RMtorq)
		# print(torqL,torqd,RMtorq)
		# print(torqR,torqvm,PM,PbloodPT)
		# print(factor1,factor2)
		
		# print(fac1,fac2)
		# print(visc,cell.vol_init[0],Vref,torqL,Radref)
		# print(torqL,torqd,Radref)
		#print(Scaletorq)
	# if cell.segment == 'SDL':
	# 	jvol = water.compute_water_fluxes(cell)
	# 	jsol = np.zeros([NS,NC,NC])

	# 	pk = np.zeros(10)
	# 	cint = np.zeros(10)
	# 	cext = np.zeros(10)
	# 	jk = np.zeros(10)
	# 	factor = np.zeros(10)

	# 	pk[0:5] = cell.area[0,1]*cell.h[0:5,0,1]+cell.area[0,4]*cell.h[0:5,0,4]
	# 	pk[6:7] = cell.area[0,1]*cell.h[8:9,0,1]+cell.area[0,4]*cell.h[8:9,0,4]
	# 	pk[8:9] = cell.area[0,1]*cell.h[10:11,0,1]+cell.area[0,4]*cell.h[10:11,0,4]

	# 	cint[0:5] = cell.conc[0:5,0]
	# 	cint[6:7] = cell.conc[8:9,0]
	# 	cint[8:9] = cell.conc[10:11,0]

	# 	cext[0:5] = cell.conc[0:5,5]
	# 	cext[6:7] = cell.conc[8:9,5]
	# 	cext[8:9] = cell.conc[10:11,5]

	# 	z = np.array([1,1,-1,-1,0,0,0,0,1,1])

	# 	factor = -1*z*F*EPref/RT
	# 	jk[0:3] = pk[0:3]*factor[0:3]*((cint[0:3]-cext[0:3]*np.exp(-1*factor[0:3]))/(1-np.exp(-1*factor[0:3])))
	# 	jk[4:7] = pk[4:7]*(cint[4:7]-cext[4:7])
	# 	jk[8:9] = pk[8:9]*factor[8:9]*((cint[8:9]-cext[8:9]*np.exp(-1*factor[8:9]))/(1-np.exp(-1*factor[8:9])))

	# 	jsol[0:5,0,5] = jk[0:5]
	# 	jsol[8:9,0,5] = jk[6:7]
	# 	jsol[10:11,0,5] = jk[8:9]
	return jvol,jsol

	
