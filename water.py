#---------------------------------------------------------------------72
#---------------------------------------------------------------------72
#    COMPUTE WATER FLUXES
#---------------------------------------------------------------------72
#---------------------------------------------------------------------72
#     The hydraulic and oncotic pressures are made non-dimensional
#     by dividing by RT*Cref

from values import *
from defs import *
import numpy as np

def compute_water_fluxes (cell):

	# need to take care of these constants
	if cell.segment == 'CNT' or cell.segment == 'CCD':
		compl = 0.3
	else:
		compl = 0.1
	# also, should define range of relevant compartments for each cell

	# local variables
	PRES = np.zeros(NC)
	ONC  = np.zeros(NC)
	jvol = np.zeros(NC*NC).reshape((NC,NC))

	# TAKE CARE OF THESE
	PM=cell.pres[0]
	if cell.segment == 'PT' or cell.segment == 'S3':
		PB = 9.0
	else:
		PB=0
	# if cell.segment == 'SDL':
	# 	PRES[0] = PM/(RTosm*Cref)
	# 	PRES[5] = 0.0
	# else:
	PRES[0] = PM/(RTosm*Cref)
	PRES[1] = PRES[0]
	PRES[2] = PRES[0]
	PRES[3] = PRES[0]
	PRES[5] = PB/(RTosm*Cref)
	PRES[4] = PRES[5]+(cell.vol[4]/cell.volref[4]-1)/compl/(RTosm*Cref)
	#print(PRES[5],cell.vol[4],cell.volref[4],compl,RTosm,Cref)
	#input()

	if cell.segment == 'SDL' or cell.segment == 'LDL' or cell.segment == 'LAL':
		PRES[1] = PRES[5]
		#PRES[4] = 0.1*PRES[4]

	# if cell.segment == 'SDL':
	# 	ONC[0] = cell.cimpref[0]
	# 	ONC[5] = cell.cimpref[5]
	# else:
	ONC[0]=cell.cimpref[0]  
	ONC[1]=cell.cimpref[1]*cell.volref[1]/cell.vol[1]
	ONC[4]=cell.cimpref[4]
	ONC[5]=cell.cimpref[5]
	if cell.segment == 'CNT' or cell.segment == 'CCD' or cell.segment == 'OMCD':
		ONC[2]=cell.cimpref[2]*cell.volref[1]/cell.vol[2]
		ONC[3]=cell.cimpref[3]*cell.volref[1]/cell.vol[3]
	#elif cell.segment == 'SDL' or cell.segment == 'LDL' or cell.segment == 'LAL':
		#ONC[1] = ONC[5]
		#ONC[4] = ONC[5]
	# if cell.segment == 'SDL':
	# 	osm = 0.0
	# 	for j in range(NS):
	# 		osm += cell.sig[j,0,5]*(cell.conc[j,0]-cell.conc[j,5])
	# 	jvol[0,5] = (cell.area[0,1]*cell.dLPV[0,1]+cell.area[0,4]*cell.dLPV[0,4])*(PRES[0]-PRES[5]-ONC[0]+ONC[5]-osm)
	# 	jvol[5,0] = -jvol[0,5]
		#print(cell.area[0,1],cell.dLPV[0,1],cell.area[0,4],cell.dLPV[0,4],PRES[0],PRES[5],ONC[0],ONC[5],osm)
	# else:
	for k in range(NC-1):
		for l in range(k+1,NC):
			osm = sum(cell.sig[:,k,l]*(cell.conc.T[k,:]-cell.conc.T[l,:]))
			# print(cell.conc.T[k,:])
			# print(cell.conc.T[l,:])
			# osm = 0.0
			# for i in range(NS):
			# 	osm += cell.sig[i,k,l]*(cell.conc[i,k]-cell.conc[i,l])
			# #check osm
			# if k == 4:
			# 	if l ==5:
			# 		print(osm)
			# 		input()
			
			jvol[k][l] = cell.area[k][l] * cell.dLPV[k][l] * (PRES[k]-PRES[l]-ONC[k]+ONC[l]-osm)

			jvol[l][k] = -jvol[k][l]

			if (cell.flag==2) and (i==0) and (j==1):
				print(cell.area[i][j] * cell.dLPV[i][j],'\t',PRES[i]-PRES[j],'\t',-ONC[i]+ONC[j],'\t',-osm)
				print('\n')
				
			if (cell.flag==2) and (i==0) and (j==4):
				print(cell.area[i][j] * cell.dLPV[i][j],'\t',PRES[i]-PRES[j],'\t',-ONC[i]+ONC[j],'\t',-osm)
				print('\n')
	#print(cell.area[4,5],cell.dLPV[4,5],PRES[4],PRES[5],ONC[4],ONC[5])
	#print(cell.area[0,1],cell.dLPV[0,1],cell.area[0,4],cell.dLPV[0,4],PRES[0],PRES[5],ONC[0],ONC[5])
	#input()
	return jvol


	
