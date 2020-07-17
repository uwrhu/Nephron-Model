# This is used to simulate the whole superficial nephron. Type 'python simulation.py' in terminal to run. If you want to run individual segment, use computation.py.

from driver import compute
from values import *
from defs import *
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
import flux
import os
import argparse

solute = ['Na','K','Cl','HCO3','H2CO3','CO2','HPO4','H2PO4','urea','NH3','NH4','H','HCO2','H2CO2','glu']
compart = ['Lumen','Cell','ICA','ICB','LIS','Bath']
cw=Vref*60e6

#gender = input('Which gender do you want to simulate? (Male/Female) ')
#diabete = input('Does this kidney have diabete? (Y/N) ')
#inhib = input('Any transporter inhibition? (nhe3 50%/ nhe3 80%/ NKCC2 70%/ NKCC2 100%) ')

parser = argparse.ArgumentParser()
parser.add_argument('--sex',choices=['Male','Female'],required = True,type = str,help = 'sex of rat')
parser.add_argument('--diabete',choices=['Y','N'],required = True,type = str,help = 'diabete status')
parser.add_argument('--inhibition',choices=['NHE3','NKCC2','NCC','ENaC','SNB'],default = None,type = str,help = 'any transporter inhibition')
parser.add_argument('--percentage',default = 0,type = float,help = 'percentage of inhibition')
parser.add_argument('--salt_load',choices = ['Y','N'],default = 'N',type = str,help = 'whether to simulate saline load simulation')

args = parser.parse_args()
gender = args.sex
diabete = args.diabete
inhib = args.inhibition
perc = args.percentage
salt = args.salt_load

file_to_save = input('Input a file name to save data: ')
if os.path.isdir(file_to_save) == False:
	os.makedirs(file_to_save)

parts = ['sup','jux1','jux2','jux3','jux4','jux5']

for sup_or_jux in parts:
	print(sup_or_jux+' begin.')
	#========================================================
	# Proximal convolute tubule
	#========================================================
	NPT = 176
	if gender == 'Male':
		filename = 'PTparams_M.dat'
	elif gender == 'Female':
		filename = 'PTparams_F.dat'
	else:
		filename ='PTparams_F.dat'
	pt=compute(NPT,filename,'Broyden',sup_or_jux,diabete,inhib,perc,salt)
	#========================================================
	# output PCT Concentrations in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+pt[0].sex+'pt_con_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NPT):
			file.write(str(pt[j].conc[i,0])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+pt[0].sex+'pt_con_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NPT):
			file.write(str(pt[j].conc[i,1])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+pt[0].sex+'pt_con_of_'+solute[i]+'_in_Bath_'+sup_or_jux+'.txt','w')
		for j in range(NPT):
			file.write(str(pt[j].conc[i,5])+'\n')
		file.close()
	#========================================================
	# output PCT Water volume in Lumen and Cell
	#========================================================
	file=open('./'+file_to_save+'/'+pt[0].sex+'pt_water_volume_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NPT):
		file.write(str(pt[j].vol[0]*cw)+'\n')
	file.close()
	file=open('./'+file_to_save+'/'+pt[0].sex+'pt_water_volume_in_Cell_'+sup_or_jux+'.txt','w')
	for j in range(NPT):
		file.write(str(pt[j].vol[1]*cw)+'\n')
	file.close()
	#========================================================
	# output PT solute flows in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+pt[0].sex+'pt_flow_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NPT):
			file.write(str(pt[j].conc[i,0]*pt[j].vol[0]*cw)+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+pt[0].sex+'pt_flow_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NPT):
			file.write(str(pt[j].conc[i,1]*pt[j].vol[1]*cw)+'\n')
		file.close()
	#========================================================
	# output PCT osmolality in Lumen, Cell, LIS, Bath
	#========================================================
	file_lumen = open('./'+file_to_save+'/'+pt[0].sex+pt[0].segment+'_osmolality_in_Lumen_'+sup_or_jux+'.txt','w')
	file_cell = open('./'+file_to_save+'/'+pt[0].sex+pt[0].segment+'_osmolality_in_Cell_'+sup_or_jux+'.txt','w')
	file_lis = open('./'+file_to_save+'/'+pt[0].sex+pt[0].segment+'_osmolality_in_LIS_'+sup_or_jux+'.txt','w')
	file_bath = open('./'+file_to_save+'/'+pt[0].sex+pt[0].segment+'_osmolality_in_Bath_'+sup_or_jux+'.txt','w')
	for j in range(NPT):
		osm_l = 0
		osm_c = 0
		osm_lis = 0
		osm_b = 0
		for i in range(NS):
			osm_l = osm_l +pt[j].conc[i,0]
			osm_c = osm_c +pt[j].conc[i,1]
			osm_lis = osm_lis+pt[j].conc[i,4]
			osm_b = osm_b +pt[j].conc[i,5]

		file_lumen.write(str(osm_l)+'\n')
		file_cell.write(str(osm_c)+'\n')
		file_lis.write(str(osm_lis)+'\n')
		file_bath.write(str(osm_b)+'\n')
	file_lumen.close()
	file_cell.close()
	file_lis.close()
	file_bath.close()
	#========================================================
	# output luminal pressure
	#========================================================
	file=open('./'+file_to_save+'/'+pt[0].sex+pt[0].segment+'_pressure_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NPT):
		file.write(str(pt[j].pres[0])+'\n')
	file.close()
	#========================================================
	# output fluxes through transporters
	#========================================================
	for j in range(NPT):
		for i in range(len(pt[j].trans)):
			transporter_type = pt[j].trans[i].type
			memb_id = pt[j].trans[i].membrane_id

			if transporter_type == 'SGLT1':
				solute_id,fluxs = glucose.sglt1(pt[j],pt[j].ep,memb_id,pt[j].trans[i].act,pt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'SGLT2':
				solute_id,fluxs = glucose.sglt2(pt[j],pt[j].ep,memb_id,pt[j].trans[i].act,pt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'GLUT1':
				solute_id,fluxs=glucose.glut1(pt[j],pt[j].ep,memb_id,pt[j].trans[i].act,pt[j].area)
				for k in range(len([solute_id])):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs)+'\n')
			elif transporter_type == 'GLUT2':
				solute_id,fluxs=glucose.glut2(pt[j],pt[j].ep,memb_id,pt[j].trans[i].act,pt[j].area)
				for k in range(len([solute_id])):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs)+'\n')			
			elif transporter_type == 'NHE3':
				solute_id,fluxs=NHE3.nhe3(pt[j],pt[j].ep,memb_id,pt[j].trans[i].act,pt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NaKATPase':
				solute_id,fluxs=ATPase.nakatpase(pt[j],pt[j].ep,memb_id,pt[j].trans[i].act,pt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')

			elif transporter_type == 'HATPase':
				solute_id,fluxs=ATPase.hatpase(pt[j],pt[j].ep,memb_id,pt[j].trans[i].act,pt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2A':
				solute_id,fluxs=NKCC.nkcc2(pt[j],memb_id,pt[j].trans[i].act,pt[j].area,'A')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2B':
				solute_id,fluxs=NKCC.nkcc2(pt[j],memb_id,pt[j].trans[i].act,pt[j].area,'B')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2F':
				solute_id,fluxs=NKCC.nkcc2(pt[j],memb_id,pt[j].trans[i].act,pt[j].area,'F')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')       
			elif transporter_type == 'KCC4':
				solute_id,fluxs=KCC.kcc4(pt[j].conc,memb_id,pt[j].trans[i].act,pt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'ENaC':
				solute_id,fluxs=ENaC.ENaC(pt[j],j,memb_id,pt[j].trans[i].act,pt[j].area,jvol)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NCC':
				solute_id,fluxs=NCC.NCC(pt[j],j,memb_id,pt[j].trans[i].act,pt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'Pendrin':
				solute_id,fluxs=Pendrin.Pendrin(pt[j],memb_id,pt[j].trans[i].act,pt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type =='AE1':
				solute_id,fluxs=AE1.AE1(pt[j],memb_id,pt[j].trans[i].act,pt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'HKATPase':
				solute_id,fluxs=ATPase.hkatpase(pt[j],memb_id,pt[j].trans[i].act,pt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NHE1':
				solute_id,fluxs=NHE1.NHE1(pt[j],memb_id,pt[j].trans[i].act,pt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC1':
				solute_id,fluxs=NKCC.nkcc1(pt[j],memb_id,pt[j].trans[i].act,delmu)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+pt[j].sex+'_'+pt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			else:
				print('What is this?',transporter_type)		
	print('PCT finished.')
	#========================================================
	# S3
	#========================================================
	NS3 = 25
	if gender == 'Male':
		filename = 'S3params_M.dat'
	elif gender == 'Female':
		filename = 'S3params_F.dat'
	else:
		filename ='S3params_F.dat'
	s3=compute(NS3,filename,'Newton',sup_or_jux,diabete,inhib,perc)
	#========================================================
	# output S3 Concentrations in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+s3[0].sex+'s3_con_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NS3):
			file.write(str(s3[j].conc[i,0])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+s3[0].sex+'s3_con_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NS3):
			file.write(str(s3[j].conc[i,1])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+s3[0].sex+'s3_con_of_'+solute[i]+'_in_Bath_'+sup_or_jux+'.txt','w')
		for j in range(NS3):
			file.write(str(s3[j].conc[i,5])+'\n')
		file.close()
	#========================================================
	# output S3 Water volume in Lumen and Cell
	#========================================================
	file=open('./'+file_to_save+'/'+s3[0].sex+'s3_water_volume_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NS3):
		file.write(str(s3[j].vol[0]*cw)+'\n')
	file.close()
	file=open('./'+file_to_save+'/'+s3[0].sex+'s3_water_volume_in_Cell_'+sup_or_jux+'.txt','w')
	for j in range(NS3):
		file.write(str(s3[j].vol[1]*cw)+'\n')
	file.close()
	#========================================================
	# output S3 solute flows in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+s3[0].sex+'s3_flow_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NS3):
			file.write(str(s3[j].conc[i,0]*s3[j].vol[0]*cw)+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+s3[0].sex+'s3_flow_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NS3):
			file.write(str(s3[j].conc[i,1]*s3[j].vol[1]*cw)+'\n')
		file.close()
	#========================================================
	# output S3 osmolality in Lumen, Cell, LIS, Bath
	#========================================================
	file_lumen = open('./'+file_to_save+'/'+s3[0].sex+s3[0].segment+'_osmolality_in_Lumen_'+sup_or_jux+'.txt','w')
	file_cell = open('./'+file_to_save+'/'+s3[0].sex+s3[0].segment+'_osmolality_in_Cell_'+sup_or_jux+'.txt','w')
	file_lis = open('./'+file_to_save+'/'+s3[0].sex+s3[0].segment+'_osmolality_in_LIS_'+sup_or_jux+'.txt','w')
	file_bath = open('./'+file_to_save+'/'+s3[0].sex+s3[0].segment+'_osmolality_in_Bath_'+sup_or_jux+'.txt','w')
	for j in range(NS3):
		osm_l = 0
		osm_c = 0
		osm_lis = 0
		osm_b = 0
		for i in range(NS):
			osm_l = osm_l +s3[j].conc[i,0]
			osm_c = osm_c +s3[j].conc[i,1]
			osm_lis = osm_lis+s3[j].conc[i,4]
			osm_b = osm_b +s3[j].conc[i,5]

		file_lumen.write(str(osm_l)+'\n')
		file_cell.write(str(osm_c)+'\n')
		file_lis.write(str(osm_lis)+'\n')
		file_bath.write(str(osm_b)+'\n')
	file_lumen.close()
	file_cell.close()
	file_lis.close()
	file_bath.close()
	#========================================================
	# output luminal pressure
	#========================================================
	file=open('./'+file_to_save+'/'+s3[0].sex+s3[0].segment+'_pressure_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NS3):
		file.write(str(s3[j].pres[0])+'\n')
	file.close()
	
	print('S3 finished.')
	#========================================================
	# Short descending limb
	#========================================================
	NSDL = 200
	if gender == 'Male':
		filename = 'SDLparams_M.dat'
	elif gender == 'Female':
		filename = 'SDLparams_F.dat'
	else:
		filename ='SDLparams_F.dat'
	sdl=compute(NSDL,filename,'Broyden',sup_or_jux,diabete,inhib,perc)
	#========================================================
	# output SDL Concentrations in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+sdl[0].sex+'sdl_con_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NSDL):
			file.write(str(sdl[j].conc[i,0])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+sdl[0].sex+'sdl_con_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NSDL):
			file.write(str(sdl[j].conc[i,1])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+sdl[0].sex+'sdl_con_of_'+solute[i]+'_in_Bath_'+sup_or_jux+'.txt','w')
		for j in range(NSDL):
			file.write(str(sdl[j].conc[i,5])+'\n')
		file.close()
	#========================================================
	# output SDL Water volume in Lumen and Cell
	#========================================================
	file=open('./'+file_to_save+'/'+sdl[0].sex+'sdl_water_volume_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NSDL):
		file.write(str(sdl[j].vol[0]*cw)+'\n')
	file.close()
	file=open('./'+file_to_save+'/'+sdl[0].sex+'sdl_water_volume_in_Cell_'+sup_or_jux+'.txt','w')
	for j in range(NSDL):
		file.write(str(sdl[j].vol[1]*cw)+'\n')
	file.close()
	#========================================================
	# output SDL solute flows in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+sdl[0].sex+'sdl_flow_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NSDL):
			file.write(str(sdl[j].conc[i,0]*sdl[j].vol[0]*cw)+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+sdl[0].sex+'sdl_flow_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NSDL):
			file.write(str(sdl[j].conc[i,1]*sdl[j].vol[1]*cw)+'\n')
		file.close()
	#========================================================
	# output SDL osmolality in Lumen, Cell, LIS, Bath
	#========================================================
	file_lumen = open('./'+file_to_save+'/'+sdl[0].sex+sdl[0].segment+'_osmolality_in_Lumen_'+sup_or_jux+'.txt','w')
	file_cell = open('./'+file_to_save+'/'+sdl[0].sex+sdl[0].segment+'_osmolality_in_Cell_'+sup_or_jux+'.txt','w')
	file_lis = open('./'+file_to_save+'/'+sdl[0].sex+sdl[0].segment+'_osmolality_in_LIS_'+sup_or_jux+'.txt','w')
	file_bath = open('./'+file_to_save+'/'+sdl[0].sex+sdl[0].segment+'_osmolality_in_Bath_'+sup_or_jux+'.txt','w')
	for j in range(NSDL):
		osm_l = 0
		osm_c = 0
		osm_lis = 0
		osm_b = 0
		for i in range(NS):
			osm_l = osm_l +sdl[j].conc[i,0]
			osm_c = osm_c +sdl[j].conc[i,1]
			osm_lis = osm_lis+sdl[j].conc[i,4]
			osm_b = osm_b +sdl[j].conc[i,5]

		file_lumen.write(str(osm_l)+'\n')
		file_cell.write(str(osm_c)+'\n')
		file_lis.write(str(osm_lis)+'\n')
		file_bath.write(str(osm_b)+'\n')
	file_lumen.close()
	file_cell.close()
	file_lis.close()
	file_bath.close()
	#========================================================
	# output luminal pressure
	#========================================================
	file=open('./'+file_to_save+'/'+sdl[0].sex+sdl[0].segment+'_pressure_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NSDL):
		file.write(str(sdl[j].pres[0])+'\n')
	file.close()

	print('SDL finished.')

	#========================================================
	# Long descending limb
	#========================================================
	if sup_or_jux != 'sup':
		NLDL = 200
		if gender == 'Male':
			filename = 'LDLparams_M.dat'
		elif gender == 'Female':
			filename = 'LDLparams_F.dat'
		else:
			filename ='LDLparams_F.dat'
		ldl=compute(NLDL,filename,'Newton',sup_or_jux,diabete,inhib,perc)
	#========================================================
	# output LDL Concentrations in Lumen and Cell
	#========================================================
		for i in range(NS):
			file=open('./'+file_to_save+'/'+ldl[0].sex+'ldl_con_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
			for j in range(NLDL):
				file.write(str(ldl[j].conc[i,0])+'\n')
			file.close()
		for i in range(NS):
			file=open('./'+file_to_save+'/'+ldl[0].sex+'ldl_con_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
			for j in range(NLDL):
				file.write(str(ldl[j].conc[i,1])+'\n')
			file.close()
		for i in range(NS):
			file=open('./'+file_to_save+'/'+ldl[0].sex+'ldl_con_of_'+solute[i]+'_in_Bath_'+sup_or_jux+'.txt','w')
			for j in range(NLDL):
				file.write(str(ldl[j].conc[i,5])+'\n')
			file.close()
	#========================================================
	# output LDL Water volume in Lumen and Cell
	#========================================================
		file=open('./'+file_to_save+'/'+ldl[0].sex+'ldl_water_volume_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NLDL):
			file.write(str(ldl[j].vol[0]*cw)+'\n')
		file.close()
		file=open('./'+file_to_save+'/'+ldl[0].sex+'ldl_water_volume_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NLDL):
			file.write(str(ldl[j].vol[1]*cw)+'\n')
		file.close()
	#========================================================
	# output LDL solute flows in Lumen and Cell
	#========================================================
		for i in range(NS):
			file=open('./'+file_to_save+'/'+ldl[0].sex+'ldl_flow_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
			for j in range(NLDL):
				file.write(str(ldl[j].conc[i,0]*ldl[j].vol[0]*cw)+'\n')
			file.close()
		for i in range(NS):
			file=open('./'+file_to_save+'/'+ldl[0].sex+'ldl_flow_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
			for j in range(NLDL):
				file.write(str(ldl[j].conc[i,1]*ldl[j].vol[1]*cw)+'\n')
			file.close()
	#========================================================
	# output LDL osmolality in Lumen, Cell, LIS, Bath
	#========================================================
		file_lumen = open('./'+file_to_save+'/'+ldl[0].sex+ldl[0].segment+'_osmolality_in_Lumen_'+sup_or_jux+'.txt','w')
		file_cell = open('./'+file_to_save+'/'+ldl[0].sex+ldl[0].segment+'_osmolality_in_Cell_'+sup_or_jux+'.txt','w')
		file_lis = open('./'+file_to_save+'/'+ldl[0].sex+ldl[0].segment+'_osmolality_in_LIS_'+sup_or_jux+'.txt','w')
		file_bath = open('./'+file_to_save+'/'+ldl[0].sex+ldl[0].segment+'_osmolality_in_Bath_'+sup_or_jux+'.txt','w')
		for j in range(NLDL):
			osm_l = 0
			osm_c = 0
			osm_lis = 0
			osm_b = 0
			for i in range(NS):
				osm_l = osm_l +ldl[j].conc[i,0]
				osm_c = osm_c +ldl[j].conc[i,1]
				osm_lis = osm_lis+ldl[j].conc[i,4]
				osm_b = osm_b +ldl[j].conc[i,5]

			file_lumen.write(str(osm_l)+'\n')
			file_cell.write(str(osm_c)+'\n')
			file_lis.write(str(osm_lis)+'\n')
			file_bath.write(str(osm_b)+'\n')
		file_lumen.close()
		file_cell.close()
		file_lis.close()
		file_bath.close()
	#========================================================
	# output luminal pressure
	#========================================================
		file=open('./'+file_to_save+'/'+ldl[0].sex+ldl[0].segment+'_pressure_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NLDL):
			file.write(str(ldl[j].pres[0])+'\n')
		file.close()
		print('LDL finished.')
	#========================================================
	# Long ascending limb
	#========================================================
		NLAL = 200
		if gender == 'Male':
			filename = 'LALparams_M.dat'
		elif gender == 'Female':
			filename = 'LALparams_F.dat'
		else:
			filename ='LALparams_F.dat'
		lal=compute(NLAL,filename,'Newton',sup_or_jux,diabete,inhib,perc)
	#========================================================
	# output LAL Concentrations in Lumen and Cell
	#========================================================
		for i in range(NS):
			file=open('./'+file_to_save+'/'+lal[0].sex+'lal_con_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
			for j in range(NLAL):
				file.write(str(lal[j].conc[i,0])+'\n')
			file.close()
		for i in range(NS):
			file=open('./'+file_to_save+'/'+lal[0].sex+'lal_con_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
			for j in range(NLAL):
				file.write(str(lal[j].conc[i,1])+'\n')
			file.close()
		for i in range(NS):
			file=open('./'+file_to_save+'/'+lal[0].sex+'lal_con_of_'+solute[i]+'_in_Bath_'+sup_or_jux+'.txt','w')
			for j in range(NLAL):
				file.write(str(lal[j].conc[i,5])+'\n')
			file.close()
	#========================================================
	# output LAL Water volume in Lumen and Cell
	#========================================================
		file=open('./'+file_to_save+'/'+lal[0].sex+'lal_water_volume_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NLAL):
			file.write(str(lal[j].vol[0]*cw)+'\n')
		file.close()
		file=open('./'+file_to_save+'/'+lal[0].sex+'lal_water_volume_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NLAL):
			file.write(str(lal[j].vol[1]*cw)+'\n')
		file.close()
	#========================================================
	# output LAL solute flows in Lumen and Cell
	#========================================================
		for i in range(NS):
			file=open('./'+file_to_save+'/'+lal[0].sex+'lal_flow_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
			for j in range(NLAL):
				file.write(str(lal[j].conc[i,0]*lal[j].vol[0]*cw)+'\n')
			file.close()
		for i in range(NS):
			file=open('./'+file_to_save+'/'+lal[0].sex+'lal_flow_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
			for j in range(NLAL):
				file.write(str(lal[j].conc[i,1]*lal[j].vol[1]*cw)+'\n')
			file.close()
	#========================================================
	# output LAL osmolality in Lumen, Cell, LIS, Bath
	#========================================================
		file_lumen = open('./'+file_to_save+'/'+lal[0].sex+lal[0].segment+'_osmolality_in_Lumen_'+sup_or_jux+'.txt','w')
		file_cell = open('./'+file_to_save+'/'+lal[0].sex+lal[0].segment+'_osmolality_in_Cell_'+sup_or_jux+'.txt','w')
		file_lis = open('./'+file_to_save+'/'+lal[0].sex+lal[0].segment+'_osmolality_in_LIS_'+sup_or_jux+'.txt','w')
		file_bath = open('./'+file_to_save+'/'+lal[0].sex+lal[0].segment+'_osmolality_in_Bath_'+sup_or_jux+'.txt','w')
		for j in range(NLAL):
			osm_l = 0
			osm_c = 0
			osm_lis = 0
			osm_b = 0
			for i in range(NS):
				osm_l = osm_l +lal[j].conc[i,0]
				osm_c = osm_c +lal[j].conc[i,1]
				osm_lis = osm_lis+lal[j].conc[i,4]
				osm_b = osm_b +lal[j].conc[i,5]

			file_lumen.write(str(osm_l)+'\n')
			file_cell.write(str(osm_c)+'\n')
			file_lis.write(str(osm_lis)+'\n')
			file_bath.write(str(osm_b)+'\n')
		file_lumen.close()
		file_cell.close()
		file_lis.close()
		file_bath.close()
	#========================================================
	# output luminal pressure
	#========================================================
		file=open('./'+file_to_save+'/'+lal[0].sex+lal[0].segment+'_pressure_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NLAL):
			file.write(str(lal[j].pres[0])+'\n')
		file.close()
		print('LAL finished.')
	#========================================================
	# Medulla thick ascending limb
	#========================================================
	NmTAL = 200
	if gender == 'Male':
		filename = 'mTALparams_M.dat'
		method = 'Newton'
	elif gender == 'Female':
		filename = 'mTALparams_F.dat'
		method = 'Newton'
	else:
		filename ='mTALparams_F.dat'
	mtal=compute(NmTAL,filename,method,sup_or_jux,diabete,inhib,perc,salt)
	#========================================================
	# output mTAL Concentrations in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+mtal[0].sex+'mtal_con_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NmTAL):
			file.write(str(mtal[j].conc[i,0])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+mtal[0].sex+'mtal_con_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NmTAL):
			file.write(str(mtal[j].conc[i,1])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+mtal[0].sex+'mtal_con_of_'+solute[i]+'_in_Bath_'+sup_or_jux+'.txt','w')
		for j in range(NmTAL):
			file.write(str(mtal[j].conc[i,5])+'\n')
		file.close()
	#========================================================
	# output mTAL Water volume in Lumen and Cell
	#========================================================
	file=open('./'+file_to_save+'/'+mtal[0].sex+'mtal_water_volume_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NmTAL):
		file.write(str(mtal[j].vol[0]*cw)+'\n')
	file.close()
	file=open('./'+file_to_save+'/'+mtal[0].sex+'mtal_water_volume_in_Cell_'+sup_or_jux+'.txt','w')
	for j in range(NmTAL):
		file.write(str(mtal[j].vol[1]*cw)+'\n')
	file.close()
	#========================================================
	# output mTAL solute flows in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+mtal[0].sex+'mtal_flow_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NmTAL):
			file.write(str(mtal[j].conc[i,0]*mtal[j].vol[0]*cw)+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+mtal[0].sex+'mtal_flow_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NmTAL):
			file.write(str(mtal[j].conc[i,1]*mtal[j].vol[1]*cw)+'\n')
		file.close()
	#========================================================
	# output mTAL osmolality in Lumen, Cell, LIS, Bath
	#========================================================
	file_lumen = open('./'+file_to_save+'/'+mtal[0].sex+mtal[0].segment+'_osmolality_in_Lumen_'+sup_or_jux+'.txt','w')
	file_cell = open('./'+file_to_save+'/'+mtal[0].sex+mtal[0].segment+'_osmolality_in_Cell_'+sup_or_jux+'.txt','w')
	file_lis = open('./'+file_to_save+'/'+mtal[0].sex+mtal[0].segment+'_osmolality_in_LIS_'+sup_or_jux+'.txt','w')
	file_bath = open('./'+file_to_save+'/'+mtal[0].sex+mtal[0].segment+'_osmolality_in_Bath_'+sup_or_jux+'.txt','w')
	for j in range(NmTAL):
		osm_l = 0
		osm_c = 0
		osm_lis = 0
		osm_b = 0
		for i in range(NS):
			osm_l = osm_l +mtal[j].conc[i,0]
			osm_c = osm_c +mtal[j].conc[i,1]
			osm_lis = osm_lis+mtal[j].conc[i,4]
			osm_b = osm_b +mtal[j].conc[i,5]

		file_lumen.write(str(osm_l)+'\n')
		file_cell.write(str(osm_c)+'\n')
		file_lis.write(str(osm_lis)+'\n')
		file_bath.write(str(osm_b)+'\n')
	file_lumen.close()
	file_cell.close()
	file_lis.close()
	file_bath.close()
	#========================================================
	# output luminal pressure
	#========================================================
	file=open('./'+file_to_save+'/'+mtal[0].sex+mtal[0].segment+'_pressure_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NmTAL):
		file.write(str(mtal[j].pres[0])+'\n')
	file.close()
	#========================================================
	# output fluxes through transporters
	#========================================================
	for j in range(NmTAL):
		for i in range(len(mtal[j].trans)):
			transporter_type = mtal[j].trans[i].type
			memb_id = mtal[j].trans[i].membrane_id

			if transporter_type == 'SGLT1':
				solute_id,fluxs = glucose.sglt1(mtal[j],mtal[j].ep,memb_id,mtal[j].trans[i].act,mtal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'SGLT2':
				solute_id,fluxs = glucose.sglt2(mtal[j],mtal[j].ep,memb_id,mtal[j].trans[i].act,mtal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'GLUT1':
				solute_id,fluxs=glucose.glut1(mtal[j],mtal[j].ep,memb_id,mtal[j].trans[i].act,mtal[j].area)
				for k in range(len([solute_id])):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs)+'\n')
			elif transporter_type == 'GLUT2':
				solute_id,fluxs=glucose.glut2(mtal[j],mtal[j].ep,memb_id,mtal[j].trans[i].act,mtal[j].area)
				for k in range(len([solute_id])):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs)+'\n')			
			elif transporter_type == 'NHE3':
				solute_id,fluxs=NHE3.nhe3(mtal[j],mtal[j].ep,memb_id,mtal[j].trans[i].act,mtal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NaKATPase':
				solute_id,fluxs=ATPase.nakatpase(mtal[j],mtal[j].ep,memb_id,mtal[j].trans[i].act,mtal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')

			elif transporter_type == 'HATPase':
				solute_id,fluxs=ATPase.hatpase(mtal[j],mtal[j].ep,memb_id,mtal[j].trans[i].act,mtal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2A':
				solute_id,fluxs=NKCC.nkcc2(mtal[j],memb_id,mtal[j].trans[i].act,mtal[j].area,'A')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2B':
				solute_id,fluxs=NKCC.nkcc2(mtal[j],memb_id,mtal[j].trans[i].act,mtal[j].area,'B')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2F':
				solute_id,fluxs=NKCC.nkcc2(mtal[j],memb_id,mtal[j].trans[i].act,mtal[j].area,'F')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')       
			elif transporter_type == 'KCC4':
				solute_id,fluxs=KCC.kcc4(mtal[j].conc,memb_id,mtal[j].trans[i].act,mtal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'ENaC':
				solute_id,fluxs=ENaC.ENaC(mtal[j],j,memb_id,mtal[j].trans[i].act,mtal[j].area,jvol)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NCC':
				solute_id,fluxs=NCC.NCC(mtal[j],j,memb_id,mtal[j].trans[i].act,mtal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'Pendrin':
				solute_id,fluxs=Pendrin.Pendrin(mtal[j],memb_id,mtal[j].trans[i].act,mtal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type =='AE1':
				solute_id,fluxs=AE1.AE1(mtal[j],memb_id,mtal[j].trans[i].act,mtal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'HKATPase':
				solute_id,fluxs=ATPase.hkatpase(mtal[j],memb_id,mtal[j].trans[i].act,mtal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NHE1':
				solute_id,fluxs=NHE1.NHE1(mtal[j],memb_id,mtal[j].trans[i].act,mtal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC1':
				solute_id,fluxs=NKCC.nkcc1(mtal[j],memb_id,mtal[j].trans[i].act,delmu)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+mtal[j].sex+'_'+mtal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			else:
				print('What is this?',transporter_type)
	print('mTAL finished.')
	#========================================================
	# Cortex thick ascending limb
	#========================================================
	NcTAL = 200
	if gender == 'Male':
		filename = 'cTALparams_M.dat'
	elif gender == 'Female':
		filename = 'cTALparams_F.dat'
	else:
		filename ='cTALparams_F.dat'

	if gender == 'Male':
		method = 'Newton'#'Broyden'
	elif gender == 'Female':
		method = 'Newton'
	ctal=compute(NcTAL,filename,method,sup_or_jux,diabete,inhib,perc,salt)
	#========================================================
	# output cTAL Concentrations in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+ctal[0].sex+'ctal_con_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NcTAL):
			file.write(str(ctal[j].conc[i,0])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+ctal[0].sex+'ctal_con_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NcTAL):
			file.write(str(ctal[j].conc[i,1])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+ctal[0].sex+'ctal_con_of_'+solute[i]+'_in_Bath_'+sup_or_jux+'.txt','w')
		for j in range(NcTAL):
			file.write(str(ctal[j].conc[i,5])+'\n')
		file.close()
	#========================================================
	# output cTAL Water volume in Lumen and Cell
	#========================================================
	file=open('./'+file_to_save+'/'+ctal[0].sex+'ctal_water_volume_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NcTAL):
		file.write(str(ctal[j].vol[0]*cw)+'\n')
	file.close()
	file=open('./'+file_to_save+'/'+ctal[0].sex+'ctal_water_volume_in_Cell_'+sup_or_jux+'.txt','w')
	for j in range(NcTAL):
		file.write(str(ctal[j].vol[1]*cw)+'\n')
	file.close()
	#========================================================
	# output cTAL solute flows in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+ctal[0].sex+'ctal_flow_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NcTAL):
			file.write(str(ctal[j].conc[i,0]*ctal[j].vol[0]*cw)+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+ctal[0].sex+'ctal_flow_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NcTAL):
			file.write(str(ctal[j].conc[i,1]*ctal[j].vol[1]*cw)+'\n')
		file.close()
	#========================================================
	# output cTAL osmolality in Lumen, Cell, LIS, Bath
	#========================================================
	file_lumen = open('./'+file_to_save+'/'+ctal[0].sex+ctal[0].segment+'_osmolality_in_Lumen_'+sup_or_jux+'.txt','w')
	file_cell = open('./'+file_to_save+'/'+ctal[0].sex+ctal[0].segment+'_osmolality_in_Cell_'+sup_or_jux+'.txt','w')
	file_lis = open('./'+file_to_save+'/'+ctal[0].sex+ctal[0].segment+'_osmolality_in_LIS_'+sup_or_jux+'.txt','w')
	file_bath = open('./'+file_to_save+'/'+ctal[0].sex+ctal[0].segment+'_osmolality_in_Bath_'+sup_or_jux+'.txt','w')
	for j in range(NcTAL):
		osm_l = 0
		osm_c = 0
		osm_lis = 0
		osm_b = 0
		for i in range(NS):
			osm_l = osm_l +ctal[j].conc[i,0]
			osm_c = osm_c +ctal[j].conc[i,1]
			osm_lis = osm_lis+ctal[j].conc[i,4]
			osm_b = osm_b +ctal[j].conc[i,5]

		file_lumen.write(str(osm_l)+'\n')
		file_cell.write(str(osm_c)+'\n')
		file_lis.write(str(osm_lis)+'\n')
		file_bath.write(str(osm_b)+'\n')
	file_lumen.close()
	file_cell.close()
	file_lis.close()
	file_bath.close()
	#========================================================
	# output luminal pressure
	#========================================================
	file=open('./'+file_to_save+'/'+ctal[0].sex+ctal[0].segment+'_pressure_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NcTAL):
		file.write(str(ctal[j].pres[0])+'\n')
	file.close()
	#========================================================
	# output fluxes through transporters
	#========================================================
	for j in range(NcTAL):
		for i in range(len(ctal[j].trans)):
			transporter_type = ctal[j].trans[i].type
			memb_id = ctal[j].trans[i].membrane_id

			if transporter_type == 'SGLT1':
				solute_id,fluxs = glucose.sglt1(ctal[j],ctal[j].ep,memb_id,ctal[j].trans[i].act,ctal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'SGLT2':
				solute_id,fluxs = glucose.sglt2(ctal[j],ctal[j].ep,memb_id,ctal[j].trans[i].act,ctal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'GLUT1':
				solute_id,fluxs=glucose.glut1(ctal[j],ctal[j].ep,memb_id,ctal[j].trans[i].act,ctal[j].area)
				for k in range(len([solute_id])):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs)+'\n')
			elif transporter_type == 'GLUT2':
				solute_id,fluxs=glucose.glut2(ctal[j],ctal[j].ep,memb_id,ctal[j].trans[i].act,ctal[j].area)
				for k in range(len([solute_id])):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs)+'\n')			
			elif transporter_type == 'NHE3':
				solute_id,fluxs=NHE3.nhe3(ctal[j],ctal[j].ep,memb_id,ctal[j].trans[i].act,ctal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NaKATPase':
				solute_id,fluxs=ATPase.nakatpase(ctal[j],ctal[j].ep,memb_id,ctal[j].trans[i].act,ctal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')

			elif transporter_type == 'HATPase':
				solute_id,fluxs=ATPase.hatpase(ctal[j],ctal[j].ep,memb_id,ctal[j].trans[i].act,ctal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2A':
				solute_id,fluxs=NKCC.nkcc2(ctal[j],memb_id,ctal[j].trans[i].act,ctal[j].area,'A')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2B':
				solute_id,fluxs=NKCC.nkcc2(ctal[j],memb_id,ctal[j].trans[i].act,ctal[j].area,'B')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2F':
				solute_id,fluxs=NKCC.nkcc2(ctal[j],memb_id,ctal[j].trans[i].act,ctal[j].area,'F')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')       
			elif transporter_type == 'KCC4':
				solute_id,fluxs=KCC.kcc4(ctal[j].conc,memb_id,ctal[j].trans[i].act,ctal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'ENaC':
				solute_id,fluxs=ENaC.ENaC(ctal[j],j,memb_id,ctal[j].trans[i].act,ctal[j].area,jvol)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NCC':
				solute_id,fluxs=NCC.NCC(ctal[j],j,memb_id,ctal[j].trans[i].act,ctal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'Pendrin':
				solute_id,fluxs=Pendrin.Pendrin(ctal[j],memb_id,ctal[j].trans[i].act,ctal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type =='AE1':
				solute_id,fluxs=AE1.AE1(ctal[j],memb_id,ctal[j].trans[i].act,ctal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'HKATPase':
				solute_id,fluxs=ATPase.hkatpase(ctal[j],memb_id,ctal[j].trans[i].act,ctal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NHE1':
				solute_id,fluxs=NHE1.NHE1(ctal[j],memb_id,ctal[j].trans[i].act,ctal[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC1':
				solute_id,fluxs=NKCC.nkcc1(ctal[j],memb_id,ctal[j].trans[i].act,delmu)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+ctal[j].sex+'_'+ctal[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			else:
				print('What is this?',transporter_type)
	print('cTAL finished.')
	#========================================================
	# Macula densa
	#========================================================
	NMD = 2
	if gender == 'Male':
		filename = 'MDparams_M.dat'
	elif gender == 'Female':
		filename = 'MDparams_F.dat'
	else:
		filename ='MDparams_F.dat'
	md=compute(NMD,filename,'Newton',sup_or_jux,diabete,inhib,perc)
	#========================================================
	# output Concentrations of Macula densa in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+md[0].sex+'md_con_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NMD):
			file.write(str(md[j].conc[i,0])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+md[0].sex+'md_con_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NMD):
			file.write(str(md[j].conc[i,1])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+md[0].sex+'md_con_of_'+solute[i]+'_in_Bath_'+sup_or_jux+'.txt','w')
		for j in range(NMD):
			file.write(str(md[j].conc[i,5])+'\n')
		file.close()
	#========================================================
	# output Water volume of Macula densa in Lumen and Cell
	#========================================================
	file=open('./'+file_to_save+'/'+md[0].sex+'md_water_volume_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NMD):
		file.write(str(md[j].vol[0]*cw)+'\n')
	file.close()
	file=open('./'+file_to_save+'/'+md[0].sex+'md_water_volume_in_Cell_'+sup_or_jux+'.txt','w')
	for j in range(NMD):
		file.write(str(md[j].vol[1]*cw)+'\n')
	file.close()
	#========================================================
	# output solute flows of Macula densa in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+md[0].sex+'md_flow_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NMD):
			file.write(str(md[j].conc[i,0]*md[j].vol[0]*cw)+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+md[0].sex+'md_flow_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NMD):
			file.write(str(md[j].conc[i,1]*md[j].vol[1]*cw)+'\n')
		file.close()

	print('Macula densa finished.')
	#========================================================
	# Distal convoluted tubule
	#========================================================
	NDCT = 200
	if gender == 'Male':
		filename = 'DCTparams_M.dat'
	elif gender == 'Female':
		filename = 'DCTparams_F.dat'
	else:
		filename ='DCTparams_F.dat'

	if gender == 'Male':
		method = 'Newton'#'Broyden'
	elif gender == 'Female' and sup_or_jux =='sup':
		method = 'Newton'
	elif gender == 'Female' and sup_or_jux == 'jux1':
		method = 'Newton'#'Broyden'
	else:
		method = 'Newton'#'Broyden'
	dct=compute(NDCT,filename,method,sup_or_jux,diabete,inhib,perc) #female: Broyden
	#========================================================
	# output DCT Concentrations in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+dct[0].sex+'dct_con_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NDCT):
			file.write(str(dct[j].conc[i,0])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+dct[0].sex+'dct_con_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NDCT):
			file.write(str(dct[j].conc[i,1])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+dct[0].sex+'dct_con_of_'+solute[i]+'_in_Bath_'+sup_or_jux+'.txt','w')
		for j in range(NDCT):
			file.write(str(dct[j].conc[i,5])+'\n')
		file.close()
	#========================================================
	# output DCT Water volume in Lumen and Cell
	#========================================================
	file=open('./'+file_to_save+'/'+dct[0].sex+'dct_water_volume_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NDCT):
		file.write(str(dct[j].vol[0]*cw)+'\n')
	file.close()
	file=open('./'+file_to_save+'/'+dct[0].sex+'dct_water_volume_in_Cell_'+sup_or_jux+'.txt','w')
	for j in range(NDCT):
		file.write(str(dct[j].vol[1]*cw)+'\n')
	file.close()
	#========================================================
	# output DCT solute flows in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+dct[0].sex+'dct_flow_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NDCT):
			file.write(str(dct[j].conc[i,0]*dct[j].vol[0]*cw)+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+dct[0].sex+'dct_flow_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NDCT):
			file.write(str(dct[j].conc[i,1]*dct[j].vol[1]*cw)+'\n')
		file.close()
	#========================================================
	# output DCT osmolality in Lumen, Cell, LIS, Bath
	#========================================================
	file_lumen = open('./'+file_to_save+'/'+dct[0].sex+dct[0].segment+'_osmolality_in_Lumen_'+sup_or_jux+'.txt','w')
	file_cell = open('./'+file_to_save+'/'+dct[0].sex+dct[0].segment+'_osmolality_in_Cell_'+sup_or_jux+'.txt','w')
	file_lis = open('./'+file_to_save+'/'+dct[0].sex+dct[0].segment+'_osmolality_in_LIS_'+sup_or_jux+'.txt','w')
	file_bath = open('./'+file_to_save+'/'+dct[0].sex+dct[0].segment+'_osmolality_in_Bath_'+sup_or_jux+'.txt','w')
	for j in range(NDCT):
		osm_l = 0
		osm_c = 0
		osm_lis = 0
		osm_b = 0
		for i in range(NS):
			osm_l = osm_l +dct[j].conc[i,0]
			osm_c = osm_c +dct[j].conc[i,1]
			osm_lis = osm_lis+dct[j].conc[i,4]
			osm_b = osm_b +dct[j].conc[i,5]

		file_lumen.write(str(osm_l)+'\n')
		file_cell.write(str(osm_c)+'\n')
		file_lis.write(str(osm_lis)+'\n')
		file_bath.write(str(osm_b)+'\n')
	file_lumen.close()
	file_cell.close()
	file_lis.close()
	file_bath.close()
	#========================================================
	# output luminal pressure
	#========================================================
	file=open('./'+file_to_save+'/'+dct[0].sex+dct[0].segment+'_pressure_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NDCT):
		file.write(str(dct[j].pres[0])+'\n')
	file.close()
	#========================================================
	# output fluxes through transporters
	#========================================================
	for j in range(NDCT):
		jvol = water.compute_water_fluxes(dct[j])
		jsol,delmu = electrochemical.compute_ecd_fluxes(dct[j],jvol)
		for i in range(len(dct[j].trans)):
			transporter_type = dct[j].trans[i].type
			memb_id = dct[j].trans[i].membrane_id

			if transporter_type == 'SGLT1':
				solute_id,fluxs = glucose.sglt1(dct[j],dct[j].ep,memb_id,dct[j].trans[i].act,dct[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'SGLT2':
				solute_id,fluxs = glucose.sglt2(dct[j],dct[j].ep,memb_id,dct[j].trans[i].act,dct[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'GLUT1':
				solute_id,fluxs=glucose.glut1(dct[j],dct[j].ep,memb_id,dct[j].trans[i].act,dct[j].area)
				for k in range(len([solute_id])):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs)+'\n')
			elif transporter_type == 'GLUT2':
				solute_id,fluxs=glucose.glut2(dct[j],dct[j].ep,memb_id,dct[j].trans[i].act,dct[j].area)
				for k in range(len([solute_id])):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs)+'\n')			
			elif transporter_type == 'NHE3':
				solute_id,fluxs=NHE3.nhe3(dct[j],dct[j].ep,memb_id,dct[j].trans[i].act,dct[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NaKATPase':
				solute_id,fluxs=ATPase.nakatpase(dct[j],dct[j].ep,memb_id,dct[j].trans[i].act,dct[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')

			elif transporter_type == 'HATPase':
				solute_id,fluxs=ATPase.hatpase(dct[j],dct[j].ep,memb_id,dct[j].trans[i].act,dct[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2A':
				solute_id,fluxs=NKCC.nkcc2(dct[j],memb_id,dct[j].trans[i].act,dct[j].area,'A')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2B':
				solute_id,fluxs=NKCC.nkcc2(dct[j],memb_id,dct[j].trans[i].act,dct[j].area,'B')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2F':
				solute_id,fluxs=NKCC.nkcc2(dct[j],memb_id,dct[j].trans[i].act,dct[j].area,'F')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')       
			elif transporter_type == 'KCC4':
				solute_id,fluxs=KCC.kcc4(dct[j].conc,memb_id,dct[j].trans[i].act,dct[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'ENaC':
				solute_id,fluxs=ENaC.ENaC(dct[j],j,memb_id,dct[j].trans[i].act,dct[j].area,jvol)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NCC':
				solute_id,fluxs=NCC.NCC(dct[j],j,memb_id,dct[j].trans[i].act,dct[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'Pendrin':
				solute_id,fluxs=Pendrin.Pendrin(dct[j],memb_id,dct[j].trans[i].act,dct[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type =='AE1':
				solute_id,fluxs=AE1.AE1(dct[j],memb_id,dct[j].trans[i].act,dct[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'HKATPase':
				solute_id,fluxs=ATPase.hkatpase(dct[j],memb_id,dct[j].trans[i].act,dct[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NHE1':
				solute_id,fluxs=NHE1.NHE1(dct[j],memb_id,dct[j].trans[i].act,dct[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC1':
				solute_id,fluxs=NKCC.nkcc1(dct[j],memb_id,dct[j].trans[i].act,delmu)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+dct[j].sex+'_'+dct[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			else:
				print('What is this?',transporter_type)
	print('DCT finished.')
	#========================================================
	# Connecting tubule
	#========================================================
	NCNT = 200
	if gender == 'Male':
		filename = 'CNTparams_M.dat'
	elif gender == 'Female':
		filename = 'CNTparams_F.dat'
	else:
		filename ='CNTparams_F.dat'
	if gender == 'Female': #and sup_or_jux == 'jux5':
		method = 'Newton'
	else:
		method = 'Newton'

	cnt=compute(NCNT,filename,method,sup_or_jux,diabete,inhib,perc,salt)
	#========================================================
	# output CNT Concentrations in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+cnt[0].sex+'cnt_con_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NCNT):
			file.write(str(cnt[j].conc[i,0])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+cnt[0].sex+'cnt_con_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NCNT):
			file.write(str(cnt[j].conc[i,1])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+cnt[0].sex+'cnt_con_of_'+solute[i]+'_in_Bath_'+sup_or_jux+'.txt','w')
		for j in range(NCNT):
			file.write(str(cnt[j].conc[i,5])+'\n')
		file.close()
	#========================================================
	# output CNT Water volume in Lumen and Cell
	#========================================================
	file=open('./'+file_to_save+'/'+cnt[0].sex+'cnt_water_volume_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NCNT):
		file.write(str(cnt[j].vol[0]*cw)+'\n')
	file.close()
	file=open('./'+file_to_save+'/'+cnt[0].sex+'cnt_water_volume_in_Cell_'+sup_or_jux+'.txt','w')
	for j in range(NCNT):
		file.write(str(cnt[j].vol[1]*cw)+'\n')
	file.close()
	#========================================================
	# output CNT solute flows in Lumen and Cell
	#========================================================
	for i in range(NS):
		file=open('./'+file_to_save+'/'+cnt[0].sex+'cnt_flow_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(NCNT):
			file.write(str(cnt[j].conc[i,0]*cnt[j].vol[0]*cw)+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+cnt[0].sex+'cnt_flow_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(NCNT):
			file.write(str(cnt[j].conc[i,1]*cnt[j].vol[1]*cw)+'\n')
		file.close()
	#========================================================
	# output CNT osmolality in Lumen, Cell, LIS, Bath
	#========================================================
	file_lumen = open('./'+file_to_save+'/'+cnt[0].sex+cnt[0].segment+'_osmolality_in_Lumen_'+sup_or_jux+'.txt','w')
	file_cell = open('./'+file_to_save+'/'+cnt[0].sex+cnt[0].segment+'_osmolality_in_Cell_'+sup_or_jux+'.txt','w')
	file_lis = open('./'+file_to_save+'/'+cnt[0].sex+cnt[0].segment+'_osmolality_in_LIS_'+sup_or_jux+'.txt','w')
	file_bath = open('./'+file_to_save+'/'+cnt[0].sex+cnt[0].segment+'_osmolality_in_Bath_'+sup_or_jux+'.txt','w')
	for j in range(NCNT):
		osm_l = 0
		osm_c = 0
		osm_lis = 0
		osm_b = 0
		for i in range(NS):
			osm_l = osm_l +cnt[j].conc[i,0]
			osm_c = osm_c +cnt[j].conc[i,1]
			osm_lis = osm_lis+cnt[j].conc[i,4]
			osm_b = osm_b +cnt[j].conc[i,5]

		file_lumen.write(str(osm_l)+'\n')
		file_cell.write(str(osm_c)+'\n')
		file_lis.write(str(osm_lis)+'\n')
		file_bath.write(str(osm_b)+'\n')
	file_lumen.close()
	file_cell.close()
	file_lis.close()
	file_bath.close()
	#========================================================
	# output luminal pressure
	#========================================================
	file=open('./'+file_to_save+'/'+cnt[0].sex+cnt[0].segment+'_pressure_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(NCNT):
		file.write(str(cnt[j].pres[0])+'\n')
	file.close()
	print('CNT finished.')
	#========================================================
	# output fluxes through transporters
	#========================================================
	for j in range(NCNT):
		jvol = water.compute_water_fluxes(cnt[j])
		jsol,delmu = electrochemical.compute_ecd_fluxes(cnt[j],jvol)
		for i in range(len(cnt[j].trans)):
			transporter_type = cnt[j].trans[i].type
			memb_id = cnt[j].trans[i].membrane_id

			if transporter_type == 'SGLT1':
				solute_id,fluxs = glucose.sglt1(cnt[j],cnt[j].ep,memb_id,cnt[j].trans[i].act,cnt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'SGLT2':
				solute_id,fluxs = glucose.sglt2(cnt[j],cnt[j].ep,memb_id,cnt[j].trans[i].act,cnt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'GLUT1':
				solute_id,fluxs=glucose.glut1(cnt[j],cnt[j].ep,memb_id,cnt[j].trans[i].act,cnt[j].area)
				for k in range(len([solute_id])):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs)+'\n')
			elif transporter_type == 'GLUT2':
				solute_id,fluxs=glucose.glut2(cnt[j],cnt[j].ep,memb_id,cnt[j].trans[i].act,cnt[j].area)
				for k in range(len([solute_id])):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs)+'\n')			
			elif transporter_type == 'NHE3':
				solute_id,fluxs=NHE3.nhe3(cnt[j],cnt[j].ep,memb_id,cnt[j].trans[i].act,cnt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NaKATPase':
				solute_id,fluxs=ATPase.nakatpase(cnt[j],cnt[j].ep,memb_id,cnt[j].trans[i].act,cnt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+str(memb_id[0])+str(memb_id[1])+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')

			elif transporter_type == 'HATPase':
				solute_id,fluxs=ATPase.hatpase(cnt[j],cnt[j].ep,memb_id,cnt[j].trans[i].act,cnt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2A':
				solute_id,fluxs=NKCC.nkcc2(cnt[j],memb_id,cnt[j].trans[i].act,cnt[j].area,'A')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2B':
				solute_id,fluxs=NKCC.nkcc2(cnt[j],memb_id,cnt[j].trans[i].act,cnt[j].area,'B')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC2F':
				solute_id,fluxs=NKCC.nkcc2(cnt[j],memb_id,cnt[j].trans[i].act,cnt[j].area,'F')
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')       
			elif transporter_type == 'KCC4':
				solute_id,fluxs=KCC.kcc4(cnt[j].conc,memb_id,cnt[j].trans[i].act,cnt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'ENaC':
				solute_id,fluxs=ENaC.ENaC(cnt[j],j,memb_id,cnt[j].trans[i].act,cnt[j].area,jvol)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NCC':
				solute_id,fluxs=NCC.NCC(cnt[j],j,memb_id,cnt[j].trans[i].act,cnt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'Pendrin':
				solute_id,fluxs=Pendrin.Pendrin(cnt[j],memb_id,cnt[j].trans[i].act,cnt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type =='AE1':
				solute_id,fluxs=AE1.AE1(cnt[j],memb_id,cnt[j].trans[i].act,cnt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'HKATPase':
				solute_id,fluxs=ATPase.hkatpase(cnt[j],memb_id,cnt[j].trans[i].act,cnt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NHE1':
				solute_id,fluxs=NHE1.NHE1(cnt[j],memb_id,cnt[j].trans[i].act,cnt[j].area)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			elif transporter_type == 'NKCC1':
				solute_id,fluxs=NKCC.nkcc1(cnt[j],memb_id,cnt[j].trans[i].act,delmu)
				for k in range(len(solute_id)):
					file = open('./'+file_to_save+'/'+cnt[j].sex+'_'+cnt[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'_'+sup_or_jux+'.txt','a')
					file.write(str(fluxs[k])+'\n')
			else:
				print('What is this?',transporter_type)
	print(sup_or_jux+' finished')
#========================================================
# Cortical collecting duct
#========================================================
NCCD = 200
if gender == 'Male':
	filename = 'CCDparams_M.dat'
elif gender == 'Female':
	filename = 'CCDparams_F.dat'
else:
	filename ='CCDparams_F.dat'
ccd=compute(NCCD,filename,'Newton',None,diabete,inhib,perc)
#========================================================
# output CCD Concentrations in Lumen and Cell
#========================================================
for i in range(NS):
	file=open('./'+file_to_save+'/'+ccd[0].sex+'ccd_con_of_'+solute[i]+'_in_Lumen.txt','w')
	for j in range(NCCD):
		file.write(str(ccd[j].conc[i,0])+'\n')
	file.close()
for i in range(NS):
	file=open('./'+file_to_save+'/'+ccd[0].sex+'ccd_con_of_'+solute[i]+'_in_Cell.txt','w')
	for j in range(NCCD):
		file.write(str(ccd[j].conc[i,1])+'\n')
	file.close()
for i in range(NS):
	file=open('./'+file_to_save+'/'+ccd[0].sex+'ccd_con_of_'+solute[i]+'_in_Bath.txt','w')
	for j in range(NCCD):
		file.write(str(ccd[j].conc[i,5])+'\n')
	file.close()
#========================================================
# output CCD Water volume in Lumen and Cell
#========================================================
file=open('./'+file_to_save+'/'+ccd[0].sex+'ccd_water_volume_in_Lumen.txt','w')
for j in range(NCCD):
	file.write(str(ccd[j].vol[0]*cw)+'\n')
file.close()
file=open('./'+file_to_save+'/'+ccd[0].sex+'ccd_water_volume_in_Cell.txt','w')
for j in range(NCCD):
	file.write(str(ccd[j].vol[1]*cw)+'\n')
file.close()
#========================================================
# output CCD solute flows in Lumen and Cell
#========================================================
for i in range(NS):
	file=open('./'+file_to_save+'/'+ccd[0].sex+'ccd_flow_of_'+solute[i]+'_in_Lumen.txt','w')
	for j in range(NCCD):
		file.write(str(ccd[j].conc[i,0]*ccd[j].vol[0]*cw)+'\n')
	file.close()
for i in range(NS):
	file=open('./'+file_to_save+'/'+ccd[0].sex+'ccd_flow_of_'+solute[i]+'_in_Cell.txt','w')
	for j in range(NCCD):
		file.write(str(ccd[j].conc[i,1]*ccd[j].vol[1]*cw)+'\n')
	file.close()
#========================================================
# output CCD osmolality in Lumen, Cell, LIS, Bath
#========================================================
file_lumen = open('./'+file_to_save+'/'+ccd[0].sex+ccd[0].segment+'_osmolality_in_Lumen.txt','w')
file_cell = open('./'+file_to_save+'/'+ccd[0].sex+ccd[0].segment+'_osmolality_in_Cell.txt','w')
file_lis = open('./'+file_to_save+'/'+ccd[0].sex+ccd[0].segment+'_osmolality_in_LIS.txt','w')
file_bath = open('./'+file_to_save+'/'+ccd[0].sex+ccd[0].segment+'_osmolality_in_Bath.txt','w')
for j in range(NCCD):
	osm_l = 0
	osm_c = 0
	osm_lis = 0
	osm_b = 0
	for i in range(NS):
		osm_l = osm_l +ccd[j].conc[i,0]
		osm_c = osm_c +ccd[j].conc[i,1]
		osm_lis = osm_lis+ccd[j].conc[i,4]
		osm_b = osm_b +ccd[j].conc[i,5]

	file_lumen.write(str(osm_l)+'\n')
	file_cell.write(str(osm_c)+'\n')
	file_lis.write(str(osm_lis)+'\n')
	file_bath.write(str(osm_b)+'\n')
file_lumen.close()
file_cell.close()
file_lis.close()
file_bath.close()
#========================================================
# output luminal pressure
#========================================================
file=open('./'+file_to_save+'/'+ccd[0].sex+ccd[0].segment+'_pressure_in_Lumen_'+sup_or_jux+'.txt','w')
for j in range(NCCD):
	file.write(str(ccd[j].pres[0])+'\n')
file.close()
#========================================================
# output flux through transporters
#========================================================
for j in range(NCCD):
	jvol = water.compute_water_fluxes(ccd[j])
	jsol,delmu = electrochemical.compute_ecd_fluxes(ccd[j],jvol)
	for i in range(len(ccd[j].trans)):
		transporter_type = ccd[j].trans[i].type
		memb_id = ccd[j].trans[i].membrane_id

		if transporter_type == 'SGLT1':
			solute_id,fluxs = glucose.sglt1(ccd[j],ccd[j].ep,memb_id,ccd[j].trans[i].act,ccd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'SGLT2':
			solute_id,fluxs = glucose.sglt2(ccd[j],ccd[j].ep,memb_id,ccd[j].trans[i].act,ccd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'GLUT1':
			solute_id,fluxs=glucose.glut1(ccd[j],ccd[j].ep,memb_id,ccd[j].trans[i].act,ccd[j].area)
			for k in range(len([solute_id])):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'.txt','a')
				file.write(str(fluxs)+'\n')
		elif transporter_type == 'GLUT2':
			solute_id,fluxs=glucose.glut2(ccd[j],ccd[j].ep,memb_id,ccd[j].trans[i].act,ccd[j].area)
			for k in range(len([solute_id])):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'.txt','a')
				file.write(str(fluxs)+'\n')			
		elif transporter_type == 'NHE3':
			solute_id,fluxs=NHE3.nhe3(ccd[j],ccd[j].ep,memb_id,ccd[j].trans[i].act,ccd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NaKATPase':
			solute_id,fluxs=ATPase.nakatpase(ccd[j],ccd[j].ep,memb_id,ccd[j].trans[i].act,ccd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+str(memb_id[0])+str(memb_id[1])+'.txt','a')
				file.write(str(fluxs[k])+'\n')

		elif transporter_type == 'HATPase':
			solute_id,fluxs=ATPase.hatpase(ccd[j],ccd[j].ep,memb_id,ccd[j].trans[i].act,ccd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC2A':
			solute_id,fluxs=NKCC.nkcc2(ccd[j],memb_id,ccd[j].trans[i].act,ccd[j].area,'A')
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC2B':
			solute_id,fluxs=NKCC.nkcc2(ccd[j],memb_id,ccd[j].trans[i].act,ccd[j].area,'B')
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC2F':
			solute_id,fluxs=NKCC.nkcc2(ccd[j],memb_id,ccd[j].trans[i].act,ccd[j].area,'F')
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')       
		elif transporter_type == 'KCC4':
			solute_id,fluxs=KCC.kcc4(ccd[j].conc,memb_id,ccd[j].trans[i].act,ccd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'ENaC':
			solute_id,fluxs=ENaC.ENaC(ccd[j],j,memb_id,ccd[j].trans[i].act,ccd[j].area,jvol)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NCC':
			solute_id,fluxs=NCC.NCC(ccd[j],j,memb_id,ccd[j].trans[i].act,ccd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'Pendrin':
			solute_id,fluxs=Pendrin.Pendrin(ccd[j],memb_id,ccd[j].trans[i].act,ccd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type =='AE1':
			solute_id,fluxs=AE1.AE1(ccd[j],memb_id,ccd[j].trans[i].act,ccd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'HKATPase':
			solute_id,fluxs=ATPase.hkatpase(ccd[j],memb_id,ccd[j].trans[i].act,ccd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NHE1':
			solute_id,fluxs=NHE1.NHE1(ccd[j],memb_id,ccd[j].trans[i].act,ccd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC1':
			solute_id,fluxs=NKCC.nkcc1(ccd[j],memb_id,ccd[j].trans[i].act,delmu)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+ccd[j].sex+'_'+ccd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		else:
			print('What is this?',transporter_type)
print('CCD finished.')
#========================================================
# Outer medullary collecting duct
#========================================================
NOMCD = 200
if gender == 'Male':
	filename = 'OMCDparams_M.dat'
elif gender == 'Female':
	filename = 'OMCDparams_F.dat'
else:
	filename ='OMCDparams_F.dat'
if ccd[0].sex == 'male':
	omcd=compute(NOMCD,filename,'Newton',None,diabete,inhib,perc)
elif ccd[0].sex == 'female':
	omcd=compute(NOMCD,filename,'Newton',None,diabete,inhib,perc)
#========================================================
# output OMCD Concentrations in Lumen and Cell
#========================================================
for i in range(NS):
	file=open('./'+file_to_save+'/'+omcd[0].sex+'omcd_con_of_'+solute[i]+'_in_Lumen.txt','w')
	for j in range(NOMCD):
		file.write(str(omcd[j].conc[i,0])+'\n')
	file.close()
for i in range(NS):
	file=open('./'+file_to_save+'/'+omcd[0].sex+'omcd_con_of_'+solute[i]+'_in_Cell.txt','w')
	for j in range(NOMCD):
		file.write(str(omcd[j].conc[i,1])+'\n')
	file.close()
for i in range(NS):
	file=open('./'+file_to_save+'/'+omcd[0].sex+'omcd_con_of_'+solute[i]+'_in_Bath.txt','w')
	for j in range(NOMCD):
		file.write(str(omcd[j].conc[i,5])+'\n')
	file.close()
#========================================================
# output OMCD Water volume in Lumen and Cell
#========================================================
file=open('./'+file_to_save+'/'+omcd[0].sex+'omcd_water_volume_in_Lumen.txt','w')
for j in range(NOMCD):
	file.write(str(omcd[j].vol[0]*cw)+'\n')
file.close()
file=open('./'+file_to_save+'/'+omcd[0].sex+'omcd_water_volume_in_Cell.txt','w')
for j in range(NOMCD):
	file.write(str(omcd[j].vol[1]*cw)+'\n')
file.close()
#========================================================
# output OMCD solute flows in Lumen and Cell
#========================================================
for i in range(NS):
	file=open('./'+file_to_save+'/'+omcd[0].sex+'omcd_flow_of_'+solute[i]+'_in_Lumen.txt','w')
	for j in range(NOMCD):
		file.write(str(omcd[j].conc[i,0]*omcd[j].vol[0]*cw)+'\n')
	file.close()
for i in range(NS):
	file=open('./'+file_to_save+'/'+omcd[0].sex+'omcd_flow_of_'+solute[i]+'_in_Cell.txt','w')
	for j in range(NOMCD):
		file.write(str(omcd[j].conc[i,1]*omcd[j].vol[1]*cw)+'\n')
	file.close()
#========================================================
# output OMCD osmolality in Lumen, Cell, LIS, Bath
#========================================================
file_lumen = open('./'+file_to_save+'/'+omcd[0].sex+omcd[0].segment+'_osmolality_in_Lumen.txt','w')
file_cell = open('./'+file_to_save+'/'+omcd[0].sex+omcd[0].segment+'_osmolality_in_Cell.txt','w')
file_lis = open('./'+file_to_save+'/'+omcd[0].sex+omcd[0].segment+'_osmolality_in_LIS.txt','w')
file_bath = open('./'+file_to_save+'/'+omcd[0].sex+omcd[0].segment+'_osmolality_in_Bath.txt','w')
for j in range(NOMCD):
	osm_l = 0
	osm_c = 0
	osm_lis = 0
	osm_b = 0
	for i in range(NS):
		osm_l = osm_l +omcd[j].conc[i,0]
		osm_c = osm_c +omcd[j].conc[i,1]
		osm_lis = osm_lis+omcd[j].conc[i,4]
		osm_b = osm_b +omcd[j].conc[i,5]

	file_lumen.write(str(osm_l)+'\n')
	file_cell.write(str(osm_c)+'\n')
	file_lis.write(str(osm_lis)+'\n')
	file_bath.write(str(osm_b)+'\n')
file_lumen.close()
file_cell.close()
file_lis.close()
file_bath.close()
#========================================================
# output luminal pressure
#========================================================
file=open('./'+file_to_save+'/'+omcd[0].sex+omcd[0].segment+'_pressure_in_Lumen_'+sup_or_jux+'.txt','w')
for j in range(NOMCD):
	file.write(str(omcd[j].pres[0])+'\n')
file.close()
#========================================================
# output flux through transporters
#========================================================
for j in range(NOMCD):
	jvol = water.compute_water_fluxes(omcd[j])
	jsol,delmu = electrochemical.compute_ecd_fluxes(omcd[j],jvol)
	for i in range(len(omcd[j].trans)):
		transporter_type = omcd[j].trans[i].type
		memb_id = omcd[j].trans[i].membrane_id

		if transporter_type == 'SGLT1':
			solute_id,fluxs = glucose.sglt1(omcd[j],omcd[j].ep,memb_id,omcd[j].trans[i].act,omcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'SGLT2':
			solute_id,fluxs = glucose.sglt2(omcd[j],omcd[j].ep,memb_id,omcd[j].trans[i].act,omcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'GLUT1':
			solute_id,fluxs=glucose.glut1(omcd[j],omcd[j].ep,memb_id,omcd[j].trans[i].act,omcd[j].area)
			for k in range(len([solute_id])):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'.txt','a')
				file.write(str(fluxs)+'\n')
		elif transporter_type == 'GLUT2':
			solute_id,fluxs=glucose.glut2(omcd[j],omcd[j].ep,memb_id,omcd[j].trans[i].act,omcd[j].area)
			for k in range(len([solute_id])):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'.txt','a')
				file.write(str(fluxs)+'\n')			
		elif transporter_type == 'NHE3':
			solute_id,fluxs=NHE3.nhe3(omcd[j],omcd[j].ep,memb_id,omcd[j].trans[i].act,omcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NaKATPase':
			solute_id,fluxs=ATPase.nakatpase(omcd[j],omcd[j].ep,memb_id,omcd[j].trans[i].act,omcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+str(memb_id[0])+str(memb_id[1])+'.txt','a')
				file.write(str(fluxs[k])+'\n')

		elif transporter_type == 'HATPase':
			solute_id,fluxs=ATPase.hatpase(omcd[j],omcd[j].ep,memb_id,omcd[j].trans[i].act,omcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC2A':
			solute_id,fluxs=NKCC.nkcc2(omcd[j],memb_id,omcd[j].trans[i].act,omcd[j].area,'A')
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC2B':
			solute_id,fluxs=NKCC.nkcc2(omcd[j],memb_id,omcd[j].trans[i].act,omcd[j].area,'B')
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC2F':
			solute_id,fluxs=NKCC.nkcc2(omcd[j],memb_id,omcd[j].trans[i].act,omcd[j].area,'F')
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')       
		elif transporter_type == 'KCC4':
			solute_id,fluxs=KCC.kcc4(omcd[j].conc,memb_id,omcd[j].trans[i].act,omcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'ENaC':
			solute_id,fluxs=ENaC.ENaC(omcd[j],j,memb_id,omcd[j].trans[i].act,omcd[j].area,jvol)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NCC':
			solute_id,fluxs=NCC.NCC(omcd[j],j,memb_id,omcd[j].trans[i].act,omcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'Pendrin':
			solute_id,fluxs=Pendrin.Pendrin(omcd[j],memb_id,omcd[j].trans[i].act,omcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type =='AE1':
			solute_id,fluxs=AE1.AE1(omcd[j],memb_id,omcd[j].trans[i].act,omcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'HKATPase':
			solute_id,fluxs=ATPase.hkatpase(omcd[j],memb_id,omcd[j].trans[i].act,omcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NHE1':
			solute_id,fluxs=NHE1.NHE1(omcd[j],memb_id,omcd[j].trans[i].act,omcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC1':
			solute_id,fluxs=NKCC.nkcc1(omcd[j],memb_id,omcd[j].trans[i].act,delmu)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+omcd[j].sex+'_'+omcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		else:
			print('What is this?',transporter_type)
print('OMCD finished.')
#========================================================
# Inner medullary collecting duct
#========================================================
NIMCD = 200
if gender == 'Male':
	filename = 'IMCDparams_M.dat'
elif gender == 'Female':
	filename = 'IMCDparams_F.dat'
else:
	filename ='IMCDparams_F.dat'
imcd=compute(NIMCD,filename,'Newton',None,diabete,inhib,perc)
#========================================================
# output IMCD Concentrations in Lumen and Cell
#========================================================
for i in range(NS):
	file=open('./'+file_to_save+'/'+imcd[0].sex+'imcd_con_of_'+solute[i]+'_in_Lumen.txt','w')
	for j in range(NIMCD):
		file.write(str(imcd[j].conc[i,0])+'\n')
	file.close()
for i in range(NS):
	file=open('./'+file_to_save+'/'+imcd[0].sex+'imcd_con_of_'+solute[i]+'_in_Cell.txt','w')
	for j in range(NIMCD):
		file.write(str(imcd[j].conc[i,1])+'\n')
	file.close()
for i in range(NS):
	file=open('./'+file_to_save+'/'+imcd[0].sex+'imcd_con_of_'+solute[i]+'_in_Bath.txt','w')
	for j in range(NIMCD):
		file.write(str(imcd[j].conc[i,5])+'\n')
	file.close()
#========================================================
# output IMCD Water volume in Lumen and Cell
#========================================================
file=open('./'+file_to_save+'/'+imcd[0].sex+'imcd_water_volume_in_Lumen.txt','w')
for j in range(NIMCD):
	file.write(str(imcd[j].vol[0]*cw)+'\n')
file.close()
file=open('./'+file_to_save+'/'+imcd[0].sex+'imcd_water_volume_in_Cell.txt','w')
for j in range(NIMCD):
	file.write(str(imcd[j].vol[1]*cw)+'\n')
file.close()
#========================================================
# output IMCD solute flows in Lumen and Cell
#========================================================
for i in range(NS):
	file=open('./'+file_to_save+'/'+imcd[0].sex+'imcd_flow_of_'+solute[i]+'_in_Lumen.txt','w')
	for j in range(NIMCD):
		file.write(str(imcd[j].conc[i,0]*imcd[j].vol[0]*cw)+'\n')
	file.close()
for i in range(NS):
	file=open('./'+file_to_save+'/'+imcd[0].sex+'imcd_flow_of_'+solute[i]+'_in_Cell.txt','w')
	for j in range(NIMCD):
		file.write(str(imcd[j].conc[i,1]*imcd[j].vol[1]*cw)+'\n')
	file.close()
#========================================================
# output IMCD osmolality in Lumen, Cell, LIS, Bath
#========================================================
file_lumen = open('./'+file_to_save+'/'+imcd[0].sex+imcd[0].segment+'_osmolality_in_Lumen.txt','w')
file_cell = open('./'+file_to_save+'/'+imcd[0].sex+imcd[0].segment+'_osmolality_in_Cell.txt','w')
file_lis = open('./'+file_to_save+'/'+imcd[0].sex+imcd[0].segment+'_osmolality_in_LIS.txt','w')
file_bath = open('./'+file_to_save+'/'+imcd[0].sex+imcd[0].segment+'_osmolality_in_Bath.txt','w')
for j in range(NIMCD):
	osm_l = 0
	osm_c = 0
	osm_lis = 0
	osm_b = 0
	for i in range(NS):
		osm_l = osm_l +imcd[j].conc[i,0]
		osm_c = osm_c +imcd[j].conc[i,1]
		osm_lis = osm_lis+imcd[j].conc[i,4]
		osm_b = osm_b +imcd[j].conc[i,5]

	file_lumen.write(str(osm_l)+'\n')
	file_cell.write(str(osm_c)+'\n')
	file_lis.write(str(osm_lis)+'\n')
	file_bath.write(str(osm_b)+'\n')
file_lumen.close()
file_cell.close()
file_lis.close()
file_bath.close()
#========================================================
# output luminal pressure
#========================================================
file=open('./'+file_to_save+'/'+imcd[0].sex+imcd[0].segment+'_pressure_in_Lumen_'+sup_or_jux+'.txt','w')
for j in range(NIMCD):
	file.write(str(imcd[j].pres[0])+'\n')
file.close()
#========================================================
# output flux through transporters
#========================================================
for j in range(NIMCD):
	jvol = water.compute_water_fluxes(imcd[j])
	jsol,delmu = electrochemical.compute_ecd_fluxes(imcd[j],jvol)
	for i in range(len(imcd[j].trans)):
		transporter_type = imcd[j].trans[i].type
		memb_id = imcd[j].trans[i].membrane_id

		if transporter_type == 'SGLT1':
			solute_id,fluxs = glucose.sglt1(imcd[j],imcd[j].ep,memb_id,imcd[j].trans[i].act,imcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'SGLT2':
			solute_id,fluxs = glucose.sglt2(imcd[j],imcd[j].ep,memb_id,imcd[j].trans[i].act,imcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'GLUT1':
			solute_id,fluxs=glucose.glut1(imcd[j],imcd[j].ep,memb_id,imcd[j].trans[i].act,imcd[j].area)
			for k in range(len([solute_id])):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'.txt','a')
				file.write(str(fluxs)+'\n')
		elif transporter_type == 'GLUT2':
			solute_id,fluxs=glucose.glut2(imcd[j],imcd[j].ep,memb_id,imcd[j].trans[i].act,imcd[j].area)
			for k in range(len([solute_id])):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'.txt','a')
				file.write(str(fluxs)+'\n')			
		elif transporter_type == 'NHE3':
			solute_id,fluxs=NHE3.nhe3(imcd[j],imcd[j].ep,memb_id,imcd[j].trans[i].act,imcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NaKATPase':
			solute_id,fluxs=ATPase.nakatpase(imcd[j],imcd[j].ep,memb_id,imcd[j].trans[i].act,imcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+str(memb_id[0])+str(memb_id[1])+'.txt','a')
				file.write(str(fluxs[k])+'\n')

		elif transporter_type == 'HATPase':
			solute_id,fluxs=ATPase.hatpase(imcd[j],imcd[j].ep,memb_id,imcd[j].trans[i].act,imcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC2A':
			solute_id,fluxs=NKCC.nkcc2(imcd[j],memb_id,imcd[j].trans[i].act,imcd[j].area,'A')
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC2B':
			solute_id,fluxs=NKCC.nkcc2(imcd[j],memb_id,imcd[j].trans[i].act,imcd[j].area,'B')
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC2F':
			solute_id,fluxs=NKCC.nkcc2(imcd[j],memb_id,imcd[j].trans[i].act,imcd[j].area,'F')
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')       
		elif transporter_type == 'KCC4':
			solute_id,fluxs=KCC.kcc4(imcd[j].conc,memb_id,imcd[j].trans[i].act,imcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'ENaC':
			solute_id,fluxs=ENaC.ENaC(imcd[j],j,memb_id,imcd[j].trans[i].act,imcd[j].area,jvol)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NCC':
			solute_id,fluxs=NCC.NCC(imcd[j],j,memb_id,imcd[j].trans[i].act,imcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'Pendrin':
			solute_id,fluxs=Pendrin.Pendrin(imcd[j],memb_id,imcd[j].trans[i].act,imcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type =='AE1':
			solute_id,fluxs=AE1.AE1(imcd[j],memb_id,imcd[j].trans[i].act,imcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'HKATPase':
			solute_id,fluxs=ATPase.hkatpase(imcd[j],memb_id,imcd[j].trans[i].act,imcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NHE1':
			solute_id,fluxs=NHE1.NHE1(imcd[j],memb_id,imcd[j].trans[i].act,imcd[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC1':
			solute_id,fluxs=NKCC.nkcc1(imcd[j],memb_id,imcd[j].trans[i].act,delmu)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+imcd[j].sex+'_'+imcd[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		else:
			print('What is this?',transporter_type)
print('IMCD finished.')

for i in range(NS):
	file=open('./'+file_to_save+'/'+pt[0].sex+'_con_of_'+solute[i]+'_in_Lumen.txt','w')
	for j in range(NPT):
		file.write(str(pt[j].conc[i,0])+'\n')
	for j in range(NS3):
		file.write(str(s3[j].conc[i,0])+'\n')
	for j in range(NSDL):
		file.write(str(sdl[j].conc[i,0])+'\n')
	for j in range(NmTAL):
		file.write(str(mtal[j].conc[i,0])+'\n')
	for j in range(NcTAL):
		file.write(str(ctal[j].conc[i,0])+'\n')
	for j in range(NDCT):
		file.write(str(dct[j].conc[i,0])+'\n')
	for j in range(NCNT):
		file.write(str(cnt[j].conc[i,0])+'\n')
	for j in range(NCCD):
		file.write(str(ccd[j].conc[i,0])+'\n')
	for j in range(NOMCD):
		file.write(str(omcd[j].conc[i,0])+'\n')
	for j in range(NIMCD):
		file.write(str(imcd[j].conc[i,0])+'\n')
	file.close()

file=open('./'+file_to_save+'/'+pt[0].sex+'_water_volume_in_Lumen.txt','w')
for j in range(NPT):
	file.write(str(pt[j].vol[0]*cw)+'\n')
for j in range(NS3):
	file.write(str(s3[j].vol[0]*cw)+'\n')
for j in range(NSDL):
	file.write(str(sdl[j].vol[0]*cw)+'\n')
for j in range(NmTAL):
	file.write(str(mtal[j].vol[0]*cw)+'\n')
for j in range(NcTAL):
	file.write(str(ctal[j].vol[0]*cw)+'\n')
for j in range(NDCT):
	file.write(str(dct[j].vol[0]*cw)+'\n')
for j in range(NCNT):
	file.write(str(cnt[j].vol[0]*cw)+'\n')
for j in range(NCCD):
	file.write(str(ccd[j].vol[0]*cw)+'\n')
for j in range(NOMCD):
	file.write(str(omcd[j].vol[0]*cw)+'\n')
for j in range(NIMCD):
	file.write(str(imcd[j].vol[0]*cw)+'\n')
file.close()

for i in range(NS):
	file=open('./'+file_to_save+'/'+pt[0].sex+'_flow_of_'+solute[i]+'_in_Lumen.txt','w')
	for j in range(NPT):
		file.write(str(pt[j].conc[i,0]*pt[j].vol[0]*cw)+'\n')
	for j in range(NS3):
		file.write(str(s3[j].conc[i,0]*s3[j].vol[0]*cw)+'\n')
	for j in range(NSDL):
		file.write(str(sdl[j].conc[i,0]*sdl[j].vol[0]*cw)+'\n')
	for j in range(NmTAL):
		file.write(str(mtal[j].conc[i,0]*mtal[j].vol[0]*cw)+'\n')
	for j in range(NcTAL):
		file.write(str(ctal[j].conc[i,0]*ctal[j].vol[0]*cw)+'\n')
	for j in range(NDCT):
		file.write(str(dct[j].conc[i,0]*dct[j].vol[0]*cw)+'\n')
	for j in range(NCNT):
		file.write(str(cnt[j].conc[i,0]*cnt[j].vol[0]*cw)+'\n')
	for j in range(NCCD):
		file.write(str(ccd[j].conc[i,0]*ccd[j].vol[0]*cw)+'\n')
	for j in range(NOMCD):
		file.write(str(omcd[j].conc[i,0]*omcd[j].vol[0]*cw)+'\n')
	for j in range(NIMCD):
		file.write(str(imcd[j].conc[i,0]*imcd[j].vol[0]*cw)+'\n')							
	file.close()

file=open('./'+file_to_save+'/'+pt[0].sex+'_pres_in_Lumen.txt','w')
for j in range(NPT):
	file.write(str(pt[j].pres[0])+'\n')
for j in range(NS3):
	file.write(str(s3[j].pres[0])+'\n')
for j in range(NSDL):
	file.write(str(sdl[j].pres[0])+'\n')
for j in range(NmTAL):
	file.write(str(mtal[j].pres[0])+'\n')
for j in range(NcTAL):
	file.write(str(ctal[j].pres[0])+'\n')
for j in range(NDCT):
	file.write(str(dct[j].pres[0])+'\n')
for j in range(NCNT):
	file.write(str(cnt[j].pres[0])+'\n')
for j in range(NCCD):
	file.write(str(ccd[j].pres[0])+'\n')
for j in range(NOMCD):
	file.write(str(omcd[j].pres[0])+'\n')
for j in range(NIMCD):
	file.write(str(imcd[j].pres[0])+'\n')
file.close()

file=open('./'+file_to_save+'/'+pt[0].sex+'_ph_in_Lumen.txt','w')
for j in range(NPT):
	ph = -np.log(pt[j].conc[11,0]/1000)/np.log(10)
	file.write(str(ph)+'\n')
for j in range(NS3):
	ph = -np.log(s3[j].conc[11,0]/1000)/np.log(10)
	file.write(str(ph)+'\n')
for j in range(NSDL):
	ph = -np.log(sdl[j].conc[11,0]/1000)/np.log(10)
	file.write(str(ph)+'\n')
for j in range(NmTAL):
	ph = -np.log(mtal[j].conc[11,0]/1000)/np.log(10)
	file.write(str(ph)+'\n')
for j in range(NcTAL):
	ph = -np.log(ctal[j].conc[11,0]/1000)/np.log(10)
	file.write(str(ph)+'\n')
for j in range(NDCT):
	ph = -np.log(dct[j].conc[11,0]/1000)/np.log(10)
	file.write(str(ph)+'\n')
for j in range(NCNT):
	ph = -np.log(cnt[j].conc[11,0]/1000)/np.log(10)
	file.write(str(ph)+'\n')
for j in range(NCCD):
	ph = -np.log(ccd[j].conc[11,0]/1000)/np.log(10)
	file.write(str(ph)+'\n')
for j in range(NOMCD):
	ph = -np.log(omcd[j].conc[11,0]/1000)/np.log(10)
	file.write(str(ph)+'\n')
for j in range(NIMCD):
	ph = -np.log(imcd[j].conc[11,0]/1000)/np.log(10)
	file.write(str(ph)+'\n')
file.close()