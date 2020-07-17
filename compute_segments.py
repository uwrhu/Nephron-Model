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

sup_segments=['PT','S3','SDL','mTAL','cTAL','MD','DCT','CNT']
jux_segments=['PT','S3','SDL','LDL','LAL','mTAL','cTAL','MD','DCT','CNT']
types = ['sup','jux1','jux2','jux3','jux4','jux5']

#gender=input('Which gender do you want to simulate? (Male or Female) ')

#diabete=input('Dose it have diabete? (Y/N) ')

#inhib=input('Any transport inhibition? (nhe3 50%/ nhe3 80%/ NKCC2 70%/ NKCC2 100%) ')

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

begin=input('Which segment do you want to start with? ')
end=input('Which segment do you want to stop with? ')
sup_or_jux=input('Is this superfical or juxtamedullary? (sup/jux1/jux2/jux3/jux4/jux5) ')
file_to_save=input('Input a file name to save data: ')

if os.path.isdir(file_to_save) == False:
	os.makedirs(file_to_save)

if sup_or_jux == 'sup':
	segments = sup_segments
else:
	segments = jux_segments

begin_index=segments.index(begin)
end_index=segments.index(end)

for segment in segments[begin_index:end_index+1]:
	if gender == 'Male':
		filename = segment+'params_M.dat'
	elif gender == 'Female':
		filename = segment+'params_F.dat'
	else:
		print('Did you choose correct segments?')

	if segment == 'PT':
		N = 176
	elif segment == 'S3':
		N = 25
	elif segment == 'MD':
		N = 2
	else:
		N = 200

	if gender == 'Male':
		if segment == 'PT' or segment == 'SDL':# or segment == 'mTAL' or segment == 'DCT':
			method = 'Broyden'
		else:
			method = 'Newton'
	elif gender == 'Female':
		if segment == 'PT' or segment == 'SDL' or segment == 'DCT': # for nhe3 50: mTAL can't use broyden
			if (segment == 'DCT' and sup_or_jux == 'sup') or (segment=='DCT' and sup_or_jux == 'jux5'):
				method = 'Newton'
			else:
				method = 'Broyden'
		else:
			method = 'Newton'
	cell = compute(N,filename,method,sup_or_jux,diabete,inhib,perc,salt)
	for i in range(NS):
		file=open('./'+file_to_save+'/'+cell[0].sex+cell[0].segment+'_con_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(N):
			file.write(str(cell[j].conc[i,0])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+cell[0].sex+cell[0].segment+'_con_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(N):
			file.write(str(cell[j].conc[i,1])+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+cell[0].sex+cell[0].segment+'_con_of_'+solute[i]+'_in_Bath_'+sup_or_jux+'.txt','w')
		for j in range(N):
			file.write(str(cell[j].conc[i,5])+'\n')
		file.close()

	file=open('./'+file_to_save+'/'+cell[0].sex+cell[0].segment+'_water_volume_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(N):
		file.write(str(cell[j].vol[0]*cw)+'\n')
	file.close()
	file=open('./'+file_to_save+'/'+cell[0].sex+cell[0].segment+'_water_volume_in_Cell_'+sup_or_jux+'.txt','w')
	for j in range(N):
		file.write(str(cell[j].vol[1]*cw)+'\n')
	file.close()

	for i in range(NS):
		file=open('./'+file_to_save+'/'+cell[0].sex+cell[0].segment+'_flow_of_'+solute[i]+'_in_Lumen_'+sup_or_jux+'.txt','w')
		for j in range(N):
			file.write(str(cell[j].conc[i,0]*cell[j].vol[0]*cw)+'\n')
		file.close()
	for i in range(NS):
		file=open('./'+file_to_save+'/'+cell[0].sex+cell[0].segment+'_flow_of_'+solute[i]+'_in_Cell_'+sup_or_jux+'.txt','w')
		for j in range(N):
			file.write(str(cell[j].conc[i,1]*cell[j].vol[1]*cw)+'\n')
		file.close()

	file=open('./'+file_to_save+'/'+cell[0].sex+cell[0].segment+'_pH_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(N):
		file.write(str(-np.log(cell[j].conc[11,0]/1000)/np.log(10))+'\n')
	file.close()

	file=open('./'+file_to_save+'/'+cell[0].sex+cell[0].segment+'_pressure_in_Lumen_'+sup_or_jux+'.txt','w')
	for j in range(N):
		file.write(str(cell[j].pres[0])+'\n')
	file.close()

	print(segment+' finished')