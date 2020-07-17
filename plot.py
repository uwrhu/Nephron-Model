import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

solutes = ['Na','K','Cl','HCO3','H2CO3','CO2','HPO4','H2PO4','urea','NH3','NH4','H','HCO2','H2CO2','glu']
comparts = ['Lumen','Cell','ICA','ICB','LIS','Bath']
segments_early = ['pt','s3','sdl','ldl','lal','mtal','ctal','dct','cnt']
segments_later = ['ccd','omcd','imcd']
segments_early_cap = ['PT','S3','SDL','LDL','LAL','mTAL','cTAL','DCT','CNT']
segments_later_cap = ['CCD','OMCD','IMCD']
genders = ['male','female']
types = ['sup','jux1','jux2','jux3','jux4','jux5']

parser = argparse.ArgumentParser()
parser.add_argument('--male',required = True, type = str,help='male input file')
parser.add_argument('--female',required = True,type = str,help='female input file')
parser.add_argument('--output',required = True,type = str,help='output file to save plots')

args = parser.parse_args()
male_file = args.male
female_file = args.female
output = args.output

plot_folder = output
comments = input('Any comments to help track changes? ')

os.makedirs(plot_folder)

comments_file = open('./'+plot_folder+'/'+'comments.txt','w')
comments_file.write(comments)
comments_file.close()

#====================================================
# Plot of segmental water flow
#====================================================
for seg in segments_early:
	fig,ax = plt.subplots()
	fig.set_figheight(20)
	fig.set_figwidth(20)

	if seg == 'pt' or seg == 's3' or seg == 'sdl' or seg == 'mtal' or seg == 'ctal' or seg == 'dct' or seg == 'cnt':
		male_sup = open('./'+male_file+'/'+'male'+seg+'_water_volume_in_Lumen_sup.txt','r')
		female_sup = open('./'+female_file+'/'+'female'+seg+'_water_volume_in_Lumen_sup.txt','r')
	
	male_jux1 = open('./'+male_file+'/'+'male'+seg+'_water_volume_in_Lumen_jux1.txt','r')
	female_jux1 = open('./'+female_file+'/'+'female'+seg+'_water_volume_in_Lumen_jux1.txt','r')
	male_jux2 = open('./'+male_file+'/'+'male'+seg+'_water_volume_in_Lumen_jux2.txt','r')
	female_jux2 = open('./'+female_file+'/'+'female'+seg+'_water_volume_in_Lumen_jux2.txt','r')
	male_jux3 = open('./'+male_file+'/'+'male'+seg+'_water_volume_in_Lumen_jux3.txt','r')
	female_jux3 = open('./'+female_file+'/'+'female'+seg+'_water_volume_in_Lumen_jux3.txt','r')
	male_jux4 = open('./'+male_file+'/'+'male'+seg+'_water_volume_in_Lumen_jux4.txt','r')
	female_jux4 = open('./'+female_file+'/'+'female'+seg+'_water_volume_in_Lumen_jux4.txt','r')
	male_jux5 = open('./'+male_file+'/'+'male'+seg+'_water_volume_in_Lumen_jux5.txt','r')
	female_jux5 = open('./'+female_file+'/'+'female'+seg+'_water_volume_in_Lumen_jux5.txt','r')

	if seg == 'pt' or seg == 's3' or seg == 'sdl' or seg == 'mtal' or seg == 'ctal' or seg == 'dct' or seg == 'cnt':
		male_sup_list = []
		female_sup_list = []
		
	male_jux1_list = []
	female_jux1_list = []
	male_jux2_list = []
	female_jux2_list = []
	male_jux3_list = []
	female_jux3_list = []
	male_jux4_list = []
	female_jux4_list = []
	male_jux5_list = []
	female_jux5_list = []

	if seg == 'pt' or seg == 's3' or seg == 'sdl' or seg == 'mtal' or seg == 'ctal' or seg == 'dct' or seg == 'cnt':
		for i in male_sup:
			line = i.split(' ')
			male_sup_list.append(float(line[0]))
		for i in female_sup:
			line = i.split(' ')
			female_sup_list.append(float(line[0]))

	for i in male_jux1:
		line = i.split(' ')
		male_jux1_list.append(float(line[0]))
	for i in female_jux1:
		line = i.split(' ')
		female_jux1_list.append(float(line[0]))

	for i in male_jux2:
		line = i.split(' ')
		male_jux2_list.append(float(line[0]))
	for i in female_jux2:
		line = i.split(' ')
		female_jux2_list.append(float(line[0]))

	for i in male_jux3:
		line = i.split(' ')
		male_jux3_list.append(float(line[0]))
	for i in female_jux3:
		line = i.split(' ')
		female_jux3_list.append(float(line[0]))

	for i in male_jux4:
		line = i.split(' ')
		male_jux4_list.append(float(line[0]))
	for i in female_jux4:
		line = i.split(' ')
		female_jux4_list.append(float(line[0]))

	for i in male_jux5:
		line = i.split(' ')
		male_jux5_list.append(float(line[0]))
	for i in female_jux5:
		line = i.split(' ')
		female_jux5_list.append(float(line[0]))

	if seg == 'pt':
		N = 176
	elif seg == 's3':
		N = 25
	else:
		N = 200
	pos = [i for i in range(N)]

	if seg == 'pt' or seg == 's3' or seg == 'sdl' or seg == 'mtal' or seg == 'ctal' or seg == 'dct' or seg == 'cnt':
		ax.plot(pos,male_sup_list,'-',color = 'black',label = 'Male sup',linewidth = 2)
		ax.plot(pos,female_sup_list,'--',dashes = (10,10),color = 'black',label = 'Female sup',linewidth = 2)
	
	# ax.plot(pos,male_jux1_list,'-',color = 'red',label = 'Male jux1',linewidth = 2)
	# ax.plot(pos,female_jux1_list,'--',dashes = (10,10),color = 'red',label = 'Female jux1',linewidth = 2)
	# ax.plot(pos,male_jux2_list,'-',color = 'blue',label = 'Male jux2',linewidth = 2)
	# ax.plot(pos,female_jux2_list,'--',dashes = (10,10),color = 'blue',label = 'Female jux2',linewidth = 2)
	# ax.plot(pos,male_jux3_list,'-',color = 'green',label = 'Male jux3',linewidth = 2)
	# ax.plot(pos,female_jux3_list,'--',dashes = (10,10),color = 'green',label = 'Female jux3',linewidth = 2)
	# ax.plot(pos,male_jux4_list,'-',color = 'purple',label = 'Male jux4',linewidth = 2)
	# ax.plot(pos,female_jux4_list,'--',dashes = (10,10),color = 'purple',label = 'Female jux4',linewidth = 2)
	ax.plot(pos,male_jux5_list,'-',color = 'yellow',label = 'Male jux5',linewidth = 2)
	ax.plot(pos,female_jux5_list,'--',dashes = (10,10),color = 'yellow',label = 'Female jux5',linewidth = 2)

	ax.legend(fontsize = 30,markerscale = 5)
	ax.tick_params(labelsize = 30)
	ax.set_ylabel('Volume flow (nl/min)',fontsize = 40)
	ax.get_xaxis().set_visible(False)

	plt.savefig('./'+plot_folder+'/'+seg+' volume flow male vs female', bbox_inches = 'tight')
	plt.close()

for seg in segments_later:
	fig,ax = plt.subplots()
	fig.set_figheight(20)
	fig.set_figwidth(20)

	male = open('./'+male_file+'/'+'male'+seg+'_water_volume_in_Lumen.txt','r')
	female = open('./'+female_file+'/'+'female'+seg+'_water_volume_in_Lumen.txt','r')

	male_list = []
	female_list = []

	for i in male:
		line = i.split(' ')
		male_list.append(float(line[0]))
	for i in female:
		line = i.split(' ')
		female_list.append(float(line[0]))

	if seg == 'pt':
		N = 176
	elif seg == 's3':
		N = 25
	else:
		N = 200
	pos = [i for i in range(N)]

	ax.plot(pos,male_list,'-',color = 'black',label = 'Male',linewidth = 2)
	ax.plot(pos,female_list,'--',dashes = (10,10),color = 'black',label = 'Female',linewidth = 2)

	ax.legend(fontsize = 30,markerscale = 5)
	ax.tick_params(labelsize = 30)
	ax.set_ylabel('Volume flow (nl/min)',fontsize = 40)
	ax.get_xaxis().set_visible(False)

	plt.savefig('./'+plot_folder+'/'+seg+' volume flow male vs female', bbox_inches = 'tight')
	plt.close()

#====================================================
# Plot of segmental solutes flow
#====================================================
for solute in solutes:
	for seg in segments_early:
		fig,ax = plt.subplots()
		fig.set_figheight(20)
		fig.set_figwidth(20)

		if seg == 'pt' or seg == 's3' or seg == 'sdl' or seg == 'mtal' or seg == 'ctal' or seg == 'dct' or seg == 'cnt':
			male_sup = open('./'+male_file+'/'+'male'+seg+'_flow_of_'+solute+'_in_Lumen_sup.txt','r')
			female_sup = open('./'+female_file+'/'+'female'+seg+'_flow_of_'+solute+'_in_Lumen_sup.txt','r')
		
		male_jux1 = open('./'+male_file+'/'+'male'+seg+'_flow_of_'+solute+'_in_Lumen_jux1.txt','r')
		female_jux1 = open('./'+female_file+'/'+'female'+seg+'_flow_of_'+solute+'_in_Lumen_jux1.txt','r')
		male_jux2 = open('./'+male_file+'/'+'male'+seg+'_flow_of_'+solute+'_in_Lumen_jux2.txt','r')
		female_jux2 = open('./'+female_file+'/'+'female'+seg+'_flow_of_'+solute+'_in_Lumen_jux2.txt','r')
		male_jux3 = open('./'+male_file+'/'+'male'+seg+'_flow_of_'+solute+'_in_Lumen_jux3.txt','r')
		female_jux3 = open('./'+female_file+'/'+'female'+seg+'_flow_of_'+solute+'_in_Lumen_jux3.txt','r')
		male_jux4 = open('./'+male_file+'/'+'male'+seg+'_flow_of_'+solute+'_in_Lumen_jux4.txt','r')
		female_jux4 = open('./'+female_file+'/'+'female'+seg+'_flow_of_'+solute+'_in_Lumen_jux4.txt','r')
		male_jux5 = open('./'+male_file+'/'+'male'+seg+'_flow_of_'+solute+'_in_Lumen_jux5.txt','r')
		female_jux5 = open('./'+female_file+'/'+'female'+seg+'_flow_of_'+solute+'_in_Lumen_jux5.txt','r')

		if seg == 'pt' or seg == 's3' or seg == 'sdl' or seg == 'mtal' or seg == 'ctal' or seg == 'dct' or seg == 'cnt':
			male_sup_list = []
			female_sup_list = []
	
		male_jux1_list = []
		female_jux1_list = []
		male_jux2_list = []
		female_jux2_list = []
		male_jux3_list = []
		female_jux3_list = []
		male_jux4_list = []
		female_jux4_list = []
		male_jux5_list = []
		female_jux5_list = []

		if seg == 'pt' or seg == 's3' or seg == 'sdl' or seg == 'mtal' or seg == 'ctal' or seg == 'dct' or seg == 'cnt':
			for i in male_sup:
				line = i.split(' ')
				male_sup_list.append(float(line[0]))
			for i in female_sup:
				line = i.split(' ')
				female_sup_list.append(float(line[0]))

		for i in male_jux1:
			line = i.split(' ')
			male_jux1_list.append(float(line[0]))
		for i in female_jux1:
			line = i.split(' ')
			female_jux1_list.append(float(line[0]))

		for i in male_jux2:
			line = i.split(' ')
			male_jux2_list.append(float(line[0]))
		for i in female_jux2:
			line = i.split(' ')
			female_jux2_list.append(float(line[0]))

		for i in male_jux3:
			line = i.split(' ')
			male_jux3_list.append(float(line[0]))
		for i in female_jux3:
			line = i.split(' ')
			female_jux3_list.append(float(line[0]))

		for i in male_jux4:
			line = i.split(' ')
			male_jux4_list.append(float(line[0]))
		for i in female_jux4:
			line = i.split(' ')
			female_jux4_list.append(float(line[0]))

		for i in male_jux5:
			line = i.split(' ')
			male_jux5_list.append(float(line[0]))
		for i in female_jux5:
			line = i.split(' ')
			female_jux5_list.append(float(line[0]))

		if seg == 'pt':
			N = 176
		elif seg == 's3':
			N = 25
		else:
			N = 200
		pos = [i for i in range(N)]

		if seg == 'pt' or seg == 's3' or seg == 'sdl' or seg == 'mtal' or seg == 'ctal' or seg == 'dct' or seg == 'cnt':
			ax.plot(pos,male_sup_list,'-',color = 'black',label = 'Male sup',linewidth = 2)
			ax.plot(pos,female_sup_list,'--',dashes = (10,10),color = 'black',label = 'Female sup',linewidth = 2)
	
		# ax.plot(pos,male_jux1_list,'-',color = 'red',label = 'Male jux1',linewidth = 2)
		# ax.plot(pos,female_jux1_list,'--',dashes = (10,10),color = 'red',label = 'Female jux1',linewidth = 2)
		# ax.plot(pos,male_jux2_list,'-',color = 'blue',label = 'Male jux2',linewidth = 2)
		# ax.plot(pos,female_jux2_list,'--',dashes = (10,10),color = 'blue',label = 'Female jux2',linewidth = 2)
		# ax.plot(pos,male_jux3_list,'-',color = 'green',label = 'Male jux3',linewidth = 2)
		# ax.plot(pos,female_jux3_list,'--',dashes = (10,10),color = 'green',label = 'Female jux3',linewidth = 2)
		# ax.plot(pos,male_jux4_list,'-',color = 'purple',label = 'Male jux4',linewidth = 2)
		# ax.plot(pos,female_jux4_list,'--',dashes = (10,10),color = 'purple',label = 'Female jux4',linewidth = 2)
		ax.plot(pos,male_jux5_list,'-',color = 'yellow',label = 'Male jux5',linewidth = 2)
		ax.plot(pos,female_jux5_list,'--',dashes = (10,10),color = 'yellow',label = 'Female jux5',linewidth = 2)

		ax.legend(fontsize = 30,markerscale = 5)
		ax.tick_params(labelsize = 30)
		ax.set_ylabel(solute+' flow in Lumen along '+seg+' (nl/min)',fontsize = 40)
		ax.get_xaxis().set_visible(False)

		plt.savefig('./'+plot_folder+'/'+seg+' '+solute+' flow male vs female', bbox_inches = 'tight')
		plt.close()

	for seg in segments_later:
		fig,ax = plt.subplots()
		fig.set_figheight(20)
		fig.set_figwidth(20)

		male = open('./'+male_file+'/'+'male'+seg+'_flow_of_'+solute+'_in_Lumen.txt','r')
		female = open('./'+female_file+'/'+'female'+seg+'_flow_of_'+solute+'_in_Lumen.txt','r')

		male_list = []
		female_list = []

		for i in male:
			line = i.split(' ')
			male_list.append(float(line[0]))
		for i in female:
			line = i.split(' ')
			female_list.append(float(line[0]))

		if seg == 'pt':
			N = 176
		elif seg == 's3':
			N = 25
		else:
			N = 200
		pos = [i for i in range(N)]

		ax.plot(pos,male_list,'-',color = 'black',label = 'Male',linewidth = 2)
		ax.plot(pos,female_list,'--',dashes = (10,10),color = 'black',label = 'Female',linewidth = 2)

		ax.legend(fontsize = 30,markerscale = 5)
		ax.tick_params(labelsize = 30)
		ax.set_ylabel(solute+' flow in Lumen along '+seg+' (nl/min)',fontsize = 40)
		ax.get_xaxis().set_visible(False)

		plt.savefig('./'+plot_folder+'/'+seg+' '+solute+' flow male vs female', bbox_inches = 'tight')
		plt.close()

#====================================================
# Plot of segmental solutes concentration
#====================================================
for solute in solutes:
	for seg in segments_early:
		fig,ax = plt.subplots()
		fig.set_figheight(20)
		fig.set_figwidth(20)

		if seg == 'pt' or seg == 's3' or seg == 'sdl' or seg == 'mtal' or seg == 'ctal' or seg == 'dct' or seg == 'cnt':
			male_sup = open('./'+male_file+'/'+'male'+seg+'_con_of_'+solute+'_in_Lumen_sup.txt','r')
			female_sup = open('./'+female_file+'/'+'female'+seg+'_con_of_'+solute+'_in_Lumen_sup.txt','r')

		male_jux1 = open('./'+male_file+'/'+'male'+seg+'_con_of_'+solute+'_in_Lumen_jux1.txt','r')
		female_jux1 = open('./'+female_file+'/'+'female'+seg+'_con_of_'+solute+'_in_Lumen_jux1.txt','r')
		male_jux2 = open('./'+male_file+'/'+'male'+seg+'_con_of_'+solute+'_in_Lumen_jux2.txt','r')
		female_jux2 = open('./'+female_file+'/'+'female'+seg+'_con_of_'+solute+'_in_Lumen_jux2.txt','r')
		male_jux3 = open('./'+male_file+'/'+'male'+seg+'_con_of_'+solute+'_in_Lumen_jux3.txt','r')
		female_jux3 = open('./'+female_file+'/'+'female'+seg+'_con_of_'+solute+'_in_Lumen_jux3.txt','r')
		male_jux4 = open('./'+male_file+'/'+'male'+seg+'_con_of_'+solute+'_in_Lumen_jux4.txt','r')
		female_jux4 = open('./'+female_file+'/'+'female'+seg+'_con_of_'+solute+'_in_Lumen_jux4.txt','r')
		male_jux5 = open('./'+male_file+'/'+'male'+seg+'_con_of_'+solute+'_in_Lumen_jux5.txt','r')
		female_jux5 = open('./'+female_file+'/'+'female'+seg+'_con_of_'+solute+'_in_Lumen_jux5.txt','r')

		if seg == 'pt' or seg == 's3' or seg == 'sdl' or seg == 'mtal' or seg == 'ctal' or seg == 'dct' or seg == 'cnt':
			male_sup_list = []
			female_sup_list = []

		male_jux1_list = []
		female_jux1_list = []
		male_jux2_list = []
		female_jux2_list = []
		male_jux3_list = []
		female_jux3_list = []
		male_jux4_list = []
		female_jux4_list = []
		male_jux5_list = []
		female_jux5_list = []

		if seg == 'pt' or seg == 's3' or seg == 'sdl' or seg == 'mtal' or seg == 'ctal' or seg == 'dct' or seg == 'cnt':
			for i in male_sup:
				line = i.split(' ')
				male_sup_list.append(float(line[0]))
			for i in female_sup:
				line = i.split(' ')
				female_sup_list.append(float(line[0]))

		for i in male_jux1:
			line = i.split(' ')
			male_jux1_list.append(float(line[0]))
		for i in female_jux1:
			line = i.split(' ')
			female_jux1_list.append(float(line[0]))

		for i in male_jux2:
			line = i.split(' ')
			male_jux2_list.append(float(line[0]))
		for i in female_jux2:
			line = i.split(' ')
			female_jux2_list.append(float(line[0]))

		for i in male_jux3:
			line = i.split(' ')
			male_jux3_list.append(float(line[0]))
		for i in female_jux3:
			line = i.split(' ')
			female_jux3_list.append(float(line[0]))

		for i in male_jux4:
			line = i.split(' ')
			male_jux4_list.append(float(line[0]))
		for i in female_jux4:
			line = i.split(' ')
			female_jux4_list.append(float(line[0]))

		for i in male_jux5:
			line = i.split(' ')
			male_jux5_list.append(float(line[0]))
		for i in female_jux5:
			line = i.split(' ')
			female_jux5_list.append(float(line[0]))

		if seg == 'pt':
			N = 176
		elif seg == 's3':
			N = 25
		else:
			N = 200
		pos = [i for i in range(N)]

		if seg == 'pt' or seg == 's3' or seg == 'sdl' or seg == 'mtal' or seg == 'ctal' or seg == 'dct' or seg == 'cnt':
			ax.plot(pos,male_sup_list,'-',color = 'black',label = 'Male sup',linewidth = 2)
			ax.plot(pos,female_sup_list,'--',dashes = (10,10),color = 'black',label = 'Female sup',linewidth = 2)

		# ax.plot(pos,male_jux1_list,'-',color = 'red',label = 'Male jux1',linewidth = 2)
		# ax.plot(pos,female_jux1_list,'--',dashes = (10,10),color = 'red',label = 'Female jux1',linewidth = 2)
		# ax.plot(pos,male_jux2_list,'-',color = 'blue',label = 'Male jux2',linewidth = 2)
		# ax.plot(pos,female_jux2_list,'--',dashes = (10,10),color = 'blue',label = 'Female jux2',linewidth = 2)
		# ax.plot(pos,male_jux3_list,'-',color = 'green',label = 'Male jux3',linewidth = 2)
		# ax.plot(pos,female_jux3_list,'--',dashes = (10,10),color = 'green',label = 'Female jux3',linewidth = 2)
		# ax.plot(pos,male_jux4_list,'-',color = 'purple',label = 'Male jux4',linewidth = 2)
		# ax.plot(pos,female_jux4_list,'--',dashes = (10,10),color = 'purple',label = 'Female jux4',linewidth = 2)
		ax.plot(pos,male_jux5_list,'-',color = 'yellow',label = 'Male jux5',linewidth = 2)
		ax.plot(pos,female_jux5_list,'--',dashes = (10,10),color = 'yellow',label = 'Female jux5',linewidth = 2)

		ax.legend(fontsize = 30,markerscale = 5)
		ax.tick_params(labelsize = 30)
		ax.set_ylabel('Concentration of '+solute+' in Lumen along '+seg+' (mM)',fontsize = 40)
		ax.get_xaxis().set_visible(False)

		plt.savefig('./'+plot_folder+'/'+seg+' '+solute+' concentration male vs female', bbox_inches = 'tight')
		plt.close()

	for seg in segments_later:
		fig,ax = plt.subplots()
		fig.set_figheight(20)
		fig.set_figwidth(20)

		male = open('./'+male_file+'/'+'male'+seg+'_con_of_'+solute+'_in_Lumen.txt','r')
		female = open('./'+female_file+'/'+'female'+seg+'_con_of_'+solute+'_in_Lumen.txt','r')

		male_list = []
		female_list = []

		for i in male:
			line = i.split(' ')
			male_list.append(float(line[0]))
		for i in female:
			line = i.split(' ')
			female_list.append(float(line[0]))

		if seg == 'pt':
			N = 176
		elif seg == 's3':
			N = 25
		else:
			N = 200
		pos = [i for i in range(N)]

		ax.plot(pos,male_list,'-',color = 'black',label = 'Male',linewidth = 2)
		ax.plot(pos,female_list,'--',dashes = (10,10),color = 'black',label = 'Female',linewidth = 2)

		ax.legend(fontsize = 30,markerscale = 5)
		ax.tick_params(labelsize = 30)
		ax.set_ylabel('Concentration of '+solute+' in Lumen along '+seg+' (mM)',fontsize = 40)
		ax.get_xaxis().set_visible(False)

		plt.savefig('./'+plot_folder+'/'+seg+' '+solute+' concentration male vs female', bbox_inches = 'tight')
		plt.close()

#====================================================
# Plot of segmental osmolality
#====================================================
for seg in segments_early_cap:
	if seg == 'LDL' or seg == 'LAL':
		types = ['jux1','jux2','jux3','jux4','jux5']
	else:
		types = ['sup','jux1','jux2','jux3','jux4','jux5']

	for sup_jux in types:
		fig,ax = plt.subplots()
		fig.set_figheight(20)
		fig.set_figwidth(20)

		male_lumen = open('./'+male_file+'/'+'male'+seg+'_osmolality_in_Lumen_'+sup_jux+'.txt','r')
		male_cell = open('./'+male_file+'/'+'male'+seg+'_osmolality_in_Cell_'+sup_jux+'.txt','r')
		male_lis = open('./'+male_file+'/'+'male'+seg+'_osmolality_in_LIS_'+sup_jux+'.txt','r')
		male_bath = open('./'+male_file+'/'+'male'+seg+'_osmolality_in_Bath_'+sup_jux+'.txt','r')
		female_lumen = open('./'+female_file+'/'+'female'+seg+'_osmolality_in_Lumen_'+sup_jux+'.txt','r')
		female_cell = open('./'+female_file+'/'+'female'+seg+'_osmolality_in_Cell_'+sup_jux+'.txt','r')
		female_lis = open('./'+female_file+'/'+'female'+seg+'_osmolality_in_LIS_'+sup_jux+'.txt','r')
		female_bath = open('./'+female_file+'/'+'female'+seg+'_osmolality_in_Bath_'+sup_jux+'.txt','r')


		male_lumen_list = []
		male_cell_list = []
		male_lis_list = []
		male_bath_list = []
		female_lumen_list = []
		female_cell_list = []
		female_lis_list = []
		female_bath_list = []

		for i in male_lumen:
			line = i.split(' ')
			male_lumen_list.append(float(line[0]))

		for i in male_cell:
			line = i.split(' ')
			male_cell_list.append(float(line[0]))

		for i in male_lis:
			line = i.split(' ')
			male_lis_list.append(float(line[0]))

		for i in male_bath:
			line = i.split(' ')
			male_bath_list.append(float(line[0]))

		for i in female_lumen:
			line = i.split(' ')
			female_lumen_list.append(float(line[0]))

		for i in female_cell:
			line = i.split(' ')
			female_cell_list.append(float(line[0]))

		for i in female_lis:
			line = i.split(' ')
			female_lis_list.append(float(line[0]))

		for i in female_bath:
			line = i.split(' ')
			female_bath_list.append(float(line[0]))

		if seg == 'PT':
			N = 176
		elif seg == 'S3':
			N = 25
		else:
			N = 200
		pos = [i for i in range(N)]
		ax.plot(pos,male_lumen_list,'-',color='black',label='Male lumen',linewidth=2)
		ax.plot(pos,male_cell_list,'-',color='red',label='Male cell',linewidth=2)
		ax.plot(pos,male_lis_list,'-',color='green',label='Male lis',linewidth=2)
		ax.plot(pos,male_bath_list,'-',color='blue',label='Male bath',linewidth=2)

		ax.plot(pos,female_lumen_list,'--',dashes=(10,10),color='black',label='Female lumen',linewidth=2)
		ax.plot(pos,female_cell_list,'--',dashes=(10,10),color='red',label='Female cell',linewidth=2)
		ax.plot(pos,female_lis_list,'--',dashes=(10,10),color='green',label='Female lis',linewidth=2)
		ax.plot(pos,female_bath_list,'--',dashes=(10,10),color='blue',label='Female bath',linewidth=2)

		ax.legend(fontsize=30,markerscale=5)
		ax.tick_params(labelsize=30)
		ax.set_ylabel('Osmolality (mM)',fontsize=40)
		ax.get_xaxis().set_visible(False)

		plt.savefig('./'+plot_folder+'/'+seg+' osmolality '+sup_jux,bbox_inches='tight')
		plt.close()

for seg in segments_later_cap:
	fig,ax = plt.subplots()
	fig.set_figheight(20)
	fig.set_figwidth(20)

	male_lumen = open('./'+male_file+'/'+'male'+seg+'_osmolality_in_Lumen.txt','r')
	male_cell = open('./'+male_file+'/'+'male'+seg+'_osmolality_in_Cell.txt','r')
	male_lis = open('./'+male_file+'/'+'male'+seg+'_osmolality_in_LIS.txt','r')
	male_bath = open('./'+male_file+'/'+'male'+seg+'_osmolality_in_Bath.txt','r')
	female_lumen = open('./'+female_file+'/'+'female'+seg+'_osmolality_in_Lumen.txt','r')
	female_cell = open('./'+female_file+'/'+'female'+seg+'_osmolality_in_Cell.txt','r')
	female_lis = open('./'+female_file+'/'+'female'+seg+'_osmolality_in_LIS.txt','r')
	female_bath = open('./'+female_file+'/'+'female'+seg+'_osmolality_in_Bath.txt','r')


	male_lumen_list = []
	male_cell_list = []
	male_lis_list = []
	male_bath_list = []
	female_lumen_list = []
	female_cell_list = []
	female_lis_list = []
	female_bath_list = []

	for i in male_lumen:
		line = i.split(' ')
		male_lumen_list.append(float(line[0]))

	for i in male_cell:
		line = i.split(' ')
		male_cell_list.append(float(line[0]))

	for i in male_lis:
		line = i.split(' ')
		male_lis_list.append(float(line[0]))

	for i in male_bath:
		line = i.split(' ')
		male_bath_list.append(float(line[0]))

	for i in female_lumen:
		line = i.split(' ')
		female_lumen_list.append(float(line[0]))

	for i in female_cell:
		line = i.split(' ')
		female_cell_list.append(float(line[0]))

	for i in female_lis:
		line = i.split(' ')
		female_lis_list.append(float(line[0]))

	for i in female_bath:
		line = i.split(' ')
		female_bath_list.append(float(line[0]))

	if seg == 'PT':
		N = 176
	elif seg == 'S3':
		N = 25
	else:
		N = 200
	pos = [i for i in range(N)]
	ax.plot(pos,male_lumen_list,'-',color='black',label='Male lumen',linewidth=2)
	ax.plot(pos,male_cell_list,'-',color='red',label='Male cell',linewidth=2)
	ax.plot(pos,male_lis_list,'-',color='green',label='Male lis',linewidth=2)
	ax.plot(pos,male_bath_list,'-',color='blue',label='Male bath',linewidth=2)

	ax.plot(pos,female_lumen_list,'--',dashes=(10,10),color='black',label='Female lumen',linewidth=2)
	ax.plot(pos,female_cell_list,'--',dashes=(10,10),color='red',label='Female cell',linewidth=2)
	ax.plot(pos,female_lis_list,'--',dashes=(10,10),color='green',label='Female lis',linewidth=2)
	ax.plot(pos,female_bath_list,'--',dashes=(10,10),color='blue',label='Female bath',linewidth=2)

	ax.legend(fontsize=30,markerscale=5)
	ax.tick_params(labelsize=30)
	ax.set_ylabel('Osmolality (mM)',fontsize=40)
	ax.get_xaxis().set_visible(False)

	plt.savefig('./'+plot_folder+'/'+seg+' osmolality '+sup_jux,bbox_inches='tight')
	plt.close()

#====================================================
# Plot of segmental luminal osmolality
#====================================================
for seg in segments_early_cap:
	if seg == 'LDL' or seg == 'LAL':
		types = ['jux1','jux2','jux3','jux4','jux5']
	else:
		types = ['sup','jux1','jux2','jux3','jux4','jux5']

	
	fig,ax = plt.subplots()
	fig.set_figheight(20)
	fig.set_figwidth(20)
	for sup_jux in types:
		male_lumen = open('./'+male_file+'/'+'male'+seg+'_osmolality_in_Lumen_'+sup_jux+'.txt','r')
		female_lumen = open('./'+female_file+'/'+'female'+seg+'_osmolality_in_Lumen_'+sup_jux+'.txt','r')


		male_lumen_list = []
		female_lumen_list = []


		for i in male_lumen:
			line = i.split(' ')
			male_lumen_list.append(float(line[0]))


		for i in female_lumen:
			line = i.split(' ')
			female_lumen_list.append(float(line[0]))

		if seg == 'PT':
			N = 176
		elif seg == 'S3':
			N = 25
		else:
			N = 200

		colors = {'sup':'black','jux1':'red','jux2':'blue','jux3':'green','jux4':'purple','jux5':'yellow'}

		pos = [i for i in range(N)]
		ax.plot(pos,male_lumen_list,'-',color=colors[sup_jux],label='Male lumen '+sup_jux,linewidth=2)

		ax.plot(pos,female_lumen_list,'--',dashes=(10,10),color=colors[sup_jux],label='Female lumen '+sup_jux,linewidth=2)

	ax.legend(fontsize=30,markerscale=5)
	ax.tick_params(labelsize=30)
	ax.set_ylabel('Osmolality (mM)',fontsize=40)
	ax.get_xaxis().set_visible(False)

	plt.savefig('./'+plot_folder+'/'+seg+' luminal osmolality',bbox_inches='tight')
	plt.close()

for seg in segments_later_cap:
	fig,ax = plt.subplots()
	fig.set_figheight(20)
	fig.set_figwidth(20)

	male_lumen = open('./'+male_file+'/'+'male'+seg+'_osmolality_in_Lumen.txt','r')
	female_lumen = open('./'+female_file+'/'+'female'+seg+'_osmolality_in_Lumen.txt','r')


	male_lumen_list = []
	female_lumen_list = []

	for i in male_lumen:
		line = i.split(' ')
		male_lumen_list.append(float(line[0]))

	for i in female_lumen:
		line = i.split(' ')
		female_lumen_list.append(float(line[0]))

	if seg == 'PT':
		N = 176
	elif seg == 'S3':
		N = 25
	else:
		N = 200
	pos = [i for i in range(N)]
	ax.plot(pos,male_lumen_list,'-',color='black',label='Male lumen',linewidth=2)
	ax.plot(pos,female_lumen_list,'--',dashes=(10,10),color='black',label='Female lumen',linewidth=2)

	ax.legend(fontsize=30,markerscale=5)
	ax.tick_params(labelsize=30)
	ax.set_ylabel('Osmolality (mM)',fontsize=40)
	ax.get_xaxis().set_visible(False)

	plt.savefig('./'+plot_folder+'/'+seg+' luminal osmolality',bbox_inches='tight')
	plt.close()

#====================================================
# Plot of osmolality along loop of henle
#====================================================
loop_of_henle_sup = ['sdl','mtal','ctal']
loop_of_henle_jux = ['sdl','ldl','lal','mtal','ctal']

pos_sdl_sup_male = [1.1+0.14*i/199 for i in range(1,200)]
pos_mtal_sup_male = [1.1+0.14+0.2*i/199 for i in range(1,200)]
pos_ctal_sup_male = [1.1+0.14+0.2+0.2*i/199 for i in range(1,200)]

pos_sdl_sup_female = [0.88+0.119*i/199 for i in range(1,200)]
pos_mtal_sup_female = [0.88+0.119+0.17*i/199 for i in range(1,200)]
pos_ctal_sup_female = [0.88+0.119+0.17+0.17*i/199 for i in range(1,200)]

pos_sdl_jux_male = [1.1+0.14*i/199 for i in range(1,200)]
pos_ldl_jux_male = [1.1+0.14+0.2*0.5*i/199 for i in range(1,200)]
pos_lal_jux_male = [1.1+0.14+0.2*0.5+0.2*0.5*i/199 for i in range(1,200)]
pos_mtal_jux_male = [1.1+0.14+0.2*0.5+0.2*0.5+0.2*i/199 for i in range(1,200)]
pos_ctal_jux_male = [1.1+0.14+0.2*0.5+0.2*0.5+0.2+0.05*i/199 for i in range(1,200)]

pos_sdl_jux_female = [0.88+0.119*i/199 for i in range(1,200)]
pos_ldl_jux_female = [0.88+0.119+0.2*0.425*i/199 for i in range(1,200)]
pos_lal_jux_female = [0.88+0.119+0.2*0.425+0.2*0.425*i/199 for i in range(1,200)]
pos_mtal_jux_female = [0.88+0.119+0.2*0.425+0.2*0.425+0.17*i/199 for i in range(1,200)]
pos_ctal_jux_female = [0.88+0.119+0.2*0.425+0.2*0.425+0.17+0.0425*i/199 for i in range(1,200)]

sdl_male_sup = open('./'+male_file+'/'+'maleSDL_osmolality_in_Lumen_sup.txt','r')
mtal_male_sup = open('./'+male_file+'/'+'malemTAL_osmolality_in_Lumen_sup.txt','r')
ctal_male_sup = open('./'+male_file+'/'+'malecTAL_osmolality_in_Lumen_sup.txt','r')

sdl_female_sup = open('./'+female_file+'/'+'femaleSDL_osmolality_in_Lumen_sup.txt','r')
mtal_female_sup = open('./'+female_file+'/'+'femalemTAL_osmolality_in_Lumen_sup.txt','r')
ctal_female_sup = open('./'+female_file+'/'+'femalecTAL_osmolality_in_Lumen_sup.txt','r')

sdl_male_jux = open('./'+male_file+'/'+'maleSDL_osmolality_in_Lumen_jux1.txt','r')
ldl_male_jux = open('./'+male_file+'/'+'maleLDL_osmolality_in_Lumen_jux1.txt','r')
lal_male_jux = open('./'+male_file+'/'+'maleLAL_osmolality_in_Lumen_jux1.txt','r')
mtal_male_jux = open('./'+male_file+'/'+'malemTAL_osmolality_in_Lumen_jux1.txt','r')
ctal_male_jux = open('./'+male_file+'/'+'malecTAL_osmolality_in_Lumen_jux1.txt','r')

sdl_female_jux = open('./'+female_file+'/'+'femaleSDL_osmolality_in_Lumen_jux1.txt','r')
ldl_female_jux = open('./'+female_file+'/'+'femaleLDL_osmolality_in_Lumen_jux1.txt','r')
lal_female_jux = open('./'+female_file+'/'+'femaleLAL_osmolality_in_Lumen_jux1.txt','r')
mtal_female_jux = open('./'+female_file+'/'+'femalemTAL_osmolality_in_Lumen_jux1.txt','r')
ctal_female_jux = open('./'+female_file+'/'+'femalecTAL_osmolality_in_Lumen_jux1.txt','r')

sdl_male_sup_list = []
mtal_male_sup_list = []
ctal_male_sup_list = []

sdl_female_sup_list = []
mtal_female_sup_list = []
ctal_female_sup_list = []

sdl_male_jux_list = []
ldl_male_jux_list = []
lal_male_jux_list = []
mtal_male_jux_list = []
ctal_male_jux_list = []

sdl_female_jux_list = []
ldl_female_jux_list = []
lal_female_jux_list = []
mtal_female_jux_list = []
ctal_female_jux_list = []

for i in sdl_male_sup:
	line = i.split(' ')
	sdl_male_sup_list.append(float(line[0]))

for i in mtal_male_sup:
	line = i.split(' ')
	mtal_male_sup_list.append(float(line[0]))

for i in ctal_male_sup:
	line = i.split(' ')
	ctal_male_sup_list.append(float(line[0]))

for i in sdl_female_sup:
	line = i.split(' ')
	sdl_female_sup_list.append(float(line[0]))

for i in mtal_female_sup:
	line = i.split(' ')
	mtal_female_sup_list.append(float(line[0]))

for i in ctal_female_sup:
	line = i.split(' ')
	ctal_female_sup_list.append(float(line[0]))

for i in sdl_male_jux:
	line = i.split(' ')
	sdl_male_jux_list.append(float(line[0]))

for i in ldl_male_jux:
	line = i.split(' ')
	ldl_male_jux_list.append(float(line[0]))

for i in lal_male_jux:
	line = i.split(' ')
	lal_male_jux_list.append(float(line[0]))

for i in mtal_male_jux:
	line = i.split(' ')
	mtal_male_jux_list.append(float(line[0]))

for i in ctal_male_jux:
	line = i.split(' ')
	ctal_male_jux_list.append(float(line[0]))

for i in sdl_female_jux:
	line = i.split(' ')
	sdl_female_jux_list.append(float(line[0]))

for i in ldl_female_jux:
	line = i.split(' ')
	ldl_female_jux_list.append(float(line[0]))

for i in lal_female_jux:
	line = i.split(' ')
	lal_female_jux_list.append(float(line[0]))

for i in mtal_female_jux:
	line = i.split(' ')
	mtal_female_jux_list.append(float(line[0]))

for i in ctal_female_jux:
	line = i.split(' ')
	ctal_female_jux_list.append(float(line[0]))

fig,ax = plt.subplots()
fig.set_figheight(20)
fig.set_figwidth(20)

ax.plot(pos_sdl_sup_male,sdl_male_sup_list[1:],'-',color='red',label='Male sup sdl',linewidth=2)
ax.plot(pos_mtal_sup_male,mtal_male_sup_list[1:],'-',color='blue',label='Male sup mtal',linewidth=2)
ax.plot(pos_ctal_sup_male,ctal_male_sup_list[1:],'-',color='green',label='Male sup ctal',linewidth=2)

ax.plot(pos_sdl_sup_female,sdl_female_sup_list[1:],'--',dashes=(10,10),color='red',label='Female sup sdl',linewidth=2)
ax.plot(pos_mtal_sup_female,mtal_female_sup_list[1:],'--',dashes=(10,10),color='blue',label='Female sup mtal',linewidth=2)
ax.plot(pos_ctal_sup_female,ctal_female_sup_list[1:],'--',dashes=(10,10),color='green',label='Female sup ctal',linewidth=2)

ax.legend(fontsize=30,markerscale=5)
ax.tick_params(labelsize=30)
ax.set_ylabel('Osmolality (mM)',fontsize=40)
#ax.get_xaxis().set_visible(False)

plt.savefig('./'+plot_folder+'/'+'Luminal osmolality along superficial loop of henle',bbox_inches='tight')
plt.close()

fig,ax = plt.subplots()
fig.set_figheight(20)
fig.set_figwidth(20)

#print(len(sdl_male_jux_list),len(ldl_male_jux_list),len(lal_male_jux_list),len(mtal_male_jux_list),len(ctal_male_jux_list))
#print(len(sdl_female_jux_list),len(ldl_female_jux_list),len(lal_female_jux_list),len(mtal_female_jux_list),len(ctal_female_jux_list))

ax.plot(pos_sdl_jux_male,sdl_male_jux_list[1:],'-',color='red',label='Male jux sdl',linewidth=2)
ax.plot(pos_ldl_jux_male,ldl_male_jux_list[1:],'-',color='purple',label='Male jux ldl',linewidth=2)
ax.plot(pos_lal_jux_male,lal_male_jux_list[1:],'-',color='yellow',label='Male jux ldl',linewidth=2)
ax.plot(pos_mtal_jux_male,mtal_male_jux_list[1:],'-',color='blue',label='Male jux mtal',linewidth=2)
ax.plot(pos_ctal_jux_male,ctal_male_jux_list[1:],'-',color='green',label='Male jux ctal',linewidth=2)

ax.plot(pos_sdl_jux_female,sdl_female_jux_list[1:],'--',dashes=(10,10),color='red',label='Female jux sdl',linewidth=2)
ax.plot(pos_ldl_jux_female,ldl_female_jux_list[1:],'--',dashes=(10,10),color='purple',label='Female jux ldl',linewidth=2)
ax.plot(pos_lal_jux_female,lal_female_jux_list[1:],'--',dashes=(10,10),color='yellow',label='Female jux lal',linewidth=2)
ax.plot(pos_mtal_jux_female,mtal_female_jux_list[1:],'--',dashes=(10,10),color='blue',label='Female jux mtal',linewidth=2)
ax.plot(pos_ctal_jux_female,ctal_female_jux_list[1:],'--',dashes=(10,10),color='green',label='Female jux ctal',linewidth=2)

ax.legend(fontsize=30,markerscale=5)
ax.tick_params(labelsize=30)
ax.set_ylabel('Osmolality (mM)',fontsize=40)
#ax.get_xaxis().set_visible(False)

plt.savefig('./'+plot_folder+'/'+'Luminal osmolality along juxtamedullary loop of henle',bbox_inches='tight')
plt.close()

fig,ax = plt.subplots()
fig.set_figheight(20)
fig.set_figwidth(20)

pos=[i for i in range(1,996)]

ax.plot(pos[0:199],sdl_male_jux_list[1:],'-',color='red',label='Male jux sdl',linewidth=2)
ax.plot(pos[199:398],ldl_male_jux_list[1:],'-',color='purple',label='Male jux ldl',linewidth=2)
ax.plot(pos[398:597],lal_male_jux_list[1:],'-',color='yellow',label='Male jux ldl',linewidth=2)
ax.plot(pos[597:796],mtal_male_jux_list[1:],'-',color='blue',label='Male jux mtal',linewidth=2)
ax.plot(pos[796:995],ctal_male_jux_list[1:],'-',color='green',label='Male jux ctal',linewidth=2)

ax.plot(pos[0:199],sdl_female_jux_list[1:],'--',dashes=(10,10),color='red',label='Female jux sdl',linewidth=2)
ax.plot(pos[199:398],ldl_female_jux_list[1:],'--',dashes=(10,10),color='purple',label='Female jux ldl',linewidth=2)
ax.plot(pos[398:597],lal_female_jux_list[1:],'--',dashes=(10,10),color='yellow',label='Female jux lal',linewidth=2)
ax.plot(pos[597:796],mtal_female_jux_list[1:],'--',dashes=(10,10),color='blue',label='Female jux mtal',linewidth=2)
ax.plot(pos[796:995],ctal_female_jux_list[1:],'--',dashes=(10,10),color='green',label='Female jux ctal',linewidth=2)

ax.legend(fontsize=30,markerscale=5)
ax.tick_params(labelsize=30)
ax.set_ylabel('Osmolality (mM)',fontsize=40)
ax.get_xaxis().set_visible(False)

plt.savefig('./'+plot_folder+'/'+'Luminal osmolality along juxtamedullary loop of henle index x axis',bbox_inches='tight')
plt.close()

#====================================================
# Plot of luminal flow along loop of henle
#====================================================
loop_of_henle_sup = ['sdl','mtal','ctal']
loop_of_henle_jux = ['sdl','ldl','lal','mtal','ctal']

pos_sdl_sup_male = [1.1+0.14*i/199 for i in range(1,200)]
pos_mtal_sup_male = [1.1+0.14+0.2*i/199 for i in range(1,200)]
pos_ctal_sup_male = [1.1+0.14+0.2+0.2*i/199 for i in range(1,200)]

pos_sdl_sup_female = [0.88+0.119*i/199 for i in range(1,200)]
pos_mtal_sup_female = [0.88+0.119+0.17*i/199 for i in range(1,200)]
pos_ctal_sup_female = [0.88+0.119+0.17+0.17*i/199 for i in range(1,200)]

pos_sdl_jux_male = [1.1+0.14*i/199 for i in range(1,200)]
pos_ldl_jux_male = [1.1+0.14+0.2*0.5*i/199 for i in range(1,200)]
pos_lal_jux_male = [1.1+0.14+0.2*0.5+0.2*0.5*i/199 for i in range(1,200)]
pos_mtal_jux_male = [1.1+0.14+0.2*0.5+0.2*0.5+0.2*i/199 for i in range(1,200)]
pos_ctal_jux_male = [1.1+0.14+0.2*0.5+0.2*0.5+0.2+0.05*i/199 for i in range(1,200)]

pos_sdl_jux_female = [0.88+0.119*i/199 for i in range(1,200)]
pos_ldl_jux_female = [0.88+0.119+0.2*0.425*i/199 for i in range(1,200)]
pos_lal_jux_female = [0.88+0.119+0.2*0.425+0.2*0.425*i/199 for i in range(1,200)]
pos_mtal_jux_female = [0.88+0.119+0.2*0.425+0.2*0.425+0.17*i/199 for i in range(1,200)]
pos_ctal_jux_female = [0.88+0.119+0.2*0.425+0.2*0.425+0.17+0.0425*i/199 for i in range(1,200)]

for solute in solutes:

	sdl_male_sup = open('./'+male_file+'/'+'malesdl_flow_of_'+solute+'_in_Lumen_sup.txt','r')
	mtal_male_sup = open('./'+male_file+'/'+'malemtal_flow_of_'+solute+'_in_Lumen_sup.txt','r')
	ctal_male_sup = open('./'+male_file+'/'+'malectal_flow_of_'+solute+'_in_Lumen_sup.txt','r')

	sdl_female_sup = open('./'+female_file+'/'+'femalesdl_flow_of_'+solute+'_in_Lumen_sup.txt','r')
	mtal_female_sup = open('./'+female_file+'/'+'femalemtal_flow_of_'+solute+'_in_Lumen_sup.txt','r')
	ctal_female_sup = open('./'+female_file+'/'+'femalectal_flow_of_'+solute+'_in_Lumen_sup.txt','r')

	sdl_male_jux = open('./'+male_file+'/'+'malesdl_flow_of_'+solute+'_in_Lumen_jux1.txt','r')
	ldl_male_jux = open('./'+male_file+'/'+'maleldl_flow_of_'+solute+'_in_Lumen_jux1.txt','r')
	lal_male_jux = open('./'+male_file+'/'+'malelal_flow_of_'+solute+'_in_Lumen_jux1.txt','r')
	mtal_male_jux = open('./'+male_file+'/'+'malemtal_flow_of_'+solute+'_in_Lumen_jux1.txt','r')
	ctal_male_jux = open('./'+male_file+'/'+'malectal_flow_of_'+solute+'_in_Lumen_jux1.txt','r')

	sdl_female_jux = open('./'+female_file+'/'+'femalesdl_flow_of_'+solute+'_in_Lumen_jux1.txt','r')
	ldl_female_jux = open('./'+female_file+'/'+'femaleldl_flow_of_'+solute+'_in_Lumen_jux1.txt','r')
	lal_female_jux = open('./'+female_file+'/'+'femalelal_flow_of_'+solute+'_in_Lumen_jux1.txt','r')
	mtal_female_jux = open('./'+female_file+'/'+'femalemtal_flow_of_'+solute+'_in_Lumen_jux1.txt','r')
	ctal_female_jux = open('./'+female_file+'/'+'femalectal_flow_of_'+solute+'_in_Lumen_jux1.txt','r')

	sdl_male_sup_list = []
	mtal_male_sup_list = []
	ctal_male_sup_list = []

	sdl_female_sup_list = []
	mtal_female_sup_list = []
	ctal_female_sup_list = []

	sdl_male_jux_list = []
	ldl_male_jux_list = []
	lal_male_jux_list = []
	mtal_male_jux_list = []
	ctal_male_jux_list = []

	sdl_female_jux_list = []
	ldl_female_jux_list = []
	lal_female_jux_list = []
	mtal_female_jux_list = []
	ctal_female_jux_list = []

	for i in sdl_male_sup:
		line = i.split(' ')
		sdl_male_sup_list.append(float(line[0]))

	for i in mtal_male_sup:
		line = i.split(' ')
		mtal_male_sup_list.append(float(line[0]))

	for i in ctal_male_sup:
		line = i.split(' ')
		ctal_male_sup_list.append(float(line[0]))

	for i in sdl_female_sup:
		line = i.split(' ')
		sdl_female_sup_list.append(float(line[0]))

	for i in mtal_female_sup:
		line = i.split(' ')
		mtal_female_sup_list.append(float(line[0]))

	for i in ctal_female_sup:
		line = i.split(' ')
		ctal_female_sup_list.append(float(line[0]))

	for i in sdl_male_jux:
		line = i.split(' ')
		sdl_male_jux_list.append(float(line[0]))

	for i in ldl_male_jux:
		line = i.split(' ')
		ldl_male_jux_list.append(float(line[0]))

	for i in lal_male_jux:
		line = i.split(' ')
		lal_male_jux_list.append(float(line[0]))

	for i in mtal_male_jux:
		line = i.split(' ')
		mtal_male_jux_list.append(float(line[0]))

	for i in ctal_male_jux:
		line = i.split(' ')
		ctal_male_jux_list.append(float(line[0]))

	for i in sdl_female_jux:
		line = i.split(' ')
		sdl_female_jux_list.append(float(line[0]))

	for i in ldl_female_jux:
		line = i.split(' ')
		ldl_female_jux_list.append(float(line[0]))

	for i in lal_female_jux:
		line = i.split(' ')
		lal_female_jux_list.append(float(line[0]))

	for i in mtal_female_jux:
		line = i.split(' ')
		mtal_female_jux_list.append(float(line[0]))

	for i in ctal_female_jux:
		line = i.split(' ')
		ctal_female_jux_list.append(float(line[0]))

	fig,ax = plt.subplots()
	fig.set_figheight(20)
	fig.set_figwidth(20)

	ax.plot(pos_sdl_sup_male,sdl_male_sup_list[1:],'-',color='red',label='Male sup sdl',linewidth=2)
	ax.plot(pos_mtal_sup_male,mtal_male_sup_list[1:],'-',color='blue',label='Male sup mtal',linewidth=2)
	ax.plot(pos_ctal_sup_male,ctal_male_sup_list[1:],'-',color='green',label='Male sup ctal',linewidth=2)

	ax.plot(pos_sdl_sup_female,sdl_female_sup_list[1:],'--',dashes=(10,10),color='red',label='Female sup sdl',linewidth=2)
	ax.plot(pos_mtal_sup_female,mtal_female_sup_list[1:],'--',dashes=(10,10),color='blue',label='Female sup mtal',linewidth=2)
	ax.plot(pos_ctal_sup_female,ctal_female_sup_list[1:],'--',dashes=(10,10),color='green',label='Female sup ctal',linewidth=2)

	ax.legend(fontsize=30,markerscale=5)
	ax.tick_params(labelsize=30)
	ax.set_ylabel('Flow of '+solute+' (pmol/min)',fontsize=40)
	#ax.get_xaxis().set_visible(False)

	plt.savefig('./'+plot_folder+'/'+'Luminal flow of '+solute+' along superficial loop of henle',bbox_inches='tight')
	plt.close()

	fig,ax = plt.subplots()
	fig.set_figheight(20)
	fig.set_figwidth(20)

#print(len(sdl_male_jux_list),len(ldl_male_jux_list),len(lal_male_jux_list),len(mtal_male_jux_list),len(ctal_male_jux_list))
#print(len(sdl_female_jux_list),len(ldl_female_jux_list),len(lal_female_jux_list),len(mtal_female_jux_list),len(ctal_female_jux_list))

	ax.plot(pos_sdl_jux_male,sdl_male_jux_list[1:],'-',color='red',label='Male jux sdl',linewidth=2)
	ax.plot(pos_ldl_jux_male,ldl_male_jux_list[1:],'-',color='purple',label='Male jux ldl',linewidth=2)
	ax.plot(pos_lal_jux_male,lal_male_jux_list[1:],'-',color='yellow',label='Male jux ldl',linewidth=2)
	ax.plot(pos_mtal_jux_male,mtal_male_jux_list[1:],'-',color='blue',label='Male jux mtal',linewidth=2)
	ax.plot(pos_ctal_jux_male,ctal_male_jux_list[1:],'-',color='green',label='Male jux ctal',linewidth=2)

	ax.plot(pos_sdl_jux_female,sdl_female_jux_list[1:],'--',dashes=(10,10),color='red',label='Female jux sdl',linewidth=2)
	ax.plot(pos_ldl_jux_female,ldl_female_jux_list[1:],'--',dashes=(10,10),color='purple',label='Female jux ldl',linewidth=2)
	ax.plot(pos_lal_jux_female,lal_female_jux_list[1:],'--',dashes=(10,10),color='yellow',label='Female jux lal',linewidth=2)
	ax.plot(pos_mtal_jux_female,mtal_female_jux_list[1:],'--',dashes=(10,10),color='blue',label='Female jux mtal',linewidth=2)
	ax.plot(pos_ctal_jux_female,ctal_female_jux_list[1:],'--',dashes=(10,10),color='green',label='Female jux ctal',linewidth=2)

	ax.legend(fontsize=30,markerscale=5)
	ax.tick_params(labelsize=30)
	ax.set_ylabel('Flow of '+solute+' (pmol/min)',fontsize=40)
#ax.get_xaxis().set_visible(False)

	plt.savefig('./'+plot_folder+'/'+'Luminal flow of '+solute+' along juxtamedullary loop of henle',bbox_inches='tight')
	plt.close()

	fig,ax = plt.subplots()
	fig.set_figheight(20)
	fig.set_figwidth(20)

	pos=[i for i in range(1,996)]

	ax.plot(pos[0:199],sdl_male_jux_list[1:],'-',color='red',label='Male jux sdl',linewidth=2)
	ax.plot(pos[199:398],ldl_male_jux_list[1:],'-',color='purple',label='Male jux ldl',linewidth=2)
	ax.plot(pos[398:597],lal_male_jux_list[1:],'-',color='yellow',label='Male jux ldl',linewidth=2)
	ax.plot(pos[597:796],mtal_male_jux_list[1:],'-',color='blue',label='Male jux mtal',linewidth=2)
	ax.plot(pos[796:995],ctal_male_jux_list[1:],'-',color='green',label='Male jux ctal',linewidth=2)

	ax.plot(pos[0:199],sdl_female_jux_list[1:],'--',dashes=(10,10),color='red',label='Female jux sdl',linewidth=2)
	ax.plot(pos[199:398],ldl_female_jux_list[1:],'--',dashes=(10,10),color='purple',label='Female jux ldl',linewidth=2)
	ax.plot(pos[398:597],lal_female_jux_list[1:],'--',dashes=(10,10),color='yellow',label='Female jux lal',linewidth=2)
	ax.plot(pos[597:796],mtal_female_jux_list[1:],'--',dashes=(10,10),color='blue',label='Female jux mtal',linewidth=2)
	ax.plot(pos[796:995],ctal_female_jux_list[1:],'--',dashes=(10,10),color='green',label='Female jux ctal',linewidth=2)

	ax.legend(fontsize=30,markerscale=5)
	ax.tick_params(labelsize=30)
	ax.set_ylabel('Flow of '+solute+' (pmol/min)',fontsize=40)
	ax.get_xaxis().set_visible(False)

	plt.savefig('./'+plot_folder+'/'+'Luminal flow of '+solute+' along juxtamedullary loop of henle index x axis',bbox_inches='tight')
	plt.close()
#====================================================
# Plot of luminal concentration along loop of henle
#====================================================
loop_of_henle_sup = ['sdl','mtal','ctal']
loop_of_henle_jux = ['sdl','ldl','lal','mtal','ctal']

pos_sdl_sup_male = [1.1+0.14*i/199 for i in range(1,200)]
pos_mtal_sup_male = [1.1+0.14+0.2*i/199 for i in range(1,200)]
pos_ctal_sup_male = [1.1+0.14+0.2+0.2*i/199 for i in range(1,200)]

pos_sdl_sup_female = [0.88+0.119*i/199 for i in range(1,200)]
pos_mtal_sup_female = [0.88+0.119+0.17*i/199 for i in range(1,200)]
pos_ctal_sup_female = [0.88+0.119+0.17+0.17*i/199 for i in range(1,200)]

pos_sdl_jux_male = [1.1+0.14*i/199 for i in range(1,200)]
pos_ldl_jux_male = [1.1+0.14+0.2*0.5*i/199 for i in range(1,200)]
pos_lal_jux_male = [1.1+0.14+0.2*0.5+0.2*0.5*i/199 for i in range(1,200)]
pos_mtal_jux_male = [1.1+0.14+0.2*0.5+0.2*0.5+0.2*i/199 for i in range(1,200)]
pos_ctal_jux_male = [1.1+0.14+0.2*0.5+0.2*0.5+0.2+0.05*i/199 for i in range(1,200)]

pos_sdl_jux_female = [0.88+0.119*i/199 for i in range(1,200)]
pos_ldl_jux_female = [0.88+0.119+0.2*0.425*i/199 for i in range(1,200)]
pos_lal_jux_female = [0.88+0.119+0.2*0.425+0.2*0.425*i/199 for i in range(1,200)]
pos_mtal_jux_female = [0.88+0.119+0.2*0.425+0.2*0.425+0.17*i/199 for i in range(1,200)]
pos_ctal_jux_female = [0.88+0.119+0.2*0.425+0.2*0.425+0.17+0.0425*i/199 for i in range(1,200)]

for solute in solutes:

	sdl_male_sup = open('./'+male_file+'/'+'malesdl_con_of_'+solute+'_in_Lumen_sup.txt','r')
	mtal_male_sup = open('./'+male_file+'/'+'malemtal_con_of_'+solute+'_in_Lumen_sup.txt','r')
	ctal_male_sup = open('./'+male_file+'/'+'malectal_con_of_'+solute+'_in_Lumen_sup.txt','r')

	sdl_female_sup = open('./'+female_file+'/'+'femalesdl_con_of_'+solute+'_in_Lumen_sup.txt','r')
	mtal_female_sup = open('./'+female_file+'/'+'femalemtal_con_of_'+solute+'_in_Lumen_sup.txt','r')
	ctal_female_sup = open('./'+female_file+'/'+'femalectal_con_of_'+solute+'_in_Lumen_sup.txt','r')

	sdl_male_jux = open('./'+male_file+'/'+'malesdl_con_of_'+solute+'_in_Lumen_jux1.txt','r')
	ldl_male_jux = open('./'+male_file+'/'+'maleldl_con_of_'+solute+'_in_Lumen_jux1.txt','r')
	lal_male_jux = open('./'+male_file+'/'+'malelal_con_of_'+solute+'_in_Lumen_jux1.txt','r')
	mtal_male_jux = open('./'+male_file+'/'+'malemtal_con_of_'+solute+'_in_Lumen_jux1.txt','r')
	ctal_male_jux = open('./'+male_file+'/'+'malectal_con_of_'+solute+'_in_Lumen_jux1.txt','r')

	sdl_female_jux = open('./'+female_file+'/'+'femalesdl_con_of_'+solute+'_in_Lumen_jux1.txt','r')
	ldl_female_jux = open('./'+female_file+'/'+'femaleldl_con_of_'+solute+'_in_Lumen_jux1.txt','r')
	lal_female_jux = open('./'+female_file+'/'+'femalelal_con_of_'+solute+'_in_Lumen_jux1.txt','r')
	mtal_female_jux = open('./'+female_file+'/'+'femalemtal_con_of_'+solute+'_in_Lumen_jux1.txt','r')
	ctal_female_jux = open('./'+female_file+'/'+'femalectal_con_of_'+solute+'_in_Lumen_jux1.txt','r')

	sdl_male_sup_list = []
	mtal_male_sup_list = []
	ctal_male_sup_list = []

	sdl_female_sup_list = []
	mtal_female_sup_list = []
	ctal_female_sup_list = []

	sdl_male_jux_list = []
	ldl_male_jux_list = []
	lal_male_jux_list = []
	mtal_male_jux_list = []
	ctal_male_jux_list = []

	sdl_female_jux_list = []
	ldl_female_jux_list = []
	lal_female_jux_list = []
	mtal_female_jux_list = []
	ctal_female_jux_list = []

	for i in sdl_male_sup:
		line = i.split(' ')
		sdl_male_sup_list.append(float(line[0]))

	for i in mtal_male_sup:
		line = i.split(' ')
		mtal_male_sup_list.append(float(line[0]))

	for i in ctal_male_sup:
		line = i.split(' ')
		ctal_male_sup_list.append(float(line[0]))

	for i in sdl_female_sup:
		line = i.split(' ')
		sdl_female_sup_list.append(float(line[0]))

	for i in mtal_female_sup:
		line = i.split(' ')
		mtal_female_sup_list.append(float(line[0]))

	for i in ctal_female_sup:
		line = i.split(' ')
		ctal_female_sup_list.append(float(line[0]))

	for i in sdl_male_jux:
		line = i.split(' ')
		sdl_male_jux_list.append(float(line[0]))

	for i in ldl_male_jux:
		line = i.split(' ')
		ldl_male_jux_list.append(float(line[0]))

	for i in lal_male_jux:
		line = i.split(' ')
		lal_male_jux_list.append(float(line[0]))

	for i in mtal_male_jux:
		line = i.split(' ')
		mtal_male_jux_list.append(float(line[0]))

	for i in ctal_male_jux:
		line = i.split(' ')
		ctal_male_jux_list.append(float(line[0]))

	for i in sdl_female_jux:
		line = i.split(' ')
		sdl_female_jux_list.append(float(line[0]))

	for i in ldl_female_jux:
		line = i.split(' ')
		ldl_female_jux_list.append(float(line[0]))

	for i in lal_female_jux:
		line = i.split(' ')
		lal_female_jux_list.append(float(line[0]))

	for i in mtal_female_jux:
		line = i.split(' ')
		mtal_female_jux_list.append(float(line[0]))

	for i in ctal_female_jux:
		line = i.split(' ')
		ctal_female_jux_list.append(float(line[0]))

	fig,ax = plt.subplots()
	fig.set_figheight(20)
	fig.set_figwidth(20)

	ax.plot(pos_sdl_sup_male,sdl_male_sup_list[1:],'-',color='red',label='Male sup sdl',linewidth=2)
	ax.plot(pos_mtal_sup_male,mtal_male_sup_list[1:],'-',color='blue',label='Male sup mtal',linewidth=2)
	ax.plot(pos_ctal_sup_male,ctal_male_sup_list[1:],'-',color='green',label='Male sup ctal',linewidth=2)

	ax.plot(pos_sdl_sup_female,sdl_female_sup_list[1:],'--',dashes=(10,10),color='red',label='Female sup sdl',linewidth=2)
	ax.plot(pos_mtal_sup_female,mtal_female_sup_list[1:],'--',dashes=(10,10),color='blue',label='Female sup mtal',linewidth=2)
	ax.plot(pos_ctal_sup_female,ctal_female_sup_list[1:],'--',dashes=(10,10),color='green',label='Female sup ctal',linewidth=2)

	ax.legend(fontsize=30,markerscale=5)
	ax.tick_params(labelsize=30)
	ax.set_ylabel('Concentration of '+solute+' (mM)',fontsize=40)
	#ax.get_xaxis().set_visible(False)

	plt.savefig('./'+plot_folder+'/'+'Luminal concentration of '+solute+' along superficial loop of henle',bbox_inches='tight')
	plt.close()

	fig,ax = plt.subplots()
	fig.set_figheight(20)
	fig.set_figwidth(20)

#print(len(sdl_male_jux_list),len(ldl_male_jux_list),len(lal_male_jux_list),len(mtal_male_jux_list),len(ctal_male_jux_list))
#print(len(sdl_female_jux_list),len(ldl_female_jux_list),len(lal_female_jux_list),len(mtal_female_jux_list),len(ctal_female_jux_list))

	ax.plot(pos_sdl_jux_male,sdl_male_jux_list[1:],'-',color='red',label='Male jux sdl',linewidth=2)
	ax.plot(pos_ldl_jux_male,ldl_male_jux_list[1:],'-',color='purple',label='Male jux ldl',linewidth=2)
	ax.plot(pos_lal_jux_male,lal_male_jux_list[1:],'-',color='yellow',label='Male jux ldl',linewidth=2)
	ax.plot(pos_mtal_jux_male,mtal_male_jux_list[1:],'-',color='blue',label='Male jux mtal',linewidth=2)
	ax.plot(pos_ctal_jux_male,ctal_male_jux_list[1:],'-',color='green',label='Male jux ctal',linewidth=2)

	ax.plot(pos_sdl_jux_female,sdl_female_jux_list[1:],'--',dashes=(10,10),color='red',label='Female jux sdl',linewidth=2)
	ax.plot(pos_ldl_jux_female,ldl_female_jux_list[1:],'--',dashes=(10,10),color='purple',label='Female jux ldl',linewidth=2)
	ax.plot(pos_lal_jux_female,lal_female_jux_list[1:],'--',dashes=(10,10),color='yellow',label='Female jux lal',linewidth=2)
	ax.plot(pos_mtal_jux_female,mtal_female_jux_list[1:],'--',dashes=(10,10),color='blue',label='Female jux mtal',linewidth=2)
	ax.plot(pos_ctal_jux_female,ctal_female_jux_list[1:],'--',dashes=(10,10),color='green',label='Female jux ctal',linewidth=2)

	ax.legend(fontsize=30,markerscale=5)
	ax.tick_params(labelsize=30)
	ax.set_ylabel('Concentration of '+solute+' (mM)',fontsize=40)
#ax.get_xaxis().set_visible(False)

	plt.savefig('./'+plot_folder+'/'+'Luminal concentration of '+solute+' along juxtamedullary loop of henle',bbox_inches='tight')
	plt.close()

	fig,ax = plt.subplots()
	fig.set_figheight(20)
	fig.set_figwidth(20)

	pos=[i for i in range(1,996)]

	ax.plot(pos[0:199],sdl_male_jux_list[1:],'-',color='red',label='Male jux sdl',linewidth=2)
	ax.plot(pos[199:398],ldl_male_jux_list[1:],'-',color='purple',label='Male jux ldl',linewidth=2)
	ax.plot(pos[398:597],lal_male_jux_list[1:],'-',color='yellow',label='Male jux ldl',linewidth=2)
	ax.plot(pos[597:796],mtal_male_jux_list[1:],'-',color='blue',label='Male jux mtal',linewidth=2)
	ax.plot(pos[796:995],ctal_male_jux_list[1:],'-',color='green',label='Male jux ctal',linewidth=2)

	ax.plot(pos[0:199],sdl_female_jux_list[1:],'--',dashes=(10,10),color='red',label='Female jux sdl',linewidth=2)
	ax.plot(pos[199:398],ldl_female_jux_list[1:],'--',dashes=(10,10),color='purple',label='Female jux ldl',linewidth=2)
	ax.plot(pos[398:597],lal_female_jux_list[1:],'--',dashes=(10,10),color='yellow',label='Female jux lal',linewidth=2)
	ax.plot(pos[597:796],mtal_female_jux_list[1:],'--',dashes=(10,10),color='blue',label='Female jux mtal',linewidth=2)
	ax.plot(pos[796:995],ctal_female_jux_list[1:],'--',dashes=(10,10),color='green',label='Female jux ctal',linewidth=2)

	ax.legend(fontsize=30,markerscale=5)
	ax.tick_params(labelsize=30)
	ax.set_ylabel('Concentration of '+solute+' (mM)',fontsize=40)
	ax.get_xaxis().set_visible(False)

	plt.savefig('./'+plot_folder+'/'+'Luminal Concentration of '+solute+' along juxtamedullary loop of henle index x axis',bbox_inches='tight')
	plt.close()

#====================================================
# Plot of water volume along loop of henle
#====================================================
loop_of_henle_sup = ['sdl','mtal','ctal']
loop_of_henle_jux = ['sdl','ldl','lal','mtal','ctal']

pos_sdl_sup_male = [1.1+0.14*i/199 for i in range(1,200)]
pos_mtal_sup_male = [1.1+0.14+0.2*i/199 for i in range(1,200)]
pos_ctal_sup_male = [1.1+0.14+0.2+0.2*i/199 for i in range(1,200)]

pos_sdl_sup_female = [0.88+0.119*i/199 for i in range(1,200)]
pos_mtal_sup_female = [0.88+0.119+0.17*i/199 for i in range(1,200)]
pos_ctal_sup_female = [0.88+0.119+0.17+0.17*i/199 for i in range(1,200)]

pos_sdl_jux_male = [1.1+0.14*i/199 for i in range(1,200)]
pos_ldl_jux_male = [1.1+0.14+0.2*0.5*i/199 for i in range(1,200)]
pos_lal_jux_male = [1.1+0.14+0.2*0.5+0.2*0.5*i/199 for i in range(1,200)]
pos_mtal_jux_male = [1.1+0.14+0.2*0.5+0.2*0.5+0.2*i/199 for i in range(1,200)]
pos_ctal_jux_male = [1.1+0.14+0.2*0.5+0.2*0.5+0.2+0.05*i/199 for i in range(1,200)]

pos_sdl_jux_female = [0.88+0.119*i/199 for i in range(1,200)]
pos_ldl_jux_female = [0.88+0.119+0.2*0.425*i/199 for i in range(1,200)]
pos_lal_jux_female = [0.88+0.119+0.2*0.425+0.2*0.425*i/199 for i in range(1,200)]
pos_mtal_jux_female = [0.88+0.119+0.2*0.425+0.2*0.425+0.17*i/199 for i in range(1,200)]
pos_ctal_jux_female = [0.88+0.119+0.2*0.425+0.2*0.425+0.17+0.0425*i/199 for i in range(1,200)]

sdl_male_sup = open('./'+male_file+'/'+'malesdl_water_volume_in_Lumen_sup.txt','r')
mtal_male_sup = open('./'+male_file+'/'+'malemtal_water_volume_in_Lumen_sup.txt','r')
ctal_male_sup = open('./'+male_file+'/'+'malectal_water_volume_in_Lumen_sup.txt','r')

sdl_female_sup = open('./'+female_file+'/'+'femalesdl_water_volume_in_Lumen_sup.txt','r')
mtal_female_sup = open('./'+female_file+'/'+'femalemtal_water_volume_in_Lumen_sup.txt','r')
ctal_female_sup = open('./'+female_file+'/'+'femalectal_water_volume_in_Lumen_sup.txt','r')

sdl_male_jux = open('./'+male_file+'/'+'malesdl_water_volume_in_Lumen_jux1.txt','r')
ldl_male_jux = open('./'+male_file+'/'+'maleldl_water_volume_in_Lumen_jux1.txt','r')
lal_male_jux = open('./'+male_file+'/'+'malelal_water_volume_in_Lumen_jux1.txt','r')
mtal_male_jux = open('./'+male_file+'/'+'malemtal_water_volume_in_Lumen_jux1.txt','r')
ctal_male_jux = open('./'+male_file+'/'+'malectal_water_volume_in_Lumen_jux1.txt','r')

sdl_female_jux = open('./'+female_file+'/'+'femalesdl_water_volume_in_Lumen_jux1.txt','r')
ldl_female_jux = open('./'+female_file+'/'+'femaleldl_water_volume_in_Lumen_jux1.txt','r')
lal_female_jux = open('./'+female_file+'/'+'femalelal_water_volume_in_Lumen_jux1.txt','r')
mtal_female_jux = open('./'+female_file+'/'+'femalemtal_water_volume_in_Lumen_jux1.txt','r')
ctal_female_jux = open('./'+female_file+'/'+'femalectal_water_volume_in_Lumen_jux1.txt','r')

sdl_male_sup_list = []
mtal_male_sup_list = []
ctal_male_sup_list = []

sdl_female_sup_list = []
mtal_female_sup_list = []
ctal_female_sup_list = []

sdl_male_jux_list = []
ldl_male_jux_list = []
lal_male_jux_list = []
mtal_male_jux_list = []
ctal_male_jux_list = []

sdl_female_jux_list = []
ldl_female_jux_list = []
lal_female_jux_list = []
mtal_female_jux_list = []
ctal_female_jux_list = []

for i in sdl_male_sup:
	line = i.split(' ')
	sdl_male_sup_list.append(float(line[0]))

for i in mtal_male_sup:
	line = i.split(' ')
	mtal_male_sup_list.append(float(line[0]))

for i in ctal_male_sup:
	line = i.split(' ')
	ctal_male_sup_list.append(float(line[0]))

for i in sdl_female_sup:
	line = i.split(' ')
	sdl_female_sup_list.append(float(line[0]))

for i in mtal_female_sup:
	line = i.split(' ')
	mtal_female_sup_list.append(float(line[0]))

for i in ctal_female_sup:
	line = i.split(' ')
	ctal_female_sup_list.append(float(line[0]))

for i in sdl_male_jux:
	line = i.split(' ')
	sdl_male_jux_list.append(float(line[0]))

for i in ldl_male_jux:
	line = i.split(' ')
	ldl_male_jux_list.append(float(line[0]))

for i in lal_male_jux:
	line = i.split(' ')
	lal_male_jux_list.append(float(line[0]))

for i in mtal_male_jux:
	line = i.split(' ')
	mtal_male_jux_list.append(float(line[0]))

for i in ctal_male_jux:
	line = i.split(' ')
	ctal_male_jux_list.append(float(line[0]))

for i in sdl_female_jux:
	line = i.split(' ')
	sdl_female_jux_list.append(float(line[0]))

for i in ldl_female_jux:
	line = i.split(' ')
	ldl_female_jux_list.append(float(line[0]))

for i in lal_female_jux:
	line = i.split(' ')
	lal_female_jux_list.append(float(line[0]))

for i in mtal_female_jux:
	line = i.split(' ')
	mtal_female_jux_list.append(float(line[0]))

for i in ctal_female_jux:
	line = i.split(' ')
	ctal_female_jux_list.append(float(line[0]))

fig,ax = plt.subplots()
fig.set_figheight(20)
fig.set_figwidth(20)

ax.plot(pos_sdl_sup_male,sdl_male_sup_list[1:],'-',color='red',label='Male sup sdl',linewidth=2)
ax.plot(pos_mtal_sup_male,mtal_male_sup_list[1:],'-',color='blue',label='Male sup mtal',linewidth=2)
ax.plot(pos_ctal_sup_male,ctal_male_sup_list[1:],'-',color='green',label='Male sup ctal',linewidth=2)

ax.plot(pos_sdl_sup_female,sdl_female_sup_list[1:],'--',dashes=(10,10),color='red',label='Female sup sdl',linewidth=2)
ax.plot(pos_mtal_sup_female,mtal_female_sup_list[1:],'--',dashes=(10,10),color='blue',label='Female sup mtal',linewidth=2)
ax.plot(pos_ctal_sup_female,ctal_female_sup_list[1:],'--',dashes=(10,10),color='green',label='Female sup ctal',linewidth=2)

ax.legend(fontsize=30,markerscale=5)
ax.tick_params(labelsize=30)
ax.set_ylabel('Water volume (nl/min)',fontsize=40)
#ax.get_xaxis().set_visible(False)

plt.savefig('./'+plot_folder+'/'+'Luminal water volume along superficial loop of henle',bbox_inches='tight')
plt.close()

fig,ax = plt.subplots()
fig.set_figheight(20)
fig.set_figwidth(20)

#print(len(sdl_male_jux_list),len(ldl_male_jux_list),len(lal_male_jux_list),len(mtal_male_jux_list),len(ctal_male_jux_list))
#print(len(sdl_female_jux_list),len(ldl_female_jux_list),len(lal_female_jux_list),len(mtal_female_jux_list),len(ctal_female_jux_list))

ax.plot(pos_sdl_jux_male,sdl_male_jux_list[1:],'-',color='red',label='Male jux sdl',linewidth=2)
ax.plot(pos_ldl_jux_male,ldl_male_jux_list[1:],'-',color='purple',label='Male jux ldl',linewidth=2)
ax.plot(pos_lal_jux_male,lal_male_jux_list[1:],'-',color='yellow',label='Male jux ldl',linewidth=2)
ax.plot(pos_mtal_jux_male,mtal_male_jux_list[1:],'-',color='blue',label='Male jux mtal',linewidth=2)
ax.plot(pos_ctal_jux_male,ctal_male_jux_list[1:],'-',color='green',label='Male jux ctal',linewidth=2)

ax.plot(pos_sdl_jux_female,sdl_female_jux_list[1:],'--',dashes=(10,10),color='red',label='Female jux sdl',linewidth=2)
ax.plot(pos_ldl_jux_female,ldl_female_jux_list[1:],'--',dashes=(10,10),color='purple',label='Female jux ldl',linewidth=2)
ax.plot(pos_lal_jux_female,lal_female_jux_list[1:],'--',dashes=(10,10),color='yellow',label='Female jux lal',linewidth=2)
ax.plot(pos_mtal_jux_female,mtal_female_jux_list[1:],'--',dashes=(10,10),color='blue',label='Female jux mtal',linewidth=2)
ax.plot(pos_ctal_jux_female,ctal_female_jux_list[1:],'--',dashes=(10,10),color='green',label='Female jux ctal',linewidth=2)

ax.legend(fontsize=30,markerscale=5)
ax.tick_params(labelsize=30)
ax.set_ylabel('Water volume (nl/min)',fontsize=40)
#ax.get_xaxis().set_visible(False)

plt.savefig('./'+plot_folder+'/'+'Luminal water volume along juxtamedullary loop of henle',bbox_inches='tight')
plt.close()

fig,ax = plt.subplots()
fig.set_figheight(20)
fig.set_figwidth(20)

pos=[i for i in range(1,996)]

ax.plot(pos[0:199],sdl_male_jux_list[1:],'-',color='red',label='Male jux sdl',linewidth=2)
ax.plot(pos[199:398],ldl_male_jux_list[1:],'-',color='purple',label='Male jux ldl',linewidth=2)
ax.plot(pos[398:597],lal_male_jux_list[1:],'-',color='yellow',label='Male jux ldl',linewidth=2)
ax.plot(pos[597:796],mtal_male_jux_list[1:],'-',color='blue',label='Male jux mtal',linewidth=2)
ax.plot(pos[796:995],ctal_male_jux_list[1:],'-',color='green',label='Male jux ctal',linewidth=2)

ax.plot(pos[0:199],sdl_female_jux_list[1:],'--',dashes=(10,10),color='red',label='Female jux sdl',linewidth=2)
ax.plot(pos[199:398],ldl_female_jux_list[1:],'--',dashes=(10,10),color='purple',label='Female jux ldl',linewidth=2)
ax.plot(pos[398:597],lal_female_jux_list[1:],'--',dashes=(10,10),color='yellow',label='Female jux lal',linewidth=2)
ax.plot(pos[597:796],mtal_female_jux_list[1:],'--',dashes=(10,10),color='blue',label='Female jux mtal',linewidth=2)
ax.plot(pos[796:995],ctal_female_jux_list[1:],'--',dashes=(10,10),color='green',label='Female jux ctal',linewidth=2)

ax.legend(fontsize=30,markerscale=5)
ax.tick_params(labelsize=30)
ax.set_ylabel('Water volume (nl/min)',fontsize=40)
ax.get_xaxis().set_visible(False)

plt.savefig('./'+plot_folder+'/'+'Luminal water volume along juxtamedullary loop of henle index x axis',bbox_inches='tight')
plt.close()

#===============================================
# Plot along the whole nephron
#===============================================

pos_pt_sup=[1.1*i/199 for i in range(1,176)]
pos_s3_sup=[1.1*i/199 for i in range(176,200)]
pos_sdl_sup=[1.1+0.14*i/199 for i in range(1,200)]
pos_mtal_sup=[1.1+0.14+0.2*i/199 for i in range(1,200)]
pos_ctal_sup=[1.1+0.14+0.2+0.2*i/199 for i in range(1,200)]
pos_dct_sup=[1.1+0.14+0.2+0.2+0.1*i/199 for i in range(1,200)]
pos_cnt_sup=[1.1+0.14+0.2+0.2+0.1+0.2*i/199 for i in range(1,200)]
pos_ccd=[1.1+0.14+0.2+0.2+0.1+0.2+0.2*i/199 for i in range(1,200)]
pos_omcd=[1.1+0.14+0.2+0.2+0.1+0.2+0.2+0.2*i/199 for i in range(1,200)]
pos_imcd=[1.1+0.14+0.2+0.2+0.1+0.2+0.2+0.2+0.5*i/199 for i in range(1,200)]
#pos_sup_early=pos_pt+pos_s3+pos_sdl+pos_mtal+pos_ctal+pos_dct+pos_cnt

pos_pt_jux = [1.1*i/199 for i in range(1,176)]
pos_s3_jux = [1.1*i/199 for i in range(176,200)]
pos_sdl_jux = [1.1+0.14*i/199 for i in range(1,200)]
pos_ldl_jux = [1.1+0.14+0.2*0.5*i/199 for i in range(1,200)]
pos_lal_jux = [1.1+0.14+0.2*0.5+0.2*0.5*i/199 for i in range(1,200)]
pos_mtal_jux = [1.1+0.14+0.2*0.5+0.2*0.5+0.2*i/199 for i in range(1,200)]
pos_ctal_jux = [1.1+0.14+0.2*0.5+0.2*0.5+0.2+0.05*i/199 for i in range(1,200)]
pos_dct_jux = [1.1+0.14+0.2*0.5+0.2*0.5+0.2+0.05+0.1*i/199 for i in range(1,200)]
pos_cnt_jux = [1.1+0.14+0.2*0.5+0.2*0.5+0.2+0.05+0.1+0.3*i/199 for i in range(1,200)]