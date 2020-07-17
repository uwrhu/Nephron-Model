import re
from defs import *
from values import *

# pick out compartment IDs. E.g. given label=Area_Lumen_Cell, it would return
# (0,1) for (lumen, cell)
def get_interface_id(label):
	tmp = (label).split('_')
	if len(tmp)==3:  # only have two compartment IDs
		ind1,ind2 = compart_id[tmp[1]],compart_id[tmp[2]]
		return ind1, ind2
	elif len(tmp)==4:  # solute ID, followed by two compartment IDs
		sid,ind1,ind2 = solute_id[tmp[1]],compart_id[tmp[2]],compart_id[tmp[3]]
		return sid,ind1,ind2

# for coupled transpoters
# pick out solute IDs, compartment IDs, and coefficients
def get_coupled_id(label):
	tmp = (label).split('_')
#    print(tmp)
	sid1,sid2 = solute_id[tmp[3]],solute_id[tmp[4]]
	ind1,ind2 = compart_id[tmp[1]],compart_id[tmp[2]]
	if len(tmp)==5:  # only involves 2 solutes
		return ind1,ind2,sid1,sid2
	elif len(tmp)==6:  # involves 3 solutes
		sid3 = solute_id[tmp[5]]
		return ind1,ind2,sid1,sid2,sid3
	else:
		print("Wrong label",tmp)

# is str_short at the beginning of str_long?
# case insensitive
def compare_string_prefix(str_long,str_short):
	return (str_long.lower())[:len(str_short)] == str_short.lower()

def compare_sex(sex, cell):
	if sex.lower() == 'male' or sex.lower() == 'female':
		return cell.sex.lower() == sex.lower()
	else:
		return True

def read_params(cell,filename,j):
#    cell = membrane()
	#filename=input('Choose a file')
	#filename='cTALparams.dat'
	file = open(filename,'r')
	cell.segment=filename[0:-12]
	data = []
	line = file.readline()
	while (line):
		line = line.replace('\t',' ')
		terms = line.split(' ')
		if line[0][0]!='#':
			id = terms[0] #re.findall(r'[A-Za-z_]+', line)
			sex = id.split('_')[-1].lower()
			if compare_sex(sex, cell):
				id = id.replace('_male', '')
				id = id.replace('_female', '')
			else:
				line = file.readline()
				continue;
			# skip over the label, which may contain numbers like in HCO3
			first_space_pos = line.index(' ')
			num = re.findall(r'-?\d+\.?\d*[Ee]?[+-]?\d*', line[first_space_pos:len(line)])
			if num: # if this line is numerical param
				value = float(num[0])


			if id.lower() == "Sex".lower():
				find = re.findall("Female".lower(), terms[-1].lower())
				if find:
					cell.sex = "Female".lower()
			
			elif compare_string_prefix(id,"Diameter"):
				if cell.diabete == 'No':
					cell.diam = value
				elif cell.diabete == 'Yes':
					if cell.segment == 'PT' or cell.segment == 'S3':
						cell.diam = value*32/25
					elif cell.segment == 'SDL' or cell.segment == 'mTAL' or cell.segment == 'cTAL' or cell.segment == 'DCT' or cell.segment == 'CNT' or cell.segment == 'CCD' or cell.segment == 'OMCD':
						cell.diam = value*1.42
					else:
						cell.diam = value
				else:
					print('What is the diabete status?')

				if cell.inhib == 'NHE3' and cell.inhib_perc == 0.5:
					if cell.sex == 'male':
						if cell.segment == 'CCD' or cell.segment == 'OMCD' or cell.segment == 'IMCD':
							cell.diam = value*1.05
					if cell.sex == 'female':
						if cell.segment == 'CCD' or cell.segment == 'OMCD' or cell.segment == 'IMCD':
							cell.diam = value*0.95

				if cell.inhib == 'NKCC2' and cell.inhib_perc == 1.0:
					if cell.sex == 'male':
						if cell.segment == 'CCD' or cell.segment == 'OMCD' or cell.segment == 'IMCD':
							cell.diam = value*0.95
					if cell.sex == 'female':
						if cell.segment == 'CCD' or cell.segment == 'OMCD' or cell.segment == 'IMCD':
							cell.diam = value*0.95

			elif compare_string_prefix(id,"Length"):
				if cell.segment == 'LDL' or cell.segment == 'LAL':
					if cell.type == 'jux1':
						looplen = 0.2
					elif cell.type == 'jux2':
						looplen = 0.4
					elif cell.type == 'jux3':
						looplen = 0.6
					elif cell.type == 'jux4':
						looplen = 0.8
					elif cell.type == 'jux5':
						looplen = 1.0

					cell.len = value*looplen
				else:
					cell.len = value

				if cell.type != 'sup':
					if cell.segment == 'cTAL':
						if cell.sex == 'male':
							cell.len = 0.05
						elif cell.sex == 'female':
							cell.len = 0.05*0.85
					elif cell.segment == 'CNT':
						if cell.sex == 'male':
							cell.len = 0.3
						elif cell.sex == 'female':
							cell.len = 0.3*0.85
				
				if cell.diabete == 'No':
					cell.len = cell.len
				elif cell.diabete == 'Yes':
					if cell.segment == 'PT' or cell.segment == 'S3':
						cell.len = cell.len*1.28
					elif cell.segment == 'SDL' or cell.segment == 'mTAL' or cell.segment == 'cTAL' or cell.segment == 'DCT' or cell.segment == 'CNT' or cell.segment == 'CCD' or cell.segment == 'OMCD':
						cell.len = cell.len*1.07
					else:
						cell.len = cell.len
				else:
					print('What is the diabete status?')
			#elif compare_string_prefix(id,"Total"):
			#    cell.total = value
			elif compare_string_prefix(id,"Pressure"):
				cell.pres[0] = value
				if cell.type !='sup' and cell.segment == 'PT':
					if cell.sex == 'male':
						cell.pres[0] = 12.5
					elif cell.sex == 'female':
						cell.pres[0] = 12.5

			elif compare_string_prefix(id,"pH"):
				for i in range(6):
					cell.pH[i] = num[i]
			elif compare_string_prefix(id,"total"):
				cell.total = value
			# surface area
			elif compare_string_prefix(id,"Area"):
				ind1,ind2 = get_interface_id(id)
				cell.area[ind1][ind2] = value
				cell.area[ind2][ind1] = value  # symmetry
				cell.area_init[ind1][ind2] = value
				cell.area_init[ind2][ind1] = value
				if cell.type != 'sup':
					if cell.segment == 'PT' or cell.segment == 'S3':
						cell.area[ind1][ind2] = 1.75*cell.area[ind1][ind2]
						cell.area[ind2][ind1] = cell.area[ind1][ind2]

			# water permeabilities
			elif compare_string_prefix(id,"Pf"):
				ind1,ind2 = get_interface_id(id)
				# Units of dimensional water flux (in 'value'): cm3/s/cm2 epith
				# Non-dimensional factor for water flux: (Pfref)*Vwbar*Cref
				# Calculate non-dimensional dLPV = Pf*Vwbar*Cref / (Pfref*Vwbar*Cref)
				# dLPV = Pf/Pfref
				cell.dLPV[ind1][ind2] = value/Pfref
				# symmetry
				cell.dLPV[ind2][ind1] = value/Pfref
				if cell.segment == 'SDL' and cell.type == 'sup':
					if j>=0.46*cell.total:
						cell.dLPV[0,1]=0.00*cell.dLPV[0,1]
						cell.dLPV[0,4]=0.00*cell.dLPV[0,4]
				elif cell.segment == 'LDL':
					if cell.sex == 'male':
						if j>=0.4*cell.total:
							cell.dLPV[0,1]=0.00*cell.dLPV[0,1]
							cell.dLPV[0,4]=0.00*cell.dLPV[0,4]
					elif cell.sex == 'female':
						if j>=0.5*cell.total:
							cell.dLPV[0,1]=0.00*cell.dLPV[0,1]
							cell.dLPV[0,4]=0.00*cell.dLPV[0,4]
				#elif cell.segment == 'SDL' and cell.type !='sup' and ind1 == 0 and ind2 == 4:
				#	cell.dLPV[0,4]=10*cell.dLPV[0,4]
				
				if cell.segment == 'CNT' and cell.type !='sup':
					if cell.sex == 'female':
						cell.dLPV = cell.dLPV*4/3
				
				if cell.diabete == 'Yes':
					if cell.segment == 'CCD' or cell.segment == 'IMCD':
						cell.dLPV[0,1] = cell.dLPV[0,1]*1.4
						cell.dLPV[1,5] = cell.dLPV[1,5]*1.4
			
			# reflective coefficients
			elif compare_string_prefix(id,"sig"):
				sid,ind1,ind2 = get_interface_id(id)
				cell.sig[sid][ind1][ind2] = value
#                if sid==0 and ind1==0 and ind2==1:
#                    print('got sig[0][0][1]]',cell.sig[0][0][1])
				# symmetry
				cell.sig[sid][ind2][ind1] = value

			# membrane solute permeabilities
			elif compare_string_prefix(id,"perm"):
				sid,ind1,ind2 = get_interface_id(id)
				cell.h[sid][ind1][ind2] = value*1.0e-5/href
				# symmetry
				cell.h[sid][ind2][ind1] = value*1.0e-5/href
				# same basolateral around bath or LIS
				if ind1==1 and ind2==5:
					cell.h[sid][ind1][4] = value*1.0e-5/href
					cell.h[sid][4][ind1] = value*1.0e-5/href
				elif ind1==1 and ind2==4:
					cell.h[sid][ind1][5] = value*1.0e-5/href
					cell.h[sid][5][ind1] = value*1.0e-5/href
				elif ind1==2 and ind2==4:
					cell.h[sid][ind1][5] = value*1.0e-5/href
				elif ind1==3 and ind2==4:
					cell.h[sid][ind1][5] = value*1.0e-5/href
				if cell.segment == 'OMCD':
					cell.h[sid][0][3] = 0.0
					cell.h[sid][3][4] = 0.0
					cell.h[sid][3][5] = 0.0
				if cell.segment == 'IMCD':
					if cell.sex == 'male':
						if j>3*cell.total/4-1:
							cell.h[8,0,1] = 300.0*1.0e-5/href
					elif cell.sex == 'female':
						if j>2*cell.total/3-1:
							cell.h[8,0,1] = 300.0*1.0e-5/href
				if cell.segment == 'LDL':
					if cell.sex == 'male':
						if j>=0.4*cell.total:
							cell.h[0,0,1]=80.0
							cell.h[0,0,4]=80.0
							cell.h[1,0,1]=100.0
							cell.h[1,0,4]=100.0
							cell.h[2,0,1]=80.0
							cell.h[2,0,4]=80.0
							cell.h[3,0,1]=20.0
							cell.h[3,0,4]=20.0
							cell.h[10,0,1]=20.0
							cell.h[10,0,4]=20.0
							cell.h[8,0,1]=80.0 #80
							cell.h[8,0,4]=80.0 #80
					elif cell.sex == 'female':
						if j>=0.5*cell.total:
							cell.h[0,0,1]=40.0
							cell.h[0,0,4]=80.0
							cell.h[1,0,1]=100.0
							cell.h[1,0,4]=100.0
							cell.h[2,0,1]=40.0
							cell.h[2,0,4]=80.0
							cell.h[3,0,1]=20.0
							cell.h[3,0,4]=20.0
							cell.h[10,0,1]=20.0
							cell.h[10,0,4]=20.0
							cell.h[8,0,1]=80.0
							cell.h[8,0,4]=80.0
					
				# if cell.inhib == 'NHE3':
				# 	if cell.segment == 'PT':
				# 		cell.h[8,0,1] = cell.h[8,0,1]*1.2
				# 		cell.h[8,0,4] = cell.h[8,0,4]*1.2
				# 	if cell.segment == 'SDL':
				# 		cell.h[8,0,1] = 3.0
				# 		cell.h[8,0,4] = 3.0
				# 	if cell.segment == 'LDL':
				# 		cell.h[8,0,1] = 0.0
				# 		cell.h[8,0,4] = 0.0
				# 	if cell.segment == 'LAL':
				# 		cell.h[8,0,1] = 0.0
				# 		cell.h[8,0,4] = 0.0
			# coupled transporters
			elif compare_string_prefix(id,"coupled"):
				# retrieve interface and solute id
				vals = get_coupled_id(id)
				newdLA = coupled_transport()
				if cell.segment == 'mTAL' and cell.type !='sup':
					newdLA.perm = value/(href*Cref)#*1.25
				elif cell.segment == 'cTAL' and cell.type !='sup':
					newdLA.perm = value/(href*Cref)#*1.15
				else:
					newdLA.perm = value / (href*Cref)
				newdLA.membrane_id = [vals[0],vals[1]]

				coef = []  # retrieve coupling coefficients
#                print(num)
				for i in range(1,len(num)):
					coef.append(int(num[i]))
				newdLA.coef = coef
				newdLA.solute_id = vals[2:len(vals)]

				if cell.salt == 'Y':
					if cell.segment == 'mTAL' or cell.segment == 'cTAL':
						if newdLA.solute_id[-1]=='HPO4':
							newdLA.perm = (1+1)*value/(href*Cref)

				cell.dLA.append(newdLA)
				# same basolateral around bath or LIS
				#if vals[0]==1 and vals[1]==5:
				#    newdLA2 = coupled_transport()
				#    newdLA2.perm = newdLA.perm
				#    newdLA2.membrane_id = [1,4]
				#    newdLA2.coef = coef
				#    newdLA2.solute_id = vals[2:len(vals)]
				#    cell.dLA.append(newdLA2)
				#elif vals[0]==1 and vals[1]==4: 
				#    newdLA2 = coupled_transport()
				#    newdLA2.perm = newdLA.perm
				#    newdLA2.membrane_id = [1,5]
				#    newdLA2.coef = coef
				#    newdLA2.solute_id = vals[2:len(vals)]
				#    cell.dLA.append(newdLA2)
			# specific transporters
			elif compare_string_prefix(id,"transport"):
				tmp = (id).split('_')
				ind1,ind2 = compart_id[tmp[1]],compart_id[tmp[2]]
				newTransp = transporter()
				newTransp.membrane_id = [ind1,ind2]
				newTransp.type = tmp[3]
				newTransp.act = value/(href*Cref)
				#print('transporter')
				#print(newTransp.membrane_id,newTransp.type,newTransp.act)
				if cell.segment == 'mTAL' and cell.type!='sup':
					newTransp.act = value/(href*Cref)#*1.25
				elif cell.segment == 'mTAL' and cell.type!='sup':
					newTransp.act = value/(href*Cref)#*1.15
				else:
					newTransp.act = value/(href*Cref)				
				#if ind1==1 and ind2==5:
				#    newTransp2 = transporter()
				#    newTransp2.membrane_id = [ind1,4]
				#    newTransp2.type = tmp[3]
				#    newTransp2.act = newTransp.act
				#    cell.trans.append(newTransp2)
				#elif ind1==1 and ind2==4:
				#    newTransp2 = transporter()
				#    newTransp2.membrane_id = [ind1,5]
				#    newTransp2.type = tmp[3]
				#    newTransp2.act = newTransp.act
				#    cell.trans.append(newTransp2)
				if cell.type != 'sup' and cell.sex == 'female':
					if cell.segment == 'mTAL' or cell.segment == 'cTAL':
						if newTransp.type == 'NKCC2A' or newTransp.type == 'NKCC2B' or newTransp.type == 'NKCC2F' or newTransp.type == 'NaKATPase':
							newTransp.act = 1.5*value/(href*Cref)
						elif newTransp.type == 'KCC4':
							newTransp.act = 2.0*value/(href*Cref)

				# NHE3 inhibition parameter settings
				if cell.inhib == 'NHE3':
					if cell.segment == 'PT':
						if newTransp.type == 'NHE3':
							newTransp.act = (1-cell.inhib_perc)*value/(href*Cref)
						#if newTransp.type == 'NaKATPase':
						#	newTransp.act = (1+cell.inhib_perc*0.25)*value/(href*Cref)
					elif cell.segment == 'mTAL' or cell.segment == 'cTAL':
						if newTransp.type == 'NHE3':
							newTransp.act = (1-cell.inhib_perc*0.5)*value/(href*Cref)
						if cell.segment == 'mTAL':
							if newTransp.type == 'NaKATPase':
								newTransp.act = (1+cell.inhib_perc*0.5)*value/(href*Cref)
						if cell.segment == 'cTAL':
							if newTransp.type == 'NaKATPase':
								newTransp.act = (1-cell.inhib_perc*0.5)*value/(href*Cref)
					elif cell.segment == 'DCT':
						if newTransp.type == 'NHE3':
							newTransp.act = (1-cell.inhib_perc)*value/(href*Cref)
				elif cell.inhib == 'NKCC2':
					if cell.segment == 'mTAL' or cell.segment == 'cTAL':
						if newTransp.type == 'NKCC2A' or newTransp.type == 'NKCC2B' or newTransp.type == 'NKCC2F':
							newTransp.act = (1-cell.inhib_perc)*value/(href*Cref)
				elif cell.inhib == 'NCC':
					if cell.segment == 'DCT':
						if newTransp.type == 'NCC':
							newTransp.act = (1-cell.inhib_perc)*value/(href*Cref)
				elif cell.inhib == 'ENaC':
					if newTransp.type == 'ENaC':
						newTransp.act = (1-cell.inhib_perc)*value/(href*Cref)
				elif cell.inhib == 'SNB':
					if cell.segment == 'mTAL' or cell.segment == 'cTAL':
						if newTransp.type == 'NKCC2A' or newTransp.type == 'NKCC2B' or newTransp.type == 'NKCC2F':
							newTransp.act = (1-cell.inhib_perc)*value/(href*Cref)
					if cell.segment == 'DCT':
						if newTransp.type == 'NCC':
							newTransp.act = (1-cell.inhib_perc)*value/(href*Cref)
					if newTransp.type == 'ENaC':
						newTransp.act = (1-cell.inhib_perc)*value/(href*Cref)

				if cell.salt == 'Y':
					if cell.segment == 'mTAL':
						if newTransp.type == 'NaKATPase':
							newTransp.act = (1+0.4)*value/(href*Cref)
					elif cell.segment == 'cTAL':
						if newTransp.type == 'NaKATPase':
							newTransp.act = (1-0.4)*value/(href*Cref)
				# elif cell.inhib == 'nhe3 80%':
				# 	if cell.segment == 'PT':
				# 		if newTransp.type == 'NHE3':
				# 			newTransp.act = 0.2*value/(href*Cref)
				# 	elif cell.segment == 'mTAL' or cell.segment == 'cTAL':
				# 		if newTransp.type == 'NHE3':
				# 			newTransp.act = 0.6*value/(href*Cref)
				# 	elif cell.segment == 'DCT':
				# 		if newTransp.type == 'NHE3':
				# 			newTransp.act = 0.2*value/(href*Cref)

				cell.trans.append(newTransp)
			# solute concentrations
			elif compare_string_prefix(id,"conc"):
				tmp = (id).split('_')
				cell.conc[solute_id[tmp[1]]][0] = float(num[0])
				cell.conc[solute_id[tmp[1]]][1] = float(num[1])
				cell.conc[solute_id[tmp[1]]][4] = float(num[2])
				cell.conc[solute_id[tmp[1]]][5] = float(num[3])
				if len(num) > 4:
					cell.conc[solute_id[tmp[1]]][2] = float(num[4])
					if len(num) > 5:
						cell.conc[solute_id[tmp[1]]][3] = float(num[5])
				if cell.diabete == 'Yes':
					if cell.segment == 'PT':
						cell.conc[14,0] = 25.0

			# reference impermeat concnetration (for cell)
			# or oncotic pressure for lumen/LIS/bath
			elif compare_string_prefix(id,"cimpref"):
				tmp = (id).split('_')
				cell.cimpref[compart_id[tmp[1]]] = float(num[0])

			# reference volumes
			elif compare_string_prefix(id,"volref"):
				tmp = (id).split('_')
				cell.volref[compart_id[tmp[1]]] = float(num[0])


			# actual volumes
			elif compare_string_prefix(id,"vol"):
				tmp = (id).split('_')
				if cell.segment == 'PT' and cell.type != 'sup':
					if compart_id[tmp[1]] == 0:	
						if cell.type == 'jux1':
							if cell.sex == 'male':
								cell.vol[0] = 0.0075
							elif cell.sex == 'female':
								cell.vol[0] = 0.006
							cell.vol_init[0] = cell.vol[0]
						elif cell.type == 'jux2':
							if cell.sex == 'male':
								cell.vol[0] = 0.0075
							elif cell.sex == 'female':
								cell.vol[0] = 0.006
							cell.vol_init[0] = cell.vol[0]

						elif cell.type == 'jux3':
							if cell.sex == 'male':
								cell.vol[0] = 0.0075
							elif cell.sex == 'female':
								cell.vol[0] = 0.006
							cell.vol_init[0] = cell.vol[0]
						elif cell.type == 'jux4':
							if cell.sex == 'male':
								cell.vol[0] = 0.0075
							elif cell.sex == 'female':
								cell.vol[0] = 0.006
							cell.vol_init[0] = cell.vol[0]
						elif cell.type == 'jux5':
							if cell.sex == 'male':
								cell.vol[0] = 0.0075
							elif cell.sex == 'female':
								cell.vol[0] = 0.006
							cell.vol_init[0] = cell.vol[0]
					else:
						cell.vol[compart_id[tmp[1]]] = float(num[0])
						cell.vol_init[compart_id[tmp[1]]] = float(num[0])
				else:
					cell.vol[compart_id[tmp[1]]] = float(num[0])
					cell.vol_init[compart_id[tmp[1]]] = float(num[0])

				if cell.inhib == 'NHE3' and cell.inhib_perc == 0.5:
					if cell.segment == 'PT' and compart_id[tmp[1]] == 0:
						if cell.type == 'sup':
							if cell.sex == 'male':
								cell.vol[0] = 0.004599
							elif cell.sex == 'female':
								cell.vol[0] = 0.00362
							cell.vol_init[0] = cell.vol[0]
						elif cell.type == 'jux1':
							if cell.sex == 'male':
								cell.vol[0] = 0.0065
							elif cell.sex == 'female':
								cell.vol[0] = 0.00515
							cell.vol_init[0] = cell.vol[0]
						elif cell.type == 'jux2':
							if cell.sex == 'male':
								cell.vol[0] = 0.0065
							elif cell.sex == 'female':
								cell.vol[0] = 0.0052
							cell.vol_init[0] = cell.vol[0]
						elif cell.type == 'jux3':
							if cell.sex == 'male':
								cell.vol[0] = 0.0065
							elif cell.sex == 'female':
								cell.vol[0] = 0.0052
							cell.vol_init[0] = cell.vol[0]
						elif cell.type == 'jux4':
							if cell.sex == 'male':
								cell.vol[0] = 0.0065
							elif cell.sex == 'female':
								cell.vol[0] = 0.00521
							cell.vol_init[0] = cell.vol[0]
						elif cell.type == 'jux5':
							if cell.sex == 'male':
								cell.vol[0] = 0.0065
							elif cell.sex == 'female':
								cell.vol[0] = 0.00522
							cell.vol_init[0] = cell.vol[0]
				elif cell.inhib == 'NHE3' and cell.inhib_perc == 0.8:
					if cell.segment == 'PT' and compart_id[tmp[1]] == 0:
						if cell.type == 'sup':
							if cell.sex == 'male':
								cell.vol[0] = float(num[0])*0.725
							elif cell.sex == 'female':
								cell.vol[0] = 0.0033
							cell.vol_init[0] = cell.vol[0]
						elif cell.type == 'jux1':
							if cell.sex == 'male':
								cell.vol[0] = float(num[0])*1.41*0.75
							elif cell.sex == 'female':
								cell.vol[0] = 0.0046
							cell.vol_init[0] = cell.vol[0]
						elif cell.type == 'jux2':
							if cell.sex == 'male':
								cell.vol[0] = float(num[0])*1.41*0.75
							elif cell.sex == 'female':
								cell.vol[0] = 0.0046
							cell.vol_init[0] = cell.vol[0]
						elif cell.type == 'jux3':
							if cell.sex == 'male':
								cell.vol[0] = float(num[0])*1.396*0.75
							elif cell.sex == 'female':
								cell.vol[0] = 0.0046
							cell.vol_init[0] = cell.vol[0]
						elif cell.type == 'jux4':
							if cell.sex == 'male':
								cell.vol[0] = float(num[0])*1.393*0.75
							elif cell.sex == 'female':
								cell.vol[0] = 0.0046
							cell.vol_init[0] = cell.vol[0]
						elif cell.type == 'jux5':
							if cell.sex == 'male':
								cell.vol[0] = float(num[0])*1.389*0.75
							elif cell.sex == 'female':
								cell.vol[0] = 0.0046
							cell.vol_init[0] = cell.vol[0]


			# membrane potential
			elif compare_string_prefix(id,"ep"):
				tmp = (id).split('_')
				cell.ep[compart_id[tmp[1]]] = float(num[0])
						
			# invalue keyword
#            else:
#                print("Wrong id",id)
				# if not a comment line
		
		line = file.readline()

		for i in range(NS):
			cell.sig[i,4,5] = 0
			cell.sig[i,5,4] = 0

		if cell.segment=='cTAL' or cell.segment=='mTAL':
			cell.dkd[2]=0
			cell.dkd[3]=0
			cell.dkd[5]=0
			cell.dkh[2]=0
			cell.dkh[3]=0
			cell.dkh[5]=0
		elif cell.segment=='DCT':
			cell.dkd = np.zeros(6)
			cell.dkh = np.zeros(6)
			cell.dkd[0]=496
			cell.dkh[0]=1.45
			cell.dkd[1]=49600
			cell.dkh[1]=145
			cell.dkd[4]=49600
			cell.dkh[4]=145
		elif cell.segment=='PT':
			cell.dkd = 496000 * np.ones(NC)
			cell.dkh = 1450 * np.ones(NC)
		elif cell.segment=='S3':
			cell.dkd = 4960 * np.ones(NC)
			cell.dkh = 14.5 * np.ones(NC)
		elif cell.segment=='CNT' or cell.segment == 'CCD':
			cell.dkd = np.zeros(6)
			cell.dkh = np.zeros(6)
			cell.dkd[0]=496
			cell.dkh[0]=1.45
			cell.dkd[1]=496
			cell.dkh[1]=1.45
			cell.dkd[2]=496000
			cell.dkh[2]=1450
			cell.dkd[3]=496000
			cell.dkh[3]=1450
			cell.dkd[4]=49600
			cell.dkh[4]=145
		elif cell.segment=='OMCD':
			cell.dkd = np.zeros(6)
			cell.dkh = np.zeros(6)
			cell.dkd[0]=49.6
			cell.dkh[0]=0.145
			cell.dkd[1]=496
			cell.dkh[1]=1.45
			cell.dkd[2]=496000
			cell.dkh[2]=1450
			cell.dkd[3]=496000
			cell.dkh[3]=1450
			cell.dkd[4]=49.6
			cell.dkh[4]=0.145			
		elif cell.segment == 'IMCD':
			cell.dkd = np.zeros(6)
			cell.dkh = np.zeros(6)
			cell.dkd[0]=49.6
			cell.dkh[0]=0.145
			cell.dkd[1]=4960
			cell.dkh[1]=14.5
			cell.dkd[2]=4960
			cell.dkh[2]=14.5
			cell.dkd[3]=4960
			cell.dkh[3]=14.5
			cell.dkd[4]=49.6
			cell.dkh[4]=0.145
	file.close()
