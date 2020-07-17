import equations
import math
import numpy as np

def Jac(func,x,k):

	epsfcn = 1.0e-3
	epsmch = 1.0e-3
	eps = math.sqrt(max(epsfcn,epsmch))
	
	Jfun=[[0 for i in range(len(x))] for i in range (len(x))]
	
	
	wa1=func(x,k)
	for i in range(len(x)):
		temp=x[i]
		h=eps*abs(temp)
		if (h==0):
			h=eps
		x[i]=temp+h
		fvec=func(x,k)
		x[i]=temp
		for j in range(len(x)):
			Jfun[j][i]=(-wa1[j]+fvec[j])/h
	
	return Jfun
	
	
def newton(func,x,k,type,sup_jux,gender=None,inhib=None,perc=0,salt=None):
	fun=equations.conservation_eqs
	f = np.matrix(fun(x,k))
	print(np.linalg.norm(f))
	TOLpcn = 1
	i = 1
	iter=0
	while(np.linalg.norm(f) > 0.0001) and (iter<150):
#	    print("Iteration Times: " + str(i) + " with TOL " + str(TOLpcn) + "%")
		i += 1
		J = np.matrix(Jac(fun,x,k))
		IJ = J.I
		F = np.matrix(fun(x,k))
		if salt == 'Y':
			if type=='DCT':
				amp = 1
			elif type == 'CNT':
				if gender =='female':
					if sup_jux =='sup':
						if np.linalg.norm(f)>100: # sup: 100							
							amp = 0.17# sup: 0.5
						else:
							amp=1.0
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							if k == 0:
								amp = 0.15#0.1
							else:
								amp = 0.17
						else:
							amp = 1.0
					elif sup_jux == 'jux2':
						if np.linalg.norm(f)>5000:
							amp = 0.15#0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.15#0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
				elif gender =='male':
					if sup_jux =='sup':
						if np.linalg.norm(f)>5000: # sup: 100
							amp = 0.3# 0.5;
						else:
							amp=0.8
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							amp = 0.5 # nhe50: 0.3; 
						else:
							amp = 1.0
					elif sup_jux == 'jux2': 
						if np.linalg.norm(f)>5000:
							amp = 0.17 #0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.17 #0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.17 # nhe50: 0.15
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.17 # nhe50: 0.15 nhe80:0.1
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
			elif type == 'SDL':
				amp = 0.2
			elif type == 'IMCD':
				if gender == 'female':
					if np.linalg.norm(f)>1000: # 100
						if k==0:
							amp = 0.2
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8#1.0 male:0.7 female:0.8
						else:
							amp = 0.8
				elif gender == 'male':
					if np.linalg.norm(f)>5000:
						if k==0:
							amp = 0.17 # saline: 0.17
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8
						else:
							amp = 1.0
			elif type == 'CCD':
				if np.linalg.norm(f)>1000:
					if k == 0:
						amp = 1.0#0.005 male:0.5 female:0.005
					else:
						amp = 0.1
				else:
					amp = 0.8
			elif type == 'OMCD':
				if np.linalg.norm(f)>100:
					amp = 0.8#0.8 #normal male/female: 1.0, diabetic male: 0.8
				else:
					amp = 1.0# normal male: 1.0, diabetic male: 0.8
			elif type == 'cTAL' or type == 'MD':
				if np.linalg.norm(f)>100:
					amp = 0.2
				else:
					amp = 0.8 # normal male: 1.0
			elif type == 'mTAL':
				if np.linalg.norm(f)>100:
					amp = 0.2 #0.2
				else:
					amp = 1.0
			elif type == 'LDL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			elif type == 'LAL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			else:
				amp = 1
		
		elif inhib == 'NHE3' and perc == 0.5:
			if type=='DCT':
				amp = 1
			elif type == 'CNT':
				if gender =='female':
					if sup_jux =='sup':
						if np.linalg.norm(f)>100: # sup: 100							
							amp = 0.17# sup: 0.5
						else:
							amp=1.0
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							if k == 0:
								amp = 0.05
							else:
								amp = 0.17
						else:
							amp = 1.0
					elif sup_jux == 'jux2':
						if np.linalg.norm(f)>5000:
							amp = 0.05
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.05
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.05
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.05
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
				elif gender =='male':
					if sup_jux =='sup':
						if np.linalg.norm(f)>5000:
							amp = 0.5
						else:
							amp=0.8
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							amp = 0.13  
						else:
							amp = 1.0
					elif sup_jux == 'jux2': 
						if np.linalg.norm(f)>5000:
							amp = 0.13 
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.13 
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
			elif type == 'SDL':
				amp = 0.2
			elif type == 'IMCD':
				if gender == 'female':
					if np.linalg.norm(f)>1000: 
						if k==0:
							amp = 0.2
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8#1.0 male:0.7 female:0.8
						else:
							amp = 0.8
				elif gender == 'male':
					if np.linalg.norm(f)>5000:
						if k==0:
							amp = 0.17 # saline: 0.17
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8
						else:
							amp = 1.0
			elif type == 'CCD':
				if np.linalg.norm(f)>1000:
					if k == 0:
						amp = 1.0#0.005 male:0.5 female:0.005
					else:
						amp = 0.1
				else:
					amp = 0.8
			elif type == 'OMCD':
				if np.linalg.norm(f)>100:
					amp = 0.8#0.8 #normal male/female: 1.0, diabetic male: 0.8
				else:
					amp = 1.0# normal male: 1.0, diabetic male: 0.8
			elif type == 'cTAL' or type == 'MD':
				if np.linalg.norm(f)>100:
					amp = 0.2
				else:
					amp = 0.8 # normal male: 1.0
			elif type == 'mTAL':
				if np.linalg.norm(f)>100:
					amp = 0.2 #0.2
				else:
					amp = 1.0
			elif type == 'LDL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			elif type == 'LAL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			else:
				amp = 1

		elif inhib == 'NHE3' and perc == 0.8:
			if type=='DCT':
				amp = 1
			elif type == 'CNT':
				if gender =='female':
					if sup_jux =='sup':
						if np.linalg.norm(f)>100: # sup: 100							
							amp = 0.17
						else:
							amp=1.0
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							if k == 0:
								amp = 0.17
							else:
								amp = 0.17
						else:
							amp = 1.0
					elif sup_jux == 'jux2':
						if np.linalg.norm(f)>5000:
							amp = 0.17
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.17
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
				elif gender =='male':
					if sup_jux =='sup':
						if np.linalg.norm(f)>5000: # sup: 100
							amp = 0.5# sup: 0.5; saline: 0.5;
						else:
							amp=0.8
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							amp = 0.3
						else:
							amp = 1.0
					elif sup_jux == 'jux2': 
						if np.linalg.norm(f)>5000:
							amp = 0.15 # nhe50: 0.15 ncc: 0.05
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.15 # nhe50: 0.15
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.15 # nhe50: 0.15
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.1 # nhe50: 0.15 nhe80:0.1
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
			elif type == 'SDL':
				amp = 0.2
			elif type == 'IMCD':
				if gender == 'female':
					if np.linalg.norm(f)>1000: # 100
						if k==0:
							amp = 0.2
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8#1.0 male:0.7 female:0.8
						else:
							amp = 0.8
				elif gender == 'male':
					if np.linalg.norm(f)>5000:
						if k==0:
							amp = 0.17 # saline: 0.17
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8
						else:
							amp = 1.0
			elif type == 'CCD':
				if np.linalg.norm(f)>1000:
					if k == 0:
						amp = 0.8#0.005 male:0.5 female:0.005
					else:
						amp = 0.1
				else:
					amp = 0.8
			elif type == 'OMCD':
				if np.linalg.norm(f)>100:
					amp = 0.8#0.8 #normal male/female: 1.0, diabetic male: 0.8
				else:
					amp = 1.0# normal male: 1.0, diabetic male: 0.8
			elif type == 'cTAL' or type == 'MD':
				if np.linalg.norm(f)>100:
					amp = 0.2
				else:
					amp = 0.8 # normal male: 1.0
			elif type == 'mTAL':
				if np.linalg.norm(f)>100:
					amp = 0.2 #0.2
				else:
					amp = 1.0
			elif type == 'LDL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			elif type == 'LAL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			else:
				amp = 1

		elif inhib == 'NKCC2' and perc == 0.7:
			if type=='DCT':
				amp = 1
			elif type == 'CNT':
				if gender =='female':
					if sup_jux =='sup':
						if np.linalg.norm(f)>100: # sup: 100							
							amp = 0.17# sup: 0.5
						else:
							amp=1.0
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							if k == 0:
								amp = 0.13
							else:
								amp = 0.17
						else:
							amp = 1.0
					elif sup_jux == 'jux2':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
				elif gender =='male':
					if sup_jux =='sup':
						if np.linalg.norm(f)>5000: # sup: 100
							amp = 0.5# sup: 0.5; saline: 0.5;
						else:
							amp=0.8
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux2': 
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
			elif type == 'SDL':
				amp = 0.2
			elif type == 'IMCD':
				if gender == 'female':
					if np.linalg.norm(f)>1000: # 100
						if k==0:
							amp = 0.2
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8#1.0 male:0.7 female:0.8
						else:
							amp = 0.8
				elif gender == 'male':
					if np.linalg.norm(f)>5000:
						if k==0:
							amp = 0.17 # saline: 0.17
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8
						else:
							amp = 1.0
			elif type == 'CCD':
				if np.linalg.norm(f)>1000:
					if k == 0:
						amp = 1.0#0.005 male:0.5 female:0.005
					else:
						amp = 0.1
				else:
					amp = 0.8
			elif type == 'OMCD':
				if np.linalg.norm(f)>100:
					amp = 0.8#0.8 #normal male/female: 1.0, diabetic male: 0.8
				else:
					amp = 1.0# normal male: 1.0, diabetic male: 0.8
			elif type == 'cTAL' or type == 'MD':
				if np.linalg.norm(f)>100:
					amp = 0.2
				else:
					amp = 0.8 # normal male: 1.0
			elif type == 'mTAL':
				if np.linalg.norm(f)>100:
					amp = 0.2 #0.2
				else:
					amp = 1.0
			elif type == 'LDL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			elif type == 'LAL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			else:
				amp = 1

		elif inhib == 'NKCC2' and perc == 1.0:
			if type=='DCT':
				amp = 1
			elif type == 'CNT':
				if gender =='female':
					if sup_jux =='sup':
						if np.linalg.norm(f)>100: # sup: 100							
							amp = 0.17# sup: 0.5
						else:
							amp=1.0
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							if k == 0:
								amp = 0.2
							else:
								amp = 0.17
						else:
							amp = 1.0
					elif sup_jux == 'jux2':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
				elif gender =='male':
					if sup_jux =='sup':
						if np.linalg.norm(f)>5000: # sup: 100
							amp = 0.17
						else:
							amp=0.8
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							amp = 0.3 # nhe50: 0.3; 
						else:
							amp = 1.0
					elif sup_jux == 'jux2': 
						if np.linalg.norm(f)>5000:
							amp = 0.3
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.3
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.3
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.3
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
			elif type == 'SDL':
				amp = 0.2
			elif type == 'IMCD':
				if gender == 'female':
					if np.linalg.norm(f)>1000: # 100
						if k==0:
							amp = 0.2
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8#1.0 male:0.7 female:0.8
						else:
							amp = 0.8
				elif gender == 'male':
					if np.linalg.norm(f)>5000:
						if k==0:
							amp = 0.17 # saline: 0.17
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8
						else:
							amp = 1.0
			elif type == 'CCD':
				if np.linalg.norm(f)>1000:
					if k == 0:
						amp = 1.0#0.005 male:0.5 female:0.005
					else:
						amp = 0.1
				else:
					amp = 0.8
			elif type == 'OMCD':
				if np.linalg.norm(f)>100:
					amp = 0.8#0.8 #normal male/female: 1.0, diabetic male: 0.8
				else:
					amp = 1.0# normal male: 1.0, diabetic male: 0.8
			elif type == 'cTAL' or type == 'MD':
				if np.linalg.norm(f)>100:
					amp = 0.2
				else:
					amp = 0.8 # normal male: 1.0
			elif type == 'mTAL':
				if np.linalg.norm(f)>100:
					amp = 0.2 #0.2
				else:
					amp = 1.0
			elif type == 'LDL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			elif type == 'LAL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			else:
				amp = 1
		
		elif inhib == 'NCC' and perc == 0.7:
			if type=='DCT':
				amp = 1
			elif type == 'CNT':
				if gender =='female':
					if sup_jux =='sup':
						if np.linalg.norm(f)>100: # sup: 100							
							amp = 0.17# sup: 0.5
						else:
							amp=1.0
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							if k == 0:
								amp = 0.1
							else:
								amp = 0.17
						else:
							amp = 1.0
					elif sup_jux == 'jux2':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
				elif gender =='male':
					if sup_jux =='sup':
						if np.linalg.norm(f)>5000:
							amp = 0.17
						else:
							amp=0.8
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux2': 
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
			elif type == 'SDL':
				amp = 0.2
			elif type == 'IMCD':
				if gender == 'female':
					if np.linalg.norm(f)>1000: # 100
						if k==0:
							amp = 0.2
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8#1.0 male:0.7 female:0.8
						else:
							amp = 0.8
				elif gender == 'male':
					if np.linalg.norm(f)>5000:
						if k==0:
							amp = 0.17 # saline: 0.17
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8
						else:
							amp = 1.0
			elif type == 'CCD':
				if np.linalg.norm(f)>1000:
					if k == 0:
						amp = 1.0#0.005 male:0.5 female:0.005
					else:
						amp = 0.1
				else:
					amp = 0.8
			elif type == 'OMCD':
				if np.linalg.norm(f)>100:
					amp = 0.8#0.8 #normal male/female: 1.0, diabetic male: 0.8
				else:
					amp = 1.0# normal male: 1.0, diabetic male: 0.8
			elif type == 'cTAL' or type == 'MD':
				if np.linalg.norm(f)>100:
					amp = 0.2
				else:
					amp = 0.8 # normal male: 1.0
			elif type == 'mTAL':
				if np.linalg.norm(f)>100:
					amp = 0.2 #0.2
				else:
					amp = 1.0
			elif type == 'LDL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			elif type == 'LAL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			else:
				amp = 1
		elif inhib == 'NCC' and perc == 1.0:
			if type=='DCT':
				amp = 1
			elif type == 'CNT':
				if gender =='female':
					if sup_jux =='sup':
						if np.linalg.norm(f)>100: # sup: 100							
							amp = 0.17
						else:
							amp=1.0
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							if k == 0:
								amp = 0.1
							else:
								amp = 0.17
						else:
							amp = 1.0
					elif sup_jux == 'jux2':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
				elif gender =='male':
					if sup_jux =='sup':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp=1.0
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							if k == 0:
								amp = 0.05
							else:
								amp = 1.0
						else:
							amp = 1.0
					elif sup_jux == 'jux2': 
						if np.linalg.norm(f)>5000:
							if k==0:
								amp = 0.05
							else:
								amp =0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							if k==0:
								amp = 0.05
							else:
								amp =0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							if k==0:
								amp = 0.05
							else:
								amp =0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							if k==0:
								amp = 0.05
							else:
								amp =0.1
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
			elif type == 'SDL':
				amp = 0.2
			elif type == 'IMCD':
				if gender == 'female':
					if np.linalg.norm(f)>1000: # 100
						if k==0:
							amp = 0.2
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8#1.0 male:0.7 female:0.8
						else:
							amp = 0.8
				elif gender == 'male':
					if np.linalg.norm(f)>5000:
						if k==0:
							amp = 0.17
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8
						else:
							amp = 1.0
			elif type == 'CCD':
				if np.linalg.norm(f)>1000:
					if k == 0:
						amp = 1.0#0.005 male:0.5 female:0.005
					else:
						amp = 0.1
				else:
					amp = 0.8
			elif type == 'OMCD':
				if gender == 'female':
					if np.linalg.norm(f)>100:
						amp = 0.8#0.8 #normal male/female: 1.0, diabetic male: 0.8
					else:
						amp = 1.0# normal male: 1.0, diabetic male: 0.8
				elif gender == 'male':
					if np.linalg.norm(f)>100:
						amp = 0.5
					else:
						amp = 1.0
			elif type == 'cTAL' or type == 'MD':
				if np.linalg.norm(f)>100:
					amp = 0.2
				else:
					amp = 0.8 # normal male: 1.0
			elif type == 'mTAL':
				if np.linalg.norm(f)>100:
					amp = 0.2 #0.2
				else:
					amp = 1.0
			elif type == 'LDL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			elif type == 'LAL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			else:
				amp = 1
		elif inhib == 'ENaC' and perc == 0.7:
			if type=='DCT':
				amp = 1
			elif type == 'CNT':
				if gender =='female':
					if sup_jux =='sup':
						if np.linalg.norm(f)>100: # sup: 100							
							amp = 0.17# sup: 0.5
						else:
							amp=1.0
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							if k == 0:
								amp = 0.05
							else:
								amp = 0.17
						else:
							amp = 1.0
					elif sup_jux == 'jux2':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
				elif gender =='male':
					if sup_jux =='sup':
						if np.linalg.norm(f)>5000: # sup: 100
							amp = 0.5# sup: 0.5; saline: 0.5;
						else:
							amp=0.8
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							amp = 0.13 # nhe50: 0.3; 
						else:
							amp = 1.0
					elif sup_jux == 'jux2': 
						if np.linalg.norm(f)>5000:
							amp = 0.13 # nhe50: 0.15 ncc: 0.05
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.13 # nhe50: 0.15
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.13 # nhe50: 0.15
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.13 # nhe50: 0.15 nhe80:0.1
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
			elif type == 'SDL':
				amp = 0.2
			elif type == 'IMCD':
				if gender == 'female':
					if np.linalg.norm(f)>1000: # 100
						if k==0:
							amp = 0.2
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8#1.0 male:0.7 female:0.8
						else:
							amp = 0.8
				elif gender == 'male':
					if np.linalg.norm(f)>5000:
						if k==0:
							amp = 0.17 # saline: 0.17
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8
						else:
							amp = 1.0
			elif type == 'CCD':
				if np.linalg.norm(f)>1000:
					if k == 0:
						amp = 1.0#0.005 male:0.5 female:0.005
					else:
						amp = 0.1
				else:
					amp = 0.8
			elif type == 'OMCD':
				if np.linalg.norm(f)>100:
					amp = 0.8#0.8 #normal male/female: 1.0, diabetic male: 0.8
				else:
					amp = 1.0# normal male: 1.0, diabetic male: 0.8
			elif type == 'cTAL' or type == 'MD':
				if np.linalg.norm(f)>100:
					amp = 0.2
				else:
					amp = 0.8 # normal male: 1.0
			elif type == 'mTAL':
				if np.linalg.norm(f)>100:
					amp = 0.2 #0.2
				else:
					amp = 1.0
			elif type == 'LDL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			elif type == 'LAL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			else:
				amp = 1
		elif inhib == 'ENaC' and perc == 1.0:
			if type=='DCT':
				amp = 1
			elif type == 'CNT':
				if gender =='female':
					if sup_jux =='sup':
						if np.linalg.norm(f)>100: # sup: 100							
							amp = 0.17# sup: 0.5
						else:
							amp=1.0
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							if k == 0:
								amp = 0.05
							else:
								amp = 0.17
						else:
							amp = 1.0
					elif sup_jux == 'jux2':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.05
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.05
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
				elif gender =='male':
					if sup_jux =='sup':
						if np.linalg.norm(f)>5000: # sup: 100
							amp = 0.5# sup: 0.5; saline: 0.5;
						else:
							amp=0.8
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							amp = 0.13 # nhe50: 0.3; 
						else:
							amp = 1.0
					elif sup_jux == 'jux2': 
						if np.linalg.norm(f)>5000:
							amp = 0.13 # nhe50: 0.15 ncc: 0.05
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
			elif type == 'SDL':
				amp = 0.2
			elif type == 'IMCD':
				if gender == 'female':
					if np.linalg.norm(f)>1000: # 100
						if k==0:
							amp = 0.2
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8#1.0 male:0.7 female:0.8
						else:
							amp = 0.8
				elif gender == 'male':
					if np.linalg.norm(f)>5000:
						if k==0:
							amp = 0.17 # saline: 0.17
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8
						else:
							amp = 1.0
			elif type == 'CCD':
				if gender == 'male':
					if np.linalg.norm(f)>1000:
						if k == 0:
							amp = 0.5
						else:
							amp = 0.1
					else:
						amp = 1.0
				elif gender == 'female':
					if np.linalg.norm(f)>1000:
						if k == 0:
							amp = 1.0
						else:
							amp = 0.1
					else:
						amp = 1.0
			elif type == 'OMCD':
				if np.linalg.norm(f)>100:
					if gender == 'male':
						amp = 0.5
					elif gender == 'female':
						amp = 0.5
				else:
					amp = 1.0
			elif type == 'cTAL' or type == 'MD':
				if np.linalg.norm(f)>100:
					amp = 0.2
				else:
					amp = 0.8 # normal male: 1.0
			elif type == 'mTAL':
				if np.linalg.norm(f)>100:
					amp = 0.2 #0.2
				else:
					amp = 1.0
			elif type == 'LDL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			elif type == 'LAL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			else:
				amp = 1
		elif inhib == 'SNB':
			if type=='DCT':
				amp = 1
			elif type == 'CNT':
				if gender =='female':
					if sup_jux =='sup':
						if np.linalg.norm(f)>100: # sup: 100							
							amp = 0.17# sup: 0.5
						else:
							amp=1.0
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							if k == 0:
								amp = 0.05
							else:
								amp = 0.17
						else:
							amp = 1.0
					elif sup_jux == 'jux2':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.05
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.05
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
				elif gender =='male':
					if sup_jux =='sup':
						if np.linalg.norm(f)>1000: # sup: 100
							amp = 0.1# sup: 0.5; saline: 0.5;
						else:
							amp=1.0
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>1000:
							amp = 0.13 # nhe50: 0.3; 
						else:
							amp = 1.0
					elif sup_jux == 'jux2': 
						if np.linalg.norm(f)>5000:
							amp = 0.13 # nhe50: 0.15 ncc: 0.05
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
			elif type == 'SDL':
				amp = 0.2
			elif type == 'IMCD':
				if gender == 'female':
					if np.linalg.norm(f)>1000: # 100
						if k==0:
							amp = 0.2
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8#1.0 male:0.7 female:0.8
						else:
							amp = 0.8
				elif gender == 'male':
					if np.linalg.norm(f)>5000:
						if k==0:
							amp = 0.17 # saline: 0.17
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8
						else:
							amp = 1.0
			elif type == 'CCD':
				if gender == 'male':
					if np.linalg.norm(f)>1000:
						if k == 0:
							amp = 0.5
						else:
							amp = 0.1
					else:
						amp = 1.0
				elif gender == 'female':
					if np.linalg.norm(f)>1000:
						if k == 0:
							amp = 1.0
						else:
							amp = 0.1
					else:
						amp = 1.0
			elif type == 'OMCD':
				if np.linalg.norm(f)>100:
					if gender == 'male':
						amp = 0.5
					elif gender == 'female':
						amp = 0.5
				else:
					amp = 1.0
			elif type == 'cTAL' or type == 'MD':
				if np.linalg.norm(f)>100:
					amp = 0.2
				else:
					amp = 0.8 # normal male: 1.0
			elif type == 'mTAL':
				if np.linalg.norm(f)>100:
					amp = 0.2 #0.2
				else:
					amp = 1.0
			elif type == 'LDL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			elif type == 'LAL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			else:
				amp = 1		
		else:
			if type=='DCT':
				amp = 1
			elif type == 'CNT':
				if gender =='female':
					if sup_jux =='sup':
						if np.linalg.norm(f)>100: # sup: 100							
							amp = 0.17# sup: 0.5
						else:
							amp=1.0
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							if k == 0:
								amp = 0.13
							else:
								amp = 0.17
						else:
							amp = 1.0
					elif sup_jux == 'jux2':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.13
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.1
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
				elif gender =='male':
					if sup_jux =='sup':
						if np.linalg.norm(f)>5000: # sup: 100
							amp = 0.5# sup: 0.5; saline: 0.5;
						else:
							amp=0.8
					elif sup_jux == 'jux1':
						if np.linalg.norm(f)>5000:
							amp = 0.13 # nhe50: 0.3; 
						else:
							amp = 1.0
					elif sup_jux == 'jux2': 
						if np.linalg.norm(f)>5000:
							amp = 0.13 # nhe50: 0.15 ncc: 0.05
						else:
							amp = 1.0
					elif sup_jux == 'jux3':
						if np.linalg.norm(f)>5000:
							amp = 0.13 # nhe50: 0.15
						else:
							amp = 1.0
					elif sup_jux == 'jux4':
						if np.linalg.norm(f)>5000:
							amp = 0.13 # nhe50: 0.15
						else:
							amp = 1.0
					elif sup_jux == 'jux5':
						if np.linalg.norm(f)>5000:
							amp = 0.13 # nhe50: 0.15 nhe80:0.1
						else:
							amp = 1.0
					else:
						if np.linalg.norm(f)>5000:
							amp = 0.15
						else:
							amp = 1.0
			elif type == 'SDL':
				amp = 0.2
			elif type == 'IMCD':
				if gender == 'female':
					if np.linalg.norm(f)>1000: # 100
						if k==0:
							amp = 0.2
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8#1.0 male:0.7 female:0.8
						else:
							amp = 0.8
				elif gender == 'male':
					if np.linalg.norm(f)>5000:
						if k==0:
							amp = 0.17 # saline: 0.17
						else:
							amp = 0.2
					else:
						if k==0:
							amp = 0.8
						else:
							amp = 1.0
			elif type == 'CCD':
				if np.linalg.norm(f)>1000:
					if k == 0:
						amp = 1.0#0.005 male:0.5 female:0.005
					else:
						amp = 0.1
				else:
					amp = 0.8
			elif type == 'OMCD':
				if np.linalg.norm(f)>100:
					amp = 0.8#0.8 #normal male/female: 1.0, diabetic male: 0.8
				else:
					amp = 1.0# normal male: 1.0, diabetic male: 0.8
			elif type == 'cTAL' or type == 'MD':
				if np.linalg.norm(f)>100:
					amp = 0.2
				else:
					amp = 0.8 # normal male: 1.0
			elif type == 'mTAL':
				if np.linalg.norm(f)>100:
					amp = 0.2 #0.2
				else:
					amp = 1.0
			elif type == 'LDL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			elif type == 'LAL':
				if np.linalg.norm(f)>5000:
					amp = 0.5
				else:
					amp = 1.0
			else:
				amp = 1
		#print(amp)
		delta = amp*np.array(F * IJ.T)[0]
		x -= delta
		f = np.matrix(fun(x,k))
		iter+=1
		print(iter,np.linalg.norm(f))
		TOLpcn = np.max(delta / x)
		#print(i)

#    print("Iteration Times: " + str(i) + " with TOL " + str(TOLpcn) + "%")
	
	return x

	
def broyden(func,x,k,type):
	fun=equations.conservation_eqs
	f=np.matrix(fun(x,k))
	J=np.matrix(Jac(fun,x,k))
	IJ=J.I
	dx=np.ones(x.shape)
	i=0
	iter=0
	#while(np.max(dx)>0.0001):
	while(np.linalg.norm(f)>0.0001) and (iter<500):
		
		f_previous=f
		x_previous=x
		
		x=x-np.array(f*IJ.T)[0]
		
		f=np.matrix(fun(x,k))
#                print('x',x)
#                print('f',f)

		df=f-f_previous
		dx=x-x_previous
		# #
		# #-------------------------------------------------------
		# #using good broyden
		# dx=np.array([dx])
		# df=np.array([df])
		# dx=dx.T
		# df=df.T

		# IJ = IJ+(dx-IJ*df)*(dx.T*IJ)/np.inner(dx.T*IJ,df)
		# #-------------------------------------------------------

		IJ=IJ-np.outer(IJ*f.T,dx)*IJ/np.inner(dx,dx+(IJ*f.T).T)
		#print(i)
		iter+=1
		print(iter,np.linalg.norm(f))

		#J=J+np.outer((df-dx*J.T),dx)/np.linalg.norm(x)**2
		#IJ=J.I

	return x
