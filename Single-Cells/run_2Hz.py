

from subprocess import call
import math
# S1 = 500
# S2 = 250
import sys


import numpy as np






# def run(Para1, Para2, Para3, S2_amp):
def run(Para1,kmfscale,ISO, f):
	# global mut
	#call(["./main","BCL", str(S1), "S2", str(S2), "Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./SAN', str(P1), str(P2)]) #"Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./NCZ_Model_New_Ito', '250', '14', '250', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	# call(['./NCZ_Model_New_Ito', '500', '14', '500', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	P=['./GB_Main', '500', str(kmfscale),str(ISO)]
	for i in range(len(Para1)):
		P.append(str(Para1[i]))
	# for i in range(len(Para2)):
	# 	P.append(str(Para2[i]))
	# for i in range(len(Para3)):  
	# 	P.append(str(Para3[i]))
	# P.append(str(-20))
	call(P, stdout=f)


def run_Restart(Para1, Para2, Para3, S2_amp):
	# global mut
	#call(["./main","BCL", str(S1), "S2", str(S2), "Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./SAN', str(P1), str(P2)]) #"Mutation", mut, "S1_number", "20", "ISO", ISO])  #
	# call(['./NCZ_Model_New_Ito', '250', '14', '250', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	# call(['./NCZ_Model_New_Ito', '500', '14', '500', 'WT', 'Normal', '0', '0', str(P1), str(P2), str(P3),str(P4),str(P5),str(P6)], stdout=f)
	P=['./NCZ_Model_Threshold_long_5ms', '1000', '14', '1000', 'WT', 'Restart', '0', '0']
	for i in range(len(Para1)):
		P.append(str(Para1[i]))
	for i in range(len(Para2)):
		P.append(str(Para2[i]))
	for i in range(len(Para3)):  
		P.append(str(Para3[i]))
	P.append(str(S2_amp))
	call(P, stdout=f)


# Data2 = np.loadtxt('AP.log.dat.1Hz', unpack=True)
def check_results():
	Data = np.loadtxt('S2_AP.log', unpack=False)
	i=0;
	if Data[1] > 317.41 - 43.19*3 and Data[1] < 317.41 + 43.19*3:   # APD90
		if Data[7] > -73.98 - 3.99*3 and Data[7] < -73.98 + 3.99*3:   # RMP
			if Data[8] > 30 and Data[8] < 219.44 +67.85*3:   # dvdtmax
			# if Data[8] > 219.44 -67.85*3 and Data[8] < 219.44 +67.85*3:   # dvdtmax
				# if Data[6] - Data[7] > 94.95-7.07*3 and Data[6] - Data[7] < 94.95+7.07*3:   # APA
				if Data[6] - Data[7] > 94.95-7.07*3 :   # APA
					if Data[3] > 138.09 - 45.14*3 and Data[3] < 138.09 + 45.14*3:   # APD50
						return True
	return False


# print check_results()



def run_Threshold(Para1, Para2, Para3, ID):
	f = open('workfile_lowwer', 'a')
	run(Para1, Para2, Para3, -20);
	call('mv Drug_CNZ_out.dat AP.BCL.1000.ID.'+index, shell=True);

	trials = 0;
	AMP = -5;
	Lastchange = 1

	run_Restart(Para1, Para2, Para3, AMP);
	go=check_results()
	lastgo = go

	compute = True
	Lastchange = 1
	while compute:
			if(lastgo == (not go) ):
				Lastchange=Lastchange/2.0

			if go:
				AMP =AMP+Lastchange -0.001
			else:
				AMP =AMP-Lastchange +0.001
			lastgo=go

			run_Restart(Para1, Para2, Para3, AMP);
			go=check_results()
			print  >>  f, "S1= ", BCL, "AMP= ", AMP, "Lastchange = ", Lastchange, "go = ", go, ' Cell_ID = ', ID
			if Lastchange<0.01 and go:
				compute = False
			if np.isnan(AMP):
				compute=False
			trials = trials +1
			if trials >100:
				AMP=np.nan
				Lastchange=np.nan
				break;
	return AMP,Lastchange






# run(1000, "SimDrug",100,100, 0.1,0.01)


# call(['make'])
# para_set = np.loadtxt('../pre_optim.dat')
para_set1 = np.loadtxt('para_log.dat.New')
# para_set2 = np.loadtxt('para_ina.dat')
# para_set3 = np.loadtxt('para_ina_v.dat')
# call(['rm', 'log.dat'])

kmfscale=1.0;
ISO=0.0

f = open("AP.log.dat.1Hz.kmf.%f.ISO.%f"%(kmfscale,ISO), "w+")
# f2 = open("para_log.dat", "w+")

BCL = 500
# run(BCL,"Normal", 0,0,0,0);
# Mode="SimDrug"

# sys.stdout = open('file', 'w')
for i in range(600):
	index = str(i)
	print ('processing ID ' +  index)
	run(para_set1[i],kmfscale, ISO, f);
	call('mv HAM_wrap_out.dat AP.BCL.1000.ID.'+index, shell=True);

	# call('mv HAM_wrap_out.dat AP.BCL.1000.ID.'+index, shell=True);
	call('mv APD_measure.dat APD_measure.dat.ID.'+index, shell=True);
	# call('mv betaAR_out.dat betaAR_out.dat.ID.'+index, shell=True);
	# call('mv CaM_cyto_out.dat CaM_cyto_out.dat.ID.'+index, shell=True);
	# call('mv CaM_dyad_out.dat CaM_dyad_out.dat.ID.'+index, shell=True);
	# call('mv CaMKII_out.dat CaMKII_out.dat.ID.'+index, shell=True);
	# call('mv para_out.dat para_out.dat.ID.'+index, shell=True);
	call('mv Ini.bin Ini.bin.'+index, shell=True);


# move all data to the new folder 
call('mkdir ../Initial_Conditions/kmf.%.1f.ISO.%.1f'%(kmfscale,ISO), shell=True)
call('mv Ini.bin.* ../Initial_Conditions/kmf.%.1f.ISO.%.1f'%(kmfscale,ISO), shell=True)