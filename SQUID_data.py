################################################

### Librairies  ###
from matplotlib import rc
rc('figure',figsize=(10,5)) # this size of window fits for the save fonction
rc('font',size=12)
rc('text',usetex=False)
from qubicpack import qubicpack as qp
import os,sys
import numpy as np
import pandas as pd
import matplotlib.mlab as mlab
import scipy.ndimage.filters as f
import glob
import string
import scipy.signal as scsig
from scipy import interpolate
import datetime as dt
import pickle
from qubicpack.pix2tes import tes2pix
import matplotlib.pyplot as plt
import scipy as sp
from datetime import datetime

################################################

### Fonctions ###
def save(path, ext='png', close=True, verbose=True):	
	# Extract the directory and filename from the given path
	directory = os.path.split(path)[0]
	filename = "%s.%s" % (os.path.split(path)[1], ext)
	if directory == '':
		directory = '.'

	# If the directory does not exist, create it
	if not os.path.exists(directory):
		os.makedirs(directory)

	# The final path to save to
	savepath = os.path.join(directory, filename)

	if verbose:
		print("Saving figure to '%s'..." % savepath),

	# Actually save the figure
	plt.savefig(savepath)
	
	# Close it
	if close:
		plt.close()

	if verbose:
		print("Done")
def plot_V_phi(Asic):
	# Mutual inductance
	Min=1./10.4E-6
	AsicNum=Asic
	### V Phi plot ###
	# !!!!!!!!!!!!!!!! Attention un peu long , a remanier !!!!!!!!!!!!!!
	for z in range (128):
	# the first 8 index are too low to induce curves, so we start at 7 
		for i in range(9):
			thedir = dirs[i+7]
	# reading the data from qubicpack
			b = qp()
			b.read_qubicstudio_dataset(thedir, asic=AsicNum)
			Rbias=b.Rbias
	# Amplitude peak to peak of the sinus
			amp=b.max_bias-b.min_bias
			# offset of the sinus signal
			offset=0.5*(b.max_bias+b.min_bias)
			Vbias=(offset+amp*b.bias_phase()/2.)
			# Ind : this allow to avoid the noncoherent points in raw data for the flux
			ind=np.argsort(Vbias) 
			# temp is a temporaly variable that will be used after for the filter
			if i == 0:
				Vsqoffset=np.mean(b.timeline(TES=z+1)*62.5/(70.*100)/(b.nsamples-b.n_masked()))
			temp = b.timeline(TES=z+1)[ind]*62.5/(70.*100)/(b.nsamples-b.n_masked())-Vsqoffset
		# savitzky golay filter
			filt = sp.signal.savgol_filter(temp, 51, 3) 


			### plot parameters ### 
			plt.plot(((Vbias[ind]/Rbias*Min)-(Vbias[ind]/Rbias*Min)[0])/2.,filt*-1)
			plt.grid()
			plt.xlabel('Flux (in quantum of flux)')
			plt.ylabel('Tension ($\mu$V)')
			plt.title('ASIC%i (input), SQUID number %i'%(AsicNum,(z+1)))

		save("./TEST/Analysis"+day+"/ASIC%i_Vphi/ASIC%i_%i"%(AsicNum,AsicNum,(z+1)),ext="png",close=True,verbose=True)
def plot_histo_sig(Asic):
	AsicNum=Asic
	histo=np.empty(16) # create tab for peak to peak val
	data=np.empty((128,16)) # create a tab for each squid to keep all ptp value for each 
	invdata=np.empty((16,128))
	for i in range (16) :
	    thedir = dirs[i]
	    h = qp()
	    h.read_qubicstudio_dataset(thedir, asic=AsicNum)
	    for z in range (128):
	        # Amplitude peak to peak of the sinus
	        amp=h.max_bias-h.min_bias
	        # offset of the sinus signal
	        offset=0.5*(h.max_bias+h.min_bias)
	        Vbias=(offset+amp*h.bias_phase()/2.)
	        # Ind : this allow to avoid the noncoherent points in raw data for the flux
	        ind=np.argsort(Vbias) 
	        # temp is a temporaly variable that will be used after for the filter
	        if i == 0 :
	            Vsqoffset=np.mean(h.timeline(TES=z+1)*62.5/(70.*100)/(h.nsamples-h.n_masked()))
	        temp = h.timeline(TES=z+1)[ind]*62.5/(70.*100)/(h.nsamples-h.n_masked())-Vsqoffset
	        # savitzky golay filter
	        filt = sp.signal.savgol_filter(temp, 51, 3) 
	        histo[i] = np.max(filt) - np.min(filt)
	        data[z,i]=histo[i]
	        invdata[i,z]=histo[i]

	plt.plot(data)
	plt.grid()
	plt.ylabel("PtP value")
	plt.xlabel("Number of SQUID")
	save("./TEST/Analysis"+day+"/Results/ASIC%i_data_plot1"%AsicNum,ext="png",close=True,verbose=True)

	plt.plot(invdata[:,:])
	plt.grid()
	plt.xlabel("Intensity (index of I)")
	plt.ylabel("SQUIDs")
	save("./TEST/Analysis"+day+"/Results/ASIC%i_data_plot2"%AsicNum,ext="png",close=True,verbose=True)

	# argmax take the position of the maxvalue for each squid
	plt.hist(np.argmax(data, axis=1), range=[0,16], bins=16)
	plt.grid()
	plt.ylabel("Number of SQUIDs")
	plt.xlabel("Index of current")
	plt.title("Histogram of the optimum current for the SQUID response for ASIC 1")
	save("./TEST/Analysis"+day+"/Results/ASIC%i_Histogram"%AsicNum,ext="png",close=True,verbose=True)

	plt.hist(data[:,9],range=[0,30], bins=30, alpha = 0.5, color= 'r' ,label="Isq = 25.5 µA")
	plt.hist(data[:,10],range=[0,30], bins=30, alpha = 0.5, color= 'b',label="Isq = 28 µA ")
	plt.hist(data[:,11],range=[0,30], bins=30, alpha = 0.5, color= 'g', label="Isq = 30.6 µA")
	plt.legend()
	plt.grid()
	plt.xlabel("Voltage ($\mu$V)")
	plt.xlim(0,25)
	plt.ylabel('Number of SQUID')
	plt.title("ASIC 1 histogram")
	save("./TEST/Analysis"+day+"/Results/ASIC%i_Histogram_multiindex"%AsicNum,ext="png",close=True,verbose=True)

	dat=np.empty(16)
	ind=np.empty(16)


	file = open("mean.txt", "w")
	for z in range (16):
	    file.write("for index %i"%z) 
	    a=np.shape((np.where(data[:,z]>= 10)))
	    prct=a[1]/128. *100 
	    file.write("%f working squid" %prct)
	    file.write("median = %f" %np.median(data[:,z]))
	    file.write("\n")
	    dat[z]= prct
	    ind[z]= z
	file.close()


	plt.plot(ind,dat)
	plt.grid()
	plt.ylabel("Distribution of working SQUID >10µV ($\mu$V)")
	plt.xlabel('Index')
	plt.title("Working SQUID by index")
	save("./TEST/Analysis"+day+"/Results/ASIC%i_percentage"%AsicNum,ext="png",close=True,verbose=True)
def plot_noise() :
	####################################################


	### variable declaration ###
	fmi=input("Enter frequence min value: ")
	fma=input("Enter frequence max value: ")
	fmin=int(fmi)
	fmax=int(fma)
	file = open("Noise.txt", "w")
	current_noise_test=np.zeros((2,16,128))#ASIC/index/TES
	Ibias=[]
	Noise_plot=[]
	#########################################################


	### loop for filling noise ###
	#selection courant et ASIC
	for i in range (2):
	    AsicNum = i+1
	    for j in range (16):
	        thedir = dirs[j]
	        print(thedir)
	        b = qp()
	        b.read_qubicstudio_dataset(thedir, asic=AsicNum)

	        #passe a travers tous les TES
	        Rbias=b.Rbias
	        Min=1./10.4E-6
	        for z in range (128):
	            amp=b.max_bias-b.min_bias
	            offset=0.5*(b.max_bias+b.min_bias)
	            Vbias=(offset+amp*b.bias_phase()/2.)
	            ind=np.argsort(Vbias) 

	            Vsqoffset=np.mean(b.timeline(TES=z+1)*62.5/(70.*100)/(b.nsamples-b.n_masked()))
	            temp = b.timeline(TES=z+1)[ind]*62.5/(70.*100)/(b.nsamples-b.n_masked())-Vsqoffset
	            filt = sp.signal.savgol_filter(temp, 51, 3) 
	            phi0=((Vbias[ind]/Rbias*Min)-(Vbias[ind]/Rbias*Min)[0])/2.

	            #Bruit de lecture
	            timeline_brut=b.timeline(TES=z+1)
	            Timeline_volt=timeline_brut*62.5/(70.*100)/(b.nsamples-b.n_masked())
	            fs = 156.25
	            f, t, Sxx = scsig.spectrogram(Timeline_volt, fs,nperseg=2**10,window='hann')
	            indf=(f>fmin) & (f<fmax)
	            a=np.median(Sxx[indf,0])
	            val=np.sqrt(a)
	            #Bruit du SQUID
	            DeltaX=max(np.gradient(phi0))
	            slope_plot=(np.gradient(filt)/DeltaX)
	            #slope=(max(slope_plot))
	            slope_filt = sp.signal.savgol_filter(slope_plot, 51, 3) 
	            slope=(max(slope_filt))
	            #plt.plot(slope_plot)
	            phi0_Hz=(val*1e-6)/(slope*1e-6)
	            A_Hz=phi0_Hz*0.2e-6
	            noise=A_Hz*1e12
	            Ibias.append(j)
	            Noise_plot.append(noise)
	            current_noise_test[i,j,z]= noise
	            file.write("ASIC n° %i ... Index n° %i ... TES n° %i ...  SQUID Noise = %d pA(Hz)^-1/2\n" %(i+1,j,z+1,noise))
	file.close()

	#TES=127
	for TES in range (128):
	    plt.plot(current_noise_test[0,:,TES],label="ASIC 1",color="blue")
	    plt.plot(current_noise_test[1,:,TES],label='ASIC 2',color="red")
	    plt.legend()
	    plt.title("Noise regarding of the Index number, TES : %i" %(TES+1))
	    plt.ylabel("Noise pA/(Hz)^1/2")
	    plt.grid()
	    save("./TEST/Analysis"+day+"/Noise"+day+"/Graph/Noise_Index_TES%i"%(TES+1),ext="png",close=True,verbose=True)

	#ind=4
	for ind in range (16):
	    plt.hist(current_noise_test[0,ind,:],bins=25,alpha=0.5,label="ASIC 1",color="blue")
	    plt.hist(current_noise_test[1,ind,:],bins=25,alpha=0.5,label="ASIC 2",color="red")
	    plt.grid()
	    plt.legend()
	    plt.title("Number of TES regarding the noise, index : %i" %(ind))
	    plt.xlabel("Noise pA/(Hz)^1/2")
	    plt.ylabel("Number of TESs")
	    save("./TEST/Analysis"+day+"/Noise"+day+"/Hist/TES_Noise_ind%i"%(ind),ext="png",close=True,verbose=True)
	    
	###########################################

	###3D plot ###


	current_noise=np.zeros((2,16,128))#ASIC/TES/Index
	for i in range (2):
	    current_noise[i,:,:]=i+1
	    for k in range (16):
	        current_noise[:,k,:]=k
	        for z in range (128):
	            current_noise[:,:,z]=z+1
	from mpl_toolkits import mplot3d
	fig = plt.figure()
	ax = plt.axes(projection='3d')


	# Data for a three-dimensional line
	zline = Noise_plot
	xline = current_noise
	yline = Ibias
	ax.scatter3D(yline, xline, zline,c=zline, cmap='jet')

	#ax.set_ylim(0,6)
	#ax.set_zlim(0,300)
	#ax.set_xlim(74.5,75.5)
	ax.set_ylabel('TES')
	ax.set_xlabel('Index')
	plt.xticks([0,2,4,6,8,10,12,14])
	ax.set_zlabel('Noise')
	ax.set_title('Noise distribution for all TESs')
	save("./TEST/Analysis"+day+"Noise"+day+"/3D_data_plot",ext="png",close=True,verbose=True)


	###########################################################################

	### Noise Histogram ###
	z=0
	#for ind in range (9):
	plt.hist(current_noise_test[z,9,:],bins=50,alpha=0.5,label="25.5 µA",color="red")
	plt.hist(current_noise_test[z,10,:],bins=50,alpha=0.5,label="28 µA",color="blue")
	plt.hist(current_noise_test[z,11,:],bins=50,alpha=0.5,label="30.6 µA",color="green")

	plt.grid()
	plt.legend()
	plt.title("ASIC %i" %(z+1))
	plt.xlabel("Noise pA/(Hz)^1/2")
	plt.xlim(0,500)
	plt.ylabel("Number of TESs")
	save("./TEST/Analysis"+day+"Noise"+day+"/Noise_ASIC1",ext="png",close=True,verbose=True)

	z=1
	#for ind in range (9):
	plt.hist(current_noise_test[z,10,:],bins=50,alpha=0.5,label="28 µA",color="red")
	plt.hist(current_noise_test[z,11,:],bins=50,alpha=0.5,label="30.6 µA",color="blue")
	plt.hist(current_noise_test[z,12,:],bins=50,alpha=0.5,label="33.2 µA",color="green")

	plt.grid()
	plt.legend()
	plt.title("ASIC %i" %(z+1))
	plt.xlabel("Noise pA/(Hz)^1/2")
	plt.xlim(0,500)
	plt.ylabel("Number of TESs")
	save("./TEST/Analysis"+day+"Noise"+day+"/Noise_ASIC2",ext="png",close=True,verbose=True)

##################################################

### Reading the datas  ###
day=input("Enter the directory name for the data : ")
data_dir = '/home/guillaume/Documents/Data/'+day+'/'+day+'/'
#dirs = np.sort(glob.glob(data_dir+'*test_sw*'))
dirs = np.sort(glob.glob(data_dir+'*bias*'))

##################################################

### Main ###
startTime = datetime.now()
for i in range(2):
	plot_V_phi(i+1)
	plot_histo_sig(i+1)

plot_noise()

print("##############")
print(datetime.now() - startTime)
print("##############")

##################################################









