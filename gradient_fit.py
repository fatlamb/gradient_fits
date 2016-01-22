import numpy as np
import matplotlib.pyplot as plt
import ROOT as rt
import scipy.optimize as opt
import math

filenames = [line.rstrip('\n') for line in open('fnames.txt')]
plotnames=["Bottom","Top"]

directory="/home/pinkpig/physics/labdata/lightguide_caen/processed/"

zdist=7.4 #cm
step=2.5 #cm

print len(filenames)

means=(np.zeros((len(filenames),2,17)))
sigmas=(np.zeros((len(filenames),2,17)))
nentries=(np.zeros((len(filenames),2,17)))
xe=(np.ones((len(filenames),2,17))*0.3) #3mm err estimate
ye=(np.zeros((len(filenames),2,17)))



for fname in enumerate(filenames):
	print fname	

	for dist in range(0,17):
		for direction in range(0,2):
			dirfile=directory+fname[1]+"_"+str(dist)+"_"+str(direction)+"/"
			dirfile+=fname[1]+"_"+str(dist)+"_"+str(direction)+".root"

			f=rt.TFile(dirfile)	
			ot=f.Get("OutTree")

			ot.GetEntry(0)
			means[fname[0]][direction][dist]=ot.Mean
			sigmas[fname[0]][direction][dist]=ot.Sigma
			nentries[fname[0]][direction][dist]=ot.Nentries
			#ye[fname[0]][direction][dist]=ot.Sigma/math.sqrt(float(ot.Nentries))
			ye[fname[0]][direction][dist]=ot.Mean*0.0052



class fitstuff(object):
	def __init__(self,x,y,xerr,yerr,bar):
		self.x=x
		self.y=y
		self.xerr=xerr
		self.yerr=yerr
		#self.barlen=[50.6,50.6,50.55,50.5] #(cm) #sbars
		#self.meniscus=[2.9,2.95,2.9,0.7] #(cm)
		#self.barlen=[50.8,50.25,51.1] #(cm) full dirtybar
		#self.meniscus=[0,0,2.2] #(cm)
		self.barlen=[51.0,50.85,50.50] #(cm)
		self.meniscus=[2.7,3.1,3.2] #(cm)
		self.bar=bar

	def setbar(self,bar):
		self.bar=bar
	

	def fitfuncforward(self,x,param):
		return param[0]*math.exp(-1.0*param[1]*x)*(1.0+(x-self.meniscus[self.bar])*param[2])
	def fitfuncbackward(self,x,param):
		return param[3]*param[0]*math.exp(-1.0*param[1]*x)*(1.0+(self.barlen[self.bar]-self.meniscus[self.bar]-x)*param[2])
	

	def chi2(self,param):
		chi2=0
		for direction in range(0,2):
			for motor in range(1,2):
				for j in range(0,len(self.x)):
					if direction==0:
						chi2_y=(self.y[direction,motor,j]-self.fitfuncforward(self.x[j],param))**2/self.yerr[direction,motor,j]**2
					else:
						chi2_y=(self.y[direction,motor,j]-self.fitfuncbackward(self.x[j],param))**2/self.yerr[direction,motor,j]**2
			
					chi2+=chi2_y
		return chi2	


edgecut=[2,-2]

x=np.ones(17)*zdist+np.linspace(0,16,17)*step
x=x[edgecut[0]:edgecut[1]]
param=np.zeros(4)

bounds=[]

bounds.append((10.0**4,10.0**5))
bounds.append((1.0/1000,1.0))
bounds.append((1.0/10000,2.0))
bounds.append((0.1,2.0))




for bar in range(0,len(filenames)/2):
	param[0]=means[2*bar,1,0]*0.5 #init guess
	#param[0]=4e+04 #init guess
	param[1]=1.0/200
	param[2]=0.005
	param[3]=1.0

	fit=fitstuff(x,means[2*bar:2*bar+2,:,edgecut[0]:edgecut[1]],xe[2*bar:2*bar+2,:,edgecut[0]:edgecut[1]],ye[2*bar:2*bar+2,:,edgecut[0]:edgecut[1]],bar)
	res=opt.minimize(fit.chi2,param,bounds=bounds,method='L-BFGS-B',tol=1e-10)


	print "CONST: ",res.x[0]
	print "SLOPE: ",(1/res.x[1])
	print "LINEAR GRADIENT: ",res.x[2]
	print "COATING FACTOR: ",res.x[3]

	yfuncf=np.zeros(len(x))
	yfuncb=np.zeros(len(x))
	for j in range(0,len(x)):
		yfuncf[j]=fit.fitfuncforward(x[j],np.asarray(res.x))
		yfuncb[j]=fit.fitfuncbackward(x[j],np.asarray(res.x))


	fig=plt.figure(1)
	ax = plt.subplot(211)
	#ax.errorbar(x,means[0,0,3:-3],xerr=xe[0,0,3:-3],yerr=ye[0,0,3:-3],fmt='ro',label='P(mu->mu): Diagonalized')
	#ax.errorbar(x,means[2*bar,0,edgecut[0]:edgecut[1]],xerr=xe[2*bar,0,edgecut[0]:edgecut[1]],yerr=ye[2*bar,0,edgecut[0]:edgecut[1]],fmt='ro',label='Outgoing')
	ax.errorbar(x,means[2*bar,1,edgecut[0]:edgecut[1]],xerr=xe[2*bar,1,edgecut[0]:edgecut[1]],yerr=ye[2*bar,1,edgecut[0]:edgecut[1]],fmt='ro',label='Data')
	ax.plot(x,yfuncf,'b-',label='Fit')
	#ax.plot(xodone,yfunc,'go')
	ax.set_title("KBar "+str(bar+1)+" with Linear Gradient")
	ax.set_ylabel("Charge (ADC*TDC)")
	legend = plt.legend(loc='lower right', shadow=False, fontsize='large')
	ax.set_ylim(0,30000)	
	ax2=plt.subplot(212)
	#ax2.errorbar(x,means[2*bar+1,0,edgecut[0]:edgecut[1]],xerr=xe[2*bar+1,0,edgecut[0]:edgecut[1]],yerr=ye[2*bar+1,0,edgecut[0]:edgecut[1]],fmt='ro',label='')
	ax2.errorbar(x,means[2*bar+1,1,edgecut[0]:edgecut[1]],xerr=xe[2*bar+1,1,edgecut[0]:edgecut[1]],yerr=ye[2*bar+1,1,edgecut[0]:edgecut[1]],fmt='ro',label='')
	ax2.plot(x,yfuncb,'b-')
	ax2.set_xlabel("Illumination Distance from PMT (cm)")
	ax2.set_ylabel("Charge (ADC*TDC)")
	
	ax2.set_ylim(0,30000)	
	

	ax2.text(0.65, -0.35, 'attenuation='+'{:04.1f}'.format(1/res.x[1])+' cm', style='italic',bbox={'facecolor':'cyan', 'alpha':0.5, 'pad':10},transform=ax.transAxes)
	ax2.text(0.65, -0.5, 'brightness='+'{:03.2f}'.format(100*res.x[2])+' %/cm', style='italic',bbox={'facecolor':'cyan', 'alpha':0.5, 'pad':10},transform=ax.transAxes)
	ax2.text(0.65, -0.65, 'coating factor='+'{:03.2f}'.format(res.x[3]), style='italic',bbox={'facecolor':'cyan', 'alpha':0.5, 'pad':10},transform=ax.transAxes)

	
	plt.show()
	
