import os
import numpy as np
import math
from PIL import Image
from skimage.morphology import disk
from skimage.morphology import square
from skimage.morphology import erosion
from skimage.morphology import closing
from sklearn.cluster import KMeans
from skimage.filters.rank import median
import cv2
from skimage.filters import threshold_otsu
from skimage.morphology import opening
import matplotlib.pyplot as plt

class comp():
	def __init__(self,propMat,AtomicWeight,method="scheil"):
		'''
		Is used to convert mass fraction to volume fraction in eutectic alloys,
		using either Schiel equation or lever rule

		propMat = [max solubility of solute in solvent, eutectic composition, partition coefficient]
		AtomicWeight = [Atomic weight of solute, atomic weight of solvent]
		'''
		self.Cs_w = propMat[0]
		self.Ce_w = propMat[1]
		self.k = propMat[2]
		self.AWslt = AtomicWeight[0]
		self.AWslv = AtomicWeight[1]
		self.method=method

	def Weight2Atomic(self,C):
		return	(1/(1+(((100-C)/float(C))*(self.AWslt/self.AWslv))))

	def mf2vf(self,C):
		Cs = self.Weight2Atomic(self.Cs_w)
		Ce = self.Weight2Atomic(self.Ce_w)
		C = self.Weight2Atomic(C)
		if self.method == "scheil":
			return (self.k*C/float(Cs))**(1/(1-float(self.k)))
		elif self.method == "lever":
			return ((C - Cs)/(Ce-float(Cs)))

class VolFracTomog():
	'''
	is used to calculate the area fraction of the secondary phase in a micrograph.

	User can implement K-means method, use their own threshold array or use Ostu's method to segment the image.

	The code then can be used to iterate over a number of images,
	this is specially useful if the user is inteding to determine the volume fraction of a phase in XMT images.
	'''
	def __init__(self,img=[],thresh=[],nClass=0):
		self.img = img
		self.thresh=thresh
		self.nClass=nClass

	def im2vec(img):
		try:
			h,w=img.shape
			return np.reshape(img,(h*w,1))
		except ValueError:
			h,w,z=img.shape
			return np.reshape(img,(h*w,z))

	def findAreaFrac(self):
		self.img = self.img.convert("L")
		self.img = np.asarray(self.img)
		if self.thresh==[] and self.nClass>0:
			''' Implement K-means '''
			phaseDic={}
			kmeans = KMeans(init='k-means++', n_clusters=nClass)
			kmeans.fit(self.im2vec(img))
			seg = kmeans.labels_
			sumMat = [np.sum(seg==0),np.sum(seg==1),np.sum(seg==2)]
			tmp = [np.sum(seg==0),np.sum(seg==1),np.sum(seg==2)]
			labelDic={}
			for phase in ['eut','pri','pore']:
				phaseDic[phase] = min(sumMat)
				labelDic[phase] = tmp.index(phaseDic[phase])
				sumMat.pop(sumMat.index(phaseDic[phase]))
			eut = np.reshape(seg==labelDic['eut'],img.shape)
			eut = median(eut, square(2))
			eut = erosion(eut, square(2))
			eut = median(eut, disk(2))
			eutRatio = np.sum(eut==np.max(eut))/float(phaseDic['eut']+phaseDic['pri'])
			return eutRatio
		elif len(self.thresh)==2:
			''' Segment using given threshold '''
			eut = (self.img>=self.thresh[1])*(self.img<256)*1.0
			eut = erosion(eut, disk(1))
			try:
				np.seterr(divide='ignore', invalid='ignore')
				eutRatio = np.sum(eut>0)/float(np.sum(self.img>=self.thresh[0]))
				if math.isnan(eutRatio):
					eutRatio=0
			except:
				eutRatio=0
			return eutRatio
		elif self.thresh==[] and self.nClass==0:
			''' Segment using Otsu method '''
			inThresh=30
			clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
			equ = clahe.apply(self.img)
			self.thresh=[inThresh,threshold_otsu(equ[equ>=30])]
			eut = (equ>=self.thresh[1])*(equ<256)*1.0

			eut=median(eut,square(12))
			#eut=closing(eut, disk(2))
			eut=erosion(eut, square(3))
			eut=median(eut,square(5))
			eut=closing(eut, disk(5))
			#eut=erosion(eut, square(3))
			#eut=erosion(eut, square(3))
			eut=median(eut,square(5))

			try:
				np.seterr(divide='ignore', invalid='ignore')
				eutRatio = np.sum(eut>0)/float(np.sum(self.img>=self.thresh[0]))
				if math.isnan(eutRatio):
					eutRatio=0
			except:
				eutRatio=0
			self.thresh=[]
			return eutRatio
		else:
			raise ValueError("Initizalization Error.")

	def findVolumeFrac(self,n, path, fname):
		folder = path + fname + "/rec_8bit"
		os.chdir(folder)
		areafrac=[]
		for i in range(n):
			tmp = str(i+1)
			if len(tmp)==1:
				tmp = "00"+tmp
			elif len(tmp)==2:
				tmp = "0"+tmp
			name = fname +tmp+".rec.8bit.tif"
			self.img = Image.open(name)
			areafrac.append(self.findAreaFrac())
			print(str(round(100*i/float(n),2))+" %")
		volfrac = sum(areafrac)/float(n)
		os.chdir(path)
		return volfrac

def plotData(comps,fnames=[],n=924,props=([5.6,33.6,0.17],[63.54,26.98])):
	alcu = comp(props[0],props[1])
	vol = VolFracTomog()
	act=np.array([])
	for i in range(len(comps)):
		act=np.append(act,alcu.mf2vf(comps[i]))
	volfrac = np.array([])
	for fname in fnames:
		print(fname)
		volfrac=np.append(volfrac,vol.findVolumeFrac(n,fname))
	e2=volfrac[[0,4,8]]
	e3=volfrac[[1,5,9]]
	e4=volfrac[[2,6,10]]
	e5=volfrac[[3,7,11]]
	e2ave = np.sum(e2)/len(e2)
	e3ave = np.sum(e3)/len(e3)
	e4ave = np.sum(e4)/len(e4)
	e5ave = np.sum(e5)/len(e5)
	e2std = np.std(e2)
	e3std = np.std(e3)
	e4std = np.std(e4)
	e5std = np.std(e5)

	ave = np.array([e2ave,e3ave,e4ave,e5ave])
	std = np.array([e2std,e3std,e4std,e5std])
	plt.plot(act,'ro')
	plt.errorbar(range(4),ave, yerr= std,fmt='b>')
	names=["Al-6.3%Cu", "Al-7.2%Cu", "Al-8.5%Cu", "Al-13.5%Cu"]
	plt.xticks(range(4),names)
	plt.xlabel("Alloy")
	plt.ylabel("Volume Fraction")
	plt.legend(["Actual Value Scheil","Measured Value"])
	plt.title("Volume Fraction Measurement")
	plt.xlim(-1, 4)
	plt.ylim(0.05,0.5)
	plt.show()

if __name__ == '__main__':

	fnames = ["E2-1_","E3-1_","E4-1_","ES-1_","E2-2_","E3-2_","E4-2_","ES-2_","E2-3_","E3-3_","E4-3_","ES-3_"]
	comps = [6.3,7.2,8.5,13.5]
	plotData(comps,fnames,924)
