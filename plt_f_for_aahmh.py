#Ag/Al2O3/hbn/mose2/hbn
import numpy as np
import matplotlib.pyplot as plt
from refrac_index_Ag import get_n_ag,get_k_ag
from refrac_index_al2o3_b import get_n_al2o3_b
#include the thickness of hbn into the plot function

ev=1.602176634*10**(-19)
h=6.626*10**(-34)
c=3*10**(8)
def lam2ev(lam):#lambda:nm, E:ev
    lam=lam*10**(-9)
    E=h*c/lam
    return E/ev

def delta(ni,di,lam):#lam dependent
    return 2*np.pi*ni*di/lam
def tmatrix(delta,n):#lam dependent
    tm=np.zeros((2,2),dtype=np.complex_)#dtype
    tm[0,0]=np.cos(delta)
    tm[0,1]=np.sin(delta)*1j/n
    tm[1,0]=np.sin(delta)*1j*n
    tm[1,1]=np.cos(delta)
    return tm
#print(tmatrix(np.pi/2,3+2*1j))#test
n1=1.85#HBN
n3=1.85#HBN
def chi(lam,E0,f,gamma):
    E=lam2ev(lam)
    return f/(E0**2-E**2-1j*E*gamma)
def epsilon(E,E0,f,gamma):
    chi=f/(E0**2-E**2-1j*E*gamma)#E0:resonance energy?
    return 1+chi
#print(epsilon(2.0,2.2,3.5,2.4))#test

def n2(epsilon):#mose2
    epabs=abs(epsilon)
    epreal=epsilon.real
    n=(((epabs+epreal))/2)**(1/2)
    k=(((epabs-epreal))/2)**(1/2)
    return n+1j*k 

def n4(lam):#al2o3
    return get_n_al2o3_b(lam)
def n5(lam):#ag
    return get_n_ag(lam)+1j*get_k_ag(lam)
def vec1(lam):
    v1=np.zeros((2,1),dtype=np.complex_)#dtype
    v1[0,0]=1
    v1[1,0]=n5(lam)
    return v1
#print(vec1(550))#test




# d1=10:up hbn
# d2=0.65:mose2
# d3=10:down hbn
# d4=90:al2o3
def plotfunc_aahmh(wavelengtharr,E0,f,gamma,d1,d2,d3,d4):###
    
    def tmatrix4(lam):
        return tmatrix(delta(n4(lam),d4,lam),n4(lam))
    #print(tmatrix4(550))#test
    
    def tmatrix3(lam):
        return tmatrix(delta(n3,d3,lam),n3)
    
    def tmatrix2(lam,E0,f,gamma):#########################
        E=lam2ev(lam)
        return tmatrix(delta(n2(lam),d2,lam),n2(epsilon(E,E0,f,gamma)))
    def tmatrix2s(lam):#########################
        return tmatrix(delta(1,d2,lam),1)
    def tmatrix1(lam):
        return tmatrix(delta(n1,d1,lam),n1)
    def contrast(x,xs):
        return (x-xs)/xs


    Barr=np.zeros(101,dtype=np.complex_)#whole
    Carr=np.zeros(101,dtype=np.complex_)
    Barrs=np.zeros(101,dtype=np.complex_)#substrate
    Carrs=np.zeros(101,dtype=np.complex_)
    chiarr_real=np.zeros(101,dtype=np.complex_)
    chiarr_imag=np.zeros(101,dtype=np.complex_)
    for i in range(np.size(wavelengtharr)):
        lam=wavelengtharr[i]
        
        vecBC=np.matmul(tmatrix1(lam),vec1(lam))
        vecBC=np.matmul(tmatrix2(lam,E0,f,gamma),vecBC)
        vecBC=np.matmul(tmatrix3(lam),vecBC)
        vecBC=np.matmul(tmatrix4(lam),vecBC)
        
        Barr[i]=(vecBC[0,0])
        Carr[i]=(vecBC[1,0])
        

        vecBCs=np.matmul(tmatrix1(lam),vec1(lam))
        vecBCs=np.matmul(tmatrix2s(lam),vecBCs)
        vecBCs=np.matmul(tmatrix3(lam),vecBCs)
        vecBCs=np.matmul(tmatrix4(lam),vecBCs)#substrate
        
        Barrs[i]=(vecBCs[0,0])
        Carrs[i]=(vecBCs[1,0])
        chiarr_real[i]=(chi(lam,E0,f,gamma).real)
        chiarr_imag[i]=(chi(lam,E0,f,gamma).imag)
    
    #print(vecBC)#test
    Yarr=np.zeros_like(Barr,dtype=np.complex_)
    Yarrs=np.zeros_like(Barrs,dtype=np.complex_)
    for i in range(np.size(Yarr)):
        Yarr[i]=Carr[i]/Barr[i]
        Yarrs[i]=Carrs[i]/Barrs[i]

    rarr=np.zeros_like(Yarr,dtype=np.complex_)
    rarrs=np.zeros_like(Yarrs,dtype=np.complex_)
    Rarr=np.zeros_like(Yarr,dtype=np.complex_)#reflected intensity
    Rarrs=np.zeros_like(Yarrs,dtype=np.complex_)#reflected intensity
    contrastarr=np.zeros_like(Yarrs,dtype=np.complex_)
    for i in range(np.size(Yarr)):
        rarr[i]=(1-Yarr[i])/(1+Yarr[i])
        
        Rarr[i]=abs(rarr[i])**2
        
        rarrs[i]=(1-Yarrs[i])/(1+Yarrs[i])
        Rarrs[i]=abs(rarrs[i])**2
        contrastarr[i]=contrast(Rarr[i],Rarrs[i])
#imag part is 0
#turn it into a real array
    x=[]
    y=[]
    ctr=[]
    z=[]
    a=[]
    for i in range(np.size(Rarr)):
        x.append(Rarr[i].real)
        y.append(Rarrs[i].real)
        ctr.append(contrastarr[i].real)
        z.append(chiarr_real[i].real)
        a.append(chiarr_imag[i].real)
#x:Rarr,y:Rarrs,ctr:contrastarr,z:chiarr_real,a:chiarr_imag
    return [wavelengtharr,x,y,ctr,z,a]

'''
wavelengtharr=np.linspace(700,800,101)
E0=1.67#ev
f=0.01#or A#####
gamma=1*10**(-3)#ev

wavelengtharr=plotfunc(wavelengtharr,E0,f,gamma)[0]
Rarr=plotfunc(wavelengtharr,E0,f,gamma)[1]
Rarrs=plotfunc(wavelengtharr,E0,f,gamma)[2]
contrastarr=plotfunc(wavelengtharr,E0,f,gamma)[3]
chiarr_real=plotfunc(wavelengtharr,E0,f,gamma)[4]
chiarr_imag=plotfunc(wavelengtharr,E0,f,gamma)[5]
fig,axis=plt.subplots(5)
axis[0].plot(wavelengtharr,Rarr,label="whole stack")
axis[0].legend()
axis[1].plot(wavelengtharr,Rarrs,label="substrate")
axis[1].legend()
axis[2].plot(wavelengtharr,contrastarr,label="contrast")
axis[2].legend()

axis[3].plot(wavelengtharr,chiarr_real,label="chi(real)")
axis[3].legend()

axis[4].plot(wavelengtharr,chiarr_imag,label="chi(imag)")
axis[4].legend()
axis[4].set_xlabel("wavelength(nm)")

plt.show()
'''