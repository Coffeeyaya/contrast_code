from plt_f_for_aahmh import plotfunc_aahmh
import numpy as np
import matplotlib.pyplot as plt
wavelengtharr=np.linspace(700,800,101)
E0=1.67#ev
f=0.01
gamma=1*10**(-3)#ev
d1=10#nm
d2=0.65#nm
d3=10#nm
d4=300#nm

Rarr=plotfunc_aahmh(wavelengtharr,E0,f,gamma,d1,d2,d3,d4)[1]
Rarrs=plotfunc_aahmh(wavelengtharr,E0,f,gamma,d1,d2,d3,d4)[2]
contrastarr=plotfunc_aahmh(wavelengtharr,E0,f,gamma,d1,d2,d3,d4)[3]
chiarr_real=plotfunc_aahmh(wavelengtharr,E0,f,gamma,d1,d2,d3,d4)[4]
chiarr_imag=plotfunc_aahmh(wavelengtharr,E0,f,gamma,d1,d2,d3,d4)[5]

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
