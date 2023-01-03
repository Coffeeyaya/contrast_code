import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
data = pd.read_csv('RefractiveIndexINFO_gold.csv')

#extract the columns from the data
df_wavelength=pd.DataFrame(data,columns=['Wavelength, µm'])
df_n=pd.DataFrame(data,columns=['n'])
df_k=pd.DataFrame(data,columns=['k'])

#change the datatype
df_list_wavelength=df_wavelength.to_numpy()#its shape is (46,1)
df_list_wavelength=df_list_wavelength[:,0]#extract the 1st column
df_list_wavelength=(df_list_wavelength)*1000#change µm to nm
df_list_n=df_n.to_numpy()
df_list_n=df_list_n[:,0]
df_list_k=df_k.to_numpy()
df_list_k=df_list_k[:,0]

#polynomial fit
polyn = np.polyfit(df_list_wavelength, df_list_n, deg=10)
polyk = np.polyfit(df_list_wavelength, df_list_k, deg=10)

def get_n_gold(lam):
    return np.polyval(polyn,lam)
def get_k_gold(lam):
    return np.polyval(polyk,lam)


#plot
fig, axis = plt.subplots(2)
axis[0].plot(df_list_wavelength,df_list_n, label='n')
axis[0].plot(df_list_wavelength,np.polyval(polyn,df_list_wavelength), label='fit n')
axis[0].legend()

axis[1].plot(df_list_wavelength,df_list_k, label='k')
axis[1].plot(df_list_wavelength,np.polyval(polyk,df_list_wavelength), label='fit k')
axis[1].legend()
plt.show()






