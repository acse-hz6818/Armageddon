# pylint: disable=invalid-name
"""
Extension 3:

Inversion to calculate the chern explosion
Try different values of Y and r in order to find the ones that give a result closer to the data 
-------------------------
returns:
a graphical output
df with errors and parameters choosen
"""
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import scipy.interpolate as si
import armageddon
#init class
earth = armageddon.Planet()
#read the csv of the values
protoCEAtable = pd.read_csv(r'C:\Users\gc2016\OneDrive - Imperial College London\ACSE\ACSE-4.2\acse-4-armageddon-hygiea\data\ChelyabinskEnergyAltitude.csv')
#initialise inital values
rc = range(10, 30, 2)
v0c = 19200
thetac = 18.3
rhoc = 3300
#put the csv values into arrays
CEAheight = np.array(protoCEAtable['Height (km)'])
CEAE = np.array(protoCEAtable['Energy Per Unit Length (kt Km^-1)'])


#error for Cheryabinsk
def rY_finder(r_min, r_max, Y_min, Y_max, nr, nY):
    """
    iterate to find the r and y from a range using the numerical solver
    """
    #convert the points to a continious function
    lp = si.interp1d(CEAheight, CEAE)
    #array for candidates
    rlist = np.linspace(r_min, r_max, nr)
    Ylist = np.linspace(Y_min, Y_max, nY)
    #containers for results
    maperror = []
    mapr = []
    mapY = []
    energies = []
    #loop nessesary: you have to loop over all combinations
    for i in range(nY):
        for j in range(nr):
            mapr.append(rlist[j])
            mapY.append(Ylist[i])
            #call numerical solver
            df = earth.solve_atmospheric_entry(rlist[j], v0c, 3300, Ylist[i], 18.3,
                                           1e5, dt=0.02, radians=False)
            df2 = earth.calculate_energy(df)
            #use only the curve for the error
            df_filtered = df2[(df2['altitude'] > 30) & (df2['altitude'] < 33) ]
            energies.append(df2.dedz)
            #rms error
            maperror.append(np.sqrt(np.sum((df_filtered.dedz)-lp(df_filtered.altitude))**2))
    errordf = pd.DataFrame({'Error': maperror, 'Radious': mapr, 'Strenght': mapY})
    return errordf, energies


def plot_model(list_e):
    """
    function to plot
    """
    plt.figure(figsize=(10,6))
    for i in list_e:
        plt.plot(i, np.linspace(100, 0, len(i)))

    plt.plot(CEAE, CEAheight, 'k', label='raw data')

    plt.xlabel('r gridpoints')
    plt.ylabel('Y gridpoints')
    plt.title('Squared Errors')
    plt.show()

# error, energies_list = (rY_finder(10, 12, 9e6, 1e7, 3, 3))
# print("error = ", error)
# plot_model(energies_list)


#print(CEAE)
#for initial conditions
df = earth.solve_atmospheric_entry(radius=10, velocity=21000, density=3000, strength=1e5, angle=45,
                                 init_altitude=100e3, dt=0.01, radians=False)
df2 = earth.calculate_energy(df)
print(df2)
plt.plot(df2.dedz, df2.altitude)
plt.show()
# print ("max energy", df2.dedz.max())
# print ("min energy", df2.dedz.min())

##################### Plot the initial values ########################################
# fig = plt.figure(figsize=(12, 10))
# CEA = fig.add_subplot(111)
# CEA.margins(0.1)

# lp = si.interp1d(CEAheight, CEAE)

# CEA.plot(CEAE, CEAheight, 'k', label='raw data')
# CEA.plot(lp(CEAheight), CEAheight, 'b*', label='approximation')

# CEA.set_xlabel('$dedz, kT/km$', fontsize=16)
# CEA.set_ylabel('$Height, z/m$', fontsize=16)
# CEA.grid(True)

# CEA.set_title('dE/dz-z Graph for Chelyabinsk and the interpolation to continuous', fontsize=16)
# CEA.legend(loc='upper left', fontsize=18)
