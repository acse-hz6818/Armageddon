import numpy as np
import pandas as pd
#import math
#import matplotlib.pyplot as plt

class Planet():
    """
    The class called Planet is initialised with constants appropriate
    for the given target planet, including the atmospheric density profile
    and other constants
    """

    def __init__(self, atmos_func='exponential', atmos_filename=None,
                 Cd=1., Ch=0.1, Q=1e7, Cl=1e-3, alpha=0.3, Rp=6371e3,
                 g=9.81, H=8000., rho0=1.2):
        """
        Set up the initial parameters and constants for the target planet
        Parameters
        ----------
        atmos_func : string, optional
            Function which computes atmospheric density, rho, at altitude, z.
            Default is the exponential function ``rho = rho0 exp(-z/H)``.
            Options are ``exponential``, ``tabular``, ``constant`` and ``mars``
        atmos_filename : string, optional
            If ``atmos_func`` = ``'tabular'``, then set the filename of the table
            to be read in here.
        Cd : float, optional
            The drag coefficient
        Ch : float, optional
            The heat transfer coefficient
        Q : float, optional
            The heat of ablation (J/kg)
        Cl : float, optional
            Lift coefficient
        alpha : float, optional
            Dispersion coefficient
        Rp : float, optional
            Planet radius (m)
        rho0 : float, optional
            Air density at zero altitude (kg/m^3)
        g : float, optional
            Surface gravity (m/s^2)
        H : float, optional
            Atmospheric scale height (m)
        Returns
        -------
        None
        """

        # Input constants
        self.Cd = Cd
        self.Ch = Ch
        self.Q = Q
        self.Cl = Cl
        self.alpha = alpha
        self.Rp = Rp
        self.g = g
        self.H = H
        self.rho0 = rho0

        if atmos_func == 'exponential':
            self.rhoa = lambda z: self.rho0 * np.exp(-z/self.H)

        elif atmos_func == 'tabular':
            # atmos_filename="../data/AltitudeDensityTable.csv"
            assert atmos_filename is not None
            data = pd.read_csv(atmos_filename, skiprows=6, delimiter=' ', names = ['Altitude', 'Density', 'Height'])
            
            def limitz(z):
                roundedz=int(z/10)
                if z < 0:
                    roundedz=0
                elif z> 86000:
                    roundedz=8600
                return data.Density[roundedz]*np.exp((data.Altitude[roundedz]-z)/(data.Height[roundedz]))
            self.rhoa = lambda z: limitz(z)
            #self.rhoa = lambda z: data.Density[int(z/10)]*np.exp((data.Altitude[int(z/10)]-z)/(data.Height[int(z/10)]))

        elif atmos_func == 'mars':
            def T(altitudez):
                if altitudez < 7000.:
                    return 242.1-0.000998*altitudez
                else:
                    return 249.7-0.00222*altitudez
            self.rhoa = lambda z: 0.699*np.exp(-0.00009*z)/(0.1921*T(z)) #p in kPa
            
        elif atmos_func == 'constant':
            self.rhoa = lambda x: rho0 #rho0
        else:
            print('''Choose from 'exponential', 'tabular', 'mars' or 'constant'.''')
            raise NotImplementedError


    #analytical v
    def r2A(self, radius):
        '''
        calculate area from given radius
        '''
        return np.pi*radius**2

    def rrho2m(self, radius, density):
        '''
        calculate mass from given radius & density
        '''
        return (4./3.)*np.pi*radius**3*density
    
    def vz(self, radius, velocity, density, strength, angle, z, degree=True):
        '''
        ANALYTICAL SOLUTION!
        
        Calculating velocity from altitude.
        
        ----------input-----------
        
        v0, initial speed, m/s
        
        A, cross-sectional area, m**2
        
        m, mass in kg
        
        angle, angle of injection, in radians
        
        z, altitude in metres, usually an array
        
        Cd, drag coefficient, preset as 1
        
        rho0, density of air, preset as 1.2kg/m**3
        
        H, height of atmosphere, preset as 8000m
        
        radian, boolean datum that tells whether the angle is given in radians; preset as False 
        
        ----------output-----------
        
        speed in m/s, usually an array
        '''
        if degree:
            angle_rad = angle*np.pi/180.
        else:
            angle_rad = angle
    #print(angle)
        return velocity*np.exp(-self.H*self.Cd*self.rho0*np.exp(-z/self.H)*\
                               self.r2A(radius)/(2*self.rrho2m(radius, density)*np.sin(angle_rad)))

    #analytical dvdz
    def dvdz(self, radius, velocity, density, strength, angle, z, degree=True):
        '''
        ANALYTICAL SOLUTION!
        
        Calculating velocity from altitude.
        
        ----------input-----------
        
        v0, initial speed, m/s
        
        A, cross-sectional area, m**2
        
        m, mass in kg
        
        angle, angle of injection, in radians
        
        z, altitude in metres, usually an array
        
        Cd, drag coefficient, preset as 1
        
        rho0, density of air, preset as 1.2kg/m**3
        
        H, height of atmosphere, preset as 8000m
        
        radian, boolean datum that tells whether the angle is given in radians; preset as False 
        
        ----------output-----------
        
        change in speed per unit altitude in /s, usually an array
        '''
        if degree:
            angle_rad = angle*np.pi/180.
        else:
            angle_rad = angle
    #print(angle)
        return self.vz(radius, velocity, density, strength, angle, z, degree=False)*\
                self.Cd*self.rho0*np.exp(-z/self.H)*self.r2A(radius)/(2*self.rrho2m(radius,density)*np.sin(angle_rad))

    #analytical dedz
    def dedz(self, radius, velocity, density, strength, angle, z, degree = True):
        '''
        ANALYTICAL SOLUTION!
        
        Calculating velocity from altitude.
        
        ----------input-----------
        
        v0, initial speed, m/s
        
        A, cross-sectional area, m**2
        
        m, mass in kg
        
        angle, angle of injection, in radians
        
        z, altitude in metres, usually an array
        
        Cd, drag coefficient, preset as 1
        
        rho0, density of air, preset as 1.2kg/m**3
        
        H, height of atmosphere, preset as 8000m
        
        radian, boolean datum that tells whether the angle is given in radians; preset as False 
        
        ----------output-----------
        
        change in energy per unit length in kT/km, usually an array
        '''
        if degree:
            angle_rad=angle*np.pi/180.
        else:
            angle_rad=angle
    #print(angle)
        return self.dvdz(radius, velocity, density, strength, angle_rad, z, degree=False)*\
                self.vz(radius, velocity, density, strength, angle_rad, z, degree=False)*self.rrho2m(radius, density)/4.184e9   #J/m to kT/km


    def deg2rad(self, angle):#function which converts from degres to radians
        return np.pi*angle/180.

    def fun(self, t, state, strength):
        """
        RHS function for impact system
        """
        f = np.zeros_like(state)
        # unpack the state vector, which is:
        # velocity, density*volume, angle, init_altitude, x, radius
        velo, mass, theta, my_z, my_x, radius = state
        #rhoa = self.rho0 * np.exp(-z/self.H)
        rhom = 3000
        #V = (4*np.pi*r**3)/3#calculate volume to get the mass from density
        A = np.pi*radius**2#calculate cross sectional area
        
        f[0] = ((-self.Cd*self.rhoa(my_z)*A*velo**2)/(2*mass)) + self.g*np.sin(theta)    #dv/dt
        f[1] = (-self.Ch*self.rhoa(my_z)*A*velo**3)/(2*self.Q) #dm/dt
        f[2] = (self.g*np.cos(theta))/velo-(self.Cl*self.rhoa(my_z)*A*velo)/(2*mass)-\
            (velo*np.cos(theta))/(self.Rp+my_z) #dtheta/dt
        f[3] = -velo*np.sin(theta) #dz/dt
        f[4] = (velo*np.cos(theta))/(1+ (my_z)/self.Rp) #dx/dt
        if self.rhoa(my_z)*velo**2 >= strength:
            f[5] = (7/2*self.alpha* self.rhoa(my_z)/rhom)**0.5 * velo #dr/dt if needed
        else:
            f[5] = 0
        return f

    def imp_eu(self, fun, radius, velocity, density, angle, strength,
               init_altitude, dt ,t_max,t0, x):
        """
        the improved Euler function from the lecture
        we have modified u to include all 5 of our variables
        """
        volume = 4.*np.pi*radius**3/3.
        my_u = np.array([velocity, density*volume, angle, init_altitude, x, radius])
        time = np.array(t0)
        u_all = [[velocity, density*volume, angle, init_altitude, x, radius]]
        t_all = [t0]
        while (time < t_max) and (my_u[0] > 0) and (my_u[1] > 0) and (my_u[3] > 0):
            ue = my_u + dt*fun(time, my_u, strength) 
            my_u = my_u + 0.5*dt* (fun(time, my_u, strength) + fun(time + dt, ue, strength))
            u_all.append(my_u)
            time = time + dt
            t_all.append(time)
            
        return np.array(u_all), np.array(t_all)


    def impact(self, radius, velocity, density, strength, angle,
               init_altitude = 100e3, dt = 0.1, radians = False):
        """
        Solve the system of differential equations for a given impact event.
        Also calculates the kinetic energy lost per unit altitude and
        analyses the result to determine the outcome of the impact.
        Parameters
        ----------
        radius : float
            The radius of the asteroid in meters
        velocity : float
            The entery speed of the asteroid in meters/second
        density : float
            The density of the asteroid in kg/m^3
        strength : float
            The strength of the asteroid (i.e., the ram pressure above which
            fragmentation and spreading occurs) in N/m^2 (Pa)
        angle : float
            The initial trajectory angle of the asteroid to the horizontal
            By default, input is in degrees. If 'radians' is set to True, the
            input should be in radians
        init_altitude : float, optional
            Initial altitude in m
        dt : float, optional
            The output timestep, in s
        radians : logical, optional
            Whether angles should be given in degrees or radians. Default=False
            Angles returned in the DataFrame will have the same units as the
            input
        Returns
        -------
        Result : DataFrame
            A pandas DataFrame containing the solution to the system.
            Includes the following columns:
            ``velocity``, ``mass``, ``angle``, ``altitude``,
            ``distance``, ``radius``, ``time``, ``dedz``
        outcome : Dict
            dictionary with details of airburst and/or cratering event.
            For an airburst, this will contain the following keys:
            ``burst_peak_dedz``, ``burst_altitude``, ``burst_total_ke_lost``.
            For a cratering event, this will contain the following keys:
            ``impact_time``, ``impact_mass``, ``impact_speed``.
            All events should also contain an entry with the key ``outcome``,
            which should contain one of the following strings:
            ``Airburst``, ``Cratering`` or ``Airburst and cratering``
        """
        res1 = self.solve_atmospheric_entry(radius, velocity, density, strength, angle,init_altitude, dt,radians)
        res2 = self.calculate_energy(res1)
        res3 = self.analyse_outcome(res2)

        return res2, res3

    def solve_atmospheric_entry(
            self, radius, velocity, density, strength, angle,
            init_altitude=100e3, dt=0.05, radians=False):
        """
        Solve the system of differential equations for a given impact scenario
        Parameters
        ----------
        radius : float
            The radius of the asteroid in meters
        velocity : float
            The entery speed of the asteroid in meters/second
        density : float
            The density of the asteroid in kg/m^3
        strength : float
            The strength of the asteroid (i.e., the ram pressure above which
            fragmentation and spreading occurs) in N/m^2 (Pa)
        angle : float
            The initial trajectory angle of the asteroid to the horizontal
            By default, input is in degrees. If 'radians' is set to True, the
            input should be in radians
        init_altitude : float, optional
            Initial altitude in m
        dt : float, optional
            The output timestep, in s
        radians : logical, optional
            Whether angles should be given in degrees or radians. Default=False
            Angles returned in the DataFrame will have the same units as the
            input
        Returns
        -------
        Result : DataFrame
            A pandas DataFrame containing the solution to the system.
            Includes the following columns:
            ``velocity``, ``mass``, ``angle``, ``altitude``,
            ``distance``, ``radius``, ``time``
        """
        if not radians:
            angle = self.deg2rad(angle)
        variables, times = self.imp_eu(self.fun, radius, velocity, density, angle, strength,
               init_altitude, dt, t_max=1e6, t0=0, x=0)
        velocity = variables[:, 0]
        mass = variables[:, 1]
        angle = variables[:, 2]
        init_altitude = variables[:, 3]
        x = variables[:, 4]
        radius = variables[:, 5]

        return pd.DataFrame({'velocity': velocity,
                         'mass': mass,
                         'angle': (angle*180/np.pi),
                         'altitude': init_altitude,
                         'distance': x,
                         'radius': radius,
                         'time': times})

        

    def calculate_energy(self, result):
        """
        Function to calculate the kinetic energy lost per unit altitude in
        kilotons TNT per km, for a given solution.
        Parameters
        ----------
        result : DataFrame
            A pandas DataFrame with columns for the velocity, mass, angle,
            altitude, horizontal distance and radius as a function of time
        Returns
        -------
        Result : DataFrame
            Returns the DataFrame with additional column ``dedz`` which is the
            kinetic energy lost per unit altitude
        """
        energy = 0.5 * result['mass'] * result['velocity']**2

        dedz = energy.diff(-1)/result.altitude.diff(-1)
        dedz = dedz/(4.184*1e9)

        result = result.copy()
        result.insert(len(result.columns),
                      'dedz', dedz)
        return result


    def analyse_outcome(self, result):
        """
        Inspect a prefound solution to calculate the impact and airburst stats
        Parameters
        ----------
        result : DataFrame
            pandas DataFrame with velocity, mass, angle, altitude, horizontal
            distance, radius and dedz as a function of time
        Returns
        -------
        outcome : Dict
            dictionary with details of airburst and/or cratering event.
            For an airburst, this will contain the following keys:
            ``burst_peak_dedz``, ``burst_altitude``, ``burst_total_ke_lost``.
            For a cratering event, this will contain the following keys:
            ``impact_time``, ``impact_mass``, ``impact_speed``.
            All events should also contain an entry with the key ``outcome``,
            which should contain one of the following strings:
            ``Airburst``, ``Cratering`` or ``Airburst and cratering``
        """

        # Enter your code here to process the result DataFrame and
        # populate the outcome dictionary.
        peak_dedz = result.dedz.max()
        p_of_burst = result.dedz.idxmax()
        burst_altitude = result['altitude'][p_of_burst]
        total_ke_altitude = (0.5 * result['mass'][0] * result['velocity'][0]**2 -\
                            0.5 * result['mass'][p_of_burst] * result['velocity'][p_of_burst]**2)/(4.184e12) #kT
        impact_time = result['time'][p_of_burst]
        impact_mass = result['mass'][p_of_burst]
        impact_speed = result['velocity'][p_of_burst]
        outcome = {}
        if burst_altitude > 5000:
            outcome['outcome'] = 'Airburst'
            outcome['burst_peak_dedz'] = peak_dedz
            outcome['burst_altitude'] = burst_altitude
            outcome['burst_total_ke_lost'] = total_ke_altitude
        elif burst_altitude > 0 and burst_altitude < 5000:
            outcome['outcome'] = 'Airburst and cratering'
            outcome['burst_peak_dedz'] = peak_dedz
            outcome['burst_altitude'] = burst_altitude
            outcome['burst_total_ke_lost'] = total_ke_altitude
        else:
            outcome['outcome'] = 'Cratering'
            outcome['impact_time'] = impact_time
            outcome['impact_mass'] = impact_mass
            outcome['impact_speed'] = impact_speed
        return outcome
