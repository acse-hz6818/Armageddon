from math import erf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter


def solve_ensemble(
        planet,
        fiducial_impact,
        variables,
        radians=False,
        rmin=8, rmax=12,
        ):
    """
    Run asteroid simulation for a distribution of initial conditions and
    find the burst distribution

    Parameters
    ----------

    planet : object
        The Planet class instance on which to perform the ensemble calculation

    fiducial_impact : dict
        Dictionary of the fiducial values of radius, angle, strength, velocity
        and density

    variables : list
        List of strings of all impact parameters to be varied in the ensemble
        calculation

    rmin : float, optional
        Minimum radius, in m, to use in the ensemble calculation,
        if radius is one of the parameters to be varied.

    rmax : float, optional
        Maximum radius, in m, to use in the ensemble calculation,
        if radius is one of the parameters to be varied.

    Returns
    -------

    ensemble : DataFrame
        DataFrame with columns of any parameters that are varied and the
        airburst altitude
    """
    n_p = 20#number of points in our distribution

    my_radius = np.array([fiducial_impact['radius']])
    angle = np.array([fiducial_impact['angle']])
    my_y = np.array([fiducial_impact['strength']])
    my_v = np.array([fiducial_impact['velocity']])
    my_d = np.array([fiducial_impact['density']])

    pdf_radius = np.array([1])
    pdf_angle = np.array([1])
    pdf_y = np.array([1])
    pdf_v = np.array([1])
    pdf_d = np.array([1])

    for my_var in variables:
        
        if my_var == 'radius':
            my_radius = np.linspace(rmin, rmax, n_p)
            cdf_radius = (my_radius-rmin)/(rmax-rmin)
            pdf_radius = cdf_radius[1:]-cdf_radius[:-1]
            pdf_radius = np.append(pdf_radius, 0)

        if my_var == 'angle':
            angle_min = 0.26
            angle_max = np.pi/2 - 0.1
            angle = np.linspace(angle_min, angle_max, n_p)
            cdf_angle = 1-np.cos(angle)*np.cos(angle)
            pdf_angle = cdf_angle[1:]-cdf_angle[:-1]
            pdf_angle = np.append(pdf_angle, 0)

        if my_var == 'strength':
            ymin = 10**3
            ymax = 10*10**6
            my_y = np.logspace(np.log10(ymin), np.log10(ymax), n_p)
            cdf_y = np.log10(my_y/ymin)/np.log10(ymax/ymin)
            pdf_y = cdf_y[1:]-cdf_y[:-1]
            pdf_y = np.append(pdf_y, 0)

        if my_var == 'velocity':
            vmin = 11000
            vmax = 50000
            my_a = 11000
            my_v = np.linspace(vmin, vmax, n_p)
            error = np.asarray([erf(v2/(np.sqrt(2)*my_a)) for v2 in my_v])
            cdf_v = error-(my_v/my_a)*np.exp((-my_v**2)/(2*my_a**2))*np.sqrt(2/np.pi)
            pdf_v = cdf_v[1:]-cdf_v[:-1]
            pdf_v = np.append(pdf_v, 0)

        if my_var == 'density':
            dmin = 1000
            dmax = 7000
            my_d = np.linspace(dmin, dmax, n_p)
            r_m = 3000
            sig = 1000
            cdf_d = np.asarray([1/2+1/2*erf((d2-r_m)/(sig*np.sqrt(2))) for d2 in my_d])
            pdf_d = cdf_d[1:]-cdf_d[:-1]
            pdf_d = np.append(pdf_d, 0)

    i_r_counter = 0
    i_a_counter = 0
    i_y_counter = 0
    i_v_counter = 0
    i_d_counter = 0
    results_dataframe = pd.DataFrame(columns=['radius']+['angle']+['strength']+['velocity']+['density']+['burst_altitude']+['probability'], index=range(0))#prepare an dataframe to hold the results
    # Implement your ensemble function here
    for radi in my_radius:
        i_a_counter = 0
        for ang in angle:
            i_y_counter = 0
            for stre in my_y:
                i_v_counter = 0
                for velo in my_v:
                    i_d_counter = 0
                    for den in my_d:
                        outcome, outcome_dic = planet.impact(radi, velo, den, stre, ang)
                        if 'burst_altitude' in outcome_dic.keys():#if needed otherwise we fail test for empty data
                            br_al = outcome_dic["burst_altitude"]
                            br_al_prob = pdf_radius[i_r_counter]*pdf_angle[i_a_counter]*pdf_y[i_y_counter]*pdf_v[i_v_counter]*pdf_d[i_d_counter]
                            row_in_dataframe = pd.DataFrame([[radi, ang, stre, velo, den, br_al, br_al_prob]], columns=['radius']+['angle']+['strength']+['velocity']+['density']+['burst_altitude']+['probability'])
                            results_dataframe = results_dataframe.append(row_in_dataframe)
                        i_d_counter += 1
                    i_v_counter += 1
                i_y_counter += 1
            i_a_counter += 1
        i_r_counter += 1
    data_frame_to_return = results_dataframe[[my_v for my_v in variables]+['burst_altitude']]#data frame with varied parameters and burst altitudes
    return data_frame_to_return

    
