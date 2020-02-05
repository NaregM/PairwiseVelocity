import numpy as np
import pandas as pd
from numpy import linalg as LA

def v12(data, number_of_bins, binsize, mass_cut):

    '''
    Direct calculation and Estimator

    data: Pandas dataframe
    '''
    h = 0.704
    L_box = 600./h

    boundry_condition = False

    M200c = data.M_200c.values
    ind = np.where(M200c > mass_cut)[0]
    M200c = M200c[ind]
    n_halos = M200c.size
    n_halosmax = n_halos
    
    XYZ = df[['x', 'y', 'z']].values
    V_xyz = df[['v_x', 'v_y', 'v_z']].values


    nbins = number_of_bins    #45
    binsize = binsize   #7
    rbins = binsize * (np.arange(nbins)+0.5)    
    n_of_r = np.zeros(nbins, dtype = np.int32)

    v_case1 = np.zeros(nbins)
    v_case2 = np.zeros(nbins)
    v_case3 = np.zeros(nbins)
    p_squared = np.zeros(nbins)
    q_squared = np.zeros(nbins)
    
    position = XYZ[ind]
    velocity = V_xyz
    
    print("There are {} halos".format(n_halosmax))

    for i in np.arange(n_halosmax):

        for j in np.arange(i):

            dr = position[j] - position[i]

            if boundry_condition:

                for k in range(3):

                    if dr[k] > 0.5 * L_box:

                        dr[k] -= L_box

                    elif dr[k] < -0.5 * L_box:

                        dr[k] += L_box

            dr_norm = np.sqrt(np.dot(dr, dr))
            dr_hat = dr / dr_norm

            p_AB = 0.5 * np.dot(dr_hat, r1_hat + r2_hat )
            r1_hat = position[i]/LA.norm(position[i])
            r2_hat = position[j]/LA.norm(position[j])
            q_AB = 0.5 * (2 * dr_hat - np.dot(dr_hat, r1_hat)*r1_hat - np.dot(dr_hat, r2_hat)*r2_hat)

            S1 = np.dot(velocity[i], r1_hat)
            S2 = np.dot(velocity[j], r2_hat)
            
            t1 = velocity[i] - S1*r1_hat
            t2 = velocity[j] - S2*r2_hat
            

            dv = velocity[j] - velocity[i]
            dv_case1 = np.dot(dv, dr_hat)

            ibin = np.floor(dr_norm/binsize)
            ibin = np.int32(ibin)
            
            if ibin < nbins:
                
                n_of_r[ibin] += 1

                v_case1[ibin] += dv_case1

                v_case2[ibin] += (S2 - S1) * p_AB
                
                v_case3[ibin] += np.dot(t1 - t2, q_AB)
                
                p_squared[ibin] += p_AB**2
                
                q_squared[ibin] += np.dot(q_AB, q_AB)

    ind = np.where(n_of_r > 0)[0]

    v_case1[ind] /= n_of_r[ind]
    v_case2[ind] /= p_squared[ind]
    v_case3[ind] /= q_squared[ind]
    err = v_case1[ind]/np.sqrt(n_of_r[ind])

    return rbins, v_case1, v_case2, v_case3, err
