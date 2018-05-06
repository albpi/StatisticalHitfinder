# ------------------------------------------------------------------------
# Copyright 2018,  Alberto Pietrini
# Statisticalhitfinder is distributed under the terms of the Simplified BSD License.
# -------------------------------------------------------------------------

import numpy as np
from scipy.stats import mode
from scipy.optimize import curve_fit



def GaussFit(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def find_1ph_peak(data_vector):
    # SET A MAX. VALUE FOR THE 0ph PEAK
    # max_ph0 = 80. works well for cxi73013 (probably working also for cxi06416?)
    max_ph0 = 80. # (usually is 80) modify so to be set externally (high or low gain)
    mode_of_datav = mode(data_vector)[0] # where the 0ph peak is expected to be, 
                                         # as so far we have more misses than hits

    min_dv, max_dv = data_vector.min(), data_vector.max()
    max_right_side_gauss = mode_of_datav+6

    # FIT THE 0ph PEAK
    reduced_vect = data_vector[data_vector<=max_right_side_gauss]
    y,x = np.histogram(data_vector, bins=range(min_dv, max_dv+1))
    x = x[:-1]

    z0,t0 = y[x<=max_right_side_gauss], x[x<=max_right_side_gauss]

    try:
        #print z0.max(), reduced_vect.mean(), np.abs(reduced_vect.std())
        popt0,pcov0 = curve_fit(GaussFit, t0, z0, p0=[z0.max(), reduced_vect.mean(), np.abs(reduced_vect.std())])
        #print popt0

        # height, mean and std. dev. of the fitted gaussian
        h0, m0, std0 = popt0[0], popt0[1], np.abs(popt0[2])

        # FIT THE 1ph PEAK

        # subtract the 0ph-peak
        y_ph0 = GaussFit(x,h0,m0,std0)
        y_sub = y-y_ph0

        # create histogram for 1-photon peak
        t1,z1 = x[(x>=m0+4*std0) * (x<max_ph0)], y_sub[(x>=m0+4*std0) * (x<max_ph0)]
        try:
            # define a function in which the height of the gaussian is already known
            # so that one can fit better the peak

            h1 = z1.max()
            def fit_gauss_1ph(xx,mean1,stdv1):
                return h1*np.exp(-(xx-mean1)**2/(2*stdv1**2))


            # define a reduced vector to have a better guess for mean/std in the fit
            red_vect1 = data_vector[(data_vector>=m0+4*std0) * (data_vector<max_ph0)]

            # actual fit
            popt1,pcov1 = curve_fit(fit_gauss_1ph, t1, z1, p0=[red_vect1.mean(),np.abs(red_vect1.std())])

            # height, mean and std. dev. of the fitted gaussian
            m1, std1 = popt1[0], np.abs(popt0[1])

            if m1<=0 or m1<t1.min() or m1>=max_ph0 or np.isnan(m1) or len(t1)<=20: m1, std1 = -1e309, -1309

        except Exception:
            h1, m1, std1 = -1e309, -1e309, -1e309
            pass

    except Exception:
        h0, m0, std0 = -1e309, -1e309, -1e309
        h1, m1, std1 = -1e309, -1e309, -1e309
        t1 = -1e309
        pass

    if np.isinf(m1):
        gain_value = -1e309
    else:
        gain_value = m1 - m0 # relative position of the 2 peaks

    return gain_value, m0, std0, np.mean(t0), np.std(t0), m1, std1, np.mean(t1), np.std(t1)

def baglivo(lambda0,lista,el):
    # BAGLIVO ALGORITHM
    N = el # nr. photons - here is the fitted photon value of the sample 
    nlambda = lambda0/np.sum(lambda0) # Normalized lambda values for true baglivo
    
    part_sum = lista*(np.log(lista)-np.log(nlambda)-np.log(N))
    part_sum[lista==0] = 0
    part_sum[np.isnan(part_sum) | np.isinf(part_sum) | np.isneginf(part_sum)] = 0
    logval = np.sum(part_sum)
    return logval
