#!/usr/bin/env python
# ------------------------------------------------------------------------
# Copyright 2018,  Alberto Pietrini
# Statisticalhitfinder is distributed under the terms of the Simplified BSD License.
# -------------------------------------------------------------------------

import numpy as np
from scipy import stats

def interval(mu,sigma):
    """ methods of estimators """
    a = mu - 4*sigma
    b = mu + 4*sigma
    return a, b

# FOR PR772 I also fit each bin
def fitfun(x,y,x1,y1):
    p3 = np.polyfit(x1,y1,3)
    fit_bkg = np.dot(p3,np.array([x**3,x**2,x,1]))
    return y - fit_bkg

def gauss(x, mu, sigma):
    return np.exp(-((x-mu)/sigma)**2/2.)/(sigma*np.sqrt(2*np.pi))

def CDF(lim,mu,sigma):
    return stats.norm.cdf((lim-mu)/sigma,mu,sigma)

def PDF(lim,mu,sigma):
    return stats.norm.pdf((lim-mu)/sigma,mu,sigma)

def moment1(a,b,mu,sigma):
    #return x*gauss(x,mu,sigma)
    return mu + sigma*(PDF(a,0,1) - PDF(b,0,1))/(CDF(b,0,1) - CDF(a,0,1))

def moment2(a,b,mu,sigma):
    return sigma*np.sqrt(1 + (((a-mu)/sigma)*PDF(a,0,1) - ((b-mu)/sigma)*PDF(b,0,1))/(CDF(b,0,1) - CDF(a,0,1)) \
           - (PDF(a,0,1) - PDF(b,0,1))/(CDF(b,0,1) - CDF(a,0,1)))

def standard_scores(x,mu,sigma):
    return (x-mu)/sigma

def prel_misses(x, y, mask):
    """ Here x is the photon count
        y is the gmd in the specific bin
        x is the fitted and subtracted photon count
        This is an adaptation for PR772 run 153"""
    it = 0
    mu0, sigma0, sigma1 = 0., -200., 1.
    select_misses = np.ones(x.shape[0]).astype('bool')

    while (x[select_misses].shape[0]>=100)*(sigma1!=sigma0):
        z = fitfun(y,x,y[select_misses],x[select_misses])

        if sigma0 == -200:
            #z = standard_scores(nyy,nyy.mean(),nyy.std())
            mu0 = np.mean(z[select_misses])
            sigma0 = np.std(z[select_misses])
            a,b = interval(mu0,sigma0)
        else:
            sigma1 = sigma0
            mu0, sigma0 = np.mean(z[select_misses]), np.std(z[select_misses])

            mu0 = moment1(a,b,mu0,sigma0)
            sigma0 = moment2(a,b,mu0,sigma0)
            a,b = interval(mu0,sigma0)
        
        """ 7000 and 450 for OmRV
        # 200000 and ??? for PR772
        if mask.sum()==388*370:
            print "NOMASK"
            select_misses = (z>a)*(z<b)*(x>200000)*(y>2.)
        else:
            select_misses = (z>a)*(z<b)*(y>2.)
        """
        # RNA polymerase II adaptation
        if mask.sum()==388*370:
            print "NOMASK"
            select_misses = (z>a)*(z<b)*(x>400)*(y>1.)
        else:
            select_misses = (z>a)*(z<b)*(y>1.)



        #print "RANK -", rank, "ITERATION", it, "SIGMA/MU", sigma0, mu0, "NR. MISSES", select_misses.sum()
        it +=1
    return select_misses
