# this script runs a gillespie simulation for the relaxed replication model
import numpy as np
import matplotlib.pyplot as plt
import random
import scipy.stats as stats
from scipy.interpolate import interp1d
from scipy.integrate import odeint
from math import log
import sys

# normal reactions for simple birth/death process
def ReactionsNormal(selected_reaction,w,m):
    x = w+m
    if 'h1' in selected_reaction:
        w += 1
    elif 'h2' in selected_reaction:
        if w<=0:
            w = w
        else:
            w -= 1
    elif 'h3' in selected_reaction:
        m += 1
    elif 'h4' in selected_reaction:
        if m<=0:
            m = m
        else:
            m -= 1
    return w,m

# simple gillespie
def GillespieRelaxed(t_max,w,m,tau,alp_ha,wopt,gamma):
    reactions = ['h1','h2','h3','h4']
# used for time points when reactions occurred
    times = []
# ms and ws are used for storing quantities of mutant and wild type overtime
    ms = []
    ws = []
    t = 0
    hs = []
    mu = 1/tau
# continue iterating as long as t is less than t_max and m/w populations are both greater than 0
    while t<t_max:
# define replication rate of the system
        rep_rate = (alp_ha*(wopt - w - gamma * m) + w + gamma * m)/tau*(w+m)
# check that lambda is greater than 0
        if rep_rate>0:
# calculate hazards for each reaction
            h1 = rep_rate
            h2 = mu
            h3 = rep_rate
            h4 = mu
# sum up all hazards
            h0 = h1 + h2 + h3 + h4
# calculate probabilities for each reaction occurring
            p1 = h1/h0
            p2 = h2/h0
            p3 = h3/h0
            p4 = h4/h0
# calculate time until next reaction
            t_prime = -log(np.random.uniform())/h0
            t = t_prime+t
# randomly select reaction using probability weights
            selected_reaction = random.choices(reactions, weights = (p1,p2,p3,p4))
# update w and m depending on which reaction was chosen
            w,m = ReactionsNormal(selected_reaction,w,m)
# calculate heteroplasmy
            h = m/(m+w)
# add information to lists
            hs.append(h)
            times.append(t)
            ws.append(w)
            ms.append(m)
# if lambda is less than 0, set lambda to 0 and continue
        else:
            rep_rate = 0
            h1 = rep_rate
            h2 = mu
            h3 = rep_rate
            h4 = mu
# sum up all hazards
            h0 = h1 + h2 + h3 + h4
# calculate probabilities for each reaction occurring
            p1 = h1/h0
            p2 = h2/h0
            p3 = h3/h0
            p4 = h4/h0
# calculate time until next reaction
            t_prime = -log(np.random.uniform())/h0
            t = t_prime+t
# randomly select reaction using probability weights
            selected_reaction = random.choices(reactions, weights = (p1,p2,p3,p4))
# update w and m depending on which reaction was chosen
            w,m = ReactionsNormal(selected_reaction,w,m)
# calculate heteroplasmy
            h = m/(m+w)
# add information to lists
            hs.append(h)
            times.append(t)
            ws.append(w)
            ms.append(m)
    return hs,ws,ms,times

# set parameters
TMAX = 1000
gam = 0
alpha = 2
tau = 5
w_opt = 1000
w0 = 5
m0 = 10

het_relaxed, w_relaxed, m_relaxed, time_relaxed = GillespieRelaxed(TMAX,w0,m0,tau,alpha,w_opt,gam)


plt.plot(w_relaxed, m_relaxed)
plt.show()
