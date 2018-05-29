# -*- coding: utf-8 -*-
# this script runs a gillespie simulation for the relaxed replication model
# it produces a plot showing copy number distribution for m and w across many iterations
import numpy as np
import matplotlib.pyplot as plt
import random
import sys
import scipy.stats as stats
from math import log
from pdb import set_trace
from matplotlib import colors as mcolors
from matplotlib.ticker import PercentFormatter

# reactions for birth/death process, i.e. what happens to w/m when each reaction is selected
def Reactions(selected_reaction,w,m):
    if 'h1' == selected_reaction:
        w += 1
    elif 'h2' == selected_reaction:
        if w<=0:
            w = w
        else:
            w -= 1
    elif 'h3' == selected_reaction:
        m += 1
    elif 'h4' == selected_reaction:
        if m<=0:
            m = m
        else:
            m -= 1
    else:
        raise Exception('Unrecognised reaction')
    return w,m

def RR_reprate(w,m,params):
    alpha = params['alpha']
    wopt = params['wopt']
    gamma = params['gamma']
    mu = params['mu']
    rep_rate = mu*(alpha*(wopt - w - gamma * m) + w + gamma * m)/float(w+m) # JA: use floats just to be safe and avoid integer division
    if rep_rate > 0:
        return rep_rate
    else:
        return 0

def LF_reprate(w,m,params):
    mu = params['mu']
    nss = params['nss']
    delta = params['delta']
    b = params['b']
    rep_rate = mu + b*(nss-w-delta*m)
    if rep_rate > 0:
        return rep_rate
    else:
        return 0

def Gillespie(t_max,w,m,params,model):
    reactions = ['h1','h2','h3','h4']
    mu = params['mu']
    # used for time points when reactions occurred
    times = []
    # ms and ws are used for storing quantities of mutant and wild type overtime
    ms = []
    ws = []
    t = 0
    # continue iterating as long as t is less than t_max and m/w populations are both greater than 0
    while t<t_max:
        # define replication rate of the system
        if model == 'RR':
            rep_rate = RR_reprate(w,m,params)
        elif model == 'LF':
            rep_rate = LF_reprate(w,m,params)
        else:
            raise Exception('Unrecognised model')

        # calculate hazards for each reaction
        h1 = rep_rate*w
        h2 = mu*w
        h3 = rep_rate*m
        h4 = mu*m
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
        selected_reaction = random.choices(reactions, weights = (p1,p2,p3,p4))[0]
        # update w and m depending on which reaction was chosen
        w,m = Reactions(selected_reaction,w,m)
        # add information to lists
        times.append(t)
        ws.append(w)
        ms.append(m)
    return w,m,times

TMAX = 110
x0 = [int(sys.argv[1]),int(sys.argv[2])]
w0 = x0[0]
m0 = x0[1]

params_LF = {
'mu': 0.2,
'b': 10**-4,
'nss': 1000,
'delta': 1
}

params_RR = {
'mu': 0.2,
'gamma' : 1, # JA: this parameter is analogous to delta, so choose its value to be 1
'alpha' : 2,
'wopt' : 1000
}

iterate = 500
w_relaxed = []
m_relaxed = []
w_linear = []
m_linear = []

# perform gillespie for both models as many times as value of iterate
for i in range(0, iterate):
    print(i)
    w_r,m_r,time_relaxed = Gillespie(TMAX,w0,m0,params_RR,'RR')
    w_l,m_l,time_linear = Gillespie(TMAX,w0,m0,params_LF,'LF')
    # add the number of w/m to list corresponding to the model
    w_relaxed.append(w_r)
    m_relaxed.append(m_r)
    w_linear.append(w_l)
    m_linear.append(m_l)


# plot
fig,(ax1, ax2) = plt.subplots(1, 2)

ax1.hist(w_relaxed,edgecolor='black',bins='fd',label="Wild types",alpha=0.5,color='orange')
ax1.hist(m_relaxed,edgecolor='black',bins='fd',label='Mutants',alpha=0.5,color='blue')
ax1.set_xlabel('Distribution of mtDNA')
ax1.set_ylabel('Frequency')
ax1.set_title('Relaxed Replication (t = 110)')
ax1.legend(loc='best')
ax2.hist(w_linear,edgecolor='black',bins='fd',label="Wild types",alpha=0.5,color='orange')
ax2.hist(m_linear,edgecolor='black',bins='fd',label='Mutants',alpha=0.5,color='blue')
ax2.set_xlabel('Distribution of mtDNA')
ax2.set_ylabel('Frequency')
ax2.set_title('Linear Feedback (t = 110)')
ax2.legend(loc='best')
fig.tight_layout()
fig.show()
