# -*- coding: utf-8 -*-
# this script runs a gillespie simulation for the relaxed replication model
# it produces a plot showing copy number distribution for m and w across many iterations
import numpy as np
import matplotlib.pyplot as plt
import random
import sys
import scipy.stats as stats
import seaborn as sns
sns.set(color_codes=True)
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

def Gillespie(t_max,t1,w,m,params,model):
    reactions = ['h1','h2','h3','h4']
    mu = params['mu']
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
        if t1 < t < t1+0.1:
            w_t1 = np.copy(w)
            t1_point = t
    return w_t1,w,t1_point

TMAX_1 = 100.25
TMAX_2 = 100.5
TMAX_3 = 101
TMAX_4 = 102
TMAX_5 = 105
TMAX_6 = 110
t_1 = 100
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

iterate = 100
w1_relaxed_1 = []
w2_relaxed_1 = []

w1_relaxed_3 = []
w2_relaxed_3 = []

w1_relaxed_5 = []
w2_relaxed_5 = []

#w1_relaxed_2 = []
#w2_relaxed_2 = []
#w1_relaxed_4 = []
#w2_relaxed_4 = []
#w1_relaxed_6 = []
#w2_relaxed_6 = []

w1_linear_1 = []
w2_linear_1 = []

w1_linear_3 = []
w2_linear_3 = []

w1_linear_5 = []
w2_linear_5 = []

#w1_linear_2 = []
#w2_linear_2 = []
#w1_linear_4 = []
#w2_linear_4 = []
#w1_linear_6 = []
#w2_linear_6 = []

# perform gillespie for both models as many times as value of iterate
for i in range(0, iterate):
    print(i)
    w1_r_1,w2_r_1,t1_relaxed = Gillespie(TMAX_1,t_1,w0,m0,params_RR,'RR')
    w1_relaxed_1.append(w1_r_1)
    w2_relaxed_1.append(w2_r_1)

    w1_l_1,w2_l_1,t1_linear = Gillespie(TMAX_1,t_1,w0,m0,params_LF,'LF')
    w1_linear_1.append(w1_l_1)
    w2_linear_1.append(w2_l_1)

    w1_r_3,w2_r_3,t1_relaxed = Gillespie(TMAX_3,t_1,w0,m0,params_RR,'RR')
    w1_relaxed_3.append(w1_r_3)
    w2_relaxed_3.append(w2_r_3)

    w1_l_3,w2_l_3,t1_linear = Gillespie(TMAX_3,t_1,w0,m0,params_LF,'LF')
    w1_linear_3.append(w1_l_3)
    w2_linear_3.append(w2_l_3)

    w1_r_5,w2_r_5,t1_relaxed = Gillespie(TMAX_5,t_1,w0,m0,params_RR,'RR')
    w1_relaxed_5.append(w1_r_5)
    w2_relaxed_5.append(w2_r_5)

    w1_l_5,w2_l_5,t1_linear = Gillespie(TMAX_5,t_1,w0,m0,params_LF,'LF')
    w1_linear_5.append(w1_l_5)
    w2_linear_5.append(w2_l_5)
    #w1_r_2,w2_r_2,t1_relaxed = Gillespie(TMAX_2,t_1,w0,m0,params_RR,'RR')
    #w1_relaxed_2.append(w1_r_2)
    #w2_relaxed_2.append(w2_r_2)
    #w1_l_2,w2_l_2,t1_linear = Gillespie(TMAX_2,t_1,w0,m0,params_LF,'LF')
    #w1_linear_2.append(w1_l_2)
    #w2_linear_2.append(w2_l_2)
    #w1_r_4,w2_r_4,t1_relaxed = Gillespie(TMAX_4,t_1,w0,m0,params_RR,'RR')
    #w1_relaxed_4.append(w1_r_4)
    #w2_relaxed_4.append(w2_r_4)
    #w1_l_4,w2_l_4,t1_linear = Gillespie(TMAX_4,t_1,w0,m0,params_LF,'LF')
    #w1_linear_4.append(w1_l_4)
    #w2_linear_4.append(w2_l_4)
    #w1_r_6,w2_r_6,t1_relaxed = Gillespie(TMAX_6,t_1,w0,m0,params_RR,'RR')
    #w1_relaxed_6.append(w1_r_6)
    #w2_relaxed_6.append(w2_r_6)
    #w1_l_6,w2_l_6,t1_linear = Gillespie(TMAX_6,t_1,w0,m0,params_LF,'LF')
    #w1_linear_6.append(w1_l_6)
    #w2_linear_6.append(w2_l_6)


f, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2, sharex=True)

sns.kdeplot(w1_relaxed_1,w2_relaxed_1,cbar=True, ax=ax1)
ax1.set_title('Relaxed Replication model')
ax1.set_ylabel('Wild type copy number at T2 = 100.25')
ax1.set_ylim(850,1150)
ax1.set_xlim(850,1150)

sns.kdeplot(w1_linear_1,w2_linear_1,cbar=True, ax=ax2)
ax2.set_title('Linear Feedback Control model')
ax2.set_ylabel('Wild type copy number at T2 = 100.25')
ax2.set_ylim(850,1150)
ax2.set_xlim(850,1150)

sns.kdeplot(w1_relaxed_3,w2_relaxed_3,cbar=True, ax=ax3)
ax3.set_ylabel('Wild type copy number at T2 = 101')
ax3.set_ylim(850,1150)
ax3.set_xlim(850,1150)

sns.kdeplot(w1_linear_3,w2_linear_3,cbar=True, ax=ax4)
ax4.set_ylabel('Wild type copy number at T2 = 101')
ax4.set_ylim(850,1150)
ax4.set_xlim(850,1150)

sns.kdeplot(w1_relaxed_5,w2_relaxed_5,cbar=True, ax=ax5)
ax5.set_xlabel('Wild type copy number at T1 = 100')
ax5.set_ylabel('Wild type copy number at T2 = 105')
ax5.set_ylim(850,1150)
ax5.set_xlim(850,1150)

sns.kdeplot(w1_linear_5,w2_linear_5,cbar=True, ax=ax6)
ax6.set_xlabel('Wild type copy number at T1 = 100')
ax6.set_ylabel('Wild type copy number at T2 = 105')
ax6.set_ylim(850,1150)
ax6.set_xlim(850,1150)


f.show()
#fig = plt.figure()
#ax1 = fig.add_subplot(211)
#sns.kdeplot(w1_relaxed,w2_relaxed,cbar=True)
#ax1.set_title('Relaxed Replication model')
#ax1.set_xlabel('Wild type copy number at T1 = 100')
#ax1.set_ylabel('Wild type copy number at T2 = 100.25')
#ax1.set_ylim(850,1150)
#ax1.set_xlim(850,1150)

#ax2 = plt.subplot(212)
#sns.kdeplot(w1_linear,w2_linear,cbar=True)
#ax2.set_title('Linear Feedback Control model')
#ax2.set_xlabel('Wild type copy number at T1 = 100')
#ax2.set_ylabel('Wild type copy number at T2 = 100.25')
#ax2.set_ylim(850,1150)
#ax2.set_xlim(850,1150)#

#plt.show()
