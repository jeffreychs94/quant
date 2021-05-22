'''
This is to calculate swap points

'''

# Rate Curve Assigned to # FIXME:

import math as m

lCurr1_DsRate = []
lCurr2_DsRate = []

def fCalculateSwapPoints (spot,tenure,fBaseRate,fQuoteRate):
    time_frac = tenure/365
    forward_point = spot * (1+fQuoteRate/100*time_frac)/(1+fBaseRate/100*time_frac)
    swap = forward_point - spot
    return swap * 10000


print(fCalculateSwapPoints(9.86,365,0.0323,0.065))
