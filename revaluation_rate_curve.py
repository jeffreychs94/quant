'''

This is for revaluation rate curve to calculate based on market rate input
- discount rates
- discount factors

function list
- market2discount

'''


import math as m

def testing(a):
    print(a)


testing("Print This")

def market2discount(marketrate,time):

    '''

    This only works for deposit rates, time is in days
    Returns Discount Rate and Discount Factor

    '''

    time_frac = time/365 #time fraction in years
    discountrate = m.log(1 + (marketrate/100)*time_frac)/time_frac
    discountfactor = 1/m.exp(discountrate*time_frac)
    return discountrate * 100 , discountfactor


lTenure = [1,31,61,92]
lMarketData = [1.725,1.86,1.9,1.94]
Output = []

for index,item in enumerate(lTenure):
    Output.append(market2discount(lMarketData[index] , lTenure[index]))


print(Output)
