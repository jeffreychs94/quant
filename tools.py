import math as m
import numpy as np


def linearinterpolaterange( lx : list, ly: list, x:float):

    '''Simple Linear Interpolation with 2 range input'''

    for index , value in enumerate(lx):
        if x < value:
            print(1)
            return ly[index]
        elif x > lx[len(lx)-1]:
            print(2)
            return ly[len(lx)-1]
        else:
            return ((x - lx[index]) / (lx[index+1] - lx[index])) * (ly[index+1] - ly[index]) + ly[index]

def linear_cof(lx : list, ly: list):
    a,b = np.polyfit(np.array(lx),np.array(ly),1)
    return a,b

lTenure = [1,31,61,92]
lMarketData = [1.725,1.86,1.9,1.94]

test = linearinterpolaterange(lTenure,lMarketData,100)

print(test)
print(linear_cof(lTenure,lMarketData))
