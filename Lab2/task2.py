
# coding: utf-8

# #### Imports

# In[122]:

import math
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
from pandas import Series, DataFrame


# #### File input

# In[123]:

with open("input.txt") as fin:
    lcg1Params = map(lambda el: int(el), fin.readline().split())
    lcg2Params = map(lambda el: int(el), fin.readline().split())
    lcg3Params = map(lambda el: int(el), fin.readline().split())


# ### Constants

# In[124]:

big_N = 10 ** 6
small_N = 10 ** 4


# ### Generators

# #### Linear congruential generator

# In[125]:

def LCG(module, a, c, seed):
    while True:
        seed = (seed * a + c) % module
        yield seed / float(module)


# #### MacLaren-Marsaglia - iterating option

# In[126]:

def MMG(moduleX, aX, cX, seedX, moduleY, aY, cY, seedY, k, n):
    randX, randY = LCG(seedX, aX, cX, moduleX), LCG(seedY, aY, cY, moduleY)
    X, Y = [randX.next() for i in range(n + k)], [randY.next() for i in range(n)]
    
    v, result = X[:k], []

    for i in range(n):
        result.append(v[int(Y[i] * k)])
        v[int(Y[i] * k)] = X[k + i]
        
    return result


# #### MacLaren-Marsaglia - generator option

# In[127]:

def MMG(moduleX, aX, cX, seedX, moduleY, aY, cY, seedY, k):
    randX, randY = LCG(seedX, aX, cX, moduleX), LCG(seedY, aY, cY, moduleY)
    v = [randY.next() for i in range(k)]
    
    while True:
        x, j = randY.next(), int(randX.next() * k)
        res, v[j] = v[j], x
        yield res


# #### Bernulli

# In[128]:

def bernulli(key_params, value_params, probability):
    mmg = MMG(*(key_params + value_params + [100]))
    while True:
        yield (1 if mmg.next() <= probability else 0)


# #### Test

# In[129]:

ber = bernulli(lcg3Params, lcg3Params, .7)

plt.hist([ber.next() for i in range(100)])
plt.show()


# #### Pirson test

# In[130]:

def pirson_for_bernulli(res, probability):
    count = len(res)
    
    expectation = np.sum(res) / count
    dispersion = np.sum([(i - expectation) ** 2 for i in res]) / (count - 1)    
    ones = np.count_nonzero(res)     

    return (ones - count) ** 2 / count / probability + (count - ones - count * (1 - probability)) ** 2 / count / (1 - probability)


# #### Discrete

# In[131]:

def discrete_distribution(key_params, value_params, a, b):
    mmg = MMG(*(key_params + value_params + [100]))
    probability = 1. / (b - a + 1)
    
    while True:
        yield a + int(mmg.next() / probability)


# #### Test

# In[132]:

dscr = discrete_distribution(lcg3Params, lcg3Params, 1, 100)

plt.hist([dscr.next() for i in range(10000)], bins="auto")
plt.show()


# #### Pirson test

# In[133]:

def pirson_for_discrete(res, a, b):
    m, n = b - a + 1, len(res)
    
    n_arr = np.zeros(m)
    for i in res:
        n_arr[i - a] += 1

    return np.sum([(i - n / m) ** 2 for i in n_arr]) / (n / m)


# #### Random

# In[134]:

plt.hist(np.random.randint(1, 100, size=10000))
plt.show()


# In[135]:

ber = bernulli(lcg3Params, lcg3Params, .7)
dscr = discrete_distribution(lcg3Params, lcg3Params, 1, 10)

ser = Series([pirson_for_bernulli([ber.next() for i in range(100)], .7), pirson_for_bernulli(np.random.randint(0, 2, size=100), .7), pirson_for_discrete([dscr.next() for i in range(100)], 1, 10), pirson_for_discrete(np.random.randint(1, 10, size=100), 1, 10)], index=["Bernulli MM", "Bernulli Random", "Discrete MM", "Discrete Random"])


# In[136]:

DataFrame({ "Statistic": ser })

