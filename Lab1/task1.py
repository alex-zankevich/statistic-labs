
# coding: utf-8

# #### Imports

# In[388]:

import math
import numpy as np
import matplotlib.pyplot as plt


# #### File input

# In[389]:

with open("input.txt") as fin:
    lcg1Params = map(lambda el: int(el), fin.readline().split())
    lcg2Params = map(lambda el: int(el), fin.readline().split())
    lcg3Params = map(lambda el: int(el), fin.readline().split())


# ### Constants

# In[390]:

big_N = 10 ** 6
small_N = 10 ** 4


# ### Generators

# #### Linear congruential generator

# In[391]:

def LCG(module, a, c, seed):
    while True:
        seed = (seed * a + c) % module
        yield seed / float(module)


# #### MacLaren-Marsaglia - iterating option

# In[392]:

def MMG(moduleX, aX, cX, seedX, moduleY, aY, cY, seedY, k, n):
    randX, randY = LCG(seedX, aX, cX, moduleX), LCG(seedY, aY, cY, moduleY)
    X, Y = [randX.next() for i in range(n + k)], [randY.next() for i in range(n)]
    
    v, result = X[:k], []

    for i in range(n):
        result.append(v[int(Y[i] * k)])
        v[int(Y[i] * k)] = X[k + i]
        
    return result


# #### MacLaren-Marsaglia - generator option

# In[393]:

def MMG(moduleX, aX, cX, seedX, moduleY, aY, cY, seedY, k):
    randX, randY = LCG(seedX, aX, cX, moduleX), LCG(seedY, aY, cY, moduleY)
    v = [randY.next() for i in range(k)]
    
    while True:
        x, j = randY.next(), int(randX.next() * k)
        res, v[j] = v[j], x
        yield res


# ### Tests

# #### Similar moments test

# In[394]:

def similar_moments_test(values):
    delta = 1.96
    n = len(values)
    
    m = sum(values) / float(n)
    sqS = sum([((values[i] - m) ** 2.) for i in range(n)]) / (n - 1.)
    
    ksiM, ksiS = m - .5, sqS - 1. / 12
    
    c1, c2 = math.sqrt(12. * n), ((n - 1.) / n) * ((0.0056 / n + 0.0028 / n ** 2 - 0.0038 / n ** 3) ** -0.5)
    
    return m, sqS, ksiM, ksiS, c1, c2


# #### Usage

# In[395]:

def logResult(generator, k):
    m, sqS, ksiM, ksiS, c1, c2 = similar_moments_test([generator.next() for i in range(k)])
    stat1 = c1 * abs(ksiM)
    stat2 = c2 * abs(ksiS)
    check_stat = lambda stat: "H0" if stat < 1.96 else "H1"
    
    print """Math Expectation: {} - hypothesis {} 
Dispersion: {} - hypothesis {}\n""".format(stat1, check_stat(stat1), stat2, check_stat(stat2))

def runLoggers(N):
    print "N is {}".format(N)
    [logResult(gen, N) for gen in [lcg1, lcg2, lcg3, mmg]]

lcg1, lcg2, lcg3, mmg = LCG(*lcg1Params), LCG(*lcg2Params), LCG(*lcg3Params), MMG(*(lcg2Params + lcg3Params + [100]))

runLoggers(small_N)
runLoggers(big_N)


# #### Covariation test

# In[396]:

def covariation_test(values):
    delta = 1.96
    bound = 10
    n = len(values)
    
    count_multi = lambda j: sum([values[i] * values[i + j] for i in range(n - j)])
    
    m = float(sum(values) / n)
    r = [1. / 12] + [0 for i in range(bound - 1)]
    _r = [((1. / (n - i - 1)) * count_multi(i) - ((n * m ** 2.) / (n - 1))) for i in range(0, bound)]
    
    h = [math.sqrt(2) * (delta / (12 * math.sqrt(n - 1)))] + [delta / (12 * math.sqrt(n - 1)) for i in range(bound - 1)]
    
    return m, r, _r, h


# #### Usage

# In[397]:

def logResult(generator, k):
    m, r, _r, h = covariation_test([generator.next() for i in range(k)])
    
    print "-" * 10
    
    for i in range(len(r)):
        print abs(r[i] - _r[i])
    
    print "-" * 10
    print h
    print "-" * 10
    
def runLoggers(N):
    print "N is {}".format(N)
    [logResult(gen, N) for gen in [lcg1, lcg2, lcg3, mmg]]

lcg1, lcg2, lcg3, mmg = LCG(*lcg1Params), LCG(*lcg2Params), LCG(*lcg3Params), MMG(*(lcg2Params + lcg3Params + [100]))

runLoggers(small_N)
runLoggers(big_N)


# #### Pirson test

# In[398]:

def pirson_test(values):
    delta = 1073.64
    n, bound = len(values), 1000
    
    n_arr = np.full(bound, 0)
    for i in range(n):
        n_arr[int(values[i] * 1000)] += 1
    
    s = sum([(n_arr[i] - float(n) / bound) ** 2 for i in range(bound)]) / (n / bound)
    
    return s


# #### Usage

# In[399]:

def logResult(generator, k):
    stat = pirson_test([generator.next() for i in range(k)])

    print stat

def runLoggers(N):
    print "N is {}".format(N)
    [logResult(gen, N) for gen in [lcg1, lcg2, lcg3, mmg]]

lcg1, lcg2, lcg3, mmg = LCG(*lcg1Params), LCG(*lcg2Params), LCG(*lcg3Params), MMG(*(lcg2Params + lcg3Params + [100]))

runLoggers(small_N)
runLoggers(big_N)


# ### Plots

# In[400]:

N = 1000

lcg = LCG(*lcg1Params)
dx = [lcg.next() for i in range(N)]

lcg1, lcg2, lcg3, mmg = LCG(*lcg1Params), LCG(*lcg2Params), LCG(*lcg3Params), MMG(*(lcg2Params + lcg3Params + [100]))

def show_plots (generator):      
    gen_values = [generator.next() for i in range(N + 1)]
    
    dx = gen_values[:-1]
    dy = gen_values[1:]

    plt.plot(dx, dy, 'ro')
    plt.show()
    
[show_plots(gen) for gen in [lcg1]]


# In[ ]:



