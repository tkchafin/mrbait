#!/usr/bin/python

import pandas as pd
import random
from collections import Counter
from collections import defaultdict
import time
import string

"""RESULTS:

Conclusion: Basic python dict is SO MUCH FASTER.

-------Dataframe lookup tests--------
100 lookups in dataframe size 100:
31 ms
1000 lookups in dataframe size 1000:
312 ms
10000 lookups in dataframe size 10000:
3396 ms


-------Dictionary lookup tests--------
100 lookups in dict size 100:
0 ms
1000 lookups in dict size 1000:
0 ms
10000 lookups in dict size 10000:
0 ms


-------Dictionary lookup tests (strings as keys)--------
100 lookups in dict size 100:
0 ms
1000 lookups in dict size 1000:
0 ms
10000 lookups in dict size 1000:
1 ms

"""



def time_me(method):
    def wrapper(*args, **kw):
        startTime = int(round(time.time() * 1000))
        result = method(*args, **kw)
        endTime = int(round(time.time() * 1000))

        print(endTime - startTime,'ms')
        return result

    return wrapper

@time_me
#lookup based on pandas df
def pandas_lookup(df, num):
	for i in range(num):
		new_num = df.loc[df["id"]==i]["weight"]

@time_me
#lookup based on base python dict (with INT as keys)
def dict_lookup(d, num):
	for key in d:
		new_num = d[key]

#Function to make a random df
def make_random_df(size):
	ids = list(range(size))
	rands = random.sample(range(size*10),size)
	df = pd.DataFrame({'id' : ids,'weight' : rands})
	return df

#Function to make a random df
def make_random_dict(size):
	ids = random.sample(range(size*10),size)
	rands = random.sample(range(size*10),size)
	d = {}
	counter = 0
	for i in ids:
		d[i] = rands[counter]
		counter +=1
	return d

def randomword(length):
   letters = string.ascii_lowercase
   return ''.join(random.choice(letters) for i in range(length))

#Function to make a random df
def make_random_dict_string(size):
	rands = random.sample(range(size*10),size)
	d = {}
	counter = 0
	for i in rands:
		key = randomword(10)
		d[i] = rands[counter]
		counter +=1
	return d

#Tests of functions

print("-------Dataframe lookup tests--------")
print("100 lookups in dataframe size 100:")
df = make_random_df(100)
pandas_lookup(df, 100)
print("1000 lookups in dataframe size 1000:")
df = make_random_df(1000)
pandas_lookup(df, 1000)
print("10000 lookups in dataframe size 10000:")
df = make_random_df(10000)
pandas_lookup(df, 10000)

print("\n\n-------Dictionary lookup tests--------")
print("100 lookups in dict size 100:")
d = make_random_dict(100)
dict_lookup(d, 100)
print("1000 lookups in dict size 1000:")
d = make_random_dict(1000)
dict_lookup(d, 1000)
print("10000 lookups in dict size 10000:")
d = make_random_dict(10000)
dict_lookup(d, 10000)

print("\n\n-------Dictionary lookup tests (strings as keys)--------")
print("100 lookups in dict size 100:")
d = make_random_dict_string(100)
dict_lookup(d, 100)
print("1000 lookups in dict size 1000:")
d = make_random_dict_string(1000)
dict_lookup(d, 1000)
print("10000 lookups in dict size 1000:")
d = make_random_dict_string(10000)
dict_lookup(d, 10000)
print("\n\n")
