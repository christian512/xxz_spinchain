"""
This module provides basic helper function for permutations of arrays filled with zeros and ones.
"""

import numpy as np
import math


def next_permutation(arr):
	"""Calculates next permutation from an array arr filled with zeros and ones"""
    # Find non-increasing suffix
    i = len(arr) - 1
    while i > 0 and arr[i - 1] >= arr[i]:
        i -= 1
    if i <= 0:
        return False

    # Find successor to pivot
    j = len(arr) - 1
    while arr[j] <= arr[i - 1]:
        j -= 1
    arr[i - 1], arr[j] = arr[j], arr[i - 1]

    # Reverse suffix
    arr[i:] = arr[len(arr) - 1: i - 1: -1]
    return True

def all_perm(N):
	"""Calculates permutations of arrays with length N containing 0, ..., N ones (otherwise zeros)"""
	#initialize array of needed size
    res = [None] * 2 ** N
    res = np.zeros((2 ** N, N), dtype=int)
    count = 0
    for i in range(N, -1, -1):
        A = np.zeros(N, dtype=int)
        for j in range(N - 1, i - 1, -1):
            A[j] = 1
        res[count] = np.copy(A)
        count = count + 1
        while True == next_permutation(A):
            res[count] = np.copy(A)
            count = count + 1
    return res

def special_perm(N, j):
	"""Calculates all permutations of array (filled with zeros and ones) with length N containing j ones"""
    if j > N:
        print('Error in permutations.py special_perm(N,j) : j larger than N')
        return 0
    # number of possible combinations
    num_comb = int(math.factorial(N) / (math.factorial(j) * math.factorial(N - j)))
    res = np.zeros((num_comb, N),dtype=int)
    arr = np.zeros(N)
    for i in range(N - 1, N - j - 1, -1):
        arr[i] = 1
    count = 0
    res[count] = np.copy(arr)
    count = count + 1
    while True == next_permutation(arr):
        res[count] = np.copy(arr)
        count = count + 1
    return res[::-1] #reversed array
