import numpy as np
import math

# consec512.dat is an array of 32-bit unsigned integers (little-endian) with as dimensions
# 512*256*256. This is a C-style array, so contiguous with row-major order.
#
# The first dimension is the position of the first byte of the pair/digraph
# The second dimension is the value of the first byte in the pair/digraph
# The third dimension is the value of the second byte in the pair/digraph
# -> The value at these indices represents the number of times this pair/digraph occurred.
stats = np.fromfile("consec512.dat", dtype="uint32").reshape((512,256,256))

# Number of RC4 keys used
numkeys = float(stats[0,:,:].sum())
print "Number of RC4 keys:", math.log(numkeys, 2)

# Example 1:
#
#       Probability Pr[Z_100 = 75 /\ Z_101 = 200] assuming we loaded byte1.dat.
#
# Note that in papers the position r in Z_r is 1-based, while python is 0-based.
pr = stats[100-1, 75, 200] / numkeys
print "Probability example event: 2^%f = %g" % (math.log(pr, 2), pr)

# Example 2:
#
#	Probability of Fluhrer-McGrew pair (255,255) in the initial keystream bytes
#
# We print the absolute relative bias compared to single-byte biases. Otherwise
# these results are too heavily influenced by single-byte biases.
single = stats.sum(axis=2)
for i in range(1,512):
    s = stats[i-1, 255, 255] / numkeys
    p = (single[i-1, 255] / numkeys) * (single[i, 255] / numkeys)
    q = abs(1 - s / p)
    print "%3d -> 2^%7.3f" % (i, math.log(q, 2))
