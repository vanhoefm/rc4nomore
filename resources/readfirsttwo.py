import numpy as np
import math

# byte1.dat and byte2.dat are arrays of 64-bit unsigned integers (little-endian) with as
# dimensions 256*256*256. This is a C-style array, so contiguous with row-major order.
#
# The first dimension is the value of the byte (first or second keystream byte)
# The second dimension is the position of the other keystream byte
# The third dimension is the value of the other keystream byte.
# -> The value at these indices represents the number of times this pair occurred.
stats = np.fromfile("byte1.dat", dtype="uint64").reshape((256,256,256))

# Number of RC4 keys used
numkeys = float(stats[:,2,:].sum())
print "Number of RC4 keys:", math.log(numkeys, 2)

# Example 1 (byte1.dat):
#
#       probability Pr[Z_1 = 123 /\ Z_45 = 215] assuming we loaded byte1.dat.
#
# Note that in papers the position r in Z_r is 1-based, while python is 0-based.
pr = stats[123, 45-1, 215] / numkeys
print "Probability example event: 2^%f = %g" % (math.log(pr, 2), pr)

# Example 2 (byte1.dat):
#
#	Bias 2 in the paper: Pr[Z_1 = 257 - i /\ Z_i = i]
#
# We print the absolute relative bias compared to single-byte biases. Otherwise
# these results are too heavily influenced by single-byte biases.
single = stats.sum(axis=0)
for i in range(2,256):
    s = stats[257-i, i-1, i] / numkeys
    p = (single[0, 257-i] / numkeys) * (single[i-1, i] / numkeys)
    q = abs(1 - s / p)
    print "%3d -> 2^%7.3f" % (i, math.log(q, 2))