# MAC-Capacity-Region
Capacity region of a two-user vector Multiple Access Channel

A memoryless two-user Gaussian vector multiple access channel can be expressed as:

Y=H_1 x_1+H_2 x_2+n

where H_1 and H_2 are channel matrices, and n is Gaussian noise.

Define two parameters a_1 and a_2 that satisfy a_1+a_2=1, and both of them are nonnegative. Thus the optimization problem can be formulated as:

Max a_1 R_1+a_2 R_2, subject to the respective power constraint on each sub-channel.

Use CVX in MATLAB to obtain the optimal covariance matrices and then the data rates can be computed.
