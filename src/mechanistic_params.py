import numpy as np

# Parameters from Fischer-Holzhausen & Röblitz (2022) / GynCycle repository
# Indices match Par(i) in MATLAB (1-indexed)
PARAMS = np.zeros(601)
PARAMS[1] = 1.8274790e+03
PARAMS[2] = 7.3099159e+03
PARAMS[3] = 1.9220410e+02
PARAMS[4] = 2.3708000e+00
PARAMS[5] = 1.9041525e-01
PARAMS[6] = 10.0
PARAMS[7] = 1.0
PARAMS[8] = 0.001
PARAMS[9] = 5.0
PARAMS[10] = 100.0
PARAMS[11] = 10.0
PARAMS[12] = 5.0
PARAMS[13] = 3.0
PARAMS[14] = 7.6779e+06
PARAMS[16] = 4.7602470e-03
PARAMS[17] = 5.6989440e-02
PARAMS[18] = 0.001
PARAMS[20] = 3.0
PARAMS[21] = 2.2129050e+04
PARAMS[28] = 1.2083408e+00

# Hormone clearance and binding rates
PARAMS[230] = 7.3042153e+00 # LH binding
PARAMS[231] = 5.6501581e-01 # LH clearance
PARAMS[240] = 3.1806080e-02 # FSH binding
PARAMS[241] = 8.9730072e-03 # FSH clearance
PARAMS[242] = 1.0213583e-03 # FSHR desensitization
PARAMS[243] = 1.3649544e-04 # FSH transport to ovary
PARAMS[244] = 1.2408453e-04 # FSHR activation
PARAMS[245] = 0.0           # FSH clearance in ovary
PARAMS[246] = 1.0           # FSH transport scale

# GnRH receptor binding parameters
PARAMS[300] = 4.4746733e-01 # GnRH degradation
PARAMS[301] = 5.5933417e-03 # GnRH mass/freq scale
PARAMS[302] = 3.2217648e+02 # GnRH-Rec binding
PARAMS[303] = 6.4435296e+02 # GnRH-Rec dissociation
PARAMS[304] = 8.9493467e-03 # GReci degradation
PARAMS[305] = 3.2217648e+01 # GReci -> RecGi
PARAMS[306] = 3.2217648e+00 # RecGa -> RecGi
PARAMS[307] = 3.2217648e+01 # RecGi -> RecGa
PARAMS[308] = 8.9493467e-02 # RecGi degradation
PARAMS[309] = 3.2217648e+01 # GReca -> GReci
PARAMS[310] = 3.2217648e+00 # GReci -> GReca
PARAMS[311] = 8.9493467e-05 # GnRH receptor synthesis

# Feedback parameters
PARAMS[203] = 9.28           # P4 threshold for GnRH freq
PARAMS[204] = 0.972          # P4 Hill exponent
PARAMS[205] = 1.7137104e+03  # E2 threshold for GnRH freq
PARAMS[206] = 9.4263469e-01  # E2 Hill exponent
PARAMS[208] = 1.44523        # E2 threshold for GnRH mass
PARAMS[209] = 2.2849472      # E2 Hill exponent
PARAMS[210] = 2.8211026e+01  # E2 threshold for GnRH mass surge
PARAMS[211] = 1.9406658e+02  # E2 Hill exponent surge

# Initial Values for the 54 hormone-related variables
INITIAL_HORMONES = np.array([
    13.018, 0, 4.8379e-11, 3.0044e+05, 3.0986, 78594, 7.9648, 3.0549, 21.213, 
    0.42422, 1.4e-05, 0.22367, 0, 0, 2e-06, 2.6e-05, 0.000351, 0.003079, 
    0.012677, 42.986, 2.2208, 0.99054, 55.526, 10.336, 0.21933, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 8.2753, 0.29968, 0.79694, 2.1793, 6.7756, 0.52814, 
    1.1968, 0, 0, 0, 0, 0, 0, 0, 0.037925, 0.0082586, 0.00096153, 
    0.00014978, 0.00013581, 0.0
])

# paraOde mapping
# v; gamma; xi; mu; k; rho; Folmax
PARA_ODE = {
    'v': 2,
    'gamma': 0.035 / 2,
    'xi': 25,
    'mu': 1,
    'k': 0.065 / (25**2),
    'rho': 0.01,
    'folmax': 22
}
