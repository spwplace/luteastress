pub const PARAMS_SIZE: usize = 601;

pub fn get_default_params() -> Vec<f64> {
    let mut params = vec![0.0; PARAMS_SIZE];
    params[1] = 1.8274790e+03;
    params[2] = 7.3099159e+03;
    params[3] = 1.9220410e+02;
    params[4] = 2.3708000e+00;
    params[5] = 1.9041525e-01;
    params[6] = 10.0;
    params[7] = 1.0;
    params[8] = 0.001;
    params[9] = 5.0;
    params[10] = 100.0;
    params[11] = 10.0;
    params[12] = 5.0;
    params[13] = 3.0;
    params[14] = 7.6779e+06;
    params[16] = 4.7602470e-03;
    params[17] = 5.6989440e-02;
    params[18] = 0.001;
    params[20] = 3.0;
    params[21] = 2.2129050e+04;
    params[28] = 1.2083408e+00;

    params[230] = 7.3042153e+00;
    params[231] = 5.6501581e-01;
    params[240] = 3.1806080e-02;
    params[241] = 8.9730072e-03;
    params[242] = 1.0213583e-03;
    params[243] = 1.3649544e-04;
    params[244] = 1.2408453e-04;
    params[245] = 0.0;
    params[246] = 1.0;

    params[300] = 4.4746733e-01;
    params[301] = 5.5933417e-03;
    params[302] = 3.2217648e+02;
    params[303] = 6.4435296e+02;
    params[304] = 8.9493467e-03;
    params[305] = 3.2217648e+01;
    params[306] = 3.2217648e+00;
    params[307] = 3.2217648e+01;
    params[308] = 8.9493467e-02;
    params[309] = 3.2217648e+01;
    params[310] = 3.2217648e+00;
    params[311] = 8.9493467e-05;

    params[203] = 9.28;
    params[204] = 0.972;
    params[205] = 1.7137104e+03;
    params[206] = 9.4263469e-01;
    params[208] = 1.44523;
    params[209] = 2.2849472;
    params[210] = 2.8211026e+01;
    params[211] = 1.9406658e+02;

    params
}

pub const INITIAL_HORMONES: [f64; 54] = [
    13.018, 0.0, 4.8379e-11, 3.0044e+05, 3.0986, 78594.0, 7.9648, 3.0549, 21.213,
    0.42422, 1.4e-05, 0.22367, 0.0, 0.0, 2e-06, 2.6e-05, 0.000351, 0.003079,
    0.012677, 42.986, 2.2208, 0.99054, 55.526, 10.336, 0.21933, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.2753, 0.29968, 0.79694, 2.1793, 6.7756, 0.52814,
    1.1968, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.037925, 0.0082586, 0.00096153,
    0.00014978, 0.00013581, 0.0,
];

pub struct ParaOde {
    pub v: f64,
    pub gamma: f64,
    pub xi: f64,
    pub mu: f64,
    pub k: f64,
    pub rho: f64,
    pub folmax: f64,
}

impl Default for ParaOde {
    fn default() -> Self {
        Self {
            v: 2.0,
            gamma: 0.035 / 2.0,
            xi: 25.0,
            mu: 1.0,
            k: 0.065 / (25.0 * 25.0),
            rho: 0.01,
            folmax: 22.0,
        }
    }
}
