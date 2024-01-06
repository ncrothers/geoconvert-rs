use std::f64::consts::PI;

use num::{Complex, Zero};

use crate::{latlon::LatLon, utility::{polyval, GeoMath, tauf, dms}, ThisOrThat, constants::{WGS84_A, WGS84_F, UTM_K0}};

// ================================
// Transverse Mercator Constants
// ================================

// Assuming GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER == 6
const B1_COEFF: [f64; 5] = [
    // b1*(n+1), polynomial in n2 of order 3
    1., 4., 64., 256., 256.,
];  // count = 5

#[allow(clippy::unreadable_literal)]
const ALP_COEFF: [f64; 27] = [
    // alp[1]/n^1, polynomial in n of order 5
    31564., -66675., 34440., 47250., -100800., 75600., 151200.,
    // alp[2]/n^2, polynomial in n of order 4
    -1983433., 863232., 748608., -1161216., 524160., 1935360.,
    // alp[3]/n^3, polynomial in n of order 3
    670412., 406647., -533952., 184464., 725760.,
    // alp[4]/n^4, polynomial in n of order 2
    6601661., -7732800., 2230245., 7257600.,
    // alp[5]/n^5, polynomial in n of order 1
    -13675556., 3438171., 7983360.,
    // alp[6]/n^6, polynomial in n of order 0
    212378941., 319334400.,
];  // count = 27

#[allow(clippy::unreadable_literal)]
const BET_COEFF: [f64; 27] = [
    // bet[1]/n^1, polynomial in n of order 5
    384796., -382725., -6720., 932400., -1612800., 1209600., 2419200.,
    // bet[2]/n^2, polynomial in n of order 4
    -1118711., 1695744., -1174656., 258048., 80640., 3870720.,
    // bet[3]/n^3, polynomial in n of order 3
    22276., -16929., -15984., 12852., 362880.,
    // bet[4]/n^4, polynomial in n of order 2
    -830251., -158400., 197865., 7257600.,
    // bet[5]/n^5, polynomial in n of order 1
    -435388., 453717., 15966720.,
    // bet[6]/n^6, polynomial in n of order 0
    20648693., 638668800.,
];  // count = 27

const MAXPOW: usize = 6;

const A: f64 = WGS84_A;
const F: f64 = WGS84_F;
const M: usize = MAXPOW / 2;
const N: f64 = F / (2. - F);
const E2: f64 = F * (2. - F);
const E2M: f64 = 1. - E2;

pub(crate) struct TransverseMercator {
    a: f64,
    k0: f64,
    es: f64,
    c: f64,
    a1: f64,
    b1: f64,
    alp: [f64; MAXPOW + 1],
    bet: [f64; MAXPOW + 1],
}

impl TransverseMercator {
    pub fn utm() -> TransverseMercator {

        let es = (F < 0.).ternary(-1., 1.) * E2.abs().sqrt();
        let c = E2M.sqrt() * 1_f64.eatanhe(es).exp();

        let b1 = polyval(&B1_COEFF[0..=M], N.powi(2)) / (B1_COEFF[M + 1] * (1. + N));
        // a1 is the equivalent radius for computing the circumference of
        // ellipse.
        let a1 = b1 * A;

        let mut alp = [0_f64; MAXPOW + 1];
        let mut bet = [0_f64; MAXPOW + 1];

        let mut o = 0;
        let mut d = N;
        let mut m = MAXPOW / 2;

        for l in 1..=MAXPOW {
            m = MAXPOW - l;
            alp[l] = d * polyval(&ALP_COEFF[o..=o+m], N) / ALP_COEFF[o + m + 1];
            bet[l] = d * polyval(&BET_COEFF[o..=o+m], N) / BET_COEFF[o + m + 1];
            o += m + 2;
            d *= N;
        }

        Self {
            a: WGS84_A,
            k0: UTM_K0,
            es,
            c,
            a1,
            b1,
            alp,
            bet,
        }
    }

    pub fn utm_to_latlon(&self, lon_input: f64, x: f64, y: f64) -> LatLon {
        let mut xi = y / (self.a1 * self.k0);
        let mut eta = x / (self.a1 * self.k0);

        let xi_sign = (!xi.is_sign_positive()).ternary(-1., 1.);
        let eta_sign = (!eta.is_sign_positive()).ternary(-1., 1.);

        xi *= xi_sign;
        eta *= eta_sign;

        let backside = xi > PI/2.;
        if backside {
            xi = PI - xi;
        }

        let c0 = (2. * xi).cos();
        let ch0 = (2. * eta).cosh();
        let s0 = (2. * xi).sin();
        let sh0 = (2. * eta).sinh();
        
        let mut a = Complex::new(2. * c0 * ch0, -2. * s0 * sh0);
        let mut n = MAXPOW;

        let mut y0 = Complex::new((n % 2 == 1).ternary(-self.bet[n], 0.), 0.);
        let mut y1 = Complex::default();
        let mut z0 = Complex::new((n % 2 == 1).ternary(-2. * n as f64 * self.bet[n], 0.), 0.);
        let mut z1 = Complex::default();

        if n % 2 == 1 {
            n -= 1;
        }

        while n > 0 {
            y1 = a * y0 - y1 - self.bet[n];
            z1 = a * z0 - z1 - 2.*(n as f64) * self.bet[n];
            n -= 1;

            y0 = a * y1 - y0 - self.bet[n];
            z0 = a * z1 - z0 - 2.*(n as f64) * self.bet[n];
            n -= 1;
        }

        a /= 2.;
        z1 = 1. - z1 + a * z0;
        a = Complex::new(s0 * ch0, c0 * sh0);
        y1 = Complex::new(xi, eta) + a * y0;
        // Ignoring k and gamma

        let xip = y1.re;
        let etap = y1.im;
        let s = etap.sinh();
        let c = 0_f64.max(xip.cos());
        let r = s.hypot(c);


        let (mut lat, mut lon) = if r.is_zero() {
            (dms::QD as f64, 0.)
        } else {
            let lon = s.atan2(c).to_degrees();
            let sxip = xip.sin();
            let tau = tauf(sxip / r, self.es);


            let lat = tau.atan().to_degrees();

            (lat, lon)
        };

        lat *= xi_sign;
        if backside {
            lon = dms::HD as f64 - lon;
        }
        lon *= eta_sign;
        lon = (lon + lon_input).ang_normalize();
        
        LatLon { latitude: lat, longitude: lon }
    }
}
