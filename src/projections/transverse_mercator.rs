use std::f64::consts::PI;

use num::{Complex, Integer};

use crate::{latlon::LatLon, utility::{polyval, GeoMath, dms}, ThisOrThat, constants::{WGS84_A, WGS84_F, UTM_K0}};

// ================================
// Transverse Mercator Constants
// ================================

// Assuming GEOGRAPHICLIB_TRANSVERSEMERCATOR_ORDER == 6
const B1_COEFF: [f64; 5] = [
    // b1*(n+1), polynomial in n2 of order 3
    1.0, 4.0, 64.0, 256.0, 256.0,
];  // count = 5

#[allow(clippy::unreadable_literal)]
const ALP_COEFF: [f64; 27] = [
    // alp[1]/n^1, polynomial in n of order 5
    31564.0, -66675.0, 34440.0, 47250.0, -100800.0, 75600.0, 151200.0,
    // alp[2]/n^2, polynomial in n of order 4
    -1983433.0, 863232.0, 748608.0, -1161216.0, 524160.0, 1935360.0,
    // alp[3]/n^3, polynomial in n of order 3
    670412.0, 406647.0, -533952.0, 184464.0, 725760.0,
    // alp[4]/n^4, polynomial in n of order 2
    6601661.0, -7732800.0, 2230245.0, 7257600.0,
    // alp[5]/n^5, polynomial in n of order 1
    -13675556.0, 3438171.0, 7983360.0,
    // alp[6]/n^6, polynomial in n of order 0
    212378941.0, 319334400.0,
];  // count = 27

#[allow(clippy::unreadable_literal)]
const BET_COEFF: [f64; 27] = [
    // bet[1]/n^1, polynomial in n of order 5
    384796.0, -382725.0, -6720.0, 932400.0, -1612800.0, 1209600.0, 2419200.0,
    // bet[2]/n^2, polynomial in n of order 4
    -1118711.0, 1695744.0, -1174656.0, 258048.0, 80640.0, 3870720.0,
    // bet[3]/n^3, polynomial in n of order 3
    22276.0, -16929.0, -15984.0, 12852.0, 362880.0,
    // bet[4]/n^4, polynomial in n of order 2
    -830251.0, -158400.0, 197865.0, 7257600.0,
    // bet[5]/n^5, polynomial in n of order 1
    -435388.0, 453717.0, 15966720.0,
    // bet[6]/n^6, polynomial in n of order 0
    20648693.0, 638668800.0,
];  // count = 27

const MAXPOW: usize = 6;

const A: f64 = WGS84_A;
const F: f64 = WGS84_F;
const M: usize = MAXPOW / 2;
const N: f64 = F / (2.0 - F);
const E2: f64 = F * (2.0 - F);

pub(crate) struct TransverseMercator {
    k0: f64,
    es: f64,
    a1: f64,
    alp: [f64; MAXPOW + 1],
    bet: [f64; MAXPOW + 1],
}

impl TransverseMercator {
    pub fn utm() -> TransverseMercator {

        let es = (F < 0.0).ternary(-1.0, 1.0) * E2.abs().sqrt();

        let b1 = polyval(&B1_COEFF[0..=M], N.powi(2)) / (B1_COEFF[M + 1] * (1.0 + N));
        // a1 is the equivalent radius for computing the circumference of
        // ellipse.
        let a1 = b1 * A;

        let mut alp = [0_f64; MAXPOW + 1];
        let mut bet = [0_f64; MAXPOW + 1];

        let mut o = 0;
        let mut d = N;
        let mut m;

        for l in 1..=MAXPOW {
            m = MAXPOW - l;
            alp[l] = d * polyval(&ALP_COEFF[o..=o+m], N) / ALP_COEFF[o + m + 1];
            bet[l] = d * polyval(&BET_COEFF[o..=o+m], N) / BET_COEFF[o + m + 1];
            o += m + 2;
            d *= N;
        }

        Self {
            k0: UTM_K0,
            es,
            a1,
            alp,
            bet,
        }
    }

    #[allow(clippy::wrong_self_convention)]
    #[allow(clippy::similar_names)]
    pub fn from_latlon(&self, lon0: f64, lat: f64, lon: f64) -> (f64, f64) {
        let mut lat = lat;
        let mut lon = lon0.ang_diff(lon);

        let mut lat_sign = lat.is_sign_negative().ternary(-1.0, 1.0);
        let lon_sign = lon.is_sign_negative().ternary(-1.0, 1.0);

        lat *= lat_sign;
        lon *= lon_sign;

        let backside = lon > f64::from(dms::QD);

        if backside {
            if lat.is_zero() {
                lat_sign = -1.;
            }

            lon = f64::from(dms::HD) - lon;
        }

        let (phi_sin, phi_cos) = lat.to_radians().sin_cos();
        let (lamda_sin, lambda_cos) = lon.to_radians().sin_cos();

        // Check if lat == QD
        let (etap, xip) = if lat.eps_eq(f64::from(dms::QD)) {
            (0.0, PI / 2.0)
        } else {
            let tau = phi_sin / phi_cos;
            let taup = tau.taupf(self.es);
            let xip = taup.atan2(lambda_cos);
            let etap = (lamda_sin / taup.hypot(lambda_cos)).asinh();

            (etap, xip)
        };

        let c0 = (2.0 * xip).cos();
        let ch0 = (2.0 * etap).cosh();
        let s0 = (2.0 * xip).sin();
        let sh0 = (2.0 * etap).sinh();

        let mut a = Complex::new(2.0 * c0 * ch0, -2.0 * s0 * sh0);
        let mut n = MAXPOW;

        let mut y0 = Complex::new(n.is_odd().ternary_lazy( ||self.alp[n], || 0.0), 0.0);
        let mut y1 = Complex::default();

        if n.is_odd() {
            n -= 1;
        }

        while n > 0 {
            y1 = a * y0 - y1 + self.alp[n];
            n -= 1;
            y0 = a * y1 - y0 + self.alp[n];
            n -= 1;
        }

        a /= 2.0;
        a = Complex::new(s0 * ch0, c0 * sh0);
        y1 = Complex::new(xip, etap) + a * y0;

        let xi = y1.re;
        let eta = y1.im;
        let y = self.a1 * self.k0 * backside.ternary_lazy(|| PI - xi, || xi) * lat_sign;
        let x = self.a1 * self.k0 * eta * lon_sign;

        (x, y)
    }

    #[allow(clippy::many_single_char_names)]
    #[allow(clippy::similar_names)]
    pub fn to_latlon(&self, lon_input: f64, x: f64, y: f64) -> LatLon {
        let mut xi = y / (self.a1 * self.k0);
        let mut eta = x / (self.a1 * self.k0);

        let xi_sign = xi.is_sign_negative().ternary(-1.0, 1.0);
        let eta_sign = eta.is_sign_negative().ternary(-1.0, 1.0);

        xi *= xi_sign;
        eta *= eta_sign;

        let backside = xi > PI/2.;
        if backside {
            xi = PI - xi;
        }

        let c0 = (2.0 * xi).cos();
        let ch0 = (2.0 * eta).cosh();
        let s0 = (2.0 * xi).sin();
        let sh0 = (2.0 * eta).sinh();
        
        let mut a = Complex::new(2.0 * c0 * ch0, -2.0 * s0 * sh0);
        let mut n = MAXPOW;

        let mut y0 = Complex::new(n.is_odd().ternary(-self.bet[n], 0.0), 0.0);
        let mut y1 = Complex::default();

        if n.is_odd() {
            n -= 1;
        }

        while n > 0 {
            y1 = a * y0 - y1 - self.bet[n];
            n -= 1;

            y0 = a * y1 - y0 - self.bet[n];
            n -= 1;
        }

        a /= 2.;
        a = Complex::new(s0 * ch0, c0 * sh0);
        y1 = Complex::new(xi, eta) + a * y0;
        // Ignoring k and gamma

        let xip = y1.re;
        let etap = y1.im;
        let s = etap.sinh();
        let c = 0_f64.max(xip.cos());
        let r = s.hypot(c);


        let (mut lat, mut lon) = if r.is_zero() {
            (f64::from(dms::QD), 0.0)
        } else {
            let lon = s.atan2(c).to_degrees();
            let sxip = xip.sin();
            let tau = (sxip / r).tauf(self.es);


            let lat = tau.atan().to_degrees();

            (lat, lon)
        };

        lat *= xi_sign;
        if backside {
            lon = f64::from(dms::HD) - lon;
        }
        lon *= eta_sign;
        lon = (lon + lon_input).ang_normalize();
        
        LatLon { latitude: lat, longitude: lon }
    }
}
