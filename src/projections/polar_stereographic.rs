use crate::{ThisOrThat, constants::{WGS84_A, WGS84_F, UPS_K0}, utility::{GeoMath, dms}, latlon::LatLon};

const F: f64 = WGS84_F;
const E2: f64 = F * (2.0 - F);

pub(crate) struct PolarStereographic {
    a: f64,
    k0: f64,
    es: f64,
    c: f64,
}

impl PolarStereographic {
    pub fn ups() -> PolarStereographic {

        let es = (F < 0.0).ternary(-1.0, 1.0) * E2.abs().sqrt();
        let c = (1.0 - F) * 1_f64.eatanhe(es).exp();

        Self {
            a: WGS84_A,
            k0: UPS_K0,
            es,
            c,
        }
    }

    #[allow(clippy::wrong_self_convention)]
    pub fn from_latlon(&self, northp: bool, lat: f64, lon: f64) -> (f64, f64) {
        let lat = lat * northp.ternary(1.0, -1.0);

        let tau = lat.to_radians().tan();
        let taup = tau.taupf(self.es);
        let mut rho = 1_f64.hypot(taup) + taup.abs();
        rho = (taup >= 0.0).ternary_lazy(|| (!lat.eps_eq(f64::from(dms::QD))).ternary_lazy(|| 1.0 / rho, || 0.0), || rho);
        rho *= 2.0 * self.k0 * self.a / self.c;

        let (mut x, mut y) = {
            let (x, y) = lon.to_radians().sin_cos();
            (x, y)
        };

        x *= rho;
        y *= northp.ternary(-rho, rho);

        (x, y)
    }

    pub fn to_latlon(&self, northp: bool, x: f64, y: f64) -> LatLon {
        let rho = x.hypot(y);
        let t = (rho != 0.0)
            .ternary_lazy(
                || rho / (2.0 * self.k0 * self.a / self.c),
                || f64::EPSILON.powi(2)
            );
        let taup = (1.0 / t - t) / 2.;
        let tau = taup.tauf(self.es);

        let lat = northp.ternary(1.0, -1.0) * tau.atan().to_degrees();
        let lon = x.atan2(northp.ternary(-y, y)).to_degrees();

        LatLon {
            latitude: lat,
            longitude: lon,
        }
    }
}
