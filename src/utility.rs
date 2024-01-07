use crate::ThisOrThat;

#[allow(dead_code)]
pub(crate) mod dms {
    /// Degrees per quarter turn
    pub const QD: i32 = 90;
    /// Minutes per degree
    pub const DM: i32 = 60;
    /// Seconds per minute
    pub const MS: i32 = 60;
    /// Degrees per half turn
    pub const HD: i32 = 2 * QD;
    /// Degrees per turn
    pub const TD: i32 = 2 * HD;
    /// Seconds per degree
    pub const DS: i32 = DM * MS;
}

fn special_sum(u: f64, v: f64) -> (f64, f64) {
    let s = u + v;
    let up = s - v;
    let vpp = s - up;

    let up = up - u;
    let vpp = vpp - v;

    let t = s.is_zero().ternary_lazy(||s, || -(up + vpp));

    (s, t)
}

/// Evaluate a polynomial
pub(crate) fn polyval(p: &[f64], x: f64) -> f64 {
    p
        .iter()
        .fold(0_f64, |acc, val| acc*x + val)
}

pub(crate) trait GeoMath {
    fn is_zero(&self) -> bool;
    fn eps_eq(&self, other: Self) -> bool;
    fn ang_normalize(&self) -> Self;
    fn ang_diff(&self, other: Self) -> Self;
    fn eatanhe(&self, es: Self) -> Self;
    fn remainder(&self, denom: Self) -> Self;
    fn taupf(&self, es: Self) -> Self;
    fn tauf(&self, es: Self) -> Self;
}

impl GeoMath for f64 {
    fn is_zero(&self) -> bool {
        self.abs() < f64::EPSILON
    }

    fn eps_eq(&self, other: f64) -> bool {
        (*self - other).abs() < f64::EPSILON
    }

    fn ang_normalize(&self) -> f64 {
        let value = self.remainder(f64::from(dms::TD));
        let hd = f64::from(dms::HD);

        if value.abs().eps_eq(hd) {
            hd.copysign(*self)
        }
        else {
            value
        }
    }

    fn ang_diff(&self, other: f64) -> f64 {
        let td = f64::from(dms::TD);
        // Use remainder instead of AngNormalize, since we treat boundary cases
        // later taking account of the error
        let (diff, err) = special_sum((-*self).remainder(td), other % td);
        // This second sum can only change d if abs(d) < 128, so don't need to
        // apply remainder yet again.
        let (diff, err) = special_sum(diff.remainder(td), err);
        
        let hd = f64::from(dms::HD);
        // Fix the sign if d = -180, 0, 180.
        if diff.is_zero() || diff.abs().eps_eq(hd) {
            // If e == 0, take sign from y - x
            // else (e != 0, implies d = +/-180), d and e must have opposite signs
            let sign = if err.is_zero() { other - *self } else { -err };
            diff.copysign(sign)
        }
        else {
            diff
        }
    }

    fn eatanhe(&self, es: f64) -> f64 {
        if es.is_sign_positive() {
            es * (es * *self).atanh()
        } else {
            -es * (es * *self).atanh()
        }
    }

    fn remainder(&self, denom: Self) -> Self {
        *self - (*self / denom).round() * denom
    }

    fn taupf(&self, es: f64) -> f64 {
        let tau1 = 1.0_f64.hypot(*self);
        let sig = (*self / tau1).eatanhe(es).sinh();

        1.0_f64.hypot(sig) * *self - sig * tau1
    }

    #[allow(clippy::similar_names)]
    fn tauf(&self, es: f64) -> f64 {
        let numit = 5;
        let tol = f64::EPSILON.sqrt() / 10.0;

        let e2m = 1.0 - es.powi(2);
        let mut tau = if self.abs() > 70.0 {
            self * 1_f64.eatanhe(es).exp()
        } else {
            self / e2m
        };

        let stol = tol * self.abs().max(1.0);
        for _ in 0..numit {
            let taupa = tau.taupf(es);
            let dtau = (self - taupa) * (1.0 + e2m * tau.powi(2))
                / (e2m * 1.0_f64.hypot(tau) * 1.0_f64.hypot(taupa));
            tau += dtau;
            if dtau.abs() < stol {
                break;
            }
        }
        tau
    }
}