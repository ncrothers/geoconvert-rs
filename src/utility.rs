use std::{ops::{Add, Sub, Neg, RemAssign}, f64::EPSILON};

use num::{FromPrimitive, Float};

use crate::ThisOrThat;

pub mod dms {
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

fn special_sum<T>(u: T, v: T) -> (T, T)
where
    T: Float + FromPrimitive
{
    let s = u + v;
    let up = s - v;
    let vpp = s - up;

    let up = up - u;
    let vpp = vpp - v;

    let t = if !s.is_zero() { -(up + vpp) } else { s };

    (s, t)
}

/// Evaluate a polynomial
pub(crate) fn polyval(p: &[f64], x: f64) -> f64 {
    p
        .iter()
        .fold(0_f64, |acc, val| acc*x + val)
}

///
/// Hypot of a given tau
///
/// # Arguments
///
/// * `tau: f64` - In radians
/// * `es: f64` - In radians
///
pub(crate) fn taupf(tau: f64, es: f64) -> f64 {
    let tau1: f64 = 1.0_f64.hypot(tau);
    let sig = (tau / tau1).eatanhe(es).sinh();

    1.0_f64.hypot(sig) * tau - sig * tau1
}

///
/// Tau float comparison
///
/// # Arguments
///
/// * `tau: f64` - In radians
/// * `es: f64` - In radians
///
pub(crate) fn tauf(taup: f64, es: f64) -> f64 {
    let numit = 5;
    let tol: f64 = EPSILON.sqrt() / 10.0;
    let taumax = 2. / EPSILON.sqrt();

    let e2m: f64 = 1.0 - es.powi(2);
    let mut tau: f64 = if taup.abs() > 70. {
        taup * 1_f64.eatanhe(es).exp()
    } else {
        taup / e2m
    };

    let stol: f64 = tol * taup.abs().max(1.0);
    for _ in 0..numit {
        let taupa: f64 = taupf(tau, es);
        let dtau: f64 = (taup - taupa) * (1.0 + e2m * tau.powi(2))
            / (e2m * 1.0_f64.hypot(tau) * 1.0_f64.hypot(taupa));
        tau += dtau;
        if dtau.abs() < stol {
            break;
        }
    }
    tau
}

pub(crate) trait GeoMath {
    fn ang_normalize(&self) -> Self;
    fn ang_diff(&self, other: Self) -> Self;
    fn eatanhe(&self, es: Self) -> Self;
    fn remainder(&self, denom: Self) -> Self;
}

impl<T> GeoMath for T where T: Float + FromPrimitive {
    fn ang_normalize(&self) -> T {
        let value = self.remainder(FromPrimitive::from_i32(dms::TD).expect("dms::TD must be convertible"));
        let hd = FromPrimitive::from_i32(dms::HD).expect("dms::HD must be convertible");

        if value.abs() == hd {
            hd.copysign(*self)
        }
        else {
            value
        }
    }

    fn ang_diff(&self, other: T) -> T {
        let td = FromPrimitive::from_i32(dms::TD).expect("dms::TD must be convertible");
        // Use remainder instead of AngNormalize, since we treat boundary cases
        // later taking account of the error
        let (diff, err) = special_sum((-*self).remainder(td), other % td);
        // This second sum can only change d if abs(d) < 128, so don't need to
        // apply remainder yet again.
        let (diff, err) = special_sum(diff.remainder(td), err);
        
        let hd = FromPrimitive::from_i32(dms::HD).expect("dms::HD must be convertible");
        // Fix the sign if d = -180, 0, 180.
        if diff.is_zero() || diff.abs() == hd {
            // If e == 0, take sign from y - x
            // else (e != 0, implies d = +/-180), d and e must have opposite signs
            let sign = if err.is_zero() { other - *self } else { -err };
            diff.copysign(sign)
        }
        else {
            diff
        }
    }

    fn eatanhe(&self, es: T) -> T {
        if es.is_sign_positive() {
            es * (es * *self).atanh()
        } else {
            -es * (es * *self).atanh()
        }
    }

    fn remainder(&self, denom: Self) -> Self {
        *self - (*self / denom).round() * denom
    }
}