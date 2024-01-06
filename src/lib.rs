#![warn(clippy::pedantic)]
#![allow(clippy::must_use_candidate)]

use thiserror::Error;

pub mod latlon;
pub mod mgrs;
pub mod utm;
pub mod utility;

pub(crate) mod projections {
    pub mod transverse_mercator;
    pub mod polar_stereographic;
}

pub(crate) mod constants;

#[derive(Debug, Error)]
pub enum Error {
    #[error("MGRS String is invalid: {0}")]
    InvalidMgrs(String),
    #[error("UTM coords are invalid: {0}")]
    InvalidUtmCoords(String),
    #[error("Coordinate type {coord_type} not valid for conversion to {dest_type}: {msg}")]
    InvalidRange {
        coord_type: String,
        dest_type: String,
        msg: String,
    },
}

pub trait ParseCoord {
    fn parse_coord(value: &str) -> Result<Self, Error>
    where Self: Sized;
}

pub fn from_str<S, T>(value: S) -> Result<T, Error>
where 
    S: AsRef<str>,
    T: ParseCoord
{
    T::parse_coord(value.as_ref())
}

trait ThisOrThat {
    fn ternary<T>(&self, r#true: T, r#false: T) -> T;
    fn ternary_lazy<F, E, T>(&self, r#true: F, r#false: E) -> T
    where
        F: Fn() -> T, 
        E: Fn() -> T;
}

impl ThisOrThat for bool {
    fn ternary<T>(&self, r#true: T, r#false: T) -> T {
        if *self { r#true } else { r#false }
    }

    fn ternary_lazy<F, E, T>(&self, r#true: F, r#false: E) -> T
    where
        F: Fn() -> T, 
        E: Fn() -> T, 
    {
        if *self { r#true() } else { r#false() }
    }
}