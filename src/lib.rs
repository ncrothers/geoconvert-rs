//! # geoconvert
//! 
//! `geoconvert` is a lightweight library for converting between different
//! geographic coordinate systems. Currently, there are three coordinate systems implemented:
//! 
//! * [`LatLon`]
//! * [`UtmUps`]
//! * [`Mgrs`]
//! 
//! The implementation of this library is a translation of a subset of 
//! [GeographicLib](https://geographiclib.sourceforge.io/C++/doc/index.html) from C++ to Rust. Specifically, `geoconvert`
//! implements some of the functionality of the [GeoConvert](https://geographiclib.sourceforge.io/C++/doc/GeoConvert.1.html) 
//! command line tool.
//! 
//! ## Usage
//! 
//! You can create coordinates manually using a struct's `create()` function, then convert to other
//! types using the `to_*`/`from_*` methods.
//! 
//! ```rust
//! use geoconvert::{LatLon, Mgrs, UtmUps};
//! 
//! // This returns a result. When calling `create()`, the arguments are validated to ensure only a valid
//! // coordinate gets created.
//! let coord = LatLon::create(40.748333, -73.985278).unwrap();
//! 
//! // Convert to a UTM/UPS coordinate
//! let coord_utm = coord.to_utmups();
//! let coord_utm = UtmUps::from_latlon(&coord);
//! 
//! // Convert to an MGRS coordinate
//! // Note that for MGRS you must specify the precision
//! let coord_mgrs = coord.to_mgrs(6);
//! let coord_mgrs = Mgrs::from_latlon(&coord, 6);
//! ```
//! 
//! ## Features
//! 
//! If you want `serde` compatibility with `Serialize`/`Deserialize`, activate the `serde` feature.

#![warn(clippy::pedantic)]
#![allow(
    // Don't require must_use
    clippy::must_use_candidate,
    clippy::return_self_not_must_use,
    // Lots of conversion between integer/float types in this library,
    // these suppress most of the unnecessary ones
    clippy::cast_precision_loss,
    clippy::cast_possible_wrap,
    clippy::cast_possible_truncation
)]

use thiserror::Error;

mod coords {
    pub mod latlon;
    pub mod mgrs;
    pub mod utm;
}

pub use coords::*;

pub(crate) mod utility;

pub use latlon::LatLon;
pub use mgrs::Mgrs;
pub use utm::UtmUps;

pub(crate) mod projections {
    pub mod transverse_mercator;
    pub mod polar_stereographic;
}

pub(crate) mod constants;

#[derive(Debug, Error)]
pub enum Error {
    #[error("The provided precision is outside of range [1, 11]")]
    InvalidPrecision(i32),
    #[error("The provided zone is outside the valid range [0, 60]")]
    InvalidZone(i32),
    #[error("Coordinate parameters are not valid: {0}")]
    InvalidCoord(String),
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