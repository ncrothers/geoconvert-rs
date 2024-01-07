use std::fmt::Display;

use crate::{Error, utm::{UtmUps, zonespec, self, central_meridian}, ThisOrThat, projections::{transverse_mercator::TransverseMercator, polar_stereographic::PolarStereographic}, mgrs::Mgrs, constants::WGS84_A};

#[derive(Clone, Copy, Debug)]
pub struct LatLon {
    pub(crate) latitude: f64,
    pub(crate) longitude: f64,
}

impl LatLon {
    /// Internal-only constructor that doesn't check the bounds of lat/lon
    pub(crate) fn new(lat: f64, lon: f64) -> LatLon {
        Self {
            latitude: lat,
            longitude: lon,
        }
    }

    pub fn latitude(&self) -> f64 {
        self.latitude
    }

    pub fn longitude(&self) -> f64 {
        self.longitude
    }

    /// Tries to create a `LatLon` point from a lat/lon pair.0 First checks if the
    /// values are valid:
    /// * Latitude must be in range [-90,90]
    /// * Longitude must be in range [-180,180]
    /// 
    /// # Errors
    /// 
    /// Returns [`Error::InvalidCoord`] if either latitude or longitude are invalid.
    /// 
    /// # Usage
    /// 
    /// ```
    /// use geoconvert::LatLon;
    /// 
    /// let coord = LatLon::create(40.748333, -73.985278);
    /// 
    /// assert!(coord.is_ok());
    /// 
    /// let coord = coord.unwrap();
    /// 
    /// assert_eq!(coord.latitude(), 40.748333);
    /// assert_eq!(coord.longitude(), -73.985278);
    /// 
    /// let invalid_coord_lat = LatLon::create(100.0, 0.0);
    /// assert!(invalid_coord_lat.is_err());
    /// 
    /// let invalid_coord_lon = LatLon::create(0.0, -200.0);
    /// assert!(invalid_coord_lon.is_err());
    /// ```
    pub fn create(lat: f64, lon: f64) -> Result<LatLon, Error> {
        if !(-90_f64..=90_f64).contains(&lat) {
            Err(Error::InvalidCoord(format!("Latitude {lat} outside of valid range [-90, 90].")))
        } else if !(-180_f64..180_f64).contains(&lon) {
            Err(Error::InvalidCoord(format!("Longitude {lon} outside of valid range [-180, 180].")))
        } else {
            Ok(LatLon::new(lat, lon))
        }
    }

    pub fn is_north(&self) -> bool {
        self.latitude.is_sign_positive()
    }

    pub fn haversine(&self, other: &LatLon) -> f64 {
        let lat1_r = self.latitude.to_radians();
        let lat2_r = other.latitude.to_radians();
        
        2.0 * 6_371_088.0 * (
            ((other.latitude - self.latitude).to_radians() / 2.0).sin().powi(2) + 
            lat1_r.cos() * lat2_r.cos() *
            ((other.longitude - self.longitude).to_radians() / 2.0).sin().powi(2)
        ).sqrt().asin()
    }
}

impl Display for LatLon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} {}",
            self.latitude,
            self.longitude,
        )
    }
}