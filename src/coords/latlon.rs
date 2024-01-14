use std::fmt::Display;

use crate::{Error, utm::UtmUps, mgrs::Mgrs};

/// Mean radius of Earth in meters
/// 
/// <https://en.wikipedia.org/wiki/Earth_radius#Arithmetic_mean_radius>
const EARTH_MEAN_RADIUS_M: f64 = 6371.0088 * 1000.0;

/// Representation of a WGS84 Latitude/Longitude point. Can be converted
/// to/from [`UtmUps`] and [`Mgrs`].
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LatLon {
    #[cfg_attr(feature = "serde", serde(alias = "lat"))]
    pub(crate) latitude: f64,
    #[cfg_attr(feature = "serde", serde(alias = "lon"))]
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

    /// Tries to create a latitude/longitude point from a lat/lon pair. First checks if the
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

    /// Returns the latitude value.
    /// 
    /// # Example
    /// ```
    /// use geoconvert::LatLon;
    /// 
    /// let coord = LatLon::create(40.748333, -73.985278).unwrap();
    /// assert_eq!(coord.latitude(), 40.748333);
    /// ```
    #[inline]
    pub fn latitude(&self) -> f64 {
        self.latitude
    }
    
    /// Returns the longitude value.
    /// 
    /// # Example
    /// ```
    /// use geoconvert::LatLon;
    /// 
    /// let coord = LatLon::create(40.748333, -73.985278).unwrap();
    /// assert_eq!(coord.longitude(), -73.985278);
    /// ```
    #[inline]
    pub fn longitude(&self) -> f64 {
        self.longitude
    }

    /// Returns whether the current point is in the northern hemisphere.
    /// 
    /// # Example
    /// 
    /// ```
    /// use geoconvert::LatLon;
    /// 
    /// let coord = LatLon::create(40.748333, -73.985278).unwrap();
    /// assert!(coord.is_north());
    /// 
    /// let coord = LatLon::create(-40.748333, -73.985278).unwrap();
    /// assert!(!coord.is_north());
    /// ```
    pub fn is_north(&self) -> bool {
        self.latitude.is_sign_positive()
    }

    /// Returns the distance in meters between two [`LatLon`] points
    /// using the [haversine formula](https://en.wikipedia.org/wiki/Haversine_formula).
    /// Uses the [mean radius of the Earth](https://en.wikipedia.org/wiki/Earth_radius#Arithmetic_mean_radius)
    /// in the calculation: `6371.0088`
    pub fn haversine(&self, other: &LatLon) -> f64 {
        let lat1_r = self.latitude.to_radians();
        let lat2_r = other.latitude.to_radians();
        
        2.0 * EARTH_MEAN_RADIUS_M * (
            ((other.latitude - self.latitude).to_radians() / 2.0).sin().powi(2) + 
            lat1_r.cos() * lat2_r.cos() *
            ((other.longitude - self.longitude).to_radians() / 2.0).sin().powi(2)
        ).sqrt().asin()
    }
    
    /// Converts from [`UtmUps`] to [`LatLon`]
    /// 
    /// # Usage
    /// 
    /// ```
    /// use geoconvert::{LatLon, UtmUps};
    /// 
    /// let coord = LatLon::create(40.748333, -73.985278).unwrap();
    /// let coord_utm = UtmUps::create(18, true, 585664.121, 4511315.422).unwrap();
    /// 
    /// let converted = LatLon::from_utmups(&coord_utm);
    /// 
    /// // Check if the converted coordinate is accurate to 6 decimals (same as reference)
    /// assert!((converted.latitude() - coord.latitude()).abs() < 1e-6);
    /// assert!((converted.longitude() - coord.longitude()).abs() < 1e-6);
    /// ```
    pub fn from_utmups(value: &UtmUps) -> LatLon {
        value.to_latlon()
    }
    
    /// Converts from [`LatLon`] to [`UtmUps`]
    /// 
    /// # Usage
    /// 
    /// ```
    /// use geoconvert::{LatLon, UtmUps};
    /// 
    /// let coord = LatLon::create(40.748333, -73.985278).unwrap();
    /// let coord_utm = UtmUps::create(18, true, 585664.121, 4511315.422).unwrap();
    /// 
    /// let converted = coord.to_utmups();
    /// 
    /// assert_eq!(converted.zone(), coord_utm.zone());
    /// assert_eq!(converted.is_north(), coord_utm.is_north());
    /// // Check if the converted coordinate is accurate to 3 decimals (same as reference)
    /// assert!((converted.easting() - coord_utm.easting()).abs() < 1e-3);
    /// assert!((converted.northing() - coord_utm.northing()).abs() < 1e-3);
    /// ```
    pub fn to_utmups(&self) -> UtmUps {
        UtmUps::from_latlon(self)
    }

    /// Converts from [`Mgrs`] to [`LatLon`]
    /// 
    /// # Usage
    /// 
    /// ```
    /// use geoconvert::{LatLon, Mgrs};
    /// 
    /// let coord = LatLon::create(40.748333, -73.985278).unwrap();
    /// let coord_mgrs = Mgrs::parse_str("18TWL856641113154").unwrap();
    /// 
    /// let converted = LatLon::from_mgrs(&coord_mgrs);
    /// 
    /// // Check if the converted coordinate is accurate to 6 decimals (same as reference)
    /// assert!((converted.latitude() - coord.latitude()).abs() < 1e-6);
    /// assert!((converted.longitude() - coord.longitude()).abs() < 1e-6);
    /// ```
    pub fn from_mgrs(value: &Mgrs) -> LatLon {
        value.to_latlon()
    }

    /// Converts from [`LatLon`] to [`Mgrs`]
    /// 
    /// # Usage
    /// 
    /// ```
    /// use geoconvert::{LatLon, Mgrs};
    /// 
    /// let coord = LatLon::create(40.748333, -73.985278).unwrap();
    /// 
    /// let converted = coord.to_mgrs(6);
    /// 
    /// assert_eq!(converted.to_string(), "18TWL856641113154");
    /// ```
    pub fn to_mgrs(&self, precision: i32) -> Mgrs {
        Mgrs::from_latlon(self, precision)
    }
}

impl Display for LatLon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut buf = ryu::Buffer::new();
        let lat = buf.format(self.latitude);
        let mut buf = ryu::Buffer::new();
        let lon = buf.format(self.longitude);
        write!(
            f,
            "{lat} {lon}",
        )
    }
}
