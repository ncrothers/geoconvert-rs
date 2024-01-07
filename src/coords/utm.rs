use crate::{latlon::LatLon, utility::{dms, GeoMath}, mgrs::{to_latitude_band, self, Mgrs}, Error, ThisOrThat, projections::{transverse_mercator::TransverseMercator, polar_stereographic::PolarStereographic}};

pub(crate) mod zonespec {
    pub(crate) const INVALID: i32 = -4;
    pub(crate) const _MATCH: i32 = -3;
    pub(crate) const UTM: i32 = -2;
    pub(crate) const STANDARD: i32 = -1;
    pub(crate) const UPS: i32 = 0;
    pub(crate) const MINZONE: i32 = 0;
    pub(crate) const MINUTMZONE: i32 = 1;
    pub(crate) const MAXUTMZONE: i32 = 60;
    pub(crate) const MAXZONE: i32 = 60;
}

const FALSE_EASTING: [i32; 4] = [
    mgrs::UPSEASTING * mgrs::TILE,
    mgrs::UPSEASTING * mgrs::TILE,
    mgrs::UTMEASTING * mgrs::TILE,
    mgrs::UTMEASTING * mgrs::TILE,
];

const FALSE_NORTHING: [i32; 4] = [
    mgrs::UPSEASTING * mgrs::TILE,
    mgrs::UPSEASTING * mgrs::TILE,
    mgrs::MAXUTM_S_ROW * mgrs::TILE,
    mgrs::MINUTM_N_ROW * mgrs::TILE,
];

const MIN_EASTING: [i32; 4] = [
    mgrs::MINUPS_S_IND * mgrs::TILE,
    mgrs::MINUPS_N_IND * mgrs::TILE,
    mgrs::MINUTMCOL * mgrs::TILE,
    mgrs::MINUTMCOL * mgrs::TILE,
];

const MAX_EASTING: [i32; 4] = [
    mgrs::MAXUPS_S_IND * mgrs::TILE,
    mgrs::MAXUPS_N_IND * mgrs::TILE,
    mgrs::MAXUTMCOL * mgrs::TILE,
    mgrs::MAXUTMCOL * mgrs::TILE,
];

const MIN_NORTHING: [i32; 4] = [
    mgrs::MINUPS_S_IND * mgrs::TILE,
    mgrs::MINUPS_N_IND * mgrs::TILE,
    mgrs::MINUTM_S_ROW * mgrs::TILE,
    (mgrs::MINUTM_N_ROW + mgrs::MINUTM_S_ROW - mgrs::MAXUTM_S_ROW) * mgrs::TILE,
];

const MAX_NORTHING: [i32; 4] = [
    mgrs::MAXUPS_S_IND * mgrs::TILE,
    mgrs::MAXUPS_N_IND * mgrs::TILE,
    (mgrs::MAXUTM_S_ROW + mgrs::MAXUTM_N_ROW - mgrs::MINUTM_N_ROW) * mgrs::TILE,
    mgrs::MAXUTM_N_ROW * mgrs::TILE,
];

/// Representation of a WGS84 
/// [UTM](https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system)
/// /
/// [UPS](https://en.wikipedia.org/wiki/Universal_polar_stereographic_coordinate_system) 
/// point. If converted to from lat/lon, it will
/// automatically determine whether it should be UTM/UPS. It becomes a UPS coordinate
/// if the latitude is outside the range `[-84,84]`. A zone value of `0`
/// designates UPS.
#[allow(clippy::module_name_repetitions)]
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct UtmUps {
    pub(crate) zone: i32,
    #[cfg_attr(feature = "serde", serde(alias = "north", alias = "is_north"))]
    pub(crate) northp: bool,
    pub(crate) easting: f64,
    pub(crate) northing: f64,
}

impl UtmUps {
    /// Internal-only constructor that doesn't check the coordinate
    pub(crate) fn new(zone: i32, northp: bool, easting: f64, northing: f64) -> UtmUps {
        Self {
            zone,
            northp,
            easting,
            northing,
        }
    }

    /// Tries to create a UTM or UPS point from its constituent parts. Zone
    /// of `0` designates UPS, otherwise it is UTM.
    /// 
    /// # Errors
    /// 
    /// Returns [`Error::InvalidZone`] if the zone is outside the range `[0, 60]`.
    /// Returns [`Error::InvalidCoord`] if the coordinate is otherwise invalid.
    /// 
    /// # Usage
    /// 
    /// ```
    /// use geoconvert::UtmUps;
    /// 
    /// let coord = UtmUps::create(18, true, 585664.121, 4511315.422);
    /// 
    /// assert!(coord.is_ok());
    /// 
    /// let coord = coord.unwrap();
    /// 
    /// assert_eq!(coord.zone(), 18);
    /// assert_eq!(coord.is_north(), true);
    /// assert!((coord.easting() - 585664.121).abs() < 1e-3);
    /// assert!((coord.northing() - 4511315.422).abs() < 1e-3);
    /// 
    /// let invalid_coord_zone_neg = UtmUps::create(-10, true, 585664.121, 4511315.422);
    /// assert!(invalid_coord_zone_neg.is_err());
    /// 
    /// let invalid_coord_zone_too_big = UtmUps::create(70, true, 585664.121, 4511315.422);
    /// assert!(invalid_coord_zone_too_big.is_err());
    /// ```
    pub fn create(zone: i32, northp: bool, easting: f64, northing: f64) -> Result<UtmUps, Error> {
        // Make sure zone is a valid value
        if !(zonespec::MINZONE..=zonespec::MAXZONE).contains(&zone) {
            return Err(Error::InvalidZone(zone));
        }

        let utmp = zone != zonespec::UPS;

        check_coords(utmp, northp, easting, northing, false)?;

        Ok(UtmUps::new(zone, northp, easting, northing))
    }

    /// Returns the UTM zone.
    /// 
    /// # Example
    /// ```
    /// use geoconvert::UtmUps;
    /// 
    /// let coord = UtmUps::create(18, true, 585664.121, 4511315.422).unwrap();
    /// assert_eq!(coord.zone(), 18);
    /// ```
    pub fn zone(&self) -> i32 {
        self.zone
    }

    /// Returns whether the coordinate is in the northern hemisphere.
    /// 
    /// # Example
    /// ```
    /// use geoconvert::UtmUps;
    /// 
    /// let coord = UtmUps::create(18, true, 585664.121, 4511315.422).unwrap();
    /// assert_eq!(coord.is_north(), true);
    /// ```
    pub fn is_north(&self) -> bool {
        self.northp
    }

    /// Returns the UTM easting.
    /// 
    /// # Example
    /// ```
    /// use geoconvert::UtmUps;
    /// 
    /// let coord = UtmUps::create(18, true, 585664.121, 4511315.422).unwrap();
    /// assert!((coord.easting() - 585664.121).abs() < 1e-3);
    /// ```
    pub fn easting(&self) -> f64 {
        self.easting
    }

    /// Returns the UTM northing.
    /// 
    /// # Example
    /// ```
    /// use geoconvert::UtmUps;
    /// 
    /// let coord = UtmUps::create(18, true, 585664.121, 4511315.422).unwrap();
    /// assert!((coord.northing() - 4511315.422).abs() < 1e-3);
    /// ```
    pub fn northing(&self) -> f64 {
        self.northing
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
    pub fn from_latlon(value: &LatLon) -> UtmUps {
        let northp = value.is_north();
        // STANDARD specifies, by default, interpret whether it should be UTM or UPS
        // TODO: Maybe do something if the zone is invalid?
        let zone = standard_zone(value.latitude, value.longitude, zonespec::STANDARD);

        let utmp = zone != zonespec::UPS;
        let (mut x, mut y) = if utmp {
            let lon0 = central_meridian(zone);

            TransverseMercator::utm().from_latlon(lon0, value.latitude, value.longitude)
        } else {
            PolarStereographic::ups().from_latlon(northp, value.latitude, value.longitude)
        };

        let ind = utmp.ternary(2, 0) + northp.ternary(1, 0);
        x += f64::from(FALSE_EASTING[ind]);
        y += f64::from(FALSE_NORTHING[ind]);

        UtmUps {
            zone,
            northp,
            northing: y,
            easting: x,
        }
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
    pub fn to_latlon(&self) -> LatLon {
        let utmp = self.zone != zonespec::UPS;

        let ind = utmp.ternary(2, 0) + self.northp.ternary(1, 0);

        let x = self.easting - f64::from(FALSE_EASTING[ind]);
        let y = self.northing - f64::from(FALSE_NORTHING[ind]);

        if utmp {
            TransverseMercator::utm().to_latlon(central_meridian(self.zone), x, y)
        } else {
            PolarStereographic::ups().to_latlon(self.northp, x, y)
        }
    }

    /// Converts from [`Mgrs`] to [`UtmUps`]
    /// 
    /// # Usage
    /// 
    /// ```
    /// use geoconvert::{Mgrs, UtmUps};
    /// 
    /// let coord = Mgrs::parse_str("18TWL856641113154").unwrap();
    /// let coord_utm = UtmUps::create(18, true, 585664.15, 4511315.45).unwrap();
    /// 
    /// let converted = coord.to_utmups();
    /// 
    /// // Check if the converted coordinate is accurate to 6 decimals (same as reference)
    /// assert_eq!(coord_utm.zone(), converted.zone());
    /// assert_eq!(coord_utm.is_north(), converted.is_north());
    /// assert!((coord_utm.easting() - converted.easting()).abs() < 1e-2);
    /// assert!((coord_utm.northing() - converted.northing()).abs() < 1e-2);
    /// ```
    pub fn from_mgrs(value: &Mgrs) -> UtmUps {
        value.utm
    }

    /// Converts from [`UtmUps`] to [`Mgrs`]
    /// 
    /// # Usage
    /// 
    /// ```
    /// use geoconvert::{Mgrs, UtmUps};
    /// 
    /// let coord = Mgrs::parse_str("18TWL856641113154").unwrap();
    /// let coord_utm = UtmUps::create(18, true, 585664.15, 4511315.45).unwrap();
    /// 
    /// let converted = Mgrs::from_utmups(&coord_utm, 6);
    /// 
    /// // Check if the converted coordinate is accurate to 6 decimals (same as reference)
    /// assert_eq!(coord.zone(), converted.zone());
    /// assert_eq!(coord.is_north(), converted.is_north());
    /// assert!((coord.easting() - converted.easting()).abs() < 1e-2);
    /// assert!((coord.northing() - converted.northing()).abs() < 1e-2);
    /// assert_eq!(coord.precision(), converted.precision());
    /// ```
    pub fn to_mgrs(&self, precision: i32) -> Mgrs {
        Mgrs {
            utm: *self,
            precision,
        }
    }
}

pub(crate) fn central_meridian(zone: i32) -> f64 {
    6.0 * f64::from(zone) - 183.
}

// Map lat/lon to zone in either UTM or UPS based on position.
fn standard_zone(lat: f64, lon: f64, setzone: i32) -> i32 {
    if setzone >= zonespec::MINZONE || setzone == zonespec::INVALID {
        return setzone;
    }

    if setzone == zonespec::UTM || ((-80_f64)..84.0).contains(&lat) {
        let mut lon_int = lon.ang_normalize().floor() as i32;
        if lon_int == dms::HD {
            lon_int = -dms::HD;
        }

        let mut zone = (lon_int + 186) / 6;
        let band = to_latitude_band(lat);
        // The Norway exception
        if band == 7 && zone == 31 && lon_int >= 3 {
            zone = 32;
        }
        // The Svalbard exception
        else if band == 9 && (0..42).contains(&lon_int) {
            zone = 2 * ((lon_int + 183) / 12) + 1;
        }

        zone
    } else {
        zonespec::UPS
    }
}

pub(crate) fn check_coords(utmp: bool, northp: bool, x: f64, y: f64, mgrs_limits: bool) -> Result<(), Error> {
    let slop = f64::from(mgrs::TILE);

    let ind = utmp.ternary(2, 0) + northp.ternary(1, 0);
    if x < f64::from(MIN_EASTING[ind]) - slop || x > f64::from(MAX_EASTING[ind]) + slop {
        return Err(Error::InvalidUtmCoords(
            format!(
                "Easting {:.2}km not in {}{} range for {} hemisphere [{:.2}km, {:.2}km]",
                x / 1000.0,
                mgrs_limits.ternary("MGRS/", ""),
                utmp.ternary("UTM", "UPS"),
                northp.ternary("N", "S"),
                (f64::from(MIN_EASTING[ind]) - slop) / 1000.0,
                (f64::from(MAX_EASTING[ind]) + slop) / 1000.0,
            )
        ));
    }

    if y < f64::from(MIN_NORTHING[ind]) - slop || y > f64::from(MAX_NORTHING[ind]) + slop {
        return Err(Error::InvalidUtmCoords(
            format!(
                "Northing {:.2}km not in {}{} range for {} hemisphere [{:.2}km, {:.2}km]",
                y / 1000.0,
                mgrs_limits.ternary("MGRS/", ""),
                utmp.ternary("UTM", "UPS"),
                northp.ternary("N", "S"),
                (f64::from(MIN_NORTHING[ind]) - slop) / 1000.0,
                (f64::from(MAX_NORTHING[ind]) + slop) / 1000.0,
            )
        ));
    }

    Ok(())
}

impl std::fmt::Display for UtmUps {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{} {} {}",
            self.zone,
            self.northp.ternary("n", "s"),
            self.easting,
            self.northing
        )
    }
}