use crate::{latlon::LatLon, utility::{dms, GeoMath}, mgrs::{to_latitude_band, self, Mgrs}, Error, ThisOrThat, projections::{transverse_mercator::TransverseMercator, polar_stereographic::PolarStereographic}};

pub mod zonespec {
    pub const INVALID: i32 = -4;
    pub const MATCH: i32 = -3;
    pub const UTM: i32 = -2;
    pub const STANDARD: i32 = -1;
    pub const UPS: i32 = 0;
    pub const MINZONE: i32 = 0;
    pub const MINUTMZONE: i32 = 1;
    pub const MAXUTMZONE: i32 = 60;
    pub const MAXZONE: i32 = 60;
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

#[allow(clippy::module_name_repetitions)]
#[derive(Clone, Copy, Debug)]
pub struct UtmUps {
    pub(crate) zone: i32,
    pub(crate) northp: bool,
    pub(crate) easting: f64,
    pub(crate) northing: f64,
}

impl UtmUps {
    pub(crate) fn new(zone: i32, northp: bool, easting: f64, northing: f64) -> UtmUps {
        Self {
            zone,
            northp,
            easting,
            northing,
        }
    }

    /// TODO
    /// 
    /// # Errors
    pub fn create(zone: i32, northp: bool, easting: f64, northing: f64) -> Result<UtmUps, Error> {
        // Make sure zone is a valid value
        if !(zonespec::MINZONE..=zonespec::MAXZONE).contains(&zone) {
            return Err(Error::InvalidRange {
                coord_type: "Utm".to_string(),
                dest_type: "LatLon".to_string(),
                msg: format!("Zone {zone} not in range [0, 60]")
            });
        }

        let utmp = zone != zonespec::UPS;

        check_coords(utmp, northp, easting, northing, false)?;

        Ok(UtmUps::new(zone, northp, easting, northing))
    }

    pub fn zone(&self) -> i32 {
        self.zone
    }

    pub fn is_north(&self) -> bool {
        self.northp
    }

    pub fn easting(&self) -> f64 {
        self.easting
    }

    pub fn northing(&self) -> f64 {
        self.northing
    }

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

    pub fn from_mgrs(value: &Mgrs) -> UtmUps {
        value.utm
    }
    
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
