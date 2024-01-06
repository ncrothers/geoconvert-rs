use crate::{latlon::LatLon, utility::{dms, GeoMath}, mgrs::{to_latitude_band, self}, Error, ThisOrThat};

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

pub const FALSE_EASTING: [i32; 4] = [
    mgrs::UPSEASTING * mgrs::TILE,
    mgrs::UPSEASTING * mgrs::TILE,
    mgrs::UTMEASTING * mgrs::TILE,
    mgrs::UTMEASTING * mgrs::TILE,
];

pub const FALSE_NORTHING: [i32; 4] = [
    mgrs::UPSEASTING * mgrs::TILE,
    mgrs::UPSEASTING * mgrs::TILE,
    mgrs::MAXUTM_S_ROW * mgrs::TILE,
    mgrs::MINUTM_N_ROW * mgrs::TILE,
];

pub const MIN_EASTING: [i32; 4] = [
    mgrs::MINUPS_S_IND * mgrs::TILE,
    mgrs::MINUPS_N_IND * mgrs::TILE,
    mgrs::MINUTMCOL * mgrs::TILE,
    mgrs::MINUTMCOL * mgrs::TILE,
];

pub const MAX_EASTING: [i32; 4] = [
    mgrs::MAXUPS_S_IND * mgrs::TILE,
    mgrs::MAXUPS_N_IND * mgrs::TILE,
    mgrs::MAXUTMCOL * mgrs::TILE,
    mgrs::MAXUTMCOL * mgrs::TILE,
];

pub const MIN_NORTHING: [i32; 4] = [
    mgrs::MINUPS_S_IND * mgrs::TILE,
    mgrs::MINUPS_N_IND * mgrs::TILE,
    mgrs::MINUTM_S_ROW * mgrs::TILE,
    (mgrs::MINUTM_N_ROW + mgrs::MINUTM_S_ROW - mgrs::MAXUTM_S_ROW) * mgrs::TILE,
];

pub const MAX_NORTHING: [i32; 4] = [
    mgrs::MAXUPS_S_IND * mgrs::TILE,
    mgrs::MAXUPS_N_IND * mgrs::TILE,
    (mgrs::MAXUTM_S_ROW + mgrs::MAXUTM_N_ROW - mgrs::MINUTM_N_ROW) * mgrs::TILE,
    mgrs::MAXUTM_N_ROW * mgrs::TILE,
];

#[derive(Debug)]
#[non_exhaustive]
pub struct Utm {
    pub zone: i32,
    pub northp: bool,
    pub easting: f64,
    pub northing: f64,
}

impl Utm {
    pub fn new(zone: i32, northp: bool, easting: f64, northing: f64) -> Utm {
        Self {
            zone,
            northp,
            easting,
            northing,
        }
    }
}

pub(crate) fn central_meridian(zone: i32) -> f64 {
    6. * zone as f64 - 183.
}

fn standard_zone(lat: f64, lon: f64, setzone: i32) -> i32 {
    if setzone >= zonespec::MINZONE || setzone == zonespec::INVALID {
        return setzone;
    }

    if setzone == zonespec::UTM || ((-80_f64)..84.).contains(&lat) {
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
    let slop = mgrs::TILE as f64;
    

    let ind = utmp.ternary(2, 0) + northp.ternary(1, 0);
    if x < MIN_EASTING[ind] as f64 - slop || x > MAX_EASTING[ind] as f64 + slop {
        return Err(Error::InvalidUtmCoords(
            format!(
                "Easting {:.2}km not in {}{} range for {} hemisphere [{:.2}km, {:.2}km]",
                x / 1000.,
                mgrs_limits.ternary("MGRS/", ""),
                utmp.ternary("UTM", "UPS"),
                northp.ternary("N", "S"),
                (MIN_EASTING[ind] as f64 - slop) / 1000.,
                (MAX_EASTING[ind] as f64 + slop) / 1000.,
            )
        ));
    }

    if y < MIN_NORTHING[ind] as f64 - slop || y > MAX_NORTHING[ind] as f64 + slop {
        return Err(Error::InvalidUtmCoords(
            format!(
                "Northing {:.2}km not in {}{} range for {} hemisphere [{:.2}km, {:.2}km]",
                y / 1000.,
                mgrs_limits.ternary("MGRS/", ""),
                utmp.ternary("UTM", "UPS"),
                northp.ternary("N", "S"),
                (MIN_NORTHING[ind] as f64 - slop) / 1000.,
                (MAX_NORTHING[ind] as f64 + slop) / 1000.,
            )
        ));
    }

    Ok(())
}

impl TryFrom<LatLon> for Utm {
    type Error = Error;
    fn try_from(value: LatLon) -> Result<Self, Self::Error> {
        if value.latitude.abs() > dms::QD as f64 {
            return Err(Error::InvalidRange {
                coord_type: "LatLon".to_string(),
                dest_type: "Utm".to_string(),
                msg: format!(
                    "Latitude {}d not in [-{1}d, {1}d]",
                    value.latitude,
                    dms::QD
                )
            });
        }

        let northp = value.is_north();
        // STANDARD specifies, by default, interpret whether it should be UTM or UPS
        // TODO: Maybe do something if the zone is invalid?
        let zone = standard_zone(value.latitude, value.longitude, zonespec::STANDARD);
        let mut x = 0;
        let mut y = 0;

        let utmp = zone != zonespec::UPS;
        let (mut x, mut y) = if utmp {
            let lon = central_meridian(zone);
            let lon_diff = lon.ang_diff(value.longitude);

            // UTM().Forward
            (0., 0.)
        } else {
            // UPS().Forward
            (0., 0.)
        };

        let ind = utmp.ternary(2, 0) + northp.ternary(1, 0);
        x += f64::from(FALSE_EASTING[ind]);
        y += f64::from(FALSE_EASTING[ind]);

        Ok(Utm {
            zone,
            northp,
            northing: y,
            easting: x,
        })
    }
}
