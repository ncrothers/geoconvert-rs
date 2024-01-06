use std::fmt::Display;

use crate::{Error, utm::{Utm, zonespec, self, central_meridian}, ThisOrThat, projections::{transverse_mercator::TransverseMercator, polar_stereographic::PolarStereographic}, mgrs::Mgrs, constants::WGS84_A};

#[derive(Debug)]
pub struct LatLon {
    pub latitude: f64,
    pub longitude: f64,
}

impl LatLon {
    pub fn is_north(&self) -> bool {
        self.latitude.is_sign_positive()
    }

    pub fn haversine(&self, other: &LatLon) -> f64 {
        let lat1 = self.latitude.to_radians();
        let lat2 = other.latitude.to_radians();
        let lon1 = self.longitude.to_radians();
        let lon2 = self.longitude.to_radians();
        
        2. * 6_371_088. * (
            ((lat2 - lat1) / 2.).sin().powi(2) + 
            lat1.cos() * lat2.cos() *
            ((lon2 - lon1) / 2.).sin().powi(2)
        ).sqrt().asin()
    }
}

impl TryFrom<Utm> for LatLon {
    type Error = Error;

    fn try_from(value: Utm) -> Result<Self, Self::Error> {
        if !(zonespec::MINZONE..=zonespec::MAXZONE).contains(&value.zone) {
            return Err(Error::InvalidRange {
                coord_type: "Utm".to_string(),
                dest_type: "LatLon".to_string(),
                msg: format!("Zone {} not in range [0, 60]", value.zone)
            });
        }

        let utmp = value.zone != zonespec::UPS;
        // TODO: Look into mgrs_limits and if/how to implement them
        utm::check_coords(utmp, value.northp, value.easting, value.northing, false)?;
        let ind = utmp.ternary(2, 0) + value.northp.ternary(1, 0);

        let x = value.easting - utm::FALSE_EASTING[ind] as f64;
        let y = value.northing - utm::FALSE_NORTHING[ind] as f64;

        let coord = if utmp {
            TransverseMercator::utm().utm_to_latlon(central_meridian(value.zone), x, y)
        } else {
            PolarStereographic::ups().to_latlon(value.northp, x, y)
        };

        Ok(coord)
    }
}

impl TryFrom<Mgrs> for LatLon {
    type Error = Error;

    fn try_from(value: Mgrs) -> Result<Self, Self::Error> {
        LatLon::try_from(value.utm)
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