use std::fmt::Display;

use lazy_static::lazy_static;

use crate::{ParseCoord, Error, utm::{zonespec::{self, MINUTMZONE, MAXUTMZONE, UPS}, Utm, check_coords}, utility::dms, ThisOrThat, latlon::LatLon};

const HEMISPHERES: &str = "SN";
const UTMCOLS: &[&str] = &["ABCDEFGH", "JKLMNPQR", "STUVWXYZ"];
const UTMROW: &str = "ABCDEFGHJKLMNPQRSTUV";
const UPSCOLS: &[&str] = &["JKLPQRSTUXYZ", "ABCFGHJKLPQR", "RSTUXYZ", "ABCFGHJ"];
const UPSROWS: &[&str] = &["ABCDEFGHJKLMNPQRSTUVWXYZ", "ABCDEFGHJKLMNP"];
const LATBAND: &str = "CDEFGHJKLMNPQRSTUVWX";
const UPSBAND: &str = "ABYZ";
const DIGITS: &str = "0123456789";
const ALPHA: &str = "ABCDEFGHJKLMNPQRSTUVWXYZabcdefghjklmnpqrstuvwxyz"; // Omit I+O

pub const TILE: i32= 100_000;
pub const MINUTMCOL: i32= 1;
pub const MAXUTMCOL: i32= 9;
pub const MINUTM_S_ROW: i32= 10;
pub const MAXUTM_S_ROW: i32= 100;
pub const MINUTM_N_ROW: i32= 0;
pub const MAXUTM_N_ROW: i32= 95;
pub const MINUPS_S_IND: i32= 8;
pub const MAXUPS_S_IND: i32= 32;
pub const MINUPS_N_IND: i32= 13;
pub const MAXUPS_N_IND: i32= 27;
pub const UPSEASTING: i32= 20;
pub const UTMEASTING: i32= 5;
pub const UTM_N_SHIFT: i32= (MAXUTM_S_ROW - MINUTM_N_ROW) * TILE;

pub const BASE: i32= 10;
pub const TILE_LEVEL: i32= 5;
pub const UTM_ROW_PERIOD: i32 = 20;
pub const UTM_EVEN_ROW_SHIFT: i32= 5;
pub const MAX_PRECISION: i32= 5 + 6;
pub const MULT: i32= 1_000_000;

#[derive(Debug)]
pub struct Mgrs {
    pub(crate) utm: Utm,
    precision: i32,
}

impl Mgrs {
    pub fn is_utm(&self) -> bool {
        self.utm.zone != UPS
    }
}

fn utm_row(band_idx: i32, col_idx: i32, row_idx: i32) -> i32 {
    let c = 100. * (8. * band_idx as f64 + 4.) / dms::QD as f64;
    let northp = band_idx >= 0;
    // These are safe bounds on the rows
    //  band_idx  minrow maxrow
    //   -10      -90    -81
    //    -9      -80    -72
    //    -8      -71    -63
    //    -7      -63    -54
    //    -6      -54    -45
    //    -5      -45    -36
    //    -4      -36    -27
    //    -3      -27    -18
    //    -2      -18     -9
    //    -1       -9     -1
    //     0        0      8
    //     1        8     17
    //     2       17     26
    //     3       26     35
    //     4       35     44
    //     5       44     53
    //     6       53     62
    //     7       62     70
    //     8       71     79
    //     9       80     94

    let min_row = if band_idx > -10 {
        (c - 4.3 - 0.1 * northp as u8 as f64).floor() as i32
    } else {
        -90
    };

    let max_row = if band_idx < 9 {
        (c + 4.4 - 0.1 * northp as u8 as f64).floor() as i32
    } else {
        94
    };

    let base_row = (min_row + max_row) / 2 - UTM_ROW_PERIOD as i32 / 2;
    // Offset row_idx by the multiple of UTM_ROW_PERIOD which brings it as close as
    // possible to the center of the latitude band, (min_row + max_row) / 2.
    // (Add MAXUTM_S_ROW = 5 * UTM_ROW_PERIOD to ensure operand is positive.)
    let mut row_idx = (row_idx - base_row + MAXUTM_S_ROW as i32) % UTM_ROW_PERIOD as i32 + base_row;
    
    if !(row_idx >= min_row && row_idx <= max_row) {
        // Outside the safe bounds, so need to check...
        // Northing = 71e5 and 80e5 intersect band boundaries
        //   y = 71e5 in scol = 2 (x = [3e5,4e5] and x = [6e5,7e5])
        //   y = 80e5 in scol = 1 (x = [2e5,3e5] and x = [7e5,8e5])
        // This holds for all the ellipsoids given in NGA.SIG.0012_2.0.0_UTMUPS.
        // The following deals with these special cases.

        // Fold [-10,-1] -> [9,0]
        let safe_band = (band_idx >= 0).ternary(band_idx, -band_idx - 1);
        // Fold [-90,-1] -> [89,0]
        let safe_row = (row_idx >= 0).ternary(row_idx, -row_idx - 1);
        // Fold [4,7] -> [3,0]
        let safe_col = (col_idx < 4).ternary(col_idx, -col_idx + 7);

        if !(
            (safe_row == 70 && safe_band == 8 && safe_col >= 2) ||
            (safe_row == 71 && safe_band == 7 && safe_col <= 2) ||
            (safe_row == 79 && safe_band == 9 && safe_col >= 1) ||
            (safe_row == 80 && safe_band == 8 && safe_col <= 1)
        ) {
            row_idx = MAXUTM_S_ROW as i32;
        }
    }

    row_idx
}

impl ParseCoord for Mgrs {
    #[allow(clippy::too_many_lines)]
    fn parse_coord(value: &str) -> Result<Self, Error> {
        let value = value.to_ascii_uppercase();
        let mut p = 0;
        let mut len = value.len();
        let chars = value.as_bytes();

        if len >= 3 && value.starts_with("INV") {
            return Err(Error::InvalidMgrs("Starts with 'INV'".to_string()))
        }

        let mut zone = 0i32;
        while p < len {
            if let Some(i) = DIGITS.find(chars[p] as char) {
                zone = 10 * zone + i as i32;
                p += 1;
            }
            else {
                break;
            }
        }
        // Check if zone is within valid range
        if p > 0 && !(MINUTMZONE..=MAXUTMZONE).contains(&zone) {
            return Err(Error::InvalidMgrs(format!("Zone {zone} not in [1,60]")));
        }

        if p > 2 {
            return Err(Error::InvalidMgrs(format!("More than 2 digits at start of MGRS {}", &value[..p])));
        }

        if len - p < 1 {
            return Err(Error::InvalidMgrs(format!("Too short: {value}")));
        }

        let utmp = zone != UPS;
        let zonem = zone - 1;
        let band = utmp.ternary(LATBAND, UPSBAND);
        
        let mut band_idx = match band.find(chars[p] as char) {
            Some(i) => i as i32,
            None => {
                let label = utmp.ternary("UTM", "UPS");
                return Err(Error::InvalidMgrs(format!("Band letter {} not in {label} set {band}", chars[p] as char)));
            }
        };
        
        p += 1;

        let northp = band_idx >= utmp.ternary(10, 2);

        if p == len { // Grid zone only (ignore centerp)
            // Approx length of a degree of meridian arc in units of tile
            let deg = (UTM_N_SHIFT as f64) / (dms::QD * TILE) as f64;
            let (x, y) = if utmp {
                // Pick central meridian except for 31V
                let x = TILE as f64 * (zone == 31 && band_idx == 17).ternary(4., 5.);
                // TODO: continue from here
                let y_add = northp.ternary(0., UTM_N_SHIFT as f64);
                let y = (8. * (band_idx as f64 - 9.5) * deg + 0.5).floor() * TILE as f64 + y_add;

                (x, y)
            } else {
                let x_cond = (band_idx % 2 == 1).ternary(1., -1.);
                let x = (x_cond * (4. * deg + 0.5).floor() + UPSEASTING as f64) * TILE as f64;
                let y = (UPSEASTING * TILE) as f64;
                (x, y)
            };

            return Ok(Mgrs {
                utm: Utm::new(zone, northp, x, y),
                precision: -1
            })
        } else if len - p < 2 {
            return Err(Error::InvalidMgrs(format!("Missing row letter in {value}")));
        }

        let col = utmp.ternary_lazy(|| UTMCOLS[(zonem % 3) as usize], || UPSCOLS[band_idx as usize]);
        let row = utmp.ternary_lazy(|| UTMROW, || UPSROWS[northp as usize]);
        let mut col_idx = col
            .find(chars[p] as char)
            .ok_or_else(|| {
                let label = if utmp { format!("zone {}", &value[..p-1]) } else { format!("UPS band {}", &value[p-1..p]) };
                Error::InvalidMgrs(format!("Column letter {} not in {label} set {col}", &value[p..=p]))
            })? as i32;

        p += 1;

        let mut row_idx = row
            .find(chars[p] as char)
            .ok_or_else(|| {
                let northp = northp as usize;
                let label = if utmp { "UTM".to_string() } else { format!("UPS {}", &HEMISPHERES[northp..=northp]) };
                Error::InvalidMgrs(format!("Row letter {} not in {label} set {row}", chars[p]))
            })? as i32;

        p += 1;

        if utmp {
            if zonem % 2 == 1 {
                row_idx = (row_idx + UTM_ROW_PERIOD as i32 - UTM_EVEN_ROW_SHIFT as i32) % UTM_ROW_PERIOD as i32;
            }

            band_idx -= 10;

            row_idx = utm_row(band_idx, col_idx as i32, row_idx as i32);
            if row_idx == MAXUTM_S_ROW as i32 {
                return Err(Error::InvalidMgrs(format!("Block {} not in zone/band {}", &value[p-2..p], &value[0..p-2])))
            }

            row_idx = northp.ternary(row_idx, row_idx + 100);
            col_idx += MINUTMCOL;
        }
        else {
            let eastp = band_idx % 2 == 1;
            col_idx += if eastp { UPSEASTING } else if northp { MINUPS_N_IND } else { MINUPS_S_IND };
            row_idx += if northp { MINUPS_N_IND as i32 } else { MINUPS_S_IND as i32 };
        }

        let precision = (len - p) / 2;
        let mut unit = 1;
        let mut x = col_idx;
        let mut y = row_idx;

        for i in 0..precision {
            unit *= BASE;
            let x_idx = DIGITS
                .find(chars[p + i] as char)
                .ok_or_else(|| {
                    Error::InvalidMgrs(format!("Encountered a non-digit in {}", &value[p..]))
                })?;
            let y_idx = DIGITS
                .find(chars[p + i + precision] as char)
                .ok_or_else(|| {
                    Error::InvalidMgrs(format!("Encountered a non-digit in {}", &value[p..]))
                })?;
            
            x = BASE * x + x_idx as i32;
            y = BASE * y + y_idx as i32;
        }

        if (len - p) % 2 == 1 {
            if DIGITS.find(chars[len - 1] as char).is_none() {
                return Err(Error::InvalidMgrs(format!("Encountered a non-digit in {}", &value[p..])));
            }

            return Err(Error::InvalidMgrs(format!("Not an even number of digits in {}", &value[p..])));
        }

        if precision > MAX_PRECISION as usize {
            return Err(Error::InvalidMgrs(format!("More than {} digits in {}", 2*MAX_PRECISION, &value[p..])));
        }

        // TODO: Include centerp somehow
        let centerp = true;
        if centerp {
            unit *= 2;
            x = 2 * x + 1;
            y = 2 * y + 1;
        }

        let x = (TILE as f64 * x as f64) / unit as f64;
        let y = (TILE as f64 * y as f64) / unit as f64;

        Ok(Self {
            utm: Utm::new(
                zone,
                northp,
                x,
                y,
            ),
            precision: precision as i32,
        })
    }
}

pub(crate) fn to_latitude_band(lat: f64) -> i32 {
    let lat_int = lat.floor() as i32;
    (-10).max(9.min((lat_int + 80) / 8 - 10))
}

lazy_static! {
    static ref ANG_EPS: f64 = 1_f64 * 2_f64.powi(-(f64::DIGITS as i32 - 7));
}

// TODO: Finish, requires UTM and UPS stuff though
impl Display for Mgrs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let lat = if self.utm.zone > 0 {
            // Does a rough estimate for latitude determine the latitude band?
            let y_est = if self.utm.northp { self.utm.northing } else { self.utm.northing - UTM_N_SHIFT as f64 };
            // A cheap calculation of the latitude which results in an "allowed"
            // latitude band would be
            //   lat = ApproxLatitudeBand(ys) * 8 + 4;
            //
            // Here we do a more careful job using the band letter corresponding to
            // the actual latitude.
            let y_est = y_est / TILE as f64;
            if y_est.abs() < 1. {
                0.9 * y_est
            }
            else {
                let pole_add = if y_est > 0. { 1. } else { -1. };
                let lat_poleward = 0.901 * y_est + pole_add * 0.135;
                let lat_eastward = 0.902 * y_est * (1. - 1.85e-6 * y_est.powi(2));

                if to_latitude_band(lat_poleward) == to_latitude_band(lat_eastward) {
                    lat_poleward
                } else {
                    let coord = LatLon::try_from(Utm::new(self.utm.zone as i32, self.utm.northp, self.utm.easting, self.utm.northing)).expect("Couldn't convert Utm to latlon");
                    coord.latitude
                }
            }
        } else {
            0.
        };
        
        // Other Forward call
        let utmp = self.utm.zone != 0;
        // TODO: Check coords on creation rather than in infallible methods like fmt()
        check_coords(utmp, self.utm.northp, self.utm.easting, self.utm.northing, true).expect("Invalid coords");
        // Create pre-allocated string of the correct length
        let mut mgrs_str = [0u8; 2 + 3 + 2*MAX_PRECISION as usize];
        let zone = self.utm.zone - 1;
        let mut z: usize = utmp.ternary(2, 0);
        let mlen = z as i32 + 3 + 2 * self.precision;

        let digits = DIGITS.as_bytes();

        if utmp {
            mgrs_str[0] = digits[(self.utm.zone / BASE) as usize];
            mgrs_str[1] = digits[(self.utm.zone % BASE) as usize];
        }

        let xx = self.utm.easting * MULT as f64;
        let yy = self.utm.northing * MULT as f64;

        let ix = xx.floor() as i64;
        let iy = yy.floor() as i64;
        let m = (MULT as i64) * (TILE as i64);

        let xh = (ix / m) as i32;
        let yh = (iy / m) as i32;

        if utmp {
            // Correct fuzziness in latitude near equator
            let band_idx = (lat.abs() < *ANG_EPS).ternary(self.utm.northp.ternary(0, -1), to_latitude_band(lat));
            let col_idx = xh - MINUTMCOL;
            let row_idx = utm_row(band_idx, col_idx, yh % UTM_ROW_PERIOD);

            if row_idx != yh - self.utm.northp.ternary(MINUTM_N_ROW, MAXUTM_S_ROW) {
                // TODO: Latitude is inconsistent with UTM coordinates
                todo!()
            }

            mgrs_str[z] = LATBAND.as_bytes()[(10 + band_idx) as usize];
            z += 1;
            mgrs_str[z] = UTMCOLS[(zone % 3) as usize].as_bytes()[col_idx as usize];
            z += 1;
            let idx = (yh + (zone % 2 == 1).ternary(UTM_EVEN_ROW_SHIFT, 0)) % UTM_ROW_PERIOD;
            mgrs_str[z] = UTMROW.as_bytes()[idx as usize];
            z += 1;
        } else {
            let eastp = xh >= UPSEASTING;
            let band_idx: usize = self.utm.northp.ternary(2, 0) + eastp.ternary(1, 0);
            mgrs_str[z] = UPSBAND.as_bytes()[band_idx];
            z += 1;
            let idx = xh - eastp.ternary(UPSEASTING, self.utm.northp.ternary(MINUPS_N_IND, MINUPS_S_IND));
            mgrs_str[z] = UPSCOLS[band_idx].as_bytes()[idx as usize];
            z += 1;
            let idx = yh - self.utm.northp.ternary(MINUPS_N_IND, MINUPS_S_IND);
            mgrs_str[z] = UPSROWS[usize::from(self.utm.northp)].as_bytes()[idx as usize];
            z += 1;
        }

        if self.precision > 0 {
            let mut ix = ix - m * xh as i64;
            let mut iy = iy - m * yh as i64;
            let d = (BASE as i64).pow((MAX_PRECISION - self.precision) as u32);
            ix /= d;
            iy /= d;

            for c in (0..=(self.precision as usize)-1).rev() {
                mgrs_str[z + c] = digits[(ix % BASE as i64) as usize];
                ix /= BASE as i64;
                mgrs_str[z + c + self.precision as usize] = digits[(iy % BASE as i64) as usize];
                iy /= BASE as i64;
            }
        }

        write!(f, "{}", String::from_utf8_lossy(&mgrs_str).trim_end_matches('\0'))
    }
}

