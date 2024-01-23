use std::{fmt::Display, str::FromStr};

use lazy_static::lazy_static;
use num::Integer;

use crate::{Error, utm::{zonespec::{MINUTMZONE, MAXUTMZONE, UPS, self}, UtmUps}, utility::{dms, GeoMath}, ThisOrThat, latlon::LatLon};

const HEMISPHERES: &str = "SN";
const UTMCOLS: &[&str] = &["ABCDEFGH", "JKLMNPQR", "STUVWXYZ"];
const UTMROW: &str = "ABCDEFGHJKLMNPQRSTUV";
const UPSCOLS: &[&str] = &["JKLPQRSTUXYZ", "ABCFGHJKLPQR", "RSTUXYZ", "ABCFGHJ"];
const UPSROWS: &[&str] = &["ABCDEFGHJKLMNPQRSTUVWXYZ", "ABCDEFGHJKLMNP"];
const LATBAND: &str = "CDEFGHJKLMNPQRSTUVWX";
const UPSBAND: &str = "ABYZ";
const DIGITS: &str = "0123456789";

pub(crate) const TILE: i32= 100_000;
pub(crate) const MINUTMCOL: i32= 1;
pub(crate) const MAXUTMCOL: i32= 9;
pub(crate) const MINUTM_S_ROW: i32= 10;
pub(crate) const MAXUTM_S_ROW: i32= 100;
pub(crate) const MINUTM_N_ROW: i32= 0;
pub(crate) const MAXUTM_N_ROW: i32= 95;
pub(crate) const MINUPS_S_IND: i32= 8;
pub(crate) const MAXUPS_S_IND: i32= 32;
pub(crate) const MINUPS_N_IND: i32= 13;
pub(crate) const MAXUPS_N_IND: i32= 27;
pub(crate) const UPSEASTING: i32= 20;
pub(crate) const UTMEASTING: i32= 5;
pub(crate) const UTM_N_SHIFT: i32= (MAXUTM_S_ROW - MINUTM_N_ROW) * TILE;

const MIN_EASTING: [i32; 4] = [
    MINUPS_S_IND,
    MINUPS_N_IND,
    MINUTMCOL,
    MINUTMCOL,
];

const MAX_EASTING: [i32; 4] = [
    MAXUPS_S_IND,
    MAXUPS_N_IND,
    MAXUTMCOL,
    MAXUTMCOL,
];

const MIN_NORTHING: [i32; 4] = [
    MINUPS_S_IND,
    MINUPS_N_IND,
    MINUTM_S_ROW,
    MINUTM_S_ROW - MAXUTM_S_ROW - MINUTM_N_ROW,
];

const MAX_NORTHING: [i32; 4] = [
    MAXUPS_S_IND,
    MAXUPS_N_IND,
    MAXUTM_N_ROW + MAXUTM_S_ROW - MINUTM_N_ROW,
    MAXUTM_N_ROW,
];

pub(crate) const BASE: i32= 10;
pub(crate) const UTM_ROW_PERIOD: i32 = 20;
pub(crate) const UTM_EVEN_ROW_SHIFT: i32= 5;
pub(crate) const MAX_PRECISION: i32= 5 + 6;
pub(crate) const MULT: i32= 1_000_000;

/// Representation of a WGS84 
/// [Military Grid Reference System](https://en.wikipedia.org/wiki/Military_Grid_Reference_System)
/// point. Stored internally as a [`UtmUps`] point with a precision.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Mgrs {
    #[cfg_attr(feature = "serde", serde(flatten))]
    pub(crate) utm: UtmUps,
    pub(crate) precision: i32,
}

impl Mgrs {
    /// Tries to create a MGRS point from its constituent parts. Validates the
    /// arguments to ensure a valid MGRS point can be created. You most likely
    /// want to instantiate this via [`parse_str`](#method.parse_str) or [`UtmUps`] instead
    /// of manually specifying the values.
    /// 
    /// # Errors
    /// 
    /// Returns [`Error::InvalidMgrs`] if the position is invalid.
    /// Returns [`Error::InvalidPrecision`] if the precision is not in range `[1, 11]`.
    /// 
    /// # Usage
    /// 
    /// ```
    /// use geoconvert::Mgrs;
    /// 
    /// let coord = Mgrs::create(18, true, 585664.121, 4511315.422, 6);
    /// 
    /// assert!(coord.is_ok());
    /// 
    /// let coord = coord.unwrap();
    /// 
    /// assert_eq!(coord.zone(), 18);
    /// assert_eq!(coord.is_north(), true);
    /// assert!((coord.easting() - 585664.121).abs() < 1e-3);
    /// assert!((coord.northing() - 4511315.422).abs() < 1e-3);
    /// assert_eq!(coord.precision(), 6);
    /// 
    /// let invalid_coord_zone_neg = Mgrs::create(-10, true, 585664.121, 4511315.422, 6);
    /// assert!(invalid_coord_zone_neg.is_err());
    /// 
    /// let invalid_coord_zone_too_big = Mgrs::create(70, true, 585664.121, 4511315.422, 6);
    /// assert!(invalid_coord_zone_too_big.is_err());
    /// ```
    pub fn create(zone: i32, northp: bool, easting: f64, northing: f64, precision: i32) -> Result<Mgrs, Error> {
        // Make sure zone is a valid value
        if !(zonespec::MINZONE..=zonespec::MAXZONE).contains(&zone) {
            return Err(Error::InvalidZone(zone));
        }

        let utmp = zone != zonespec::UPS;

        check_coords(utmp, northp, easting, northing)?;

        Ok(Mgrs {
            utm: UtmUps::new(zone, northp, easting, northing),
            precision,
        })
    }

    /// Returns whether the MGRS is stored as UTM or UPS.
    /// 
    /// # Example
    /// ```
    /// use geoconvert::Mgrs;
    /// 
    /// let coord = Mgrs::parse_str("18TWL856641113154").unwrap();
    /// assert_eq!(coord.is_utm(), true);
    /// ```
    #[inline]
    pub fn is_utm(&self) -> bool {
        self.utm.zone != UPS
    }

    /// Returns the UTM zone.
    /// 
    /// # Example
    /// ```
    /// use geoconvert::Mgrs;
    /// 
    /// let coord = Mgrs::parse_str("18TWL856641113154").unwrap();
    /// assert_eq!(coord.zone(), 18);
    /// ```
    #[inline]
    pub fn zone(&self) -> i32 {
        self.utm.zone
    }

    /// Returns whether the coordinate is in the northern hemisphere.
    /// 
    /// # Example
    /// ```
    /// use geoconvert::Mgrs;
    /// 
    /// let coord = Mgrs::parse_str("18TWL856641113154").unwrap();
    /// assert_eq!(coord.is_north(), true);
    /// ```
    #[inline]
    pub fn is_north(&self) -> bool {
        self.utm.northp
    }

    /// Returns the UTM easting.
    /// 
    /// # Example
    /// ```
    /// use geoconvert::Mgrs;
    /// 
    /// let coord = Mgrs::parse_str("18TWL856641113154").unwrap();
    /// assert!((coord.easting() - 585664.15).abs() < 1e-2);
    /// ```
    #[inline]
    pub fn easting(&self) -> f64 {
        self.utm.easting
    }

    /// Returns the UTM northing.
    /// 
    /// # Example
    /// ```
    /// use geoconvert::Mgrs;
    /// 
    /// let coord = Mgrs::parse_str("18TWL856641113154").unwrap();
    /// assert!((coord.northing() - 4511315.45).abs() < 1e-2);
    /// ```
    #[inline]
    pub fn northing(&self) -> f64 {
        self.utm.northing
    }

    /// Returns the current precision for outputting to a string.
    /// 
    /// # Example
    /// ```
    /// use geoconvert::Mgrs;
    /// 
    /// let coord = Mgrs::parse_str("18TWL856641113154").unwrap();
    /// assert_eq!(coord.precision(), 6);
    /// ```
    #[inline]
    pub fn precision(&self) -> i32 {
        self.precision
    }

    /// Set the precision.
    /// 
    /// Must be in range `[1, 11]`.
    /// 
    /// # Errors
    /// 
    /// * [`Error::InvalidPrecision`]: `precision` is not in the valid range
    /// 
    /// # Example
    /// ```
    /// use geoconvert::Mgrs;
    /// 
    /// let mut coord = Mgrs::parse_str("18TWL856641113154").unwrap();
    /// coord.set_precision(7);
    /// 
    /// assert_eq!(coord.precision(), 7);
    /// ```
    #[inline]
    pub fn set_precision(&mut self, precision: i32) -> Result<(), Error> {
        if !(1..=11).contains(&precision) {
            return Err(Error::InvalidPrecision(precision));
        }

        self.precision = precision;
        Ok(())
    }

    /// Parses a string as MGRS. Assumes the string is _only_ composed of
    /// the MGRS coordinate (e.g. no preceding/trailing whitespace) and there
    /// are no spaces in the string. Example valid strings:
    /// 
    /// * `27UXQ0314512982`
    /// * `YXL6143481146`
    /// 
    /// # Errors
    /// 
    /// * [`Error::InvalidMgrs`]: the string couldn't be parsed to a valid MGRS coordinate.
    pub fn parse_str(mgrs_str: &str) -> Result<Mgrs, Error> {
        Self::from_str(mgrs_str)
    }

    /// Converts from [`LatLon`] to [`Mgrs`]
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
    pub fn from_latlon(value: &LatLon, precision: i32) -> Mgrs {
        Mgrs {
            utm: UtmUps::from_latlon(value),
            precision,
        }
    }

    /// Converts from [`Mgrs`] to [`LatLon`]
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
    pub fn to_latlon(&self) -> LatLon {
        self.utm.to_latlon()
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
    pub fn from_utmups(value: &UtmUps, precision: i32) -> Mgrs {
        Mgrs {
            utm: *value,
            precision,
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
    pub fn to_utmups(&self) -> UtmUps {
        self.utm
    }
}

fn utm_row(band_idx: i32, col_idx: i32, row_idx: i32) -> i32 {
    let c = 100.0 * (8.0 * f64::from(band_idx) + 4.0) / f64::from(dms::QD);
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
        (c - 4.3 - 0.1 * f64::from(u8::from(northp))).floor() as i32
    } else {
        -90
    };

    let max_row = if band_idx < 9 {
        (c + 4.4 - 0.1 * f64::from(u8::from(northp))).floor() as i32
    } else {
        94
    };

    let base_row = (min_row + max_row) / 2 - UTM_ROW_PERIOD / 2;
    // Offset row_idx by the multiple of UTM_ROW_PERIOD which brings it as close as
    // possible to the center of the latitude band, (min_row + max_row) / 2.
    // (Add MAXUTM_S_ROW = 5 * UTM_ROW_PERIOD to ensure operand is positive.0)
    let mut row_idx = (row_idx - base_row + MAXUTM_S_ROW) % UTM_ROW_PERIOD + base_row;
    
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
            row_idx = MAXUTM_S_ROW;
        }
    }

    row_idx
}

impl FromStr for Mgrs {
    type Err = Error;

    #[allow(clippy::too_many_lines)]
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let value = s.to_ascii_uppercase();
        let mut p = 0;
        let len = value.len();
        if !value.is_ascii() {
            return Err(Error::InvalidMgrs("String contains unicode characters".to_string()))
        }
        let chars = value.as_bytes();

        if len >= 3 && value.starts_with("INV") {
            return Err(Error::InvalidMgrs("Starts with 'INV'".to_string()))
        }

        let mut zone = 0i32;
        while p < len {
            // if let Some(i) = DIGITS_MAP.get(&(chars[p] as char)) {
            if (chars[p] as char).is_ascii_digit() {
                zone = 10 * zone + i32::from(chars[p] - b'0');
                p += 1;
            } else {
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

        let cur_char = chars[p];
        #[allow(clippy::collapsible_else_if)]
        let mut band_idx = if utmp {
            // First check if it's a valid latband
            if (b'C'..=b'X').contains(&cur_char) && cur_char != b'I' && cur_char != b'O' {
                // Then convert to index
                // First make it relative to C
                let idx = cur_char - b'C';
                // Decrement if it's past H to account for missing I
                let idx = (cur_char > b'H').ternary_lazy(|| idx - 1, || idx);
                // Decrement if it's past N to account for missing O
                let idx = (cur_char > b'N').ternary_lazy(|| idx - 1, || idx);
                i32::from(idx)
            } else {
                -1
            }
        } else {
            if cur_char == b'A' {
                0
            } else if cur_char == b'B' {
                1
            } else if cur_char == b'Y' {
                2
            } else if cur_char == b'Z' {
                3
            } else {
                -1
            }
        };

        if band_idx == -1 {
            let band = utmp.ternary(LATBAND, UPSBAND);
            let label = utmp.ternary("UTM", "UPS");
            return Err(Error::InvalidMgrs(format!("Band letter {} not in {label} set {band}", chars[p] as char)));
        }

        p += 1;

        let northp = band_idx >= utmp.ternary(10, 2);

        if p == len { // Grid zone only (ignore centerp)
            // Approx length of a degree of meridian arc in units of tile
            let deg = (f64::from(UTM_N_SHIFT)) / f64::from(dms::QD * TILE);
            let (x, y) = if utmp {
                // Pick central meridian except for 31V
                let x = f64::from(TILE) * (zone == 31 && band_idx == 17).ternary(4.0, 5.0);
                // TODO: continue from here
                let y_add = northp.ternary(0.0, f64::from(UTM_N_SHIFT));
                let y = (8.0 * (f64::from(band_idx) - 9.5) * deg + 0.5).floor() * f64::from(TILE) + y_add;

                (x, y)
            } else {
                let x_cond = band_idx.is_odd().ternary(1.0, -1.0);
                let x = (x_cond * (4.0 * deg + 0.5).floor() + f64::from(UPSEASTING)) * f64::from(TILE);
                let y = f64::from(UPSEASTING * TILE);
                (x, y)
            };

            return Ok(Mgrs {
                utm: UtmUps::new(zone, northp, x, y),
                precision: -1
            })
        } else if len - p < 2 {
            return Err(Error::InvalidMgrs(format!("Missing row letter in {value}")));
        }

        let cur_char = chars[p];
        // More efficient than find()
        let mut col_idx = if utmp {
            match zonem % 3 {
                0 => {
                    if (b'A'..=b'H').contains(&cur_char) {
                        i32::from(cur_char - b'A')
                    } else {
                        -1
                    }
                }
                1 => {
                    if (b'J'..=b'R').contains(&cur_char) && cur_char != b'O' {
                        if cur_char < b'O' {
                            i32::from(cur_char - b'J')
                        } else {
                            i32::from(cur_char - b'J' - 1)
                        }
                    } else {
                        -1
                    }
                } 
                2 => {
                    if (b'S'..=b'Z').contains(&cur_char) {
                        i32::from(cur_char - b'S')
                    } else {
                        -1
                    }
                }
                _ => unreachable!()
            }
        } else {
            // &["JKLPQRSTUXYZ", "ABCFGHJKLPQR", "RSTUXYZ", "ABCFGHJ"]
            match band_idx {
                // JKLPQRSTUXYZ
                0 => {
                    if (b'J'..=b'Z').contains(&cur_char) && !(b'M'..=b'O').contains(&cur_char) && cur_char != b'V' && cur_char != b'W' {
                        let idx = cur_char - b'J';
                        let idx = (cur_char > b'L').ternary_lazy(|| idx - 3, || idx);
                        let idx = (cur_char > b'U').ternary_lazy(|| idx - 2, || idx);
                        i32::from(idx)
                    } else {
                        -1
                    }
                }
                // ABCFGHJKLPQR
                1 => {
                    if  (b'A'..=b'R').contains(&cur_char) && 
                        cur_char != b'D' &&
                        cur_char != b'E' &&
                        cur_char != b'I' &&
                        !(b'M'..=b'O').contains(&cur_char)
                    {
                        let idx = cur_char - b'A';
                        let idx = (cur_char > b'C').ternary_lazy(|| idx - 2, || idx);
                        let idx = (cur_char > b'H').ternary_lazy(|| idx - 1, || idx);
                        let idx = (cur_char > b'L').ternary_lazy(|| idx - 3, || idx);
                        i32::from(idx)
                    } else {
                        -1
                    }
                }
                // RSTUXYZ
                2 => {
                    if  (b'R'..=b'Z').contains(&cur_char) && cur_char != b'V' && cur_char != b'W' {
                        let idx = cur_char - b'R';
                        let idx = (cur_char > b'U').ternary_lazy(|| idx - 2, || idx);
                        i32::from(idx)
                    } else {
                        -1
                    }
                }
                // ABCFGHJ
                3 => {
                    if  (b'A'..=b'J').contains(&cur_char) && 
                        cur_char != b'D' &&
                        cur_char != b'E' &&
                        cur_char != b'I'
                    {
                        let idx = cur_char - b'A';
                        let idx = (cur_char > b'C').ternary_lazy(|| idx - 2, || idx);
                        let idx = (cur_char > b'H').ternary_lazy(|| idx - 1, || idx);
                        i32::from(idx)
                    } else {
                        -1
                    }
                }
                _ => unreachable!()
            }
        };

        if col_idx == -1 {
            #[allow(clippy::cast_sign_loss)]
            let col = utmp.ternary_lazy(|| UTMCOLS[(zonem % 3) as usize], || UPSCOLS[band_idx as usize]);
            let label = if utmp { format!("zone {}", &value[..p-1]) } else { format!("UPS band {}", &value[p-1..p]) };
            return Err(Error::InvalidMgrs(format!("Column letter {} not in {label} set {col}", &value[p..=p])));
        }

        p += 1;

        let cur_char = chars[p];
        // More efficient than find()
        let mut row_idx = if utmp {
            // "ABCDEFGHJKLMNPQRSTUV"
            // First check if it's a valid latband
            if (b'A'..=b'V').contains(&cur_char) && cur_char != b'I' && cur_char != b'O' {
                // Then convert to index
                // First make it relative to A
                let idx = cur_char - b'A';
                // Decrement if it's past H to account for missing I
                let idx = (cur_char > b'H').ternary_lazy(|| idx - 1, || idx);
                // Decrement if it's past N to account for missing O
                let idx = (cur_char > b'N').ternary_lazy(|| idx - 1, || idx);
                i32::from(idx)
            } else {
                -1
            }
        } else {
            // &["ABCDEFGHJKLMNPQRSTUVWXYZ", "ABCDEFGHJKLMNP"]
            #[allow(clippy::collapsible_else_if)]
            if northp {
                if (b'A'..=b'P').contains(&cur_char) && cur_char != b'I' && cur_char != b'O' {
                    // Then convert to index
                    // First make it relative to A
                    let idx = cur_char - b'A';
                    // Decrement if it's past H to account for missing I
                    let idx = (cur_char > b'H').ternary_lazy(|| idx - 1, || idx);
                    // Decrement if it's past N to account for missing O
                    let idx = (cur_char > b'N').ternary_lazy(|| idx - 1, || idx);
                    i32::from(idx)
                } else {
                    -1
                }
            } else {
                if cur_char.is_ascii_uppercase() && cur_char != b'I' && cur_char != b'O' {
                    // Then convert to index
                    // First make it relative to A
                    let idx = cur_char - b'A';
                    // Decrement if it's past H to account for missing I
                    let idx = (cur_char > b'H').ternary_lazy(|| idx - 1, || idx);
                    // Decrement if it's past N to account for missing O
                    let idx = (cur_char > b'N').ternary_lazy(|| idx - 1, || idx);
                    i32::from(idx)
                } else {
                    -1
                }
            }
        };
        
        if row_idx == -1 {
            #[allow(clippy::cast_sign_loss)]
            let row = utmp.ternary_lazy(|| UTMROW, || UPSROWS[usize::from(northp)]);
            let northp = usize::from(northp);
            let label = if utmp { "UTM".to_string() } else { format!("UPS {}", &HEMISPHERES[northp..=northp]) };
            return Err(Error::InvalidMgrs(format!("Row letter {} not in {label} set {row}", chars[p] as char)));
        }

        p += 1;

        if utmp {
            if zonem.is_odd() {
                row_idx = (row_idx + UTM_ROW_PERIOD - UTM_EVEN_ROW_SHIFT) % UTM_ROW_PERIOD;
            }

            band_idx -= 10;

            row_idx = utm_row(band_idx, col_idx, row_idx);
            if row_idx == MAXUTM_S_ROW {
                return Err(Error::InvalidMgrs(format!("Block {} not in zone/band {}", &value[p-2..p], &value[0..p-2])))
            }

            row_idx = northp.ternary_lazy(|| row_idx, || row_idx + 100);
            col_idx += MINUTMCOL;
        }
        else {
            let eastp = band_idx.is_odd();
            col_idx += if eastp { UPSEASTING } else if northp { MINUPS_N_IND } else { MINUPS_S_IND };
            row_idx += if northp { MINUPS_N_IND } else { MINUPS_S_IND };
        }

        let precision = (len - p) / 2;
        let mut unit = 1;
        let mut x = col_idx;
        let mut y = row_idx;

        for i in 0..precision {
            unit *= BASE;
            let x_char = chars[p + i];
            let x_idx = if x_char.is_ascii_digit() {
                i32::from(x_char - b'0')
            } else {
                return Err(Error::InvalidMgrs(format!("Encountered a non-digit in {}", &value[p..])));
            };

            let y_char = chars[p + i + precision];
            let y_idx = if y_char.is_ascii_digit() {
                i32::from(y_char - b'0')
            } else {
                return Err(Error::InvalidMgrs(format!("Encountered a non-digit in {}", &value[p..])));
            };
            
            x = BASE * x + x_idx;
            y = BASE * y + y_idx;
        }

        if (len - p) % 2 == 1 {
            if !(chars[len - 1] as char).is_ascii_digit() {
                return Err(Error::InvalidMgrs(format!("Encountered a non-digit in {}", &value[p..])));
            }

            return Err(Error::InvalidMgrs(format!("Not an even number of digits in {}", &value[p..])));
        }

        if precision > MAX_PRECISION as usize {
            return Err(Error::InvalidMgrs(format!("More than {} digits in {}", 2*MAX_PRECISION, &value[p..])));
        }

        let centerp = true;
        if centerp {
            unit *= 2;
            x = 2 * x + 1;
            y = 2 * y + 1;
        }

        let x = (f64::from(TILE) * f64::from(x)) / f64::from(unit);
        let y = (f64::from(TILE) * f64::from(y)) / f64::from(unit);

        Ok(Self {
            utm: UtmUps::new(
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

pub(crate) fn check_coords(utmp: bool, northp: bool, x: f64, y: f64) -> Result<(bool, f64, f64), Error> {
    lazy_static! {
        static ref ANG_EPS: f64 = 1_f64 * 2_f64.powi(-(f64::DIGITS as i32 - 25));
    }

    let x_int = (x / f64::from(TILE)).floor() as i32;
    let y_int = (y / f64::from(TILE)).floor() as i32;
    let ind = utmp.ternary(2, 0) + northp.ternary(1, 0);

    let mut x_new = x;
    let mut y_new = y;

    if !(MIN_EASTING[ind]..MAX_EASTING[ind]).contains(&x_int) {
        if x_int == MAX_EASTING[ind] && x.eps_eq(f64::from(MAX_EASTING[ind] * TILE)) {
            x_new -= *ANG_EPS;
        } else {
            return Err(Error::InvalidMgrs(
                format!(
                    "Easting {:.2}km not in MGRS/{} range for {} hemisphere [{:.2}km, {:.2}km]",
                    x / 1000.0,
                    utmp.ternary("UTM", "UPS"),
                    northp.ternary("N", "S"),
                    MIN_EASTING[ind] * (TILE / 1000),
                    MAX_EASTING[ind] * (TILE / 1000),
                )
            ));
        }
    }

    if !(MIN_NORTHING[ind]..MAX_NORTHING[ind]).contains(&y_int) {
        if y_int == MAX_NORTHING[ind] && y.eps_eq(f64::from(MAX_NORTHING[ind] * TILE)) {
            y_new -= *ANG_EPS;
        } else {
            return Err(Error::InvalidMgrs(
                format!(
                    "Northing {:.2}km not in MGRS/{} range for {} hemisphere [{:.2}km, {:.2}km]",
                    y / 1000.0,
                    utmp.ternary("UTM", "UPS"),
                    northp.ternary("N", "S"),
                    MIN_NORTHING[ind] * (TILE / 1000),
                    MAX_NORTHING[ind] * (TILE / 1000),
                )
            ));
        }
    }

    let (northp_new, y_new) = if utmp {
        if northp && y_int < MINUTM_S_ROW {
            (false, y_new + f64::from(UTM_N_SHIFT))
        } else if !northp && y_int >= MAXUTM_S_ROW {
            if y.eps_eq(f64::from(MAXUTM_S_ROW * TILE)) {
                (northp, y_new - *ANG_EPS)
            } else {
                (true, y - f64::from(UTM_N_SHIFT))
            }
        } else {
            (northp, y_new)
        }
    } else {
        (northp, y_new)
    };

    Ok((northp_new, x_new, y_new))
}

impl Display for Mgrs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        lazy_static! {
            static ref ANG_EPS: f64 = 1_f64 * 2_f64.powi(-(f64::MANTISSA_DIGITS as i32 - 7));
        }

        let lat = if self.utm.zone > 0 {
            // Does a rough estimate for latitude determine the latitude band?
            let y_est = self.utm.northp.ternary_lazy(|| self.utm.northing, || self.utm.northing - f64::from(UTM_N_SHIFT));
            // A cheap calculation of the latitude which results in an "allowed"
            // latitude band would be
            //   lat = ApproxLatitudeBand(ys) * 8 + 4;
            //
            // Here we do a more careful job using the band letter corresponding to
            // the actual latitude.
            let y_est = y_est / f64::from(TILE);
            if y_est.abs() < 1.0 {
                0.9 * y_est
            }
            else {
                let pole_add = (y_est > 0.0).ternary(1.0, -1.0);
                let lat_poleward = 0.901 * y_est + pole_add * 0.135;
                let lat_eastward = 0.902 * y_est * (1.0 - 1.85e-6 * y_est.powi(2));

                if to_latitude_band(lat_poleward) == to_latitude_band(lat_eastward) {
                    lat_poleward
                } else {
                    let coord = UtmUps::new(self.utm.zone, self.utm.northp, self.utm.easting, self.utm.northing).to_latlon();
                    coord.latitude
                }
            }
        } else {
            0.
        };
        
        // Other Forward call
        let utmp = self.utm.zone != 0;
        let (northp, easting, northing) = check_coords(utmp, self.utm.northp, self.utm.easting, self.utm.northing)
            .expect("Invalid coords; please report this to the library author");
        // Create pre-allocated string of the correct length
        let mut mgrs_str = [0u8; 2 + 3 + 2*MAX_PRECISION as usize];
        let zone = self.utm.zone - 1;
        let mut z: usize = utmp.ternary(2, 0);

        let digits = DIGITS.as_bytes();

        #[allow(clippy::cast_sign_loss)]
        if utmp {
            mgrs_str[0] = digits[(self.utm.zone / BASE) as usize];
            mgrs_str[1] = digits[(self.utm.zone % BASE) as usize];
        }

        let xx = easting * f64::from(MULT);
        let yy = northing * f64::from(MULT);

        let ix = xx.floor() as i64;
        let iy = yy.floor() as i64;
        let m = i64::from(MULT) * i64::from(TILE);

        let xh = (ix / m) as i32;
        let yh = (iy / m) as i32;

        #[allow(clippy::cast_sign_loss)]
        if utmp {
            // Correct fuzziness in latitude near equator
            let band_idx = (lat.abs() < *ANG_EPS).ternary_lazy(|| northp.ternary(0, -1), || to_latitude_band(lat));
            let col_idx = xh - MINUTMCOL;
            let row_idx = utm_row(band_idx, col_idx, yh % UTM_ROW_PERIOD);

            assert!(
                row_idx == yh - northp.ternary(MINUTM_N_ROW, MAXUTM_S_ROW),
                "Latitude is inconsistent with UTM; this should not occur."
            );

            mgrs_str[z] = LATBAND.as_bytes()[(10 + band_idx) as usize];
            z += 1;
            mgrs_str[z] = UTMCOLS[(zone % 3) as usize].as_bytes()[col_idx as usize];
            z += 1;
            let idx = (yh + zone.is_odd().ternary(UTM_EVEN_ROW_SHIFT, 0)) % UTM_ROW_PERIOD;
            mgrs_str[z] = UTMROW.as_bytes()[idx as usize];
            z += 1;
        } else {
            let eastp = xh >= UPSEASTING;
            let band_idx: usize = northp.ternary(2, 0) + eastp.ternary(1, 0);
            mgrs_str[z] = UPSBAND.as_bytes()[band_idx];
            z += 1;
            let idx = xh - eastp.ternary(UPSEASTING, northp.ternary(MINUPS_N_IND, MINUPS_S_IND));
            mgrs_str[z] = UPSCOLS[band_idx].as_bytes()[idx as usize];
            z += 1;
            let idx = yh - northp.ternary(MINUPS_N_IND, MINUPS_S_IND);
            mgrs_str[z] = UPSROWS[usize::from(northp)].as_bytes()[idx as usize];
            z += 1;
        }

        if self.precision > 0 {
            let mut ix = ix - m * i64::from(xh);
            let mut iy = iy - m * i64::from(yh);
            #[allow(clippy::cast_sign_loss)]
            let d = i64::from(BASE).pow((MAX_PRECISION - self.precision) as u32);
            ix /= d;
            iy /= d;

            #[allow(clippy::cast_sign_loss)]
            for c in (0..self.precision as usize).rev() {
                mgrs_str[z + c] = digits[(ix % i64::from(BASE)) as usize];
                ix /= i64::from(BASE);
                mgrs_str[z + c + self.precision as usize] = digits[(iy % i64::from(BASE)) as usize];
                iy /= i64::from(BASE);
            }
        }

        write!(f, "{}", String::from_utf8_lossy(&mgrs_str).trim_end_matches('\0'))
    }
}
