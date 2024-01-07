// Semi-major axis a
pub(crate) const WGS84_A: f64 = 6_378_137.;
// Flattening
#[allow(clippy::unreadable_literal)]
pub(crate) const WGS84_F: f64 = 1.0 / 298.257223563;

// UTM central scale factor
pub(crate) const UTM_K0: f64 = 9996.0 / 10_000.;
// UPS central scale factor
pub(crate) const UPS_K0: f64 = 994.0 / 1000.;