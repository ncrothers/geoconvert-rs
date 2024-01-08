# geoconvert

[![Latest Version](https://img.shields.io/crates/v/geoconvert.svg)](crates.io])
[![Docs](https://docs.rs/geoconvert/badge.svg)](docs.rs])

[crates.io]: https://crates.io/crates/geoconvert
[docs.rs]: https://docs.rs/geoconvert

`geoconvert` is a lightweight library for converting between different
geographic coordinate systems. Currently, there are three coordinate systems implemented:

* `LatLon`
* `UtmUps`
* `Mgrs`

The implementation of this library is a translation of a subset of 
[GeographicLib](https://geographiclib.sourceforge.io/C++/doc/index.html) from C++ to Rust. Specifically, `geoconvert`
implements some of the functionality of the [GeoConvert](https://geographiclib.sourceforge.io/C++/doc/GeoConvert.1.html) 
command line tool.

## Usage

You can create coordinates manually using a struct's `create()` function, then convert to other
types using the `to_*`/`from_*` methods.

```rust
use geoconvert::{LatLon, Mgrs, UtmUps};

// This returns a result. When calling `create()`, the arguments are validated to ensure only a valid
// coordinate gets created.
let coord = LatLon::create(40.748333, -73.985278).unwrap();

// Convert to a UTM/UPS coordinate
let coord_utm = coord.to_utmups();
let coord_utm = UtmUps::from_latlon(&coord);

// Convert to an MGRS coordinate
// Note that for MGRS you must specify the precision
let coord_mgrs = coord.to_mgrs(6);
let coord_mgrs = Mgrs::from_latlon(&coord, 6);
```

## Features

If you want `serde` compatibility with `Serialize`/`Deserialize`, activate the `serde` feature.

## Testing Accuracy

To test the accuracy compared to GeographicLib yourself, you'll need a dataset of lat/lon and MGRS points. I have a [gist](https://gist.github.com/ncrothers/0fc036c89cef307caa399347cda6c3f8) that contains a sample dataset of ~100K points generated using Python and converted using [GeoConvert](https://geographiclib.sourceforge.io/C++/doc/GeoConvert.1.html). If you'd like to generate your own dataset to validate the accuracy, create files named `mgrs.txt` and `latlon.txt`, where `mgrs.txt` is a list of MGRS coordinates (one per line) and `latlon.txt` is a list of latitude longitude pairs, each pair internally delimited by a space (i.e. like `<latitude> <longitude>`). You can use GeoConvert to do the conversion, or use a different source for ground truth.

Once you have these files, place them into `tests/` and run:

```bash
cargo test -- --include-ignored
```

The test will fail if the `geoconvert` conversion differs from the ground truth by a distance of 1mm or more (calculated using the haversine formula). It also prints the average distance error.