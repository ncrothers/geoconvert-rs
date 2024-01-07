use std::fs;

use geoconvert::{LatLon, Mgrs};

#[test]
#[ignore]
fn conversion_accuracy() {
    let mgrs_points = fs::read_to_string("./tests/mgrs.txt")
        .unwrap();
    let mgrs_points = mgrs_points
        .lines()
        .map(|line| line.trim());
    let latlon_points = fs::read_to_string("./tests/latlon.txt")
        .unwrap();
    let latlon_points = latlon_points
        .lines()
        .map(|line| {
            let mut pieces = line.split(' ');
            LatLon::create(
                pieces.next().unwrap().parse().unwrap(),
                pieces.next().unwrap().parse().unwrap()
            ).unwrap()
        });

    let errors = mgrs_points
        .zip(latlon_points)
        .map(|(mgrs, latlon)| {
            let val = Mgrs::parse_str(mgrs).unwrap();
            let coord = val.to_latlon();
    
            coord.haversine(&latlon)
        });

    // Check if any differences between GeographicLib and ours exceeds 1mm
    let significant_errors = errors
        .clone()
        .filter(|dist| *dist > 1e-3);

    let count = errors.clone().count();
    let sum: f64 = errors.clone().sum();

    println!("Average error: {}", sum / count as f64);

    assert_eq!(significant_errors.count(), 0);
}
