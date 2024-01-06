use std::{time::Instant, fs};

use geoconvert::{mgrs::Mgrs, ParseCoord, utm::Utm, latlon::LatLon};
use geomorph::coord::Coord;
use polars::prelude::*;

fn test_mgrs_accuracy() {
    let mgrs_points = fs::read_to_string("./helpers/output_mgrs.txt")
        .unwrap();
    let mgrs_points = mgrs_points
        .lines()
        .map(|line| line.trim());
    let latlon_points = fs::read_to_string("./helpers/output_latlon.txt")
        .unwrap();
    let latlon_points = latlon_points
        .lines()
        .map(|line| {
            let mut pieces = line.split(' ');
            LatLon {
                latitude: pieces.next().unwrap().parse().unwrap(),
                longitude: pieces.next().unwrap().parse().unwrap(),
            }
        });

    let mgrs_vec = mgrs_points.clone().collect::<Vec<_>>();
    let true_latlon = latlon_points.clone().map(|x| x.to_string()).collect::<Vec<_>>();

    let count = mgrs_points.clone().count();

    let preds = mgrs_points
        .clone()
        .map(|mgrs | {
            let val = Mgrs::parse_coord(mgrs).unwrap();
            LatLon::try_from(val).unwrap().to_string()
        })
        .collect::<Vec<_>>();

    let errors = mgrs_points
        .zip(latlon_points)
        .map(|(mgrs, latlon)| {
            let val = Mgrs::parse_coord(mgrs).unwrap();
            let coord = LatLon::try_from(val).unwrap();
    
            coord.haversine(&latlon)
        })
        .collect::<Vec<_>>();

    let source = Series::new("mgrs", &mgrs_vec);
    let pred = Series::new("pred", &preds);
    let trues = Series::new("true", &true_latlon);
    let distance = Series::new("distance", &errors);

    let df = DataFrame::new(vec![source, pred, trues, distance]).unwrap();
    let df = df.filter(&df.column("distance").unwrap().gt(-1.).unwrap()).unwrap();
    
    let df = df.sort(["distance"], false, false).unwrap();
    
    println!("{:?}", df);

    // println!("Median: {}", df.column("distance").median().unwrap());
    // println!("Min: {}, Max: {}", errors.min::<f64>().unwrap().unwrap(), errors.max::<f64>().unwrap().unwrap());
    // println!("Mean: {}", errors.mean().unwrap());
}

fn test_mgrs_parsing_accuracy() {
    let mgrs_points = fs::read_to_string("./helpers/output_mgrs.txt")
        .unwrap();
    let mgrs_points = mgrs_points
        .lines()
        .map(|line| line.trim());

    let n_incorrect = mgrs_points
        .fold(0, |sum, mgrs| {
            let val = Mgrs::parse_coord(mgrs).unwrap();
            let pred = val.to_string();
            if mgrs != pred {
                sum + 1
            }
            else {
                sum
            }
        });

    println!("Number of incorrect parsings: {}", n_incorrect);
}

fn test_mgrs() {
    let mgrs_str = "25XEN0416386465";
    let val = Mgrs::parse_coord(mgrs_str).unwrap();
    
    println!("Parsed value: {val:?}");
    println!("Original MGRS: {mgrs_str}");
    println!("Printed value: {val}");

    let utm_test = Utm::new(24, true, 294257., 4522751.);
    let geo_utm = geomorph::utm::Utm::new(294257., 4522751., true, 24, 'A', false);

    println!("Converted coords: {:?}", LatLon::try_from(utm_test).unwrap());
    println!("Converted coords: {:?}", Coord::from(geo_utm));
}

fn main() {
    test_mgrs_accuracy();
    test_mgrs_parsing_accuracy();
}
