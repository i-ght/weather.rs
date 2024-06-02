use chrono::{DateTime, Utc};
use std::{
    env::{self},
    error::Error,
};
use weather;

struct Args {
    lat: f64,
    lon: f64,
}

fn args_of_env() -> Result<Args, Box<dyn Error>> {
    let lat = env::var("LAT")?.parse()?;
    let lon = env::var("LON")?.parse()?;
    Ok(Args { lat, lon })
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = args_of_env()?;
    let date_time: DateTime<Utc> = Utc::now(); // Use current date and time as an example

    // Calculate astronomical events
    let (sunrise_time, solar_noon_time, sunset_time) =
        weather::events(args.lat, args.lon, date_time);
        
    println!("{}", sunrise_time.with_timezone(&chrono::Local));
    println!("{}", solar_noon_time.with_timezone(&chrono::Local));
    println!("{}", sunset_time.with_timezone(&chrono::Local));
    
    Ok(())
}
