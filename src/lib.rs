use chrono::{DateTime, Datelike, Duration, FixedOffset, Timelike, Utc};
use std::f64::consts::PI;

pub enum ZodiacSign {
    Aries,
    Taurus,
    Gemini,
    Cancer,
    Leo,
    Virgo,
    Libra,
    Scorpio,
    Sagittarius,
    Capricorn,
    Aquarius,
    Pisces,
}

impl ZodiacSign {
    fn value(&self) -> f64 {
        match self {
            ZodiacSign::Aries => 0.0,
            ZodiacSign::Taurus => 30.0,
            ZodiacSign::Gemini => 60.0,
            ZodiacSign::Cancer => 90.0,
            ZodiacSign::Leo => 120.0,
            ZodiacSign::Virgo => 150.0,
            ZodiacSign::Libra => 180.0,
            ZodiacSign::Scorpio => 210.0,
            ZodiacSign::Sagittarius => 240.0,
            ZodiacSign::Capricorn => 270.0,
            ZodiacSign::Aquarius => 300.0,
            ZodiacSign::Pisces => 330.0,
        }
    }
}
pub fn modulo(a: f64, b: f64) -> f64 {
    (a % b + b) % b
}

fn hour_time_to_decimal(date: DateTime<Utc>) -> f64 {
    let hour = date.hour() as f64;
    let minute = date.minute() as f64;
    hour + minute / 60.0
}

fn rad_to_deg(angle_rad: f64) -> f64 {
    180.0 * angle_rad / PI
}

fn deg_to_rad(angle_deg: f64) -> f64 {
    PI * angle_deg / 180.0
}

fn calc_geom_mean_long_sun(t: f64) -> f64 {
    let mut l0 = 280.46646 + t * (36000.76983 + t * 0.0003032);
    while l0 >= 360.0 {
        l0 -= 360.0;
    }
    while l0 < 0.0 {
        l0 += 360.0;
    }
    l0
}

fn calc_eccentricity_earth_orbit(t: f64) -> f64 {
    0.016708634 - t * (0.000042037 + 0.0000001267 * t)
}

fn calc_geom_mean_anomaly_sun(t: f64) -> f64 {
    357.52911 + t * (35999.05029 - 0.0001537 * t)
}

fn calc_mean_obliquity_of_ecliptic(t: f64) -> f64 {
    let seconds = 21.448 - t * (46.8150 + t * (0.00059 - t * 0.001813));
    23.0 + (26.0 + seconds / 60.0) / 60.0
}

fn calc_obliquity_correction(t: f64) -> f64 {
    let e0 = calc_mean_obliquity_of_ecliptic(t);
    let omega = 125.04 - 1934.136 * t;
    e0 + 0.00256 * deg_to_rad(omega).cos()
}

fn calc_sun_eq_of_center(t: f64) -> f64 {
    let m = calc_geom_mean_anomaly_sun(t);
    let mrad = deg_to_rad(m);
    let sinm = mrad.sin();
    let sin2m = (mrad + mrad).sin();
    let sin3m = (mrad + mrad + mrad).sin();
    sinm * (1.914602 - t * (0.004817 + 0.000014 * t))
        + sin2m * (0.019993 - 0.000101 * t)
        + sin3m * 0.000289
}

fn calc_sun_true_long(t: f64) -> f64 {
    calc_geom_mean_long_sun(t) + calc_sun_eq_of_center(t)
}

fn calc_sun_apparent_long(t: f64) -> f64 {
    let o = calc_sun_true_long(t);
    let omega = 125.04 - 1934.136 * t;
    o - 0.00569 - 0.00478 * deg_to_rad(omega).sin()
}

/* fn calc_sun_true_anomaly(t: f64) -> f64 {
    calc_geom_mean_anomaly_sun(t) + calc_sun_eq_of_center(t)
}
 */
fn calc_sun_declination(t: f64) -> f64 {
    let e = calc_obliquity_correction(t);
    let lambda = calc_sun_apparent_long(t);
    let sint = deg_to_rad(e).sin() * deg_to_rad(lambda).sin();
    rad_to_deg(sint.asin())
}
/* 
fn calc_sun_rad_vector(t: f64) -> f64 {
    let v = calc_sun_true_anomaly(t);
    let e = calc_eccentricity_earth_orbit(t);
    (1.000001018 * (1.0 - e * e)) / (1.0 + e * deg_to_rad(v).cos())
} */

fn calc_refraction(elev: f64) -> f64 {
    if elev > 85.0 {
        0.0
    } else {
        let te = deg_to_rad(elev).tan();
        if elev > 5.0 {
            58.1 / te - 0.07 / (te * te * te) + 0.000086 / (te * te * te * te * te)
        } else if elev > -0.575 {
            1735.0 + elev * (-518.2 + elev * (103.4 + elev * (-12.79 + elev * 0.711)))
        } else {
            -20.774 / te / 3600.0
        }
    }
}

fn calc_equation_of_time(t: f64) -> f64 {
    let epsilon = calc_obliquity_correction(t);
    let l0 = calc_geom_mean_long_sun(t);
    let e = calc_eccentricity_earth_orbit(t);
    let m = calc_geom_mean_anomaly_sun(t);

    let y = (deg_to_rad(epsilon) / 2.0).sin().powi(2);

    let sin2l0 = (2.0 * deg_to_rad(l0)).sin();
    let sinm = deg_to_rad(m).sin();
    let cos2l0 = (2.0 * deg_to_rad(l0)).cos();
    let sin4l0 = (4.0 * deg_to_rad(l0)).sin();
    let sin2m = (2.0 * deg_to_rad(m)).sin();

    let etime = y * sin2l0 - 2.0 * e * sinm + 4.0 * e * y * sinm * cos2l0
        - 0.5 * y * y * sin4l0
        - 1.25 * e * e * sin2m;
    rad_to_deg(etime) * 4.0
}

fn calc_az_el(t: f64, localtime: f64, latitude: f64, longitude: f64, zone: f64) -> (f64, f64) {
    let eq_time = calc_equation_of_time(t);
    let theta = calc_sun_declination(t);

    let solar_time_fix = eq_time + 4.0 * longitude - 60.0 * zone;
    let true_solar_time = (localtime + solar_time_fix) % 1440.0;

    let hour_angle = true_solar_time / 4.0 - 180.0;
    let ha_rad = deg_to_rad(hour_angle);

    let csz = deg_to_rad(latitude).sin() * deg_to_rad(theta).sin()
        + deg_to_rad(latitude).cos() * deg_to_rad(theta).cos() * ha_rad.cos();

    let zenith = rad_to_deg(csz.acos().min(1.0).max(-1.0));

    let az_denom = deg_to_rad(latitude).cos() * deg_to_rad(zenith).sin();

    let azimuth = if az_denom.abs() > 0.001 {
        let az_rad = ((deg_to_rad(latitude).sin() * deg_to_rad(zenith).cos())
            - deg_to_rad(theta).sin())
            / az_denom;
        let temp_azimuth = 180.0 - rad_to_deg(az_rad.min(1.0).max(-1.0).acos());
        if hour_angle > 0.0 {
            -temp_azimuth
        } else {
            temp_azimuth
        }
    } else {
        if latitude > 0.0 {
            180.0
        } else {
            0.0
        }
    };

    let azimuth = (azimuth + 360.0) % 360.0;

    let exoatm_elevation = 90.0 - zenith;

    let refraction_correction = calc_refraction(exoatm_elevation);

    let solar_zen = zenith - refraction_correction;
    let elevation = 90.0 - solar_zen;

    (azimuth, elevation)
}

fn calc_time_julian_cent(julian_day: f64) -> f64 {
    (julian_day - 2451545.0) / 36525.0
}

/* fn time_local(date_time_offset: DateTime<Utc>) -> f64 {
    let hour = date_time_offset.hour() as f64;
    let minute = date_time_offset.minute() as f64;
    let second = date_time_offset.second() as f64;
    hour * 60.0 + minute + second / 60.0
} */

fn calc_hour_angle_sunrise(lat: f64, solar_dec: f64) -> f64 {
    let lat_rad = deg_to_rad(lat);
    let sd_rad = deg_to_rad(solar_dec);
    let ha_arg =
        deg_to_rad(90.833).cos() / (lat_rad.cos() * sd_rad.cos()) - lat_rad.tan() * sd_rad.tan();
    ha_arg.acos()
}

fn calc_sunrise_set_utc(rise: bool, jd: f64, latitude: f64, longitude: f64) -> f64 {
    let t = calc_time_julian_cent(jd);
    let eq_time = calc_equation_of_time(t);
    let solar_dec = calc_sun_declination(t);
    let hour_angle = calc_hour_angle_sunrise(latitude, solar_dec);
    let hour_angle = if rise { hour_angle } else { -hour_angle };
    let delta = longitude + rad_to_deg(hour_angle);
    720.0 - (4.0 * delta) - eq_time
}

/* fn minutes_to_date_time_offset(date: DateTime<Utc>, minutes: f64) -> DateTime<Utc> {
    date + Duration::minutes(minutes as i64)
} */

fn tz_offset(date: DateTime<FixedOffset>) -> f64 {
    let utc_offset = date.offset().local_minus_utc();
    let hours = utc_offset / 3600;
    let minutes = (utc_offset % 3600) / 60;
    hours as f64 + minutes as f64 / 60.0
}

fn calc_sunrise_set(
    rise: bool,
    jd: f64,
    latitude: f64,
    longitude: f64,
    timezone: f64,
) -> (f64, f64, f64) {
    let time_utc = calc_sunrise_set_utc(rise, jd, latitude, longitude);
    let new_time_utc = calc_sunrise_set_utc(rise, jd + time_utc / 1440.0, latitude, longitude);

    let mut time_local = new_time_utc + timezone * 60.0;
    let mut jday = jd;

    while time_local < 0.0 {
        time_local += 1440.0;
        jday -= 1.0;
    }

    while time_local >= 1440.0 {
        time_local -= 1440.0;
        jday += 1.0;
    }

    let rise_t = calc_time_julian_cent(jday + time_local / 1440.0);
    let (azimuth, _elevation) = calc_az_el(rise_t, time_local, latitude, longitude, timezone);

    (azimuth, time_local, jday)
}

fn calc_sol_noon(jd: f64, longitude: f64, timezone: f64) -> f64 {
    let tnoon = calc_time_julian_cent(jd - longitude / 360.0);
    let eq_time = calc_equation_of_time(tnoon);
    let sol_noon_offset = 720.0 - (longitude * 4.0) - eq_time;
    let newt = calc_time_julian_cent(jd - 0.5 + sol_noon_offset / 1440.0);
    let eq_time = calc_equation_of_time(newt);
    let mut sol_noon_local = 720.0 - (longitude * 4.0) - eq_time + (timezone * 60.0);

    while sol_noon_local < 0.0 {
        sol_noon_local += 1440.0;
    }

    while sol_noon_local >= 1440.0 {
        sol_noon_local -= 1440.0;
    }

    sol_noon_local
}

fn jd_precise(date: DateTime<Utc>) -> f64 {
    let utc = date;
    let ut = hour_time_to_decimal(utc);
    let year = utc.year() as f64;
    let month = utc.month() as f64;
    let day = utc.day() as f64;

    367.0 * year
        - ((7.0 * (year + ((month + 9.0) / 12.0).floor())) / 4.0).floor()
        - ((3.0 * ((year + (month - 9.0) / 7.0).floor() / 100.0 + 1.0)) / 4.0).floor()
        + ((275.0 * month) / 9.0).floor()
        + day
        + 1721028.5
        + (ut / 24.0)
}

fn jd(date: DateTime<Utc>) -> f64 {
    let (year, month, day) = if date.month() <= 2 {
        (date.year() - 1, date.month() + 12, date.day())
    } else {
        (date.year(), date.month(), date.day())
    };

    let year = year as f64;
    let month = month as f64;
    let day = day as f64;

    let a = (year / 100.0).floor();
    let b = 2.0 - a + (a / 4.0).floor();

    (365.25 * (year + 4716.0)).floor() + (30.6001 * (month + 1.0)).floor() + day + b - 1524.5
}

fn sin_from_degrees(degrees: f64) -> f64 {
    deg_to_rad(degrees).sin()
}

fn cos_from_degrees(degrees: f64) -> f64 {
    deg_to_rad(degrees).cos()
}

fn tan_from_degrees(degrees: f64) -> f64 {
    deg_to_rad(degrees).tan()
}

const OBLIQUITY_ECLIPTIC: f64 = 23.4367;
const JULIAN_DAYS_JAN_1_2000: f64 = 2451545.0;

fn local_sidereal_time(date: DateTime<Utc>, longitude: f64) -> f64 {
    let jd = jd_precise(date);
    let julian_days_since_2000 = jd - JULIAN_DAYS_JAN_1_2000;
    let t_factor = julian_days_since_2000 / 36525.0; // centuries
    let degrees_rotation_in_sidereal_day = 360.98564736629;
    let lst = 280.46061837
        + degrees_rotation_in_sidereal_day * julian_days_since_2000
        + 0.000387933 * (t_factor.powi(2))
        - (t_factor.powi(3) / 38710000.0)
        + longitude;
    modulo(lst, 360.0)
}

pub fn ascendant(latitude: f64, longitude: f64, date: DateTime<Utc>) -> ZodiacSign {
    let local_sidereal_time = local_sidereal_time(date, longitude);
    let obliquity_ecliptic = OBLIQUITY_ECLIPTIC;

    let a = -cos_from_degrees(local_sidereal_time);
    let b = sin_from_degrees(obliquity_ecliptic) * tan_from_degrees(latitude);
    let c = cos_from_degrees(obliquity_ecliptic) * sin_from_degrees(local_sidereal_time);
    let d = b + c;
    let e = a / d;
    let f = e.atan();

    let mut ascendant = f.to_degrees();

    // Modulation from wikipedia
    // https://en.wikipedia.org/wiki/Ascendant
    // Citation Peter Duffett-Smith, Jonathan Zwart, Practical astronomy with your calculator or spreadsheet-4th ed., p47, 2011

    if d < 0.0 {
        ascendant += 180.0;
    } else {
        ascendant += 360.0;
    }

    if ascendant >= 180.0 {
        ascendant -= 180.0;
    } else {
        ascendant += 180.0;
    }

    ascendant = modulo(ascendant, 360.0);

    match ascendant {
        a if a >= ZodiacSign::Pisces.value() => ZodiacSign::Pisces,
        a if a >= ZodiacSign::Aquarius.value() => ZodiacSign::Aquarius,
        a if a >= ZodiacSign::Capricorn.value() => ZodiacSign::Capricorn,
        a if a >= ZodiacSign::Sagittarius.value() => ZodiacSign::Sagittarius,
        a if a >= ZodiacSign::Scorpio.value() => ZodiacSign::Scorpio,
        a if a >= ZodiacSign::Libra.value() => ZodiacSign::Libra,
        a if a >= ZodiacSign::Virgo.value() => ZodiacSign::Virgo,
        a if a >= ZodiacSign::Leo.value() => ZodiacSign::Leo,
        a if a >= ZodiacSign::Cancer.value() => ZodiacSign::Cancer,
        a if a >= ZodiacSign::Gemini.value() => ZodiacSign::Gemini,
        a if a >= ZodiacSign::Taurus.value() => ZodiacSign::Taurus,
        _ => ZodiacSign::Aries,
    }
}

#[allow(deprecated)]
fn minutes_to_datetimeoffset(date: DateTime<Utc>, minutes: f64) -> DateTime<Utc> {
    let dt = date.date();
    let duration = Duration::minutes(minutes as i64);
    dt.and_hms(0, 0, 0) + duration
}

pub fn events(
    lat: f64,
    lon: f64,
    date: DateTime<Utc>,
) -> (DateTime<Utc>, DateTime<Utc>, DateTime<Utc>) {
    let offset = tz_offset(date.fixed_offset());
    let jd = jd(date);
    let (_, rise_time_local, _rise_azimuth) = calc_sunrise_set(true, jd, lat, lon, offset);
    let (_, set_time_local, _set_azimuth) = calc_sunrise_set(false, jd, lat, lon, offset);
    let sol_noon = calc_sol_noon(jd, lon, offset);

    (
        minutes_to_datetimeoffset(date, rise_time_local),
        minutes_to_datetimeoffset(date, sol_noon),
        minutes_to_datetimeoffset(date, set_time_local),
    )
}
