use pansy_meteor_consensus::{split, Config, Point};
fn track(r: f64, v: f64, start: u64, n: usize) -> Vec<Point> {
    (0..n)
        .map(|i| Point {
            sample_idx: start + i as u64 * 8000,
            range_km: r + v * i as f64 * 0.008 + 0.03 * (i as f64).sin(),
            radial_velocity_km_s: v + 0.1 * (i as f64).cos(),
            snr_linear: 20.,
        })
        .collect()
}
#[test]
fn crossing_meteors_and_bright_outliers_remain_separate() {
    let mut points = track(101., -20., 1_746_000_000_000_000, 30);
    points.extend(track(103., -40., 1_746_000_000_000_000, 30));
    for i in 0..10 {
        points.push(Point {
            sample_idx: points[i].sample_idx,
            range_km: 140. + i as f64,
            radial_velocity_km_s: 50.,
            snr_linear: 1000.,
        });
    }
    let t = split(&points, &Config::default()).unwrap();
    assert_eq!(t.len(), 2);
    assert_eq!(t.iter().map(|x| x.members.len()).sum::<usize>(), 60);
    for x in t {
        assert!(
            x.members.iter().all(|&i| i < 30) || x.members.iter().all(|&i| (30..60).contains(&i))
        );
    }
}
#[test]
fn weak_edges_recovered_without_seeding_noise_and_low_radial_speed_allowed() {
    let mut p = track(100., -0.5, 1_746_000_000_000_000, 30);
    for (i, x) in p.iter_mut().enumerate() {
        if !(8..22).contains(&i) {
            x.snr_linear = 3.;
        }
    }
    let t = split(&p, &Config::default()).unwrap();
    assert_eq!(t.len(), 1);
    assert_eq!(t[0].members.len(), 30);
    for x in &mut p {
        x.snr_linear = 3.;
    }
    assert!(split(&p, &Config::default()).unwrap().is_empty());
}
#[test]
fn separated_events_do_not_bridge_gaps_and_no_duplicate_pulses() {
    let mut p = track(100., -20., 0, 15);
    p.extend(track(90., -20., 500_000, 15));
    p.extend(p.clone());
    let t = split(&p, &Config::default()).unwrap();
    assert_eq!(t.len(), 2);
    for x in &t {
        let first = p[x.members[0]].sample_idx;
        let last = p[*x.members.last().unwrap()].sample_idx;
        assert!(last - first < 200_000);
        assert_eq!(x.members.len(), 15);
    }
    assert!(!t.is_empty());
}
#[test]
fn nonfinite_and_overloaded_windows_are_explicit() {
    let mut p = track(100., -20., 0, 15);
    p[0].range_km = f64::NAN;
    let t = split(&p, &Config::default()).unwrap();
    assert_eq!(t[0].members.len(), 14);
    let c = Config {
        max_points: 10,
        ..Config::default()
    };
    assert!(split(&p, &c).is_err());
}

#[test]
fn changing_radial_velocity_does_not_fragment_one_meteor() {
    let mut p = track(100., -25., 1_746_000_000_000_000, 60);
    for (i, x) in p.iter_mut().enumerate() {
        let t = i as f64 * 0.008;
        x.range_km += 0.5 * 35. * t * t;
        x.radial_velocity_km_s += 35. * t;
    }
    let tracks = split(&p, &Config::default()).unwrap();
    assert_eq!(tracks.len(), 1);
    assert_eq!(tracks[0].members.len(), 60);
    assert!((tracks[0].radial_acceleration_km_s2 - 35.).abs() < 0.5);
}
