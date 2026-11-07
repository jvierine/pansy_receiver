use pansy_meteor_consensus::{phase::split_phase, Config, Point};
fn signature(angle: f64, common: f64) -> [f32; 14] {
    let mut s = [0.; 14];
    for ch in 0..7 {
        s[2 * ch] = (ch as f64 * angle + common).cos() as f32;
        s[2 * ch + 1] = (ch as f64 * angle + common).sin() as f32;
    }
    s
}
#[test]
fn identical_range_doppler_different_interferometry_separates() {
    let mut p = vec![];
    let mut s = vec![];
    for i in 0..50 {
        p.push(Point {
            sample_idx: 1_746_000_000_000_000 + i * 8000,
            range_km: 100. - 20. * i as f64 * 0.008,
            radial_velocity_km_s: -20.,
            snr_linear: 20.,
        });
        s.push(signature(
            if i % 2 == 0 { 0.1 } else { 1.5 },
            i as f64 * 2.1,
        ));
    }
    let tracks = split_phase(&p, &s, 0.65, 0.017, &Config::default()).unwrap();
    assert_eq!(tracks.len(), 2);
    for t in tracks {
        assert_eq!(t.members.len(), 25);
        assert!(t.members.iter().all(|&i| i % 2 == t.members[0] % 2));
    }
}
#[test]
fn motion_common_phase_and_missing_echoes_are_supported() {
    let mut p = vec![];
    let mut s = vec![];
    for i in 0..50 {
        if i % 9 < 3 {
            continue;
        }
        p.push(Point {
            sample_idx: 1_746_000_000_000_000 + i * 8000,
            range_km: 100. - 20. * i as f64 * 0.008,
            radial_velocity_km_s: -20.,
            snr_linear: 20.,
        });
        s.push(signature(0.1 + 0.01 * i as f64, i as f64 * 2.9));
    }
    let tracks = split_phase(&p, &s, 0.65, 0.05, &Config::default()).unwrap();
    assert_eq!(tracks.len(), 1);
    assert_eq!(tracks[0].members.len(), p.len());
    s[0][0] = f32::NAN;
    assert!(split_phase(&p, &s, 0.65, 0.017, &Config::default()).is_err());
}

#[test]
fn one_intermediate_phase_outlier_cannot_join_simultaneous_tracks() {
    let mut p = Vec::new();
    let mut s = Vec::new();
    for i in 0..60 {
        p.push(Point {
            sample_idx: 1_746_000_000_000_000 + i * 1600,
            range_km: 100. - 20. * i as f64 * 0.0016,
            radial_velocity_km_s: -20.,
            snr_linear: 20.,
        });
        s.push(signature(
            if i % 2 == 0 { 0.1 } else { 0.5 },
            i as f64 * 1.2,
        ));
    }
    p.push(Point {
        sample_idx: 1_746_000_000_047_200,
        range_km: 99.056,
        radial_velocity_km_s: -20.,
        snr_linear: 500.,
    });
    s.push(signature(0.3, 1.));
    let tracks = split_phase(&p, &s, 0.65, 0.017, &Config::default()).unwrap();
    assert_eq!(tracks.len(), 2);
    for t in tracks {
        let ids: Vec<_> = t.members.iter().copied().filter(|&i| i < 60).collect();
        assert!(ids.iter().all(|&i| i % 2 == ids[0] % 2));
    }
}
