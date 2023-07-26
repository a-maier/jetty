use std::{fs::File, time::Duration};

use anyhow::Result;
use cpu_time::ProcessTime;
use jetty::{PseudoJet, anti_kt_f, cluster_if};

fn main() -> Result<()> {
    let input = File::open("data/momenta_showered.rmp.zst")?;
    let mut events = Vec::new();
    zstd::stream::copy_decode(input, &mut events).unwrap();
    let events: Vec<Vec<[f64; 4]>> = rmp_serde::from_slice(&events)?;

    let events: Vec<Vec<PseudoJet>> = events.into_iter().map(
        |ev| ev.into_iter().map(|p| p.into()).collect()
    ).collect();
    pub const NEVENTS: usize = 10000;
    assert_eq!(NEVENTS, events.len()); // helps with optimisations

    let algo = anti_kt_f(0.4);
    let start = ProcessTime::now();
    let njets: usize = events.into_iter()
        .map(|ev| cluster_if(ev, &algo, |j| j.pt2() > 100.).len())
        .sum();
    let cpu_time: Duration = start.elapsed();
    let avg_njets = njets as f64 / NEVENTS as f64;
    println!("Found {avg_njets:1} jets per event in {cpu_time:?}");

    Ok(())
}
