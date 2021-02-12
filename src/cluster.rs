use crate::pseudojet::PseudoJet;
use crate::distance::Distance;

use std::cmp::Ord;

use log::{debug, trace};

fn minmax<T: Ord>(i: T, j: T) -> (T, T) {
    if i > j {
        (j, i)
    } else {
        (i, j)
    }
}

/// Cluster `partons` into jets using the distance measure `d`
pub fn cluster<D: Distance>(partons: Vec<PseudoJet>, d: &D) -> Vec<PseudoJet> {
    cluster_if(partons, d, |_| true)
}

/// Cluster `partons` into jets using the distance measure `d`
/// Only jets for which `accept` is true are returned
pub fn cluster_if<D, F>(
    mut partons: Vec<PseudoJet>,
    d: &D,
    mut accept: F
) -> Vec<PseudoJet>
where
    D: Distance,
    F: FnMut(PseudoJet) -> bool
{
    let n = partons.len();

    let mut dists = Vec::with_capacity((n * (n + 1))/2);
    for i in 0..n {
        for j in i+1..n {
            dists.push((
                d.distance(&partons[i], &partons[j]),
                i, j
            ));
        }
        dists.push((
            d.beam_distance(&partons[i]),
            i, i
        ))
    }
    trace!("distances: {:#?}", dists);

    // TODO: test keeping distances sorted

    let mut jets = Vec::with_capacity(partons.len());
    while let Some(&(_, i, j)) = dists.iter().min() {
        if i == j {
            // closest is beam axis, this is already a jet
            dists.retain(|(_, ii, jj)| *ii != i && *jj != i);
            let jet = partons.swap_remove(i);
            debug!("new jet: {:?}", jet);
            for (_dist, ii, jj) in &mut dists {
                if *ii == partons.len() {
                    *ii = i
                }
                if *jj == partons.len() {
                    *jj = i;
                }
            }
            trace!("distances: {:#?}", dists);
            if accept(jet) { jets.push(jet) }
        } else {
            // combine into new pseudojet
            let (i, j) = minmax(i, j);
            debug_assert!(j > i);
            dists.retain(|(_, ii, jj)| *ii != i && *jj != i && *ii != j && *jj != j);
            let p2 = partons.swap_remove(j);
            let p1 = partons.swap_remove(i);
            for (_dist, ii, jj) in &mut dists {
                if *ii == partons.len() {
                    *ii = i;
                } else if *ii == partons.len() + 1 {
                    *ii = j;
                }
                if *jj == partons.len() {
                    *jj = i;
                } else if *jj == partons.len() + 1 {
                    *jj = j;
                }
            }
            let pj = p1 + p2;
            debug!("cluster pseudojets {} {}", i, j);
            let j = partons.len();
            for (i, pi) in partons.iter().copied().enumerate() {
                dists.push((d.distance(&pi, &pj), i, j));
            }
            dists.push((d.beam_distance(&pj), j, j));
            trace!("distances: {:#?}", dists);
            partons.push(pj);
        }
    }
    jets
}
