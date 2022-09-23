use crate::distance::Distance;
use crate::pseudojet::PseudoJet;

use std::cmp::Ord;

use log::{debug, trace};
use noisy_float::prelude::*;

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
    partons: Vec<PseudoJet>,
    d: &D,
    mut accept: F,
) -> Vec<PseudoJet>
where
    D: Distance,
    F: FnMut(PseudoJet) -> bool,
{
    debug!("clustering partons: {:#?}", partons);
    let clustering = ClusterHistory::new(partons, d);

    clustering
        .filter_map(|s| match s {
            ClusterStep::Jet(jet) if accept(jet) => Some(jet),
            _ => None,
        })
        .collect()
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct ClusterHistory<D> {
    pseudojets: Vec<PseudoJet>,
    distance: D,
    distances: Vec<(N64, usize, usize)>,
}

impl<D: Distance> ClusterHistory<D> {
    pub fn new(partons: Vec<PseudoJet>, distance: D) -> Self {
        let distances = calc_distances(&partons, &distance);
        Self {
            pseudojets: partons,
            distance,
            distances,
        }
    }

    fn extract_as_jet(&mut self, i: usize) -> PseudoJet {
        self.distances.retain(|(_, ii, jj)| *ii != i && *jj != i);
        let jet = self.pseudojets.swap_remove(i);
        debug!("new jet: {:?}", jet);
        for (_dist, ii, jj) in &mut self.distances {
            if *ii == self.pseudojets.len() {
                *ii = i
            }
            if *jj == self.pseudojets.len() {
                *jj = i
            }
        }
        trace!("distances: {:#?}", self.distances);
        jet
    }

    fn combine(&mut self, i: usize, j: usize) -> [PseudoJet; 2] {
        let res = [self.pseudojets[i], self.pseudojets[j]];
        let (i, j) = minmax(i, j);
        debug_assert!(j > i);
        debug!("cluster pseudojets {} {}", i, j);
        self.distances.retain(|(_, ii, jj)| *ii != j && *jj != j);
        let p2 = self.pseudojets.swap_remove(j);
        for (_dist, ii, jj) in &mut self.distances {
            if *ii == self.pseudojets.len() {
                *ii = j
            }
            if *jj == self.pseudojets.len() {
                *jj = j
            }
        }
        self.pseudojets[i] += p2;
        // update distances
        let affected_dists = self
            .distances
            .iter_mut()
            .filter(|(_dist, ii, jj)| *ii == i || *jj == i);
        for (dist, ii, jj) in affected_dists {
            *dist = if ii != jj {
                self.distance
                    .distance(&self.pseudojets[*ii], &self.pseudojets[*jj])
            } else {
                self.distance.beam_distance(&self.pseudojets[i])
            };
        }
        trace!("distances: {:#?}", self.distances);
        res
    }
}

impl<D: Distance> Iterator for ClusterHistory<D> {
    type Item = ClusterStep;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(&(_, i, j)) = self.distances.iter().min() {
            if i == j {
                Some(self.extract_as_jet(i).into())
            } else {
                Some(self.combine(i, j).into())
            }
        } else {
            None
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum ClusterStep {
    Combine([PseudoJet; 2]),
    Jet(PseudoJet),
}

impl From<[PseudoJet; 2]> for ClusterStep {
    fn from(source: [PseudoJet; 2]) -> Self {
        Self::Combine(source)
    }
}

impl From<PseudoJet> for ClusterStep {
    fn from(jet: PseudoJet) -> Self {
        Self::Jet(jet)
    }
}

fn calc_distances<D: Distance>(
    pseudojets: &[PseudoJet],
    d: &D,
) -> Vec<(N64, usize, usize)> {
    let n = pseudojets.len();

    let mut dists = Vec::with_capacity((n * (n + 1)) / 2);
    for i in 0..n {
        for j in i + 1..n {
            dists.push((d.distance(&pseudojets[i], &pseudojets[j]), i, j));
        }
        dists.push((d.beam_distance(&pseudojets[i]), i, i))
    }
    trace!("distances: {:#?}", dists);
    dists
}
