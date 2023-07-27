use crate::{distance::Distance, cluster_naive::ClusterNaive};
use crate::pseudojet::PseudoJet;

use std::cmp::Ord;

use log::debug;

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
    let clustering = ClusterNaive::new(partons, d);

    clustering
        .filter_map(|s| match s {
            ClusterStep::Jet(jet) if accept(jet) => Some(jet),
            _ => None,
        })
        .collect()
}

/// Result of a clustering step
#[derive(Clone, Debug, Ord, PartialOrd, Hash)]
pub enum ClusterStep {
    /// Two pseudojets were combined into a new pseudojet
    Combine([PseudoJet; 2]),
    /// A jet was found
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

impl PartialEq for ClusterStep {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Combine(left), Self::Combine(right)) => left == right || (
                (&left[0] == &right[1]) && (&left[1] == &right[0])
            ),
            (Self::Jet(l0), Self::Jet(r0)) => l0 == r0,
            _ => false,
        }
    }
}

impl Eq for ClusterStep { }

pub trait ClusterHist: Iterator<Item = ClusterStep>{}

impl<T> ClusterHist for T where T: Iterator<Item = ClusterStep> { }
