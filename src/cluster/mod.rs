//!
//! In general, it is recommended to use the general-purpose
//! [ClusterHistory]. It chooses dynamically between the following
//! specialised algorithms:
//!
//! * [ClusterNaive](crate::cluster::naive::ClusterNaive): choose this
//!   algorithm if you know that the number of partons is not too big,
//!   at most about 25. Also use this algorithm for custom distances
//!   where the actual nearest neighbours are not always the nearest
//!   neighbours in ΔR.
//!
//! * [ClusterGeom](crate::cluster::geom::ClusterGeom): the fastest
//!   implemented algorithm for a number of partons roughly between 25
//!   and 50.
//!
//! * [ClusterGeomTile](crate::cluster::geom_tile::ClusterGeomTile):
//!   the fastest implemented algorithm for a large number of partons
//!   starting at about 50.
//!
/// Clustering using the geometric O(N^2) approach of [arXiv:0512210](https://arxiv.org/abs/hep-ph/0512210)
pub mod geom;
/// Clustering using the geometric O(N^2) approach of [arXiv:0512210](https://arxiv.org/abs/hep-ph/0512210) with tiling
pub mod geom_tile;
/// Naive clustering
pub mod naive;

use crate::distance::Distance;
use crate::pseudojet::PseudoJet;

use std::cmp::Ord;
use std::hash::Hash;

use log::debug;

use self::{
    geom::ClusterGeom, geom_tile::ClusterGeomTile, naive::ClusterNaive,
};

/// Cluster `partons` into jets using the distance measure `d`
#[deprecated = "Use `Cluster::cluster` instead"]
pub fn cluster<D: Distance>(partons: Vec<PseudoJet>, d: &D) -> Vec<PseudoJet> {
    partons.cluster_if(d, |_| true)
}

/// Cluster `partons` into jets using the distance measure `d`
/// Only jets for which `accept` is true are returned
#[deprecated = "Use `Cluster::cluster_if` instead"]
pub fn cluster_if<D, F>(
    partons: Vec<PseudoJet>,
    d: &D,
    accept: F,
) -> Vec<PseudoJet>
where
    D: Distance,
    F: FnMut(PseudoJet) -> bool,
{
    partons.cluster_if(d, accept)
}

/// Objects that can be clustered into jets
pub trait Cluster {
    /// Cluster into jets using the distance measure `d`
    fn cluster<D: Distance>(self, d: D) -> Vec<PseudoJet>;

    /// Cluster into jets using the distance measure `d`
    /// Only jets for which `accept` is true are returned
    fn cluster_if<D, F>(self, d: D, accept: F) -> Vec<PseudoJet>
    where
        D: Distance,
        F: FnMut(PseudoJet) -> bool;
}

impl Cluster for Vec<PseudoJet> {
    fn cluster_if<D, F>(self, d: D, mut accept: F) -> Vec<PseudoJet>
    where
        D: Distance,
        F: FnMut(PseudoJet) -> bool,
    {
        debug!("clustering partons: {self:#?}");
        let clustering = ClusterHistory::new(self, d);

        clustering
            .filter_map(|s| match s {
                ClusterStep::Jet(jet) if accept(jet) => Some(jet),
                _ => None,
            })
            .collect()
    }

    fn cluster<D: Distance>(self, d: D) -> Vec<PseudoJet> {
        self.cluster_if(d, |_| true)
    }
}

impl<'a, T> Cluster for &'a [T]
where
    &'a T: Into<PseudoJet>,
{
    fn cluster_if<D, F>(self, d: D, accept: F) -> Vec<PseudoJet>
    where
        D: Distance,
        F: FnMut(PseudoJet) -> bool,
    {
        let partons = Vec::from_iter(self.iter().map(|p| p.into()));
        partons.cluster_if(d, accept)
    }

    fn cluster<D: Distance>(self, d: D) -> Vec<PseudoJet> {
        self.cluster_if(d, |_| true)
    }
}

/// Result of a clustering step
#[derive(Clone, Debug, Ord, PartialOrd)]
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
            (Self::Combine(left), Self::Combine(right)) => {
                left == right
                    || ((left[0] == right[1]) && (left[1] == right[0]))
            }
            (Self::Jet(l0), Self::Jet(r0)) => l0 == r0,
            _ => false,
        }
    }
}

impl Eq for ClusterStep {}

impl Hash for ClusterStep {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            ClusterStep::Combine([p1, p2]) => {
                if p1 < p2 {
                    [p1, p2].hash(state)
                } else {
                    [p2, p1].hash(state)
                }
            }
            ClusterStep::Jet(p1) => p1.hash(state),
        }
    }
}

/// Trait marking a clustering algorithm
pub trait ClusterHist: Iterator<Item = ClusterStep> {}

impl<T> ClusterHist for T where T: Iterator<Item = ClusterStep> {}

/// General-purpose cluster history
pub struct ClusterHistory<'a>(Box<dyn ClusterHist + 'a>);

impl<'a> ClusterHistory<'a> {
    const START_GEOM_THRESHOLD: usize = 25;
    const END_GEOM_THRESHOLD: usize = 49;
    const START_TILE_THRESHOLD: usize = Self::END_GEOM_THRESHOLD + 1;

    /// Initialise clustering for the given `partons` and `distance`
    pub fn new<D: Distance + 'a>(partons: Vec<PseudoJet>, distance: D) -> Self {
        let hist: Box<dyn ClusterHist> = match partons.len() {
            Self::START_TILE_THRESHOLD.. => {
                Box::new(ClusterGeomTile::new(partons, distance))
            }
            Self::START_GEOM_THRESHOLD..=Self::END_GEOM_THRESHOLD => {
                Box::new(ClusterGeom::new(partons, distance))
            }
            _ => Box::new(ClusterNaive::new(partons, distance)),
        };
        Self(hist)
    }
}

impl<'a> Iterator for ClusterHistory<'a> {
    type Item = ClusterStep;

    /// Perform the next clustering step
    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }
}
