// TODO: huge amount of code duplication with `geom`
use std::{cmp::min, f64::consts::PI};

use indexmap::IndexSet;
use itertools::Itertools;
use log::{debug, trace};
use noisy_float::prelude::*;
use num_traits::cast::ToPrimitive;

use crate::{PseudoJet, distance::Distance, ClusterStep};

const MAX_RAP: f64 = 5.;
const N_RAP_BINS: usize = 10;
const N_PHI_BINS: usize = 6;

/// Cluster history using the geometric O(N^2) approach of [arXiv:0512210](https://arxiv.org/abs/hep-ph/0512210) with tiling
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ClusterGeomTile<D> {
    pseudojets: Vec<PseudoJetWithDist>,
    distance: D,
    tiles: [[IndexSet<usize>; N_PHI_BINS]; N_RAP_BINS],
}

impl<D: Distance> ClusterGeomTile<D> {
    /// Initialise clustering for the given `partons` and `distance`
    pub fn new(partons: Vec<PseudoJet>, distance: D) -> Self {
        let pseudojets = partons.into_iter().map(
            |pseudojet| PseudoJetWithDist::new(pseudojet, &distance)
        ).collect();
        let mut res = Self {
            pseudojets,
            distance,
            tiles: Default::default(),
        };
        res.init_tiles();
        res.init_nearest();
        res
    }

    fn min_idx(&self) -> Option<usize> {
        self.pseudojets.iter().enumerate()
            .min_by(|(_, a), (_, b)| a.cmp(b))
            .map(|(n, _)| n)
    }

    // Exchange two pseudojets
    fn swap(&mut self, i: usize, j: usize) {
        assert!(i < self.pseudojets.len());
        assert!(j < self.pseudojets.len());
        if i != j {
            let i_is_nearest_for = self.pseudojets[i].nearest_neighbour_for.clone();
            let nearest_i = self.pseudojets[i].nearest_neighbour_idx;
            let [rap_idx_i, phi_idx_i] = self.tile_coord(&self.pseudojets[i].pseudojet);
            let j_is_nearest_for = self.pseudojets[j].nearest_neighbour_for.clone();
            let nearest_j = self.pseudojets[j].nearest_neighbour_idx;
            let [rap_idx_j, phi_idx_j] = self.tile_coord(&self.pseudojets[j].pseudojet);

            for idx in i_is_nearest_for {
                debug_assert_eq!(self.pseudojets[idx].nearest_neighbour_idx, i);
                self.pseudojets[idx].nearest_neighbour_idx = j;
            }
            for idx in j_is_nearest_for {
                debug_assert_eq!(self.pseudojets[idx].nearest_neighbour_idx, j);
                self.pseudojets[idx].nearest_neighbour_idx = i;
            }

            if nearest_i < self.pseudojets.len() {
                let to_update_idx  = self.pseudojets[nearest_i]
                    .nearest_neighbour_for
                    .iter()
                    .position(|&k| k == i)
                    .unwrap();
                self.pseudojets[nearest_i]
                    .nearest_neighbour_for[to_update_idx] = j;
            }

            if nearest_j < self.pseudojets.len() {
                let to_update_idx  = self.pseudojets[nearest_j]
                    .nearest_neighbour_for
                    .iter()
                    .position(|&k| k == j)
                    .unwrap();
                self.pseudojets[nearest_j]
                    .nearest_neighbour_for[to_update_idx] = i;
            }

            self.tiles[rap_idx_i][phi_idx_i].swap_remove(&i);
            self.tiles[rap_idx_i][phi_idx_i].insert(j);
            self.tiles[rap_idx_j][phi_idx_j].swap_remove(&j);
            self.tiles[rap_idx_j][phi_idx_j].insert(i);

            self.pseudojets.swap(i, j);
        }
    }

    // Remove pseudojet at `idx`, updating the nearest-neighbour indices
    fn remove(&mut self, idx: usize) -> PseudoJetWithDist {
        assert!(idx < self.pseudojets.len());
        trace!("Before removing {idx}: {:#?}", self.pseudojets);
        self.swap(idx, self.pseudojets.len() - 1);
        trace!("After swap: {:#?}", self.pseudojets);

        self.remove_nearest_link(self.pseudojets.len() - 1);
        let [rap_idx, phi_idx] = self.tile_coord(
            &self.pseudojets[self.pseudojets.len() - 1].pseudojet
        );
        self.tiles[rap_idx][phi_idx].swap_remove(&(self.pseudojets.len() - 1));
        let pseudojet = self.pseudojets.pop().unwrap();
        // TODO: maybe don't recalculate nearest neighbours yet
        self.update_nearest(&pseudojet.nearest_neighbour_for);
        trace!("After removal: {:#?}", self.pseudojets);
        pseudojet
    }

    fn update_nearest(&mut self, pos: &[usize]) {
        for idx in pos {
            self.update_nearest_at_idx(*idx);
        }
    }

    fn update_nearest_at_idx(&mut self, pos: usize) {
        assert!(pos < self.pseudojets.len());
        self.remove_nearest_link(pos);

        let [rap_idx, phi_idx] = self.tile_coord(
            &self.pseudojets[pos].pseudojet
        );
        let nearest_idx = self.tile_neighbours(rap_idx, phi_idx)
            .filter(|&j| j != pos)
            .map(|idx| {
                let gdist = self.pseudojets[pos].delta_r2(&self.pseudojets[idx]);
                (gdist, idx)
            }).min_by_key(|(d, _)| *d)
            .map(|(_d, idx)| idx)
            .unwrap_or(usize::MAX);
        debug_assert_ne!(nearest_idx, pos);
        self.pseudojets[pos].nearest_neighbour_idx = nearest_idx;
        if nearest_idx < usize::MAX {
            assert!(nearest_idx < self.pseudojets.len());
            self.pseudojets[nearest_idx].nearest_neighbour_for.push(pos);
            self.pseudojets[pos].nearest_dist = self.distance(
                &self.pseudojets[pos],
                &self.pseudojets[nearest_idx]
            );
        } else {
            self.pseudojets[pos].nearest_dist = N64::max_value()
        }
    }

    fn push(&mut self, pseudojet: PseudoJet) {
        trace!("before push: {:#?}", self.pseudojets);
        let [rap_idx, phi_idx] = self.tile_coord(&pseudojet);
        let mut pseudojet = PseudoJetWithDist::new(pseudojet, &self.distance);
        let mut nearest_dist = N64::max_value();
        let mut nearest_idx = usize::MAX;
        let neighbours = Vec::from_iter(self.tile_neighbours(rap_idx, phi_idx));
        for n in neighbours {
            let d = self.distance(&pseudojet, &self.pseudojets[n]);
            if d < nearest_dist {
                nearest_dist = d;
                nearest_idx = n;
            }
            if d < self.pseudojets[n].nearest_dist {
                self.remove_nearest_link(n);
                self.pseudojets[n].nearest_neighbour_idx = self.pseudojets.len();
                pseudojet.nearest_neighbour_for.push(n);
            }
        }
        pseudojet.nearest_neighbour_idx = nearest_idx;
        if nearest_idx < usize::MAX {
            let len = self.pseudojets.len();
            assert!(nearest_idx < len);
            self.pseudojets[nearest_idx].nearest_neighbour_for.push(len);
            pseudojet.nearest_dist = self.distance(
                &pseudojet,
                &self.pseudojets[nearest_idx]
            )
        }
        self.tiles[rap_idx][phi_idx].insert(self.pseudojets.len());
        self.pseudojets.push(pseudojet);
        trace!("after push: {:#?}", self.pseudojets);
    }

    // update such that no other pseudojet considers itself the
    // nearest neighbour for the one at `pos`
    fn remove_nearest_link(&mut self, pos: usize) {
        assert!(pos < self.pseudojets.len());
        let nearest_idx = self.pseudojets[pos].nearest_neighbour_idx;
        if nearest_idx < self.pseudojets.len() {
            let to_remove_idx  = self.pseudojets[nearest_idx]
                .nearest_neighbour_for
                .iter()
                .position(|&j| j == pos)
                .unwrap();
            self.pseudojets[nearest_idx]
                .nearest_neighbour_for
                .swap_remove(to_remove_idx);
        }
    }

    fn distance(&self, p1: &PseudoJetWithDist, p2: &PseudoJetWithDist) -> N64 {
        self.distance.distance(&p1.pseudojet, &p2.pseudojet)
    }

    fn init_tiles(&mut self) {
        for (n, p) in self.pseudojets.iter().enumerate() {
            let [rap, phi] = self.tile_coord(&p.pseudojet);
            self.tiles[rap][phi].insert(n);
        }
    }

    fn tile_coord(&self, pseudojet: &PseudoJet) -> [usize; 2] {
        let rap_coord = (f64::from(pseudojet.rap()) + MAX_RAP).floor() as i32;
        let rap_coord = rap_coord.clamp(0, N_RAP_BINS as i32 - 1) as usize;
        let phi_coord = pseudojet.phi() * (N_PHI_BINS as f64 / (2. * PI));
        assert!(phi_coord >= 0.);
        assert!(phi_coord < N_PHI_BINS as f64);
        let phi_coord = phi_coord.to_usize().unwrap();
        [rap_coord, phi_coord]
    }

    fn tile_neighbours(
        &self,
        rap_idx: usize,
        phi_idx: usize
    ) -> impl Iterator<Item = usize> + '_ {
        let rap_idx_range = match rap_idx + 1 {
            1 => 0..2,
            N_RAP_BINS => (N_RAP_BINS - 2)..N_RAP_BINS,
            _ => (rap_idx - 1)..(rap_idx + 2),
        };
        let phi_idx_range = match phi_idx + 1 {
            1 => [N_PHI_BINS - 1, 0 , 1],
            N_PHI_BINS => [N_PHI_BINS - 2, N_PHI_BINS - 1 , 0],
            _ => [phi_idx - 1, phi_idx, phi_idx + 1]
        };
        rap_idx_range.cartesian_product(phi_idx_range)
            .flat_map(|(rap, phi)| self.tiles[rap][phi].iter().copied())
    }

    fn init_nearest(&mut self)  {
        for i in 0..self.pseudojets.len() {
            let [rap_idx, phi_idx] = self.tile_coord(&self.pseudojets[i].pseudojet);
            let nearest = self.tile_neighbours(rap_idx, phi_idx)
                .filter(|&j| i != j)
                .map(|j| {
                    let gdist = self.pseudojets[i].delta_r2(&self.pseudojets[j]);
                    (gdist, j)
                }).min();
            if let Some((_, nearest_idx)) = nearest {
                assert!(nearest_idx < self.pseudojets.len());
                self.pseudojets[i].nearest_neighbour_idx = nearest_idx;
                self.pseudojets[i].nearest_dist = self.distance(
                    &self.pseudojets[i],
                    &self.pseudojets[nearest_idx]
                );
                self.pseudojets[nearest_idx].nearest_neighbour_for.push(i);
            }
        }
    }
}

impl<D: Distance> Iterator for ClusterGeomTile<D> {
    type Item = ClusterStep;

    /// Perform the next clustering step
    fn next(&mut self) -> Option<Self::Item> {
        trace!("pseudojets: {:#?}", self.pseudojets);
        let i = self.min_idx()?;
        let pi = self.remove(i);
        if pi.beam_dist < pi.nearest_dist {
            let pi = pi.pseudojet;
            debug!("new jet: {pi:?}");
            Some(pi.into())
        } else {
            let j = pi.nearest_neighbour_idx;
            debug!("cluster pseudojets {i} {j}");
            let pj = self.remove(j);
            let pi = pi.pseudojet;
            let pj = pj.pseudojet;
            self.push(pi + pj);
            Some([pi, pj].into())
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
struct PseudoJetWithDist {
    pseudojet: PseudoJet,
    beam_dist: N64,
    nearest_dist: N64,
    nearest_neighbour_idx: usize,
    nearest_neighbour_for: Vec<usize>, // TODO: An IndexSet might be better
}

impl Default for PseudoJetWithDist {
    fn default() -> Self {
        Self {
            pseudojet: Default::default(),
            beam_dist: N64::max_value(),
            nearest_dist: N64::max_value(),
            nearest_neighbour_idx: usize::MAX,
            nearest_neighbour_for: Default::default(),
        }
    }
}

impl PseudoJetWithDist {
    fn new<D: Distance>(pseudojet: PseudoJet, distance: D) -> Self {
        Self {
            beam_dist: distance.beam_distance(&pseudojet),
            pseudojet,
            ..Default::default()
        }
    }

    fn min_dist(&self) -> N64 {
        min(self.nearest_dist, self.beam_dist)
    }

    fn delta_r2(&self, p: &PseudoJetWithDist) -> N64 {
        self.pseudojet.delta_r2(&p.pseudojet)
    }

}

impl PartialOrd for PseudoJetWithDist {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PseudoJetWithDist {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.min_dist().cmp(&other.min_dist())
    }
}

#[cfg(test)]
mod tests {
    use crate::{test_data::*, cluster::naive::ClusterNaive, anti_kt_f};

    use super::*;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn cmp_2_to_1() {
        log_init();

        let partons = partons_2_to_1();
        let naive = ClusterNaive::new(partons.clone(), anti_kt_f(0.4));
        let tree = ClusterGeomTile::new(partons, anti_kt_f(0.4));

        for (naive, tree) in naive.zip(tree) {
            assert_eq!(naive, tree)
        }
    }

    #[test]
    fn cmp_3_to_2() {
        log_init();

        let partons = partons_3_to_2();
        let naive = ClusterNaive::new(partons.clone(), anti_kt_f(0.4));
        let tree = ClusterGeomTile::new(partons, anti_kt_f(0.4));

        for (naive, tree) in naive.zip(tree) {
            assert_eq!(naive, tree)
        }
    }

    #[test]
    fn cmp_4_to_4() {
        log_init();

        let partons = partons_4_to_4();
        let naive = ClusterNaive::new(partons.clone(), anti_kt_f(0.4));
        let tree = ClusterGeomTile::new(partons, anti_kt_f(0.4));

        for (naive, tree) in naive.zip(tree) {
            assert_eq!(naive, tree)
        }
    }

    #[test]
    fn cmp_8_to_7() {
        log_init();

        let partons = partons_8_to_7();
        let naive = ClusterNaive::new(partons.clone(), anti_kt_f(0.4));
        let tree = ClusterGeomTile::new(partons, anti_kt_f(0.4));

        for (naive, tree) in naive.zip(tree) {
            assert_eq!(naive, tree)
        }
    }

    #[test]
    fn cmp_9_to_7() {
        log_init();

        let partons = partons_9_to_7();
        let naive = ClusterNaive::new(partons.clone(), anti_kt_f(0.4));
        let tree = ClusterGeomTile::new(partons, anti_kt_f(0.4));

        for (naive, tree) in naive.zip(tree) {
            assert_eq!(naive, tree)
        }
    }
}
