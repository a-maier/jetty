use std::cmp::min;

use log::{debug, trace};
use noisy_float::{types::N64, prelude::Float};

use crate::{PseudoJet, distance::Distance, ClusterStep};

/// Cluster history using the geometric O(N^2) approach of [arXiv:0512210](https://arxiv.org/abs/hep-ph/0512210)
#[derive(Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct ClusterGeom<D> {
    pseudojets: Vec<PseudoJetWithDist>,
    distance: D,
}

impl<D: Distance> ClusterGeom<D> {
    /// Initialise clustering for the given `partons` and `distance`
    pub fn new(partons: Vec<PseudoJet>, distance: D) -> Self {
        let mut pseudojets = Vec::from_iter(
            partons.into_iter().map(
                |pseudojet| PseudoJetWithDist { pseudojet, ..Default::default()}
            )
        );
        for i in 0..pseudojets.len() {
            pseudojets[i].beam_dist = distance.beam_distance(&pseudojets[i].pseudojet);
            let mut nearest_gdist = N64::max_value();
            let mut nearest_idx = usize::MAX;
            for j in (0..i).chain((i + 1)..pseudojets.len()) {
                let gdist = pseudojets[i].delta_r2(&pseudojets[j]);
                if gdist < nearest_gdist {
                    nearest_gdist = gdist;
                    nearest_idx = j;
                }
            }
            pseudojets[i].nearest_neighbour_idx = nearest_idx;
            if nearest_idx < usize::MAX {
                assert!(nearest_idx < pseudojets.len());
                pseudojets[i].nearest_dist = distance.distance(
                    &pseudojets[i].pseudojet,
                    &pseudojets[nearest_idx].pseudojet
                );
                pseudojets[nearest_idx].nearest_neighbour_for.push(i);
            } else {
                pseudojets[i].nearest_dist = N64::max_value();
            }
        }
        Self {
            pseudojets,
            distance,
        }
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
            let j_is_nearest_for = self.pseudojets[j].nearest_neighbour_for.clone();
            let nearest_j = self.pseudojets[j].nearest_neighbour_idx;

            for idx in i_is_nearest_for {
                debug_assert_eq!(self.pseudojets[idx].nearest_neighbour_idx, i);
                self.pseudojets[idx].nearest_neighbour_idx = j;
            }
            for idx in j_is_nearest_for {
                debug_assert_eq!(self.pseudojets[idx].nearest_neighbour_idx, j);
                self.pseudojets[idx].nearest_neighbour_idx = i;
            }

            // for particles very close to the beam axis
            // the distance to all others can be infinity
            // within floating point precision
            // in that case `nearest_i` will be `usize::MAX`
            // (see `new`)
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

        let others = (0..pos).chain((pos + 1)..self.pseudojets.len());
        let nearest_idx = others.map(|idx| {
            let gdist = self.pseudojets[pos].delta_r2(&self.pseudojets[idx]);
            (gdist, idx)
        }).min_by_key(|(d, _)| *d)
            .map(|(_d, idx)| idx)
            .unwrap_or(usize::MAX);
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
        let beam_dist = self.distance.beam_distance(&pseudojet);
        let mut pseudojet = PseudoJetWithDist {
            pseudojet,
            beam_dist,
            nearest_dist: N64::max_value(),
            ..Default::default()
        };
        let mut nearest_dist = N64::max_value();
        let mut nearest_idx = usize::MAX;
        for n in 0..self.pseudojets.len() {
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
}

impl<D: Distance> Iterator for ClusterGeom<D> {
    type Item = ClusterStep;

    /// Perform the next clustering step
    fn next(&mut self) -> Option<Self::Item> {
        trace!("pseudojets: {:#?}", self.pseudojets);
        let Some(i) = self.min_idx() else {
            return None
        };
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

#[derive(Clone, Debug, Default, Eq, PartialEq, Hash)]
struct PseudoJetWithDist {
    pseudojet: PseudoJet,
    beam_dist: N64,
    nearest_dist: N64,
    nearest_neighbour_idx: usize,
    nearest_neighbour_for: Vec<usize>, // TODO: An IndexSet might be better
}
impl PseudoJetWithDist {
    fn min_dist(&self) -> N64 {
        min(self.nearest_dist, self.beam_dist)
    }

    fn delta_r2(&self, p: &PseudoJetWithDist) -> N64 {
        self.pseudojet.delta_r2(&p.pseudojet)
    }
}

impl PartialOrd for PseudoJetWithDist {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.min_dist().partial_cmp(&other.min_dist())
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
        let tree = ClusterGeom::new(partons, anti_kt_f(0.4));

        for (naive, tree) in naive.zip(tree) {
            assert_eq!(naive, tree)
        }
    }

    #[test]
    fn cmp_3_to_2() {
        log_init();

        let partons = partons_3_to_2();
        let naive = ClusterNaive::new(partons.clone(), anti_kt_f(0.4));
        let tree = ClusterGeom::new(partons, anti_kt_f(0.4));

        for (naive, tree) in naive.zip(tree) {
            assert_eq!(naive, tree)
        }
    }

    #[test]
    fn cmp_4_to_4() {
        log_init();

        let partons = partons_4_to_4();
        let naive = ClusterNaive::new(partons.clone(), anti_kt_f(0.4));
        let tree = ClusterGeom::new(partons, anti_kt_f(0.4));

        for (naive, tree) in naive.zip(tree) {
            assert_eq!(naive, tree)
        }
    }

    #[test]
    fn cmp_8_to_7() {
        log_init();

        let partons = partons_8_to_7();
        let naive = ClusterNaive::new(partons.clone(), anti_kt_f(0.4));
        let tree = ClusterGeom::new(partons, anti_kt_f(0.4));

        for (naive, tree) in naive.zip(tree) {
            assert_eq!(naive, tree)
        }
    }

    #[test]
    fn cmp_9_to_7() {
        log_init();

        let partons = partons_9_to_7();
        let naive = ClusterNaive::new(partons.clone(), anti_kt_f(0.4));
        let tree = ClusterGeom::new(partons, anti_kt_f(0.4));

        for (naive, tree) in naive.zip(tree) {
            assert_eq!(naive, tree)
        }
    }
}
