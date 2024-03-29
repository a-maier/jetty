//! Implementations of common inclusive jet clustering algorithms.
//!
//! In the current version, the following distance measures are implemented:
//!
//! - [anti-kt](https://arxiv.org/abs/0802.1189)
//! - [Cambridge](https://arxiv.org/abs/hep-ph/9707323)/[Aachen](https://arxiv.org/abs/hep-ph/9907280)
//! - [kt](https://arxiv.org/abs/hep-ph/9305266)
//! - Generalised kt.
//!
//! For state-of-the-art implementations of many more jet algorithms,
//! have a look at the excellent [fastjet](http://fastjet.fr/)
//! library.
//!
//! # Examples
//!
//! Cluster a number of partons into jets using the anti-kt algorithm with radius 0.4:
//!
//! ```rust
//! use jetty::{anti_kt_f, pseudojet_f, Cluster, ClusterHistory, ClusterStep};
//!
//! let partons = vec![
//!     pseudojet_f(0.2626773221934335, -0.08809521946454194, -0.1141608706693822, -0.2195584284654444),
//!     pseudojet_f(2.21902459329915, -0.7529973704809976, -0.9658189214109036, -1.850475321845671)
//! ];
//!
//! // get all jets
//! let all_jets = partons.clone().cluster(anti_kt_f(0.4));
//! assert_eq!(all_jets.len(), 1);
//!
//! // get all jets with at least 40 GeV
//! let jets_40gev = partons.clone().cluster_if(
//!    anti_kt_f(0.4),
//!    |jet| jet.pt2() > 40. * 40.
//! );
//! assert_eq!(jets_40gev.len(), 0);
//!
//! // go through the cluster history step-by-step
//! let history = ClusterHistory::new(partons, anti_kt_f(0.4));
//! for step in history {
//!    match step {
//!       ClusterStep::Jet(j) => println!("Found a jet: {j:?}"),
//!       ClusterStep::Combine([_j1, _j2]) => println!("Combined two pseudojets"),
//!    }
//! }
//! ```
/// Jet clustering algorithms
pub mod cluster;
/// Distances and jet definitions
pub mod distance;
/// Pseudojets
pub mod pseudojet;

#[cfg(test)]
mod test_data;

#[allow(deprecated)]
pub use cluster::{cluster, cluster_if};

pub use cluster::{Cluster, ClusterHistory, ClusterStep};
pub use distance::{anti_kt, cambridge_aachen, gen_kt, kt};
pub use distance::{anti_kt_f, cambridge_aachen_f, gen_kt_f, kt_f};
pub use pseudojet::{pseudojet, pseudojet_f, PseudoJet};

#[cfg(test)]
mod tests {
    use super::{anti_kt_f, Cluster};
    use crate::test_data::*;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn tst_cluster() {
        log_init();

        let partons = partons_9_to_7();

        let jets = partons.cluster(anti_kt_f(0.4));
        assert_eq!(jets.len(), 7);
    }

    #[test]
    fn tst_cluster_none() {
        log_init();

        let partons = partons_4_to_4();

        let jets = partons.cluster(anti_kt_f(0.4));
        assert_eq!(jets.len(), 4);
    }

    #[test]
    fn tst_cluster_both() {
        log_init();

        let partons = partons_2_to_1();

        // get all jets
        let all_jets = partons.clone().cluster(anti_kt_f(0.4));
        assert_eq!(all_jets.len(), 1);

        // get all jets with at least 40 GeV
        let jets_40gev = partons
            .clone()
            .cluster_if(anti_kt_f(0.4), |jet| jet.pt2() > 20.);
        assert_eq!(jets_40gev.len(), 0);
    }

    #[test]
    fn tst_cluster_3_to_2() {
        log_init();

        let partons = partons_3_to_2();

        let jets = partons.cluster(anti_kt_f(0.4));
        assert_eq!(jets.len(), 2);
    }

    #[test]
    fn tst_cluster_8_to_7() {
        log_init();

        let partons = partons_8_to_7();

        let jets = partons.cluster(anti_kt_f(0.4));
        assert_eq!(jets.len(), 7);
    }
}
