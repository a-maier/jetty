//! Implementations of common inclusive jet clustering algorithms.
//!
//! In the current version, the following distance measures are implemented:
//!
//! - [anti-kt](https://arxiv.org/abs/0802.1189)
//! - [Cambridge](https://arxiv.org/abs/hep-ph/9707323)/[Aachen](https://arxiv.org/abs/hep-ph/9907280)
//! - [kt](https://arxiv.org/abs/hep-ph/9305266)
//! - Generalised kt.
//!
//! This crate uses a naive clustering implementation. For
//! state-of-the-art implementations of many more jet algorithms, have a
//! look at the excellent [fastjet](http://fastjet.fr/) library.
//!
//! # Examples
//!
//! Cluster a number of partons into jets using the anti-kt algorithm with radius 0.4:
//!
//! ```rust
//! use jetty::{anti_kt_f, cluster, cluster_if, pseudojet_f};
//!
//! let partons = vec![
//!     pseudojet_f(0.2626773221934335, -0.08809521946454194, -0.1141608706693822, -0.2195584284654444),
//!     pseudojet_f(2.21902459329915, -0.7529973704809976, -0.9658189214109036, -1.850475321845671)
//! ];
//!
//! // get all jets
//! let all_jets = cluster(partons.clone(), &anti_kt_f(0.4));
//! assert_eq!(all_jets.len(), 1);
//!
//! // get all jets with at least 40 GeV
//! let jets_40gev = cluster_if(
//!    partons.clone(),
//!    &anti_kt_f(0.4),
//!    |jet| jet.pt2() > 40. * 40.
//! );
//! assert_eq!(jets_40gev.len(), 0);
//! ```

pub mod cluster;
pub mod distance;
pub mod pseudojet;

#[cfg(test)]
mod test_data;

pub use cluster::{cluster, cluster_if, ClusterNaive, ClusterStep};
pub use distance::{anti_kt, cambridge_aachen, gen_kt, kt};
pub use distance::{anti_kt_f, cambridge_aachen_f, gen_kt_f, kt_f};
pub use pseudojet::{pseudojet, pseudojet_f, PseudoJet};

#[cfg(test)]
mod tests {
    use crate::test_data::*;

    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn tst_cluster() {
        log_init();

        use super::{anti_kt_f, cluster};
        let partons = partons_9_to_7();

        let jets = cluster(partons, &anti_kt_f(0.4));
        assert_eq!(jets.len(), 7);
    }

    #[test]
    fn tst_cluster_none() {
        log_init();

        use super::{anti_kt_f, cluster};
        let partons = partons_4_to_4();

        let jets = cluster(partons, &anti_kt_f(0.4));
        assert_eq!(jets.len(), 4);
    }

    #[test]
    fn tst_cluster_both() {
        log_init();

        use super::{anti_kt_f, cluster, cluster_if};
        let partons = partons_2_to_1();

        // get all jets
        let all_jets = cluster(partons.clone(), &anti_kt_f(0.4));
        assert_eq!(all_jets.len(), 1);

        // get all jets with at least 40 GeV
        let jets_40gev =
            cluster_if(partons.clone(), &anti_kt_f(0.4), |jet| jet.pt2() > 20.);
        assert_eq!(jets_40gev.len(), 0);
    }

    #[test]
    fn tst_cluster_3_to_2() {
        log_init();

        use super::{anti_kt_f, cluster};
        let partons = partons_3_to_2();

        let jets = cluster(partons, &anti_kt_f(0.4));
        assert_eq!(jets.len(), 2);
    }

    #[test]
    fn tst_cluster_8_to_7() {
        log_init();

        use super::{anti_kt_f, cluster};
        let partons = partons_8_to_7();

        let jets = cluster(partons, &anti_kt_f(0.4));
        assert_eq!(jets.len(), 7);
    }
}
