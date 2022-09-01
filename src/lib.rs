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

pub mod pseudojet;
pub mod distance;
pub mod cluster;

pub use pseudojet::{PseudoJet, pseudojet, pseudojet_f};
pub use distance::{anti_kt, kt, gen_kt, cambridge_aachen};
pub use distance::{anti_kt_f, kt_f, gen_kt_f, cambridge_aachen_f};
pub use cluster::{cluster, cluster_if, ClusterHistory, ClusterStep};

#[cfg(test)]
mod tests {
    fn log_init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn tst_cluster() {
        log_init();

        use super::{anti_kt_f, cluster, pseudojet_f};
        let partons = vec![
            pseudojet_f(69.26839536506921, 2.523788521334493, 3.311656952663986, -69.14314149775808),
            pseudojet_f(2.292439531535948, 1.678009926288044, -0.01258571588949442, 1.561858922116857),
            pseudojet_f(3.289626669238658, 0.7015436852248482, -1.673170474700709, -2.744081019808675),
            pseudojet_f(211.8446186434034, -1.500613625861779, -1.145111454447386, -211.8362087242677),
            pseudojet_f(36.19723806801562, 1.307180544524667, 0.1299048201403605, -36.17339419399588),
            pseudojet_f(1321.648751509191, -0.9834513292440703, 2.690421552193925, -1321.645647222113),
            pseudojet_f(614.6463548568801, -0.9800167457297735, 1.332258289303498, -614.6441297156254),
            pseudojet_f(84.89137294427485, -0.6700938500184943, -0.3903008019701275, -84.88783092929236),
            pseudojet_f(323.1911878589112, 3.631300879389308, -0.9682072466926734, -323.1693366306664),
        ];

        let jets = cluster(partons, &anti_kt_f(0.4));
        assert_eq!(jets.len(), 7);

    }

    #[test]
    fn tst_cluster_none() {
        log_init();

        use super::{anti_kt_f, cluster, pseudojet_f};
        let partons = vec![
            pseudojet_f(6.918281417330659, 0.0921329982846809, -1.37262399277452, -6.78012040117026),
            pseudojet_f(14.08869844306916, -2.416716165822407, -3.598378480403583, 13.40531906018502),
            pseudojet_f(10.58285213260104, -0.1240036336102471, -0.7792325830056485, 10.55339655944626),
            pseudojet_f(1.512949203734659, 0.4907951919308299, -0.3455630061912586, 1.388784209807401),
        ];

        let jets = cluster(partons, &anti_kt_f(0.4));
        assert_eq!(jets.len(), 4);
    }

    #[test]
    fn tst_cluster_both() {
        log_init();

        use super::{anti_kt_f, cluster, cluster_if, pseudojet_f};
        let partons = vec![
            pseudojet_f(0.2626773221934335, -0.08809521946454194, -0.1141608706693822, -0.2195584284654444),
            pseudojet_f(2.21902459329915, -0.7529973704809976, -0.9658189214109036, -1.850475321845671)
        ];

        // get all jets
        let all_jets = cluster(partons.clone(), &anti_kt_f(0.4));
        assert_eq!(all_jets.len(), 1);

        // get all jets with at least 40 GeV
        let jets_40gev = cluster_if(partons.clone(), &anti_kt_f(0.4), |jet| jet.pt2() > 20.);
        assert_eq!(jets_40gev.len(), 0);
    }

    #[test]
    fn tst_cluster_3_to_2() {
        log_init();

        use super::{anti_kt_f, cluster, pseudojet_f};
        let partons = vec![
            pseudojet_f(48.32406329129799, -3.576937946768497, 0.1029621819467338, -48.1913893418257),
            pseudojet_f(90.45021831804598, -6.668149504968421, 2.750224246879194, -90.16215415879022),
            pseudojet_f(8.331785929751781, 1.154124013760492, 1.428371653463948, -8.126894176723383)
        ];

        let jets = cluster(partons, &anti_kt_f(0.4));
        assert_eq!(jets.len(), 2);
    }

    #[test]
    fn tst_cluster_8_to_7() {
        log_init();

        use super::{anti_kt_f, cluster, pseudojet_f};
        let partons = vec![
            pseudojet_f(55.0566721275858, 2.853555315376817, -0.2177619033434216, -54.98224211124915),
            pseudojet_f(6.22300156243039, 1.897786377727834, 3.878240628357652, -4.481451209047378),
            pseudojet_f(34.72698340289853, -1.789702417265187, -0.8624446512551973, 34.67011004808532),
            pseudojet_f(0.2482353087299138, 0.1606423482788706, 0.09154657065594995, 0.1656322125698495),
            pseudojet_f(4.49341017269029, 1.186315084968943, 2.123673434575869, -3.778015701726003),
            pseudojet_f(12.33218634443971, 0.8398185974883793, -0.8370639741138476, -12.27504984350526),
            pseudojet_f(11.04475618445089, -1.016498045091359, -0.9891379652313798, -10.95330895136394),
            pseudojet_f(3.689996497465931, 2.126537512149597, 0.3956733858199553, -2.989540923366796),
        ];

        let jets = cluster(partons, &anti_kt_f(0.4));
        assert_eq!(jets.len(), 7);
    }

}
