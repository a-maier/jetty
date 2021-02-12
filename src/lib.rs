pub mod pseudojet;
pub mod distance;
pub mod cluster;

pub use pseudojet::{PseudoJet, pseudojet, pseudojet_f};
pub use distance::{anti_kt, kt, gen_kt, cambridge_aachen};
pub use distance::{anti_kt_f, kt_f, gen_kt_f, cambridge_aachen_f};
pub use cluster::{cluster, cluster_if};

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
}
