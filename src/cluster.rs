use crate::pseudojet::PseudoJet;
use crate::distance::Distance;

// Cluster `partons` into jets using the distance measure `d`
pub fn cluster<D: Distance>(mut partons: Vec<PseudoJet>, d: &D) -> Vec<PseudoJet> {
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

    // TODO: test keeping distances sorted

    let mut jets = Vec::with_capacity(partons.len());
    while let Some(&(_, i, j)) = dists.iter().min() {
        if i == j {
            // closest is beam axis, this is already a jet
            dists.retain(|(_, ii, jj)| *ii != i && *jj != i);
            let jet = partons.swap_remove(i);
            for (_dist, ii, jj) in &mut dists {
                if *ii == partons.len() {
                    *ii = i;
                }
                if *jj == partons.len() {
                    *jj = i;
                }
            }
            jets.push(jet)
        } else {
            // combine into new pseudojet
            debug_assert!(j > i);
            dists.retain(|(_, ii, jj)| *ii != i && *jj != i && *ii != j && *jj != j);
            let p1 = partons.swap_remove(i);
            let p2 = partons.swap_remove(j);
            for (_dist, ii, jj) in &mut dists {
                if *ii == partons.len() + 1 {
                    *ii = i;
                } else if *ii == partons.len() {
                    *ii = j;
                }
                if *jj == partons.len() + 1 {
                    *jj = i;
                } else if *jj == partons.len() {
                    *jj = j;
                }
            }
            let pj = p1 + p2;
            let j = partons.len();
            for (i, pi) in partons.iter().copied().enumerate() {
                dists.push((d.distance(&pi, &pj), i, j));
            }
            dists.push((d.beam_distance(&pj), j, j));
            partons.push(pj);
        }
    }
    jets
}
