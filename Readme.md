# jetty

Implementations of common inclusive jet clustering algorithms.

In the current version, the following distance measures are implemented:

- [anti-kt](https://arxiv.org/abs/0802.1189)
- [Cambridge](https://arxiv.org/abs/hep-ph/9707323)/[Aachen](https://arxiv.org/abs/hep-ph/9907280)
- [kt](https://arxiv.org/abs/hep-ph/9305266)
- Generalised kt.

For state-of-the-art implementations of many more jet algorithms,
have a look at the excellent [fastjet](http://fastjet.fr/)
library.

## Examples

Cluster a number of partons into jets using the anti-kt algorithm with radius 0.4:

```rust
use jetty::{anti_kt_f, pseudojet_f, Cluster, ClusterHistory, ClusterStep};

let partons = vec![
    pseudojet_f(0.2626773221934335, -0.08809521946454194, -0.1141608706693822, -0.2195584284654444),
    pseudojet_f(2.21902459329915, -0.7529973704809976, -0.9658189214109036, -1.850475321845671)
];

// get all jets
let all_jets = partons.clone().cluster(anti_kt_f(0.4));
assert_eq!(all_jets.len(), 1);

// get all jets with at least 40 GeV
let jets_40gev = partons.clone().cluster_if(
   anti_kt_f(0.4),
   |jet| jet.pt2() > 40. * 40.
);
assert_eq!(jets_40gev.len(), 0);

// go through the cluster history step-by-step
let history = ClusterHistory::new(partons, anti_kt_f(0.4));
for step in history {
   match step {
      ClusterStep::Jet(j) => println!("Found a jet: {j:?}"),
      ClusterStep::Combine([_j1, _j2]) => println!("Combined two pseudojets"),
   }
}
```

License: GPL-3.0-or-later
