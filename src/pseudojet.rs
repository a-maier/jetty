use std::convert::From;
use std::default::Default;
use std::f64::consts::PI;
use std::ops::{Add, AddAssign, Index, Sub, SubAssign};

use noisy_float::prelude::*;

pub const D: usize = 4;

/// A pseudojet is a particle momentum or a sum of momenta of clustered particles
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct PseudoJet {
    comp: [N64; D],
    inv_pt2: N64,
    phi: N64,
    rap: N64,
}

impl PseudoJet {
    /// Create pseudojet with vanishing four-momentum
    pub fn new() -> Self {
        Self::default()
    }

    /// Energy
    pub fn e(&self) -> N64 {
        self[0]
    }

    /// Momentum in x direction
    pub fn px(&self) -> N64 {
        self[1]
    }

    /// Momentum in y direction
    pub fn py(&self) -> N64 {
        self[2]
    }

    /// Momentum in z direction, i.e. along the beam axis
    pub fn pz(&self) -> N64 {
        self[3]
    }

    /// Azimuthal angle φ
    pub fn phi(&self) -> N64 {
        self.phi
    }

    /// Rapidity η
    pub fn rap(&self) -> N64 {
        self.rap
    }

    /// Inverse square of transverse momentum `inv_pt2 = 1/pt2`
    pub fn inv_pt2(&self) -> N64 {
        self.inv_pt2
    }

    /// Square of transverse momentum `pt2 = px*px + py*py`
    pub fn pt2(&self) -> N64 {
        n64(1.) / self.inv_pt2
    }

    /// Calculate ΔR^2 = Δφ^2 + Δη^2
    pub fn delta_r2(&self, p: &PseudoJet) -> N64 {
        self.delta_phi2(p) + self.delta_rap2(p)
    }

    /// Calculate ΔR = (Δφ^2 + Δη^2)^(1/2)
    pub fn delta_r(&self, p: &PseudoJet) -> N64 {
        self.delta_r2(p).sqrt()
    }

    /// Square Δφ^2 of azimuthal angle difference
    pub fn delta_phi2(&self, p: &PseudoJet) -> N64 {
        let dphi = self.delta_phi_abs(p);
        dphi * dphi
    }

    /// Difference Δφ in azimuthal angle
    ///
    /// The difference is normalised such that -π < Δφ <= π
    pub fn delta_phi(&self, p: &PseudoJet) -> N64 {
        let mut dphi = self.phi() - p.phi();
        if dphi > PI {
            dphi -= 2. * PI;
        } else if dphi <= -PI {
            dphi += 2. * PI;
        }
        debug_assert!(dphi > -PI);
        debug_assert!(dphi <= PI);
        dphi
    }

    /// Absolute difference |Δφ| in azimuthal angle
    ///
    /// The difference is normalised such that 0 <= |Δφ| <= π
    pub fn delta_phi_abs(&self, p: &PseudoJet) -> N64 {
        let mut abs_dphi = (self.phi() - p.phi()).abs();
        if abs_dphi > PI {
            abs_dphi = n64(2. * PI - f64::from(abs_dphi));
        }
        debug_assert!(abs_dphi >= 0.);
        debug_assert!(abs_dphi <= PI);
        abs_dphi
    }

    /// Square Δφ^2 of azimuthal angle difference
    pub fn delta_rap2(&self, p: &PseudoJet) -> N64 {
        let drap = self.delta_rap(p);
        drap * drap
    }

    /// Difference Δη in rapidity
    pub fn delta_rap(&self, p: &PseudoJet) -> N64 {
        self.rap() - p.rap()
    }

    fn init_pt2_phi_rap(&mut self) {
        let e = self[0];
        let px = self[1];
        let py = self[2];
        let pz = self[3];

        // initialisation taken from fastjet
        let pt2 = px * px + py * py;
        self.inv_pt2 = n64(1.) / pt2;

        self.phi = if pt2 > 0. { py.atan2(px) } else { n64(0.) };
        if self.phi < 0. {
            self.phi += n64(2.) * PI;
        }
        if self.phi > n64(2.) * PI {
            self.phi -= n64(2.) * PI;
        }

        self.rap = if e == 0. && pz == 0. {
            n64(0.)
        } else {
            ((e + pz) / (e - pz)).ln() / 2.
        }
    }

}

/// Create a pseudojet from the four-momentum components
impl From<[N64; D]> for PseudoJet {
    fn from(arr: [N64; D]) -> Self {
        let mut res = Self::new();
        res.comp = arr;
        res.init_pt2_phi_rap();
        res
    }
}

/// Create a pseudojet from the four-momentum components
impl From<[f64; D]> for PseudoJet {
    fn from(arr: [f64; D]) -> Self {
        let mut arr_n64 = [n64(0.); D];
        for i in 0..D {
            arr_n64[i] = n64(arr[i])
        }
        arr_n64.into()
    }
}

/// Create a pseudojet from the four-momentum components
impl From<(N64, N64, N64, N64)> for PseudoJet {
    fn from(p: (N64, N64, N64, N64)) -> Self {
        let mut res = Self::new();
        let (e, px, py, pz) = p;
        res.comp = [e, px, py, pz];
        res.init_pt2_phi_rap();
        res
    }
}

/// Create a pseudojet from the four-momentum components
impl From<(f64, f64, f64, f64)> for PseudoJet {
    fn from(p: (f64, f64, f64, f64)) -> Self {
        let (e, px, py, pz) = p;
        [n64(e), n64(px), n64(py), n64(pz)].into()
    }
}

/// Create a pseudojet from the four-momentum components
pub fn pseudojet(e: N64, px: N64, py: N64, pz: N64) -> PseudoJet {
    [e, px, py, pz].into()
}

/// Create a pseudojet from the four-momentum components
pub fn pseudojet_f(e: f64, px: f64, py: f64, pz: f64) -> PseudoJet {
    pseudojet(n64(e), n64(px), n64(py), n64(pz))
}

impl Default for PseudoJet {
    fn default() -> Self {
        PseudoJet {
            comp: Default::default(),
            inv_pt2: n64(f64::INFINITY),
            phi: Default::default(),
            rap: Default::default(),
        }
    }
}

impl Index<usize> for PseudoJet {
    type Output = N64;

    fn index(&self, i: usize) -> &Self::Output {
        &self.comp[i]
    }
}

impl AddAssign for PseudoJet {
    fn add_assign(&mut self, other: PseudoJet) {
        for i in 0..D {
            self.comp[i] += other.comp[i]
        }
        self.init_pt2_phi_rap()
    }
}

impl Add for PseudoJet {
    type Output = Self;

    fn add(mut self, other: PseudoJet) -> Self::Output {
        self += other;
        self
    }
}

impl SubAssign for PseudoJet {
    fn sub_assign(&mut self, other: PseudoJet) {
        for i in 0..D {
            self.comp[i] -= other.comp[i]
        }
        self.init_pt2_phi_rap()
    }
}

impl Sub for PseudoJet {
    type Output = Self;

    fn sub(mut self, other: PseudoJet) -> Self::Output {
        self -= other;
        self
    }
}
