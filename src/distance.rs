use crate::pseudojet::PseudoJet;

use std::cmp::min;

use noisy_float::prelude::*;

/// Distance measure between pseudojets used for clustering
pub trait Distance {
    /// Distance between pseudojets
    fn distance(&self, p1: &PseudoJet, p2: &PseudoJet) -> N64;
    /// Distance to the beam axis
    fn beam_distance(&self, p1: &PseudoJet) -> N64;
}

pub struct AntiKt {
    r2: N64,
}

/// anti-kt distance measure with radius parameter `r`
pub fn anti_kt(r: N64) -> AntiKt {
    AntiKt { r2: r * r }
}

/// anti-kt distance measure with radius parameter `r`
pub fn anti_kt_f(r: f64) -> AntiKt {
    anti_kt(n64(r))
}

impl Distance for AntiKt {
    fn distance(&self, p1: &PseudoJet, p2: &PseudoJet) -> N64 {
        min(p1.inv_pt2(), p2.inv_pt2()) * p1.delta_r2(p2) / self.r2
    }

    fn beam_distance(&self, p1: &PseudoJet) -> N64 {
        p1.inv_pt2()
    }
}

pub struct Kt {
    r2: N64,
}

/// kt distance measure with radius parameter `r`
pub fn kt(r: N64) -> Kt {
    Kt { r2: r * r }
}

/// kt distance measure with radius parameter `r`
pub fn kt_f(r: f64) -> Kt {
    kt(n64(r))
}

impl Distance for Kt {
    fn distance(&self, p1: &PseudoJet, p2: &PseudoJet) -> N64 {
        min(p1.pt2(), p2.pt2()) * p1.delta_r2(p2) / self.r2
    }

    fn beam_distance(&self, p1: &PseudoJet) -> N64 {
        p1.pt2()
    }
}

pub struct CambridgeAachen {
    r2: N64,
}

/// Cambridge/Aachen distance measure with radius parameter `r`
pub fn cambridge_aachen(r: N64) -> CambridgeAachen {
    CambridgeAachen { r2: r * r }
}

/// Cambridge/Aachen distance measure with radius parameter `r`
pub fn cambridge_aachen_f(r: f64) -> CambridgeAachen {
    cambridge_aachen(n64(r))
}

impl Distance for CambridgeAachen {
    fn distance(&self, p1: &PseudoJet, p2: &PseudoJet) -> N64 {
        p1.delta_r2(p2) / self.r2
    }

    fn beam_distance(&self, _p1: &PseudoJet) -> N64 {
        n64(1.)
    }
}

pub struct GenKt {
    r2: N64,
    p: N64,
}

/// Generalised kt distance measure with radius parameter `r` and exponent `p`
pub fn gen_kt(r: N64, p: N64) -> GenKt {
    GenKt { r2: r * r, p }
}

/// Generalised kt distance measure with radius parameter `r` and exponent `p`
pub fn gen_kt_f(r: f64, p: f64) -> GenKt {
    gen_kt(n64(r), n64(p))
}

impl Distance for GenKt {
    fn distance(&self, p1: &PseudoJet, p2: &PseudoJet) -> N64 {
        min(p1.pt2().powf(self.p), p2.pt2().powf(self.p)) * p1.delta_r2(p2)
            / self.r2
    }

    fn beam_distance(&self, p1: &PseudoJet) -> N64 {
        p1.pt2().powf(self.p)
    }
}

impl<T: Distance> Distance for &T {
    fn distance(&self, p1: &PseudoJet, p2: &PseudoJet) -> N64 {
        (*self).distance(p1, p2)
    }

    fn beam_distance(&self, p1: &PseudoJet) -> N64 {
        (*self).beam_distance(p1)
    }
}
