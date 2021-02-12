use crate::pseudojet::PseudoJet;

use std::cmp::min;

use noisy_float::prelude::*;

pub trait Distance {
    fn distance(&self, p1: &PseudoJet, p2: &PseudoJet) -> N64;
    fn beam_distance(&self, p1: &PseudoJet) -> N64;
}

pub struct AntiKt {
    r2: N64
}

pub fn anti_kt(r: N64) -> AntiKt {
    AntiKt{r2: r*r}
}

impl Distance for AntiKt {
    fn distance(&self, p1: &PseudoJet, p2: &PseudoJet) -> N64 {
        let dy = p1.rap() - p2.rap();
        let dphi = p1.phi() - p2.phi();

        let delta_sq = dy*dy + dphi*dphi;

        min(p1.inv_pt2(), p2.inv_pt2()) * delta_sq/self.r2
    }

    fn beam_distance(&self, p1: &PseudoJet) -> N64 {
        p1.inv_pt2()
    }
}

pub struct Kt {
    r2: N64
}

pub fn kt(r: N64) -> Kt {
    Kt{r2: r*r}
}

impl Distance for Kt {
    fn distance(&self, p1: &PseudoJet, p2: &PseudoJet) -> N64 {
        let dy = p1.rap() - p2.rap();
        let dphi = p1.phi() - p2.phi();

        let delta_sq = dy*dy + dphi*dphi;

        min(p1.pt2(), p2.pt2()) * delta_sq/self.r2
    }

    fn beam_distance(&self, p1: &PseudoJet) -> N64 {
        p1.pt2()
    }
}

pub struct CambridgeAachen {
    r2: N64
}

pub fn cambridge_aachen(r: N64) -> CambridgeAachen {
    CambridgeAachen{r2: r*r}
}

impl Distance for CambridgeAachen {
    fn distance(&self, p1: &PseudoJet, p2: &PseudoJet) -> N64 {
        let dy = p1.rap() - p2.rap();
        let dphi = p1.phi() - p2.phi();

        let delta_sq = dy*dy + dphi*dphi;

        delta_sq/self.r2
    }

    fn beam_distance(&self, _p1: &PseudoJet) -> N64 {
        n64(1.)
    }
}

pub struct GenKt {
    r2: N64,
    p: N64,
}

pub fn gen_kt(r: N64, p: N64) -> GenKt {
    GenKt{r2: r*r, p}
}

impl Distance for GenKt {
    fn distance(&self, p1: &PseudoJet, p2: &PseudoJet) -> N64 {
        let dy = p1.rap() - p2.rap();
        let dphi = p1.phi() - p2.phi();

        let delta_sq = dy*dy + dphi*dphi;

        min(p1.pt2().powf(self.p), p2.pt2().powf(self.p)) * delta_sq/self.r2
    }

    fn beam_distance(&self, p1: &PseudoJet) -> N64 {
        p1.pt2().powf(self.p)
    }
}
