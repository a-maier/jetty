use std::convert::From;
use std::default::Default;
use std::ops::{Add, AddAssign, Sub, SubAssign, Index};
use std::f64::consts::PI;

use noisy_float::prelude::*;

pub const D: usize = 4;

#[derive(Copy,Clone,Eq,PartialEq,Ord,PartialOrd,Hash,Debug)]
pub struct PseudoJet {
    comp: [N64; D],
    inv_pt2: N64,
    phi: N64,
    rap: N64,
}

impl PseudoJet {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn e(&self) -> N64 {
        self[0]
    }

    pub fn px(&self) -> N64  {
        self[1]
    }

    pub fn py(&self) -> N64  {
        self[2]
    }

    pub fn pz(&self) -> N64  {
        self[3]
    }

    pub fn phi(&self) -> N64 {
        self.phi
    }

    pub fn rap(&self) -> N64  {
        self.rap
    }

    pub fn inv_pt2(&self) -> N64  {
        self.inv_pt2
    }

    fn init_pt2_phi_rap(&mut self) {
        let e = self[0];
        let px = self[1];
        let py = self[2];
        let pz = self[3];

        // initialisation taken from fastjet
        let pt2 = px*px + py*py;
        self.inv_pt2 = n64(1.) / pt2;

        self.phi = if pt2 > 0. {
            py.atan2(px)
        } else {
            n64(0.)
        };
        if self.phi < 0. {
            self.phi += n64(2.)*PI;
        }
        if self.phi > n64(2.)*PI {
            self.phi -= n64(2.)*PI;
        }

        self.rap = if e == 0. {
            if pz == 0. {
                n64(0.)
            } else {
                let rat = e/pz;
                ((rat + 1.)/(rat - 1.)).ln()/2.
            }
        } else {
            let rat = pz/e;
            ((rat + 1.)/(n64(1.) - rat)).ln()/2.
        };
    }
}

impl From<[N64; D]> for PseudoJet {
    fn from(arr: [N64; D]) -> Self {
        let mut res = Self::new();
        res.comp = arr;
        res.init_pt2_phi_rap();
        res
    }
}

impl From<(N64, N64, N64, N64)> for PseudoJet {
    fn from(p: (N64, N64, N64, N64)) -> Self {
        let mut res = Self::new();
        let (e, px, py, pz) = p;
        res.comp = [e, px, py, pz];
        res.init_pt2_phi_rap();
        res
    }
}

pub fn pseudojet(e: N64, px: N64, py: N64, pz: N64) -> PseudoJet {
    [e, px, py, pz].into()
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
