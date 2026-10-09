//! Where the k-mers of the input indexes of a set operation are in its result.

use simple_sds_sbwt::ops::{BitVec, Rank, Select};

use crate::util::u64_bitvec_to_simple_sds;

/// Maps the colex positions of the k-mers of one input index to their colex positions in the merged
/// index, through the merge interleaving: input position i is at interleaving position
/// select(s, i), where s is the input's bit vector of the interleaving (s1 or s2), and that
/// interleaving position is at merged position select(s, i) minus the number of removed
/// interleaving positions before it. Only dummy nodes are removed, so this is valid for k-mers.
pub(crate) struct ToMerged<'a> {
    in_input: simple_sds_sbwt::bit_vector::BitVector, // s1 or s2 of the interleaving, with select support
    removed: Option<&'a simple_sds_sbwt::bit_vector::BitVector>, // With rank support. None if nothing was removed
}

impl<'a> ToMerged<'a> {
    pub(crate) fn new(in_input: &bitvec::vec::BitVec<u64, bitvec::order::Lsb0>, removed: Option<&'a simple_sds_sbwt::bit_vector::BitVector>) -> Self {
        // TODO: This copies s1 or s2 (n bits for n interleaving positions) only to get select
        // support from simple-sds. A select structure directly on the bitvec of the interleaving
        // would save that copy, at the cost of our own select code (sbwt had one, ForwardSelect,
        // before it switched to simple-sds).
        let mut in_input = u64_bitvec_to_simple_sds(in_input);
        in_input.enable_select();
        Self { in_input, removed }
    }

    pub(crate) fn merged_colex(&self, input_colex: usize) -> usize {
        let pos = self.in_input.select(input_colex).unwrap();
        pos - self.removed.map_or(0, |r| r.rank(pos))
    }
}

/// Maps the colex positions of the k-mers of one input index to their colex positions in the
/// intersection, through the interleaving of the inputs: input position i is at interleaving
/// position p = select(s, i), where s is the input's bit vector of the interleaving (s1 or s2).
/// The k-mer is in the intersection iff p is in both inputs and is not a dummy, and the
/// intersection k-mers are in the same order as those interleaving positions. The colex positions
/// themselves differ, because the intersection can add and remove dummy nodes.
pub(crate) struct ToIntersection<'a> {
    in_input: simple_sds_sbwt::bit_vector::BitVector, // s1 or s2 of the interleaving, with select support
    in_result: &'a simple_sds_sbwt::bit_vector::BitVector, // Interleaving positions of the result k-mers, with rank support
    result_kmers: &'a simple_sds_sbwt::bit_vector::BitVector, // Colex positions of the result k-mers, with select support
}

impl<'a> ToIntersection<'a> {
    pub(crate) fn new(in_input: &bitvec::vec::BitVec<u64, bitvec::order::Lsb0>, in_result: &'a simple_sds_sbwt::bit_vector::BitVector, result_kmers: &'a simple_sds_sbwt::bit_vector::BitVector) -> Self {
        // TODO: This copies s1 or s2 only to get select support, as in ToMerged.
        let mut in_input = u64_bitvec_to_simple_sds(in_input);
        in_input.enable_select();
        Self { in_input, in_result, result_kmers }
    }

    pub(crate) fn result_colex(&self, input_colex: usize) -> Option<usize> {
        let pos = self.in_input.select(input_colex).unwrap();
        if !self.in_result.get(pos) {
            return None;
        }
        Some(self.result_kmers.select(self.in_result.rank(pos)).unwrap())
    }
}
