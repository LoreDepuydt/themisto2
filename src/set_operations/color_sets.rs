//! Generating the color set elements of the key k-mers of the result of a set operation, from the
//! color sets of the two input indexes, for building the color set storage of the result.

use std::{ops::Range, sync::Arc};

use rayon::iter::{IntoParallelIterator, ParallelIterator};
use sbwt::MergeInterleaving;
use simple_sds_sbwt::ops::{BitVec, Rank, Select};

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetStorage, ColorSetView}, set_of_sets_construction::{ParallelElementGenerator, SetElement}};

pub struct ElementGeneratorFromMergeInterleaving<'a, CSS: ColorSetStorage + Sync + Send> {
    pub interleaving: &'a MergeInterleaving,
    pub coloring1: &'a CompactColexKmers<CSS>,
    pub coloring2: &'a CompactColexKmers<CSS>,
    pub merged_key_kmer_marks: &'a bitvec::vec::BitVec, // Only reporting set elements for these
    // Marks the positions of the interleaving that were removed from the merged SBWT (redundant dummy nodes).
    // The remaining positions keep their order, so the j-th unmarked position is merged colex position j.
    // None if nothing was removed, in which case interleaving positions are merged colex positions.
    pub removed_positions: Option<&'a bitvec::vec::BitVec<u64, bitvec::order::Lsb0>>,
    pub filter: Option<Arc<simple_sds_sbwt::bit_vector::BitVector>>, // With rank support
    // Maps each color of coloring2 to its color id in the merged coloring. Ids smaller than
    // the number of colors of coloring1 are colors shared with coloring1.
    pub color2_to_merged: &'a [usize],
}

struct ThreadInput {
    interleaving_range: Range<usize>,
    merged_start_colex: usize, // Merged colex position of the first position of interleaving_range that was not removed
    s1_start_rank: usize,
    s2_start_rank: usize,
} 

impl<'a, CSS: ColorSetStorage + Sync + Send> ParallelElementGenerator for ElementGeneratorFromMergeInterleaving<'a, CSS> {


    fn run(&mut self, callback: impl Fn(SetElement) + Send + Sync, n_threads: usize) {
        assert!(self.interleaving.s1.len() == self.interleaving.s2.len());
        let n = self.interleaving.s1.len();

        let s1 = &self.interleaving.s1;
        let s2 = &self.interleaving.s2;

        let thread_ranges = crate::util::segment_range(0..n, n_threads);
        let mut thread_inputs = Vec::<ThreadInput>::with_capacity(n_threads);

        let removed = self.removed_positions;
        let n_removed = removed.map_or(0, |r| r.count_ones());
        assert_eq!(self.merged_key_kmer_marks.len(), n - n_removed);

        let mut n_bits_s1 = 0_usize;
        let mut n_bits_s2 = 0_usize;
        let mut n_kept = 0_usize;
        for range in thread_ranges.iter() {
            let input = ThreadInput {
                interleaving_range: range.clone(),
                merged_start_colex: n_kept,
                s1_start_rank: n_bits_s1,
                s2_start_rank: n_bits_s2,
            };
            thread_inputs.push(input);
            n_bits_s1 += s1[range.clone()].count_ones();
            n_bits_s2 += s2[range.clone()].count_ones();
            n_kept += range.len() - removed.map_or(0, |r| r[range.clone()].count_ones());
        }

        let n_colors1 = self.coloring1.get_set_storage().n_colors();
        assert_eq!(self.color2_to_merged.len(), self.coloring2.get_set_storage().n_colors());
        // Marks the colors of coloring1 that coloring2 also has. Merged ids below n_colors1 are
        // exactly the colors of coloring1.
        let mut is_shared1 = bitvec::bitvec![0; n_colors1];
        for &color in self.color2_to_merged.iter().filter(|&&c| c < n_colors1) {
            is_shared1.set(color, true);
        }
        let has_shared_colors = is_shared1.any();

        let bar = indicatif::ProgressBar::new(n as u64);
        thread_inputs.into_par_iter().for_each(|input| { //TODO: not sure if some check are overkill
            let mut s1_colex = input.s1_start_rank;
            let mut s2_colex = input.s2_start_rank;
            let mut merged_colex = input.merged_start_colex;
            // The shared colors of coloring1 at the current k-mer, as a bitmap for O(1) lookups
            // and as a list for clearing the bitmap afterwards. Not allocated without shared colors.
            let mut in_set1 = bitvec::bitvec![0; if has_shared_colors { n_colors1 } else { 0 }];
            let mut in_set1_list = Vec::<usize>::new();
            for pos in input.interleaving_range {
                if pos > 0 && pos % 10000 == 0 {
                    bar.inc(10000);
                }
                let is_removed = removed.is_some_and(|r| r[pos]);
                if !is_removed && self.merged_key_kmer_marks[merged_colex] {
                    if let Some(new_set_id) = maybe_apply_filter(self.filter.as_deref(), merged_colex) {
                        let in_1 = self.interleaving.s1[pos];
                        let in_2 = self.interleaving.s2[pos];
                        // A shared color may be in both sets. It must be reported only once.
                        let dedup = has_shared_colors && in_1 && in_2;

                        if in_1 {
                            for color in self.coloring1.colex_to_set(s1_colex).iter() {
                                callback(SetElement{set_id: new_set_id, color});
                                if dedup && is_shared1[color] {
                                    in_set1.set(color, true);
                                    in_set1_list.push(color);
                                }
                            }
                        }

                        if in_2 {
                            for color in self.coloring2.colex_to_set(s2_colex).iter() {
                                let color = self.color2_to_merged[color];
                                if dedup && color < n_colors1 && in_set1[color] {
                                    continue; // Already reported from coloring1
                                }
                                callback(SetElement{set_id: new_set_id, color});
                            }
                        }

                        for color in in_set1_list.drain(..) { //TODO: not sure about this
                            in_set1.set(color, false);
                        }
                    }
                }
                s1_colex += self.interleaving.s1[pos] as usize;
                s2_colex += self.interleaving.s2[pos] as usize;
                merged_colex += !is_removed as usize; // Removed positions do not exist in the merged SBWT
            }
        });
        bar.finish();
    }

    fn set_filter(&mut self, filter: Arc<simple_sds_sbwt::bit_vector::BitVector>) {
        self.filter = Some(filter.clone());
    }

    fn rewind(&mut self) {
        // Nothing needs to done, calling run() again already works
    }
}

/// How the color sets of a k-mer in the two input indexes are combined into its color set in the
/// intersection of the indexes.
#[derive(Clone, Copy)]
pub enum ColorCombination<'a> {
    /// The union of the two sets. Maps each color of coloring2 to its color id in the result, as
    /// color2_to_merged in ElementGeneratorFromMergeInterleaving: the colors of coloring1 keep
    /// their ids, and ids smaller than the number of colors of coloring1 are shared colors.
    Union { color2_to_result: &'a [usize] },
    /// The intersection of the two sets. Maps each color of coloring1 and coloring2 to its color
    /// id in the result, or to None if the color is not in the result.
    Intersection { color1_to_result: &'a [Option<usize>], color2_to_result: &'a [Option<usize>] },
}

/// Generates the color set elements of the key k-mers of the intersection of the SBWTs of
/// coloring1 and coloring2, computed with the given interleaving.
pub struct ElementGeneratorFromIntersectionInterleaving<'a, CSS: ColorSetStorage + Sync + Send> {
    pub interleaving: &'a MergeInterleaving,
    pub coloring1: &'a CompactColexKmers<CSS>,
    pub coloring2: &'a CompactColexKmers<CSS>,
    // Marks the colex positions of the result SBWT that are k-mers, with select support. The
    // intersection can add and remove dummy nodes, but the k-mers are those of the interleaving
    // positions in both inputs that are not dummies, in the same order, so the j-th such position
    // is the result k-mer at select(j).
    pub result_kmers: &'a simple_sds_sbwt::bit_vector::BitVector,
    pub result_key_kmer_marks: &'a bitvec::vec::BitVec, // Only reporting set elements for these
    pub filter: Option<Arc<simple_sds_sbwt::bit_vector::BitVector>>, // With rank support
    pub colors: ColorCombination<'a>,
    pub n_result_colors: usize,
}

struct IntersectionThreadInput {
    interleaving_range: Range<usize>,
    s1_start_rank: usize,
    s2_start_rank: usize,
    result_kmer_start_rank: usize,
}

impl<'a, CSS: ColorSetStorage + Sync + Send> ParallelElementGenerator for ElementGeneratorFromIntersectionInterleaving<'a, CSS> {

    fn run(&mut self, callback: impl Fn(SetElement) + Send + Sync, n_threads: usize) {
        let interleaving = self.interleaving;
        assert!(interleaving.s1.len() == interleaving.s2.len());
        let n = interleaving.s1.len();
        assert_eq!(self.result_key_kmer_marks.len(), self.result_kmers.len());

        let is_result_kmer = |pos: usize| interleaving.s1[pos] && interleaving.s2[pos] && !interleaving.is_dummy[pos];

        // The number of bits of s1 and s2 and of result k-mers in each piece, counted in parallel
        let ranges = crate::util::segment_range(0..n, n_threads);
        let counts: Vec<(usize, usize, usize)> = ranges.clone().into_par_iter().map(|range| {
            (interleaving.s1[range.clone()].count_ones(), interleaving.s2[range.clone()].count_ones(), count_result_kmers(interleaving, range))
        }).collect();
        let mut thread_inputs = Vec::<IntersectionThreadInput>::with_capacity(n_threads);
        let (mut n_bits_s1, mut n_bits_s2, mut n_result_kmers) = (0_usize, 0_usize, 0_usize);
        for (range, (c1, c2, c3)) in ranges.into_iter().zip(counts) {
            thread_inputs.push(IntersectionThreadInput {
                interleaving_range: range,
                s1_start_rank: n_bits_s1,
                s2_start_rank: n_bits_s2,
                result_kmer_start_rank: n_result_kmers,
            });
            n_bits_s1 += c1;
            n_bits_s2 += c2;
            n_result_kmers += c3;
        }
        assert_eq!(n_result_kmers, self.result_kmers.count_ones());

        let n_colors1 = self.coloring1.get_set_storage().n_colors();
        let n_colors2 = self.coloring2.get_set_storage().n_colors();
        // The colors of the result that come from coloring1 and may also come from coloring2, so
        // that they must be reported only once (union) or only if they come from both (intersection).
        let mut is_shared = bitvec::bitvec![0; self.n_result_colors];
        match self.colors {
            ColorCombination::Union { color2_to_result } => {
                assert_eq!(color2_to_result.len(), n_colors2);
                for &color in color2_to_result.iter().filter(|&&c| c < n_colors1) {
                    is_shared.set(color, true);
                }
            },
            ColorCombination::Intersection { color1_to_result, color2_to_result } => {
                assert_eq!(color1_to_result.len(), n_colors1);
                assert_eq!(color2_to_result.len(), n_colors2);
                is_shared.fill(true);
            },
        }
        let has_shared_colors = is_shared.any();

        let bar = indicatif::ProgressBar::new(n as u64);
        thread_inputs.into_par_iter().for_each(|input| {
            let mut s1_colex = input.s1_start_rank;
            let mut s2_colex = input.s2_start_rank;
            // The colex positions of the result k-mers, in order, from the first one in this piece
            let mut result_kmer_positions = self.result_kmers.select_iter(input.result_kmer_start_rank);
            // The shared colors of coloring1 at the current k-mer, as a bitmap for O(1) lookups
            // and as a list for clearing the bitmap afterwards. Not allocated without shared colors.
            let mut in_set1 = bitvec::bitvec![0; if has_shared_colors { self.n_result_colors } else { 0 }];
            let mut in_set1_list = Vec::<usize>::new();
            for pos in input.interleaving_range {
                if pos > 0 && pos % 10000 == 0 {
                    bar.inc(10000);
                }
                if is_result_kmer(pos) {
                    let (_, result_colex) = result_kmer_positions.next().unwrap();
                    if self.result_key_kmer_marks[result_colex] {
                        if let Some(new_set_id) = maybe_apply_filter(self.filter.as_deref(), result_colex) {
                            let set1 = self.coloring1.colex_to_set(s1_colex);
                            let set2 = self.coloring2.colex_to_set(s2_colex);
                            match self.colors {
                                ColorCombination::Union { color2_to_result } => {
                                    for color in set1.iter() {
                                        callback(SetElement{set_id: new_set_id, color});
                                        if has_shared_colors && is_shared[color] {
                                            in_set1.set(color, true);
                                            in_set1_list.push(color);
                                        }
                                    }
                                    for color in set2.iter() {
                                        let color = color2_to_result[color];
                                        if !(color < n_colors1 && has_shared_colors && in_set1[color]) {
                                            callback(SetElement{set_id: new_set_id, color});
                                        } // Else already reported from coloring1
                                    }
                                },
                                ColorCombination::Intersection { color1_to_result, color2_to_result } => {
                                    for color in set1.iter().filter_map(|c| color1_to_result[c]) {
                                        in_set1.set(color, true);
                                        in_set1_list.push(color);
                                    }
                                    for color in set2.iter().filter_map(|c| color2_to_result[c]) {
                                        if in_set1[color] {
                                            callback(SetElement{set_id: new_set_id, color});
                                        }
                                    }
                                },
                            }
                            for color in in_set1_list.drain(..) {
                                in_set1.set(color, false);
                            }
                        }
                    }
                }
                s1_colex += interleaving.s1[pos] as usize;
                s2_colex += interleaving.s2[pos] as usize;
            }
        });
        bar.finish();
    }

    fn set_filter(&mut self, filter: Arc<simple_sds_sbwt::bit_vector::BitVector>) {
        self.filter = Some(filter.clone());
    }

    fn rewind(&mut self) {
        // Nothing needs to done, calling run() again already works
    }
}

/// The number of positions in `range` of the interleaving that are k-mers of the intersection: in
/// both inputs and not dummies. Counted a word at a time.
fn count_result_kmers(interleaving: &MergeInterleaving, range: Range<usize>) -> usize {
    if range.is_empty() {
        return 0;
    }
    let (s1, s2, is_dummy) = (interleaving.s1.as_raw_slice(), interleaving.s2.as_raw_slice(), interleaving.is_dummy.as_raw_slice());
    let (first_word, last_word) = (range.start / 64, (range.end - 1) / 64);
    (first_word..=last_word).map(|w| {
        let mut word = s1[w] & s2[w] & !is_dummy[w];
        if w == first_word {
            word &= !0_u64 << (range.start % 64); // Only bits from range.start on
        }
        if w == last_word && range.end % 64 != 0 {
            word &= (1_u64 << (range.end % 64)) - 1; // Only bits before range.end
        }
        word.count_ones() as usize
    }).sum()
}

/// Returns the id of the set of the given colex position after filtering: its rank among the
/// positions marked in the filter, or None if it is not marked. Without a filter, the colex
/// position itself.
fn maybe_apply_filter(filter: Option<&simple_sds_sbwt::bit_vector::BitVector>, colex: usize) -> Option<usize> {
    match filter {
        Some(filter) if !filter.get(colex) => None, // Do not report this
        Some(filter) => Some(filter.rank(colex)), // Assign new id
        None => Some(colex), // No filter
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn count_result_kmers_matches_counting_bit_by_bit() {
        use rand_chacha::rand_core::{RngCore, SeedableRng};
        let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(1);
        let n = 1000;
        let mut random_bits = || -> bitvec::vec::BitVec<u64, bitvec::order::Lsb0> { (0..n).map(|_| rng.next_u64() % 3 != 0).collect() };
        let interleaving = sbwt::MergeInterleaving { s1: random_bits(), s2: random_bits(), is_dummy: random_bits(), is_leader: random_bits() };
        let expected = |range: std::ops::Range<usize>| range.filter(|&p| interleaving.s1[p] && interleaving.s2[p] && !interleaving.is_dummy[p]).count();
        // Ranges within one word, across word boundaries, ending at a word boundary, and empty
        for range in [0..0, 0..1, 3..60, 0..64, 64..128, 60..70, 1..1000, 63..65, 500..500, 999..1000, 0..1000] {
            assert_eq!(super::count_result_kmers(&interleaving, range.clone()), expected(range.clone()), "range {:?}", range);
        }
    }
}
