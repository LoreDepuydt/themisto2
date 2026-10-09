use std::{collections::HashSet, ops::Range, sync::Arc};

use rayon::iter::{IntoParallelIterator, ParallelIterator};
use sbwt::{LcsArray, MergeInterleaving, SbwtIndex, SeqStream, StreamingIndex, SubsetMatrix, reverse_complement_in_place};
use simple_sds_sbwt::ops::{BitVec, Rank, Select};

use crate::{colex_colored_kmers::CompactColexKmers, coloring_interface::{ColorSetStorage, ColorSetView}, io::{self, RewindableSeqStreamGenerator}, set_of_sets_construction::{ParallelElementGenerator, SetElement}};

pub struct MsElementGenerator<'a> {
    color_stream_generator: Box<dyn RewindableSeqStreamGenerator + Sync + Send>,
    streaming_index: StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>,
    filter: Option<Arc<simple_sds_sbwt::bit_vector::BitVector>>,
    include_rev_comp: bool,
    n_parser_threads: usize,
}

impl<'a> MsElementGenerator<'a> {
    pub fn new(color_stream_generator: Box<dyn RewindableSeqStreamGenerator + Sync + Send>, streaming_index: StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>, include_rev_comp: bool, n_parser_threads: usize) -> Self {
        Self {
            color_stream_generator,
            streaming_index,
            filter: None,
            include_rev_comp,
            n_parser_threads,
        }
    }
}

pub struct MsWorker<'a, CB: Fn(crate::set_of_sets_construction::SetElement) + Send + Sync> {
    streaming_index: &'a StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>,
    filter: &'a Option<Arc<simple_sds_sbwt::bit_vector::BitVector>>,
    include_rev_comp: bool,
    rev_comp_buf: Vec<u8>,
    callback: &'a CB, 
}

impl<'a, CB: Fn(crate::set_of_sets_construction::SetElement) + Send + Sync> crate::work_dispatcher::Worker for MsWorker<'a, CB> {
    fn process(&mut self, seq: &[u8], color: usize) {
        self.process_internal(seq, color);
        if self.include_rev_comp {
            self.rev_comp_buf.clear();
            self.rev_comp_buf.extend_from_slice(seq);
            reverse_complement_in_place(&mut self.rev_comp_buf);
            self.process_internal(&self.rev_comp_buf, color);
        }
    }
}

impl<'a, CB: Fn(crate::set_of_sets_construction::SetElement) + Send + Sync> MsWorker<'a, CB> {

    // Don't call this directly. This is called from the process-function of crate::work_dispatcher::Worker.
    fn process_internal(&self, seq: &[u8], color: usize) {
        let k = self.streaming_index.k();
        let ms_iter = self.streaming_index.matching_statistics_iter(seq);
        let kmer_iter = ms_iter.skip(k-1).filter(|(len, _colex)| *len == k);
        let filtered_iter = kmer_iter.filter_map(|(_, colex)| {
            assert!(colex.len() == 1);
            let set_id = colex.start;
            if let Some(filter) = &self.filter {
                if !filter.get(set_id) {
                    None // Do not report this
                } else {
                    // Assign new id
                    let new_id = filter.rank(set_id);
                    Some(new_id)
                }
            } else {
                Some(set_id) // No filter
            }
        });

        for id in filtered_iter {
            (self.callback)(SetElement{
                set_id: id,
                color,
            });
        }
    }
}

impl<'a> crate::set_of_sets_construction::ParallelElementGenerator for MsElementGenerator<'a> {
    fn run(&mut self, callback: impl Fn(crate::set_of_sets_construction::SetElement) + Send + Sync, n_threads: usize) {

        // Here we need to get a bit tricky to avoid mutable aliasing of self. The issue is that
        // the work dispatcher needs a mutable reference to the ReWindableSeqStreamGenerator at self,
        // while the workers need non-mutable access to the rest of self. This is not possible
        // at the same time because borrowing self borrows everything. The workaround is that we
        // swap in a dummy generator into self, so that we can get separate ownership of the generator
        // and pass it into the producer. In the end we swap it back in.
        let mut dummy_color_stream_generator: Box<dyn RewindableSeqStreamGenerator + Sync + Send> = Box::new(io::EmptyRewindableSeqStreamGenerator{});
        std::mem::swap(&mut self.color_stream_generator, &mut dummy_color_stream_generator);
        let mut color_stream_generator = dummy_color_stream_generator; // Now this function owns this

        let workers: Vec<MsWorker<_>> = (0..n_threads).map(|_| {
            MsWorker {
                streaming_index: &self.streaming_index,
                filter: &self.filter,
                include_rev_comp: self.include_rev_comp,
                rev_comp_buf: vec![],
                callback: &callback,
            }
        }).collect();

        crate::work_dispatcher::dispatch_work(&mut color_stream_generator, workers, self.n_parser_threads, 1 << 23); 

        self.color_stream_generator = color_stream_generator; // Put this back in (see comment at the start of the function)
    }
    
    fn set_filter(&mut self, filter: Arc<simple_sds_sbwt::bit_vector::BitVector>) {
        self.filter = Some(filter);
    }

    fn rewind(&mut self) {
        self.color_stream_generator.rewind();
    }
}

struct DeduplicatingBuffer {
    universe_size: usize,
    hashset: Option<HashSet<usize>>,
    bitmap: Option<bitvec::vec::BitVec>, // Switch to this when the set is large enough
    empty: bool,
} 

impl DeduplicatingBuffer {
    fn new(universe_size: usize) -> Self {
        Self {
            universe_size,
            hashset: Some(HashSet::new()),
            bitmap: None,
            empty: true,
        }
    }

    fn insert(&mut self, id: usize) {
        self.empty = false;
        if let Some(hs) = &mut self.hashset {
            hs.insert(id);
            if hs.capacity() > self.universe_size / 64 {
                // Switch to bitmap
                let mut bv = bitvec::bitvec![0; self.universe_size];
                for &elem in hs.iter() {
                    bv.set(elem, true);
                }
                self.bitmap = Some(bv);
                self.hashset = None;
            }
        } else if let Some(bv) = &mut self.bitmap {
            bv.set(id, true);
        } else {
            panic!("Both hashset and bitmap are None");
        }
    }

    fn into_iter(self) -> DeduplicatingBufferIter {
        if let Some(hs) = self.hashset {
            DeduplicatingBufferIter {
                hashset_iter: Some(hs.into_iter()),
                bitmap_iter: None,
            }
        } else if let Some(bv) = self.bitmap {
            DeduplicatingBufferIter {
                hashset_iter: None,
                bitmap_iter: Some(OwningBitVecOnesIterator::new(bv)),
            }
        } else {
            panic!("Both hashset and bitmap are None");
        }
    }
}

struct OwningBitVecOnesIterator {
    bv: bitvec::vec::BitVec,
    pos: usize,
}

impl OwningBitVecOnesIterator {
    fn new(bv: bitvec::vec::BitVec) -> Self {
        Self { bv, pos: 0 }
    }
}

impl Iterator for OwningBitVecOnesIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(offset) = self.bv[self.pos..].first_one() {
            let ret = self.pos + offset;
            self.pos += offset + 1; // Starting point of the next iteration
            Some(ret)
        } else {
            None
        }
    }
}   

struct DeduplicatingBufferIter {
    hashset_iter: Option<std::collections::hash_set::IntoIter<usize>>,
    bitmap_iter: Option<OwningBitVecOnesIterator>,
}

impl Iterator for DeduplicatingBufferIter {
    type Item = usize;
    
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(hs_iter) = &mut self.hashset_iter {
            hs_iter.next()
        } else if let Some(bv_iter) = &mut self.bitmap_iter {
            bv_iter.next()
        } else {
            panic!("Both hashset_iter and bitmap_iter are None");
        }
    }
}

pub struct DistinctColexComputation<'a> {
    streaming_index: sbwt::StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>,
    input: Box<dyn SeqStream + Sync + Send>,
    set_ids: DeduplicatingBuffer,
}

impl<'a> DistinctColexComputation<'a> {
    pub fn new(sbwt: &'a sbwt::SbwtIndex<SubsetMatrix>, lcs: &'a LcsArray, input: Box<dyn SeqStream + Sync + Send>) -> Self {
        let streaming_index = StreamingIndex::new(sbwt, lcs); 
        Self {
            streaming_index,
            input,
            set_ids: DeduplicatingBuffer::new(sbwt.n_sets()),
        }
    }

    fn process_seq(&mut self, cur_seq: &[u8]) {
        let k = self.streaming_index.k();

        let ms_iter = self.streaming_index.matching_statistics_iter(cur_seq);
        for (_, colex) in ms_iter.skip(k-1).filter(|(len, _colex)| *len == k) {
            assert!(colex.len() == 1);
            self.set_ids.insert(colex.start);
        }
    }

    fn run(mut self, include_rev_comp: bool) -> DeduplicatingBuffer {
        let mut buf = Vec::<u8>::new();
        while let Some(seq) = self.input.stream_next() {
            buf.clear();
            buf.extend_from_slice(seq);
            self.process_seq(&buf);
            if include_rev_comp {
                reverse_complement_in_place(&mut buf);
                self.process_seq(&buf);
            }
        }
        self.set_ids
    }
}

pub struct DeduplicatingColorElementGenerator<'a> {
    color_stream_generator: Box<dyn RewindableSeqStreamGenerator + Sync + Send>,
    sbwt: &'a SbwtIndex<SubsetMatrix>,
    lcs: &'a LcsArray,
    filter: Option<Arc<simple_sds_sbwt::bit_vector::BitVector>>,
    include_rev_comp: bool,
}

impl<'a> DeduplicatingColorElementGenerator<'a> {
    pub fn new( sbwt: &'a SbwtIndex<SubsetMatrix>, lcs: &'a LcsArray, color_stream_generator: Box<dyn RewindableSeqStreamGenerator + Sync + Send>, include_rev_comp: bool) -> Self {
        Self { color_stream_generator, sbwt, lcs, filter: None, include_rev_comp }
    }
}

impl<'a> crate::set_of_sets_construction::ParallelElementGenerator for DeduplicatingColorElementGenerator<'a> {
    fn run(&mut self, callback: impl Fn(crate::set_of_sets_construction::SetElement) + Send + Sync, n_threads: usize) {
        // TODO: a lot of this code is duplicated with  MsElementGenerator.

        let (sender, receiver) = crossbeam::channel::bounded::<(usize, Box<dyn SeqStream + Send + Sync>)>(2*n_threads);
        let receiver_ref = &receiver; // To capture a reference

        // Here we need to get a bit tricky to avoid mutable aliasing of self. The issue is that
        // the producer thread needs a mutable reference to the ReWindableSeqStreamGenerator at self,
        // while the consumers need non-mutable access to the rest of self. This is not possible
        // at the same time because borrowing self borrows everything. The workaround is that we
        // swap in a dummy generator into self, so that we can get separate ownership of the generator
        // and pass it into the producer. In the end we swap it back in.
        let mut dummy_color_stream_generator: Box<dyn RewindableSeqStreamGenerator + Sync + Send> = Box::new(io::EmptyRewindableSeqStreamGenerator{});
        std::mem::swap(&mut self.color_stream_generator, &mut dummy_color_stream_generator);
        let mut color_stream_generator = dummy_color_stream_generator; // Now this function owns this

        std::thread::scope(|scope| {
            // Channel of pairs (color id, seq stream)
            let producer_handle = scope.spawn(|| {
                let mut color = 0_usize;
                while let Some((color_stream, _stream_idx)) = color_stream_generator.next() {
                    sender.send((color, color_stream)).unwrap();
                    color += 1;
                }
                drop(sender); // Finished
            });

            let consumer_handles: Vec<_> = (0..n_threads).map(|_| {
                scope.spawn(|| {
                    while let Ok((color, color_stream)) = receiver_ref.recv() {
                        log::info!("Processing color {}", color);
                        let gen = DistinctColexComputation::new(self.sbwt, self.lcs, color_stream);

                        let distinct_colex_positions = gen.run(self.include_rev_comp);
                        for colex in distinct_colex_positions.into_iter() {
                            let set_id = if let Some(filter) = &self.filter {
                                if !filter.get(colex) {
                                    continue; // Do not report this
                                } else {
                                    // Assign new id
                                    filter.rank(colex)
                                }
                            } else {
                                colex // No filter
                            };

                            callback(SetElement{
                                set_id,
                                color,
                            });
                        }
                    }
                })
            }).collect();

            // Wait for threads to finish
            producer_handle.join().unwrap();
            for h in consumer_handles { h.join().unwrap(); }
        });

        self.color_stream_generator = color_stream_generator; // Put this back in (see comment at the start of the function)
    }

    fn set_filter(&mut self, filter: Arc<simple_sds_sbwt::bit_vector::BitVector>) {
        self.filter = Some(filter.clone())
    }

    fn rewind(&mut self) {
        self.color_stream_generator.rewind();
    }
}

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

        let mut thread_inputs = Vec::<IntersectionThreadInput>::with_capacity(n_threads);
        let (mut n_bits_s1, mut n_bits_s2, mut n_result_kmers) = (0_usize, 0_usize, 0_usize);
        for range in crate::util::segment_range(0..n, n_threads) {
            thread_inputs.push(IntersectionThreadInput {
                interleaving_range: range.clone(),
                s1_start_rank: n_bits_s1,
                s2_start_rank: n_bits_s2,
                result_kmer_start_rank: n_result_kmers,
            });
            n_bits_s1 += interleaving.s1[range.clone()].count_ones();
            n_bits_s2 += interleaving.s2[range.clone()].count_ones();
            n_result_kmers += range.filter(|&pos| is_result_kmer(pos)).count();
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
            let mut result_kmer_rank = input.result_kmer_start_rank;
            // The shared colors of coloring1 at the current k-mer, as a bitmap for O(1) lookups
            // and as a list for clearing the bitmap afterwards. Not allocated without shared colors.
            let mut in_set1 = bitvec::bitvec![0; if has_shared_colors { self.n_result_colors } else { 0 }];
            let mut in_set1_list = Vec::<usize>::new();
            for pos in input.interleaving_range {
                if pos > 0 && pos % 10000 == 0 {
                    bar.inc(10000);
                }
                if is_result_kmer(pos) {
                    let result_colex = self.result_kmers.select(result_kmer_rank).unwrap();
                    result_kmer_rank += 1;
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
    use sbwt::{BitPackedKmerSortingMem, SeqStream, StreamingIndex, reverse_complement_in_place};
    use crate::{io::RewindableSeqStreamGenerator, set_of_sets_construction::{ParallelElementGenerator, SetElement}, util::VecVecSeqStream};
    use super::MsElementGenerator;

    struct VecColorStream {
        colors: Vec<Vec<Vec<u8>>>,
        color_idx: usize,
    }

    impl VecColorStream {
        fn new(colors: Vec<Vec<Vec<u8>>>) -> Self {
            Self { colors, color_idx: 0 }
        }
    }

    impl RewindableSeqStreamGenerator for VecColorStream {
        fn next(&mut self) -> Option<(Box<dyn SeqStream + Send + Sync>, usize)> {
            if self.color_idx == self.colors.len() {
                return None;
            }
            let seqs = self.colors[self.color_idx].clone();
            let color_idx = self.color_idx;
            self.color_idx += 1;
            Some((Box::new(VecVecSeqStream::new(seqs)), color_idx))
        }
        fn rewind(&mut self) {
            self.color_idx = 0;
        }
    }

    // Compute expected SetElements by running matching statistics directly on forward
    // and reverse-complement sequences, mirroring the WorkBatch logic.
    fn expected_elements(
        si: &StreamingIndex<'_, sbwt::SbwtIndex<sbwt::SubsetMatrix>, sbwt::LcsArray>,
        k: usize,
        color_seqs: &[Vec<Vec<u8>>],
    ) -> Vec<SetElement> {
        let mut out = Vec::new();
        for (color, seqs) in color_seqs.iter().enumerate() {
            for seq in seqs {
                let mut buf = seq.clone();
                // forward
                for (len, range) in si.matching_statistics_iter(&buf).skip(k - 1) {
                    if len == k {
                        assert_eq!(range.len(), 1);
                        out.push(SetElement { set_id: range.start, color });
                    }
                }
                // reverse complement
                reverse_complement_in_place(&mut buf);
                for (len, range) in si.matching_statistics_iter(&buf).skip(k - 1) {
                    if len == k {
                        assert_eq!(range.len(), 1);
                        out.push(SetElement { set_id: range.start, color });
                    }
                }
            }
        }
        out.sort();
        out
    }

    #[test]
    fn ms_element_generator_emits_correct_elements_including_duplicates() {
        let k = 3_usize;

        // seq0: "ACGCG"
        //   forward k-mers:  ACG, CGC, GCG
        //   rev-comp = CGCGT, k-mers: CGC, GCG, CGT
        //   → CGC and GCG each appear twice for color 0 (non-deduplication is observable)
        let seq0: &[u8] = b"ACGCG";
        let seq1: &[u8] = b"TTTGGG";

        let (sbwt, lcs) = BitPackedKmerSortingMem::new_from_slices(&[seq0, seq1], k)
            .add_rev_comp(true)
            .build_lcs(true)
            .run();
        let lcs = lcs.unwrap();

        let color_seqs: Vec<Vec<Vec<u8>>> = vec![
            vec![seq0.to_vec()], // color 0
            vec![seq1.to_vec()], // color 1
        ];

        let si_ref = StreamingIndex::new(&sbwt, &lcs);
        let expected = expected_elements(&si_ref, k, &color_seqs);

        // Run MsElementGenerator with include_rev_comp = true.
        let gen: Box<dyn RewindableSeqStreamGenerator + Sync + Send> =
            Box::new(VecColorStream::new(color_seqs));
        let si = StreamingIndex::new(&sbwt, &lcs);
        let mut ms_gen = MsElementGenerator::new(gen, si, true, 2);
        let got_mutex = std::sync::Mutex::new(Vec::<SetElement>::new());
        ms_gen.run(|e| got_mutex.lock().unwrap().push(e), 1);
        let mut got = got_mutex.into_inner().unwrap();
        got.sort();

        assert_eq!(got, expected);

        // CGC appears in both the forward and the reverse complement of seq0,
        // so it must be reported twice for color 0 (the non-deduplication property).
        #[allow(clippy::iter_skip_next)] // I don't want to touch it
        let colex_cgc = si_ref
            .matching_statistics_iter(b"CGC")
            .skip(k - 1)
            .next()
            .map(|(_, r)| r.start)
            .expect("CGC should be in the SBWT");
        let cgc_color0 = got.iter().filter(|e| e.set_id == colex_cgc && e.color == 0).count();
        assert_eq!(cgc_color0, 2, "CGC should be reported twice for color 0 (once per strand)");
    }
}
