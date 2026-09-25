use std::{cmp::max, sync::Arc};

use sbwt::{dbg::{Dbg, Node}, LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};
use simple_sds_sbwt::ops::{BitVec, Rank};

use crate::{atomic_bitmap::AtomicBitmap, colex_colored_kmers::{ColexToColorSetMap, CompactColexKmers}, coloring_interface::ColorSetStorage, parallel_ms_iteration::ElementGeneratorFromMergeInterleaving, set_of_sets_construction::{build_color_set_storage, find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates}, unitig_export::break_to_colored_subunitigs};

fn mark_kmer(colex: usize, marks: &AtomicBitmap) {
    marks.set(colex, true);
}

fn search_and_mark_kmer(kmer: &[u8], sbwt: &SbwtIndex<SubsetMatrix>, marks: &AtomicBitmap) {
    let colex = sbwt.search(kmer);
    let colex = colex.unwrap_or_else(|| panic!("k-mer not found in merged SBWT: {:?}", String::from_utf8_lossy(kmer)));
    assert!(colex.len() == 1);
    let colex = colex.start;

    mark_kmer(colex, marks);
}

fn mark_in_neighbors<'a>(colex: usize, dbg: &Dbg<'a, SubsetMatrix>, marks: &AtomicBitmap) {
    let mut in_neighbor_buf = Vec::<(Node, u8)>::new(); // TODO: avoid this allocation
    dbg.push_in_neighbors(Node{id: colex}, &mut in_neighbor_buf);
    for (in_node, _) in in_neighbor_buf.iter() {
        marks.set(in_node.id, true);
    }
}

fn search_and_mark_in_neighbors<'a>(kmer: &[u8], sbwt: &SbwtIndex<SubsetMatrix>, dbg: &Dbg<'a, SubsetMatrix>, marks: &AtomicBitmap) {
    let colex = sbwt.search(kmer);
    let colex = colex.unwrap_or_else(|| panic!("k-mer not found in merged SBWT: {:?}", String::from_utf8_lossy(kmer)));
    assert!(colex.len() == 1);
    let colex = colex.start;

    mark_in_neighbors(colex, dbg, marks);
}

fn mark_key_kmers_for<'a, CSS: ColorSetStorage + Send + Sync>(coloring: &'a CompactColexKmers<CSS>, merged_sbwt: &SbwtIndex<SubsetMatrix>, merged_lcs: &LcsArray, merged_dbg: &Dbg<'_, SubsetMatrix>, key_kmer_marks: &AtomicBitmap, visited_marks: &AtomicBitmap, n_threads: usize) -> Dbg<'a, SubsetMatrix> {

    let merged_si = StreamingIndex::new(merged_sbwt, merged_lcs);

    let k = merged_sbwt.k();
    assert_eq!(k, coloring.get_k());

    log::info!("Initializing DBG");
    let dbg = Dbg::new(coloring.sbwt(), Some(coloring.lcs()), n_threads);

    log::info!("Iterating unitigs");
    let bar = indicatif::ProgressBar::new(coloring.sbwt().n_kmers() as u64);
    dbg.iter_unitigs_with_callback(|nodes|{
        let mut unitig = Vec::<u8>::with_capacity(nodes.len());
        dbg.push_unitig_string(nodes, &mut unitig);
        assert!(unitig.len() >= k);

        let unitig_colex_ranks = nodes.iter().map(|v| v.id).collect::<Vec<usize>>(); // TODO: avoid this allocation
        let (_, subunitig_ranges) = break_to_colored_subunitigs(&unitig_colex_ranks, &unitig, coloring.get_map(), coloring.sbwt());

        // Mark last k-mer of each colored subunitig, and the in-neighbors of
        // the first k-mer of each colored subunitig.
        for subunitig_range in subunitig_ranges {
            // (s,e) = (start of first k-mer, start of the k-mer after the last k-mer)
            let (s,e) = (subunitig_range.start, subunitig_range.end); 
            assert!(s < e);

            // TODO: If we had access to the merge plan here we would not have to search
            // these kmers.

            let last_kmer = &unitig[e-1..e-1+k];
            search_and_mark_kmer(last_kmer, merged_sbwt, key_kmer_marks);

            let first_kmer = &unitig[s..s+k];
            search_and_mark_in_neighbors(first_kmer, merged_sbwt, merged_dbg, key_kmer_marks);
        }

        // Mark last k-mer of every run of k-mers that were visited before, and
        // the in-neighbors of the first k-mer of every run of k-mers that were
        // visited before.
        let mut prev_was_visited = false;
        for (kmer_start, (match_len, colex_range)) in merged_si.matching_statistics_iter(&unitig).skip(k-1).enumerate() {
            assert!(match_len == k);
            assert!(colex_range.len() == 1);

            let kmer_colex = colex_range.start;
            let visited = visited_marks.get(kmer_colex);
            // Visited is true iff the k-mer was visited while processing some earlier coloring.

            if visited & !prev_was_visited {
                // Start of a new colored subunitig
                // -> Mark all in-neighbors for sampling
                // TODO: here we don't have to search again since we get the colex ranks for the MS iterator
                search_and_mark_in_neighbors(&unitig[kmer_start..kmer_start+k], merged_sbwt, merged_dbg, key_kmer_marks);
            } else if !visited && prev_was_visited {
                // One past the end of a colored subunitig
                // -> mark previous node for sampling
                // TODO: here we don't have to search again since we get the colex ranks for the MS iterator
                assert!(kmer_start > 0);
                search_and_mark_kmer(&unitig[kmer_start-1..kmer_start-1+k], merged_sbwt, key_kmer_marks);
            }
            prev_was_visited = visited;

            // Mark this k-mer as visited. Here we are modifying the same bitmap that we are accessing,
            // but that is alright because we are iterating unitigs and hence this k-mer will not be
            // encountered a second time while processing this coloring. So in the future if we find
            // a visited-bit that is set to 1, then it must have been set during the processing of
            // some previous coloring.
            visited_marks.set(kmer_colex, true);
        }
        bar.inc(nodes.len() as u64);
    }, n_threads);
    bar.finish();

    dbg
}

fn mark_new_key_kmers<'a, 'b, CSS: ColorSetStorage + Send + Sync>(coloring1: &'a CompactColexKmers<CSS>, coloring2: &'b CompactColexKmers<CSS>, merged_sbwt: &SbwtIndex<SubsetMatrix>, merged_lcs: &LcsArray, merged_dbg: &Dbg<'_, SubsetMatrix>, sample_distance: usize, n_threads: usize) -> (bitvec::vec::BitVec, Dbg<'a, SubsetMatrix>, Dbg<'b, SubsetMatrix>) {
    let k = merged_sbwt.k();
    assert_eq!(k, coloring1.get_k());
    assert_eq!(k, coloring2.get_k());

    let key_kmer_marks = AtomicBitmap::new(merged_sbwt.n_sets());

    // Mark kmers around branches in the DBG. This part is independent of coloring
    let bar = indicatif::ProgressBar::new(merged_sbwt.n_kmers() as u64);
    merged_dbg.iter_unitigs_with_callback(|nodes|{
        mark_in_neighbors(nodes.first().unwrap().id, merged_dbg, &key_kmer_marks);
        mark_kmer(nodes.last().unwrap().id, &key_kmer_marks);

        for v in nodes.iter().rev().step_by(sample_distance) {
            key_kmer_marks.set(v.id,  true);
        }
        bar.inc(nodes.len() as u64);
    }, n_threads);
    bar.finish();
    
    // Mark around starts and ends of colored subunitigs
    let visited_marks = AtomicBitmap::new(merged_sbwt.n_sets());
    let dbg1 = mark_key_kmers_for(coloring1, merged_sbwt, merged_lcs, merged_dbg, &key_kmer_marks, &visited_marks, n_threads);
    let dbg2 = mark_key_kmers_for(coloring2, merged_sbwt, merged_lcs, merged_dbg, &key_kmer_marks, &visited_marks, n_threads);

    assert_eq!(visited_marks.into_bitvec().count_ones(), merged_sbwt.n_kmers());

    (key_kmer_marks.into_bitvec(), dbg1, dbg2)
}

pub fn merge_compact_colex_kmers<CSS: ColorSetStorage + Send + Sync>(coloring1: CompactColexKmers<CSS>, coloring2: CompactColexKmers<CSS>, optimize_peak_ram: bool, sample_distance: usize, n_threads: usize) -> CompactColexKmers<CSS> {

    log::info!("Computing the sbwt merge plan");
    let merge_plan = sbwt::MergeInterleaving::new(coloring1.sbwt(), coloring2.sbwt(), optimize_peak_ram, n_threads);
    assert_eq!(merge_plan.s1.len(), merge_plan.s2.len());

    log::info!("Merging SBWTs");
    let precalc_len = max(coloring1.sbwt().get_lookup_table().prefix_length, coloring2.sbwt().get_lookup_table().prefix_length);
    // Temporarily destructure the colorings into parts in order to be able to put the
    // SBWT into an Arc to pass to sbwt::merge. The function sbwt::merge takes an Arc
    // for good reasons by design (read the comment at sbwt::merge for an explanation).
    let (sbwt1, lcs1, map1, sets1, color_names_1) = coloring1.into_parts();
    let (sbwt2, lcs2, map2, sets2, color_names_2) = coloring2.into_parts();
    let sbwt1 = Arc::new(sbwt1);
    let sbwt2 = Arc::new(sbwt2);

    let merge_plan = Arc::new(merge_plan);

    // The clones here close just the Arcs.
    let mut merged_sbwt = sbwt::merge(sbwt1.clone(), sbwt2.clone(), merge_plan.clone(), precalc_len, n_threads); 
    merged_sbwt.build_select();

    // Put the coloring structs back together
    let sbwt1 = Arc::try_unwrap(sbwt1).unwrap();
    let sbwt2 = Arc::try_unwrap(sbwt2).unwrap();
    let coloring1 = CompactColexKmers::<CSS>::new(sbwt1, lcs1, map1, sets1, Some(&color_names_1)); 
    let coloring2 = CompactColexKmers::<CSS>::new(sbwt2, lcs2, map2, sets2, Some(&color_names_2)); 

    log::info!("Building the LCS array for the merged SBWT");
    let merged_sbwt_lcs = LcsArray::from_sbwt(&merged_sbwt, n_threads, true);

    log::info!("Initializing DBG for the merged SBWT");
    let merged_dbg = Dbg::new(&merged_sbwt, Some(&merged_sbwt_lcs), n_threads);

    log::info!("=== Phase 1/3: marking new key k-mers ===");
    let (new_key_kmer_marks, _, _) = mark_new_key_kmers(&coloring1, &coloring2, &merged_sbwt, &merged_sbwt_lcs, &merged_dbg, sample_distance, n_threads);
    log::info!("Marked {:.2} % of all k-mers", new_key_kmer_marks.count_ones() as f64 / merged_sbwt.n_kmers() as f64 * 100.0);

    log::info!("=== PHASE 2/3: Building color set finperprints for key k-mers ===");
    let random_seed = 123123; // Todo: be more random
    let gen = ElementGeneratorFromMergeInterleaving {
        interleaving: &merge_plan,
        coloring1: &coloring1,
        coloring2: &coloring2,
        merged_key_kmer_marks: &new_key_kmer_marks,
        filter: None,
    } ;

    let n_colors = coloring1.get_set_storage().n_colors() + coloring2.get_set_storage().n_colors();
    let n_colors = u32::try_from(n_colors).unwrap_or_else( |_| {
        log::error!("Maximum number of colors 2^32 exceeded");
        panic!();
    });

    // TODO: do not clone new_key_kmer_marks
    let (repr_kmer_marks, distinct_set_sizes, key_kmer_idx_to_set_id, _new_key_kmer_marks) = find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(gen, new_key_kmer_marks.clone(), n_colors, n_threads, random_seed);

    let gen = ElementGeneratorFromMergeInterleaving {
        interleaving: &merge_plan,
        coloring1: &coloring1,
        coloring2: &coloring2,
        merged_key_kmer_marks: &new_key_kmer_marks,
        filter: None,
    } ;
    log::info!("=== PHASE 3/3: Build the distinct color set storage ===");
        
    let css = build_color_set_storage(n_colors as usize, repr_kmer_marks, distinct_set_sizes, gen, n_threads);

    log::info!("Building rank support for key k-mer marks");
    let mut key_kmer_marks = crate::util::bitvec_to_simple_sds_bitvec(new_key_kmer_marks);
    key_kmer_marks.enable_rank();
    assert!(key_kmer_idx_to_set_id.len() == key_kmer_marks.rank(key_kmer_marks.len()));

    let colex_map = ColexToColorSetMap {
        sampling: key_kmer_marks, 
        color_set_ids: key_kmer_idx_to_set_id,
    };

    let mut color_names = Vec::<String>::new();
    color_names.extend(coloring1.get_color_names().iter().cloned());
    color_names.extend(coloring2.get_color_names().iter().cloned());

    CompactColexKmers::<CSS>::new(merged_sbwt, merged_sbwt_lcs, colex_map, css, Some(&color_names))
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, hash::BuildHasherDefault};

    use jseqio::seq_db::SeqDB;
    use rustc_hash::FxHasher;
    use sbwt::{BitPackedKmerSortingMem, LcsArray, SbwtIndex, SubsetMatrix};
    use simple_sds_sbwt::ops::{BitVec, Rank};

    use crate::{colex_colored_kmers::{ColexToColorSetMap, mark_key_kmers}, coloring_interface::{ColorSetStorage, ColorSetView}, int_vec::CompactIntVec, io::RewindableSeqStreamGenerator, sparse_dense_storage::SparseDenseStorage, util::VecVecRewindableGen};

    use super::CompactColexKmers;

    // Annoying plumbing to get a RewindableSeqStreamGenerator from Vec<SeqDB>
    /*
    struct SeqDBsColorStream{
        dbs: Vec<SeqDB>,
        db_idx: usize,
    }
    impl RewindableSeqStreamGenerator for SeqDBsColorStream {
        fn next(&mut self) -> Option<Box<dyn sbwt::SeqStream + Send + Sync>> {
            if self.db_idx == self.dbs.len() { None }
            else {
                let seqs: Vec<Vec<u8>> = self.dbs[self.db_idx].iter().map(|rec| rec.seq.to_owned()).collect();
                let ss = VecVecSeqStream::new(seqs);
                let it: Box<dyn sbwt::SeqStream + Send + Sync> = Box::new(ss);
                Some(it)
            }
        }
    
        fn rewind(&mut self) {
            self.db_idx = 0;
        }
    }
    */

    /// Output:
    /// - Distinct color sets encoded as ColorSetStorage
    /// - HashMap from color set to its index in ColorSets
    pub fn hash_and_encode_distinct_sets<'a, CSS: ColorSetStorage>(colex_to_set: &'a CSS, n_colors: usize) -> (CSS, HashMap::<CSS::SetView<'a>, usize, BuildHasherDefault::<FxHasher>>) {
        let n_sets = colex_to_set.n_sets();

        log::info!("Hashing distinct color sets");

        let mut distinct_sets = HashMap::<CSS::SetView<'a>, usize, BuildHasherDefault::<FxHasher>>::default(); // Set -> id
        let mut distinct_set_colex_ranks = Vec::<usize>::new();
        let bar = indicatif::ProgressBar::new(n_sets as u64);
        for colex in 0..n_sets {
            let key = colex_to_set.get_set_view(colex);
            if !distinct_sets.contains_key(&key) {
                distinct_sets.insert(key, distinct_sets.len());
                distinct_set_colex_ranks.push(colex);
            }
            if colex % 1000 == 0 {
                bar.inc(1000);
            }
        }
        bar.finish();

        log::info!("{} distinct color sets found", distinct_sets.len());

        // Create an iterator of iterators, each inner iterator iterating over one color set
        let color_sets_iterator = distinct_set_colex_ranks.into_iter().map(|colex| {
            colex_to_set.get_set_view(colex).iter()
        });

        let colorsets = CSS::new_from_iter_of_iters(color_sets_iterator, n_colors);

        (*colorsets, distinct_sets)

    }

    #[cfg(test)]
    pub(crate) fn gen_random_dna_string(len: usize, seed: u64) -> Vec<u8> {
        use rand_chacha::rand_core::{RngCore, SeedableRng};

        let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(seed);
        (0..len).map(|_| { 
            match rng.next_u64() % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => panic!("Impossible")
            }
        }).collect()
    }

    fn build_color_sets<CSS: ColorSetStorage>(sbwt1: &SbwtIndex<SubsetMatrix>, lcs1: &LcsArray, dbs1: Vec<SeqDB>, n_threads: usize) 
    -> (Vec<usize>, CSS){
        let n_colors_1 = dbs1.len();
        let bms1 = crate::bitmap_storage::build_from_seq_dbs(dbs1, sbwt1, lcs1, n_threads);

        let iter_of_iters_1 = (0..sbwt1.n_sets()).map(|colex| bms1.get_set_view(colex).iter());
        let colex_to_css_1 = *CSS::new_from_iter_of_iters(iter_of_iters_1, n_colors_1);

        let (distinct_css_1, set_to_id_1) = hash_and_encode_distinct_sets(&colex_to_css_1, n_colors_1);
        let colex_to_id: Vec<usize> = (0..sbwt1.n_sets()).map(|colex| {
            set_to_id_1[&colex_to_css_1.get_set_view(colex)]
        }).collect(); 

        (colex_to_id, distinct_css_1)
    }

    fn seqs_to_dbs(seqs: &[Vec<u8>]) -> Vec<SeqDB> {
        seqs.iter().map(|seq| {
            let mut db = SeqDB::new();
            db.push_seq(seq);
            db
        }).collect()
    }

    fn build_sbwt(k: usize, seqs: &[Vec<u8>]) -> (SbwtIndex<SubsetMatrix>, LcsArray) {
        let (mut sbwt, lcs) = sbwt::SbwtIndexBuilder::new()
            .add_rev_comp(false)
            .k(k)
            .build_lcs(true)
            .n_threads(3)
            .precalc_length(5)
            .algorithm(BitPackedKmerSortingMem::new().dedup_batches(true))
        .run_from_vecs(seqs);
        sbwt.build_select();
        (sbwt, lcs.unwrap())
    }

    /// Builds a CompactColexKmers where color i consists of the single sequence seqs[i].
    fn build_coloring(k: usize, seqs: &[Vec<u8>], sample_distance: usize, n_threads: usize) -> CompactColexKmers<SparseDenseStorage> {
        let (sbwt, lcs) = build_sbwt(k, seqs);
        let (colex_to_id, storage) = build_color_sets::<SparseDenseStorage>(&sbwt, &lcs, seqs_to_dbs(seqs), n_threads);

        let mut gen: Box<dyn RewindableSeqStreamGenerator + Sync + Send> = Box::new(VecVecRewindableGen::new(seqs.to_vec()));
        let key_kmers = mark_key_kmers(&sbwt, &lcs, sample_distance, &mut gen, n_threads, 1, false);
        let sampled_ids: Vec<usize> = colex_to_id.iter().enumerate().filter(|(i, _)| key_kmers[*i]).map(|(_,x)| *x).collect();
        assert!(key_kmers.count_ones() == sampled_ids.len());

        let mut key_kmers = crate::util::bitvec_to_simple_sds_bitvec(key_kmers);
        key_kmers.enable_rank();

        let colex_map = ColexToColorSetMap{
            sampling: key_kmers,
            color_set_ids: CompactIntVec::from_vec(sampled_ids),
        };

        CompactColexKmers::new(sbwt, lcs, colex_map, storage, None)
    }

    /// Merges the colorings of input_seqs_1 and input_seqs_2 (one color per sequence) and checks
    /// that every k-mer gets the same color set as in a coloring built directly from all sequences.
    fn check_merge(k: usize, input_seqs_1: &[Vec<u8>], input_seqs_2: &[Vec<u8>], input_sample_distance: usize, merge_sample_distance: usize, n_threads: usize) {
        let mut all_input_seqs = input_seqs_1.to_vec();
        all_input_seqs.extend(input_seqs_2.iter().cloned());

        let ccc1 = build_coloring(k, input_seqs_1, input_sample_distance, n_threads);
        let ccc2 = build_coloring(k, input_seqs_2, input_sample_distance, n_threads);
        let ccc_both = build_coloring(k, &all_input_seqs, input_sample_distance, n_threads);

        let ccc_merged = super::merge_compact_colex_kmers(ccc1, ccc2, true, merge_sample_distance, n_threads);
        let sbwt_merged = &ccc_merged.sbwt();

        for colex in 0..ccc_both.sbwt().n_sets() {
            let kmer = ccc_both.sbwt().access_kmer(colex);

            if kmer.iter().all(|c| *c != b'$') { // Not a dummy k-mer
                let true_colors: Vec<usize> = ccc_both.colex_to_set(colex).iter().collect();
                let range = sbwt_merged.search(&kmer).unwrap();
                assert_eq!(range.len(), 1);
                let colex_merged = range.start;
                let merged_colors: Vec<usize> = ccc_merged.colex_to_set(colex_merged).iter().collect();

                assert_eq!(true_colors, merged_colors, "k = {}, k-mer {}", k, String::from_utf8_lossy(&kmer));
            }
        }
    }

    #[test]
    fn test_merge() {

        // Opt-in logging: set RUST_LOG=info to see output; silent by default.
        let _ = env_logger::try_init();

        let n_threads = 3;

        for k in 3_usize..10_usize { // k < 3 does not work because construction uses 3-mer binning.

            let input_seqs_1: Vec<Vec<u8>> = (0..10).map(|i| gen_random_dna_string(20, (i + k.pow(4)) as u64)).collect();
            let input_seqs_2: Vec<Vec<u8>> = (0..10).map(|i| gen_random_dna_string(20, (123456 + i + k.pow(4)) as u64)).collect();

            check_merge(k, &input_seqs_1, &input_seqs_2, 3, 5, n_threads);
        }
    }

    // The sample distances in the tests below are larger than any unitig, so that
    // only the structurally required key k-mers get marked. Any k-mer that the
    // merge forgets to mark then results in a wrong color set.
    const NO_SAMPLING: usize = 1000;

    #[test]
    fn test_merge_multiple_colored_subunitigs_in_one_unitig() {
        let _ = env_logger::try_init();
        let k = 11;
        let genome = gen_random_dna_string(80, 1);

        // Color 1 covers the middle of the unitig of color 0, so the single unitig
        // of the first coloring breaks into three colored subunitigs.
        let input_seqs_1 = vec![genome.clone(), genome[20..50].to_vec()];
        let input_seqs_2 = vec![gen_random_dna_string(80, 2)];

        check_merge(k, &input_seqs_1, &input_seqs_2, NO_SAMPLING, NO_SAMPLING, 3);
        check_merge(k, &input_seqs_2, &input_seqs_1, NO_SAMPLING, NO_SAMPLING, 3);
    }

    #[test]
    fn test_merge_shared_kmers_in_middle_of_unitig() {
        let _ = env_logger::try_init();
        let k = 11;
        let genome = gen_random_dna_string(80, 3);

        // The k-mers of the first coloring form a run in the middle of the unitig
        // of the second coloring, so processing the second coloring enters and
        // leaves a run of previously visited k-mers.
        let inner = vec![genome[20..50].to_vec()];
        let outer = vec![genome.clone()];
        check_merge(k, &inner, &outer, NO_SAMPLING, NO_SAMPLING, 3);
        check_merge(k, &outer, &inner, NO_SAMPLING, NO_SAMPLING, 3);

        // Partial overlap: a shared run at the start of one unitig and at the end of the other.
        let left = vec![genome[0..50].to_vec()];
        let right = vec![genome[30..80].to_vec()];
        check_merge(k, &left, &right, NO_SAMPLING, NO_SAMPLING, 3);
        check_merge(k, &right, &left, NO_SAMPLING, NO_SAMPLING, 3);
    }

    #[test]
    fn test_merge_overlapping_random_substrings() {
        use rand_chacha::rand_core::{RngCore, SeedableRng};

        let _ = env_logger::try_init();
        let n_threads = 3;

        // Both colorings consist of random substrings of the same genome, so the
        // colorings share many k-mers and unitigs of the colorings partially overlap
        // in all kinds of ways.
        for k in [5_usize, 8, 11, 15] {
            for (sample_distance, seed) in [(NO_SAMPLING, 0_u64), (3, 1)] {
                let genome = gen_random_dna_string(300, 1000 + k as u64 + seed);
                let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(2000 + k as u64 * 100 + seed);
                let mut random_substrings = |n: usize| -> Vec<Vec<u8>> {
                    (0..n).map(|_| {
                        let start = rng.next_u64() as usize % (genome.len() - k);
                        let len = k + rng.next_u64() as usize % 80;
                        genome[start..(start + len).min(genome.len())].to_vec()
                    }).collect()
                };
                let input_seqs_1 = random_substrings(6);
                let input_seqs_2 = random_substrings(6);
                check_merge(k, &input_seqs_1, &input_seqs_2, sample_distance, sample_distance, n_threads);
            }
        }
    }
}
