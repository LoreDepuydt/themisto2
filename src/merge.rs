use std::{cmp::max, collections::{HashMap, HashSet}, sync::Arc};

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

fn mark_key_kmers_for<'a, CSS: ColorSetStorage + Send + Sync>(coloring: &'a CompactColexKmers<CSS>, merged_sbwt: &SbwtIndex<SubsetMatrix>, merged_lcs: &LcsArray, merged_dbg: &Dbg<'_, SubsetMatrix>, key_kmer_marks: &AtomicBitmap, visited_marks: Option<&AtomicBitmap>, n_threads: usize) -> Dbg<'a, SubsetMatrix> {

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

        // Debug-only coverage check: record every k-mer of this unitig as visited, so
        // mark_new_key_kmers can assert afterwards that every merged k-mer was visited while
        // processing coloring1 or coloring2. This costs a full streaming-index pass over the
        // unitig, so it is skipped in release builds. It is otherwise unneeded for marking:
        // every k-mer that would be marked here is already marked by the loop above or by
        // mark_new_key_kmers (proof: only the second coloring can have visited k-mers, and
        // visited = contained in the first coloring. Let x be the first k-mer of a visited run.
        // - If x is the first k-mer of this unitig, its in-neighbors were marked by the loop above.
        // - Otherwise its predecessor w in this unitig is not in the first coloring. If x has no
        //   in-neighbors in the first coloring, x starts a unitig of the first coloring, so the
        //   loop above marked its in-neighbors when processing the first coloring. Otherwise x
        //   has an in-neighbor v != w there, so x has in-degree >= 2 in the merged DBG and
        //   mark_new_key_kmers marked its in-neighbors.
        // Let x be the last k-mer of a visited run, followed by y in this unitig. y is not in the
        // first coloring, so either x has out-degree 0 in the first coloring and ends a unitig
        // there (marked by the loop above), or x has out-degree >= 2 in the merged DBG and ends
        // a merged unitig (marked by mark_new_key_kmers). Checked with an assertion on ~3700
        // random adversarial merges.)
        if let Some(visited_marks) = visited_marks {
            for (match_len, colex_range) in merged_si.matching_statistics_iter(&unitig).skip(k-1) {
                assert!(match_len == k);
                assert!(colex_range.len() == 1);

                // Modifying the same bitmap we are accessing is alright because we are iterating
                // unitigs and hence this k-mer will not be encountered a second time while
                // processing this coloring. So in the future if we find a visited-bit that is set
                // to 1, then it must have been set during the processing of some previous coloring.
                visited_marks.set(colex_range.start, true);
            }
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
    
    // Debug-only sanity check that every merged k-mer gets visited while processing coloring1
    // or coloring2 (see the coverage check in mark_key_kmers_for). Skipped in release builds:
    // it costs a full extra streaming-index pass over the merged graph and is not needed for
    // marking key k-mers, only for catching merge bugs during development.
    let visited_marks = cfg!(debug_assertions).then(|| AtomicBitmap::new(merged_sbwt.n_sets()));
    let dbg1 = mark_key_kmers_for(coloring1, merged_sbwt, merged_lcs, merged_dbg, &key_kmer_marks, visited_marks.as_ref(), n_threads);
    let dbg2 = mark_key_kmers_for(coloring2, merged_sbwt, merged_lcs, merged_dbg, &key_kmer_marks, visited_marks.as_ref(), n_threads);

    if let Some(visited_marks) = visited_marks {
        assert_eq!(visited_marks.into_bitvec().count_ones(), merged_sbwt.n_kmers());
    }

    (key_kmer_marks.into_bitvec(), dbg1, dbg2)
}

fn check_unique_names(names: &[String], which: &str) -> Result<(), String> {
    let mut seen = HashSet::<&str>::with_capacity(names.len());
    for name in names {
        if !seen.insert(name.as_str()) {
            return Err(format!("Color name {:?} appears more than once in the {} index. Color names must be unique to merge shared colors.", name, which));
        }
    }
    Ok(())
}

/// Returns the merged color id of each color of the second coloring, and the color names of
/// the merged coloring. The colors of the first coloring keep their ids. If merge_shared is
/// true, a color of the second coloring with the same name as a color of the first coloring is
/// mapped to that color. All other colors of the second coloring get new ids after the colors
/// of the first coloring, in their original order. Returns an error if merge_shared is true and
/// a name appears twice within one coloring.
fn merged_color_mapping(names1: &[String], names2: &[String], merge_shared: bool) -> Result<(Vec<usize>, Vec<String>), String> {
    let n1 = names1.len();
    let mut merged_names = names1.to_vec();

    if !merge_shared {
        merged_names.extend(names2.iter().cloned());
        return Ok(((n1..n1 + names2.len()).collect(), merged_names));
    }

    check_unique_names(names1, "first")?;
    check_unique_names(names2, "second")?;

    let name_to_id1: HashMap<&str, usize> = names1.iter().enumerate().map(|(i, name)| (name.as_str(), i)).collect();
    let color2_to_merged: Vec<usize> = names2.iter().map(|name| {
        name_to_id1.get(name.as_str()).copied().unwrap_or_else(|| {
            merged_names.push(name.clone());
            merged_names.len() - 1
        })
    }).collect();

    log::info!("{} colors are shared between the two indexes", n1 + names2.len() - merged_names.len());
    Ok((color2_to_merged, merged_names))
}

pub fn merge_compact_colex_kmers<CSS: ColorSetStorage + Send + Sync>(coloring1: CompactColexKmers<CSS>, coloring2: CompactColexKmers<CSS>, merge_shared_colors: bool, optimize_peak_ram: bool, sample_distance: usize, n_threads: usize) -> CompactColexKmers<CSS> {

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
    let (color2_to_merged, merged_color_names) = merged_color_mapping(&color_names_1, &color_names_2, merge_shared_colors).unwrap_or_else(|e| {
        log::error!("{}", e);
        panic!("{}", e);
    });
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
        color2_to_merged: &color2_to_merged,
    } ;

    let n_colors = u32::try_from(merged_color_names.len()).unwrap_or_else( |_| {
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
        color2_to_merged: &color2_to_merged,
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

    CompactColexKmers::<CSS>::new(merged_sbwt, merged_sbwt_lcs, colex_map, css, Some(&merged_color_names))
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

    // =====================================================================
    // Test helpers
    // =====================================================================

    // The sample distances in the tests below are larger than any unitig, so that
    // only the structurally required key k-mers get marked. Any k-mer that the
    // merge forgets to mark then results in a wrong color set.
    const NO_SAMPLING: usize = 1000;

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

    /// Turns (color name, sequence) pairs into (color id, sequence) pairs, where the color id
    /// is the index of the name in `names`. Names that are not yet in `names` are appended, so
    /// a name that appears several times is one color with several sequences.
    fn assign_color_ids(names: &mut Vec<String>, colors: &[(&str, Vec<u8>)]) -> Vec<(usize, Vec<u8>)> {
        colors.iter().map(|(name, seq)| {
            let color = names.iter().position(|x| x == name).unwrap_or_else(|| {
                names.push(name.to_string());
                names.len() - 1
            });
            (color, seq.clone())
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

    /// Builds a CompactColexKmers where color i is named names[i] and consists of the sequences
    /// s with (i, s) in seqs.
    fn build_named_coloring(k: usize, names: &[String], seqs: &[(usize, Vec<u8>)], sample_distance: usize, n_threads: usize) -> CompactColexKmers<SparseDenseStorage> {
        let mut dbs: Vec<SeqDB> = names.iter().map(|_| SeqDB::new()).collect();
        for (color, seq) in seqs {
            dbs[*color].push_seq(seq);
        }
        let all_seqs: Vec<Vec<u8>> = seqs.iter().map(|(_, seq)| seq.clone()).collect();
        let (sbwt, lcs) = build_sbwt(k, &all_seqs);
        let (colex_to_id, storage) = build_color_sets::<SparseDenseStorage>(&sbwt, &lcs, dbs, n_threads);

        let mut gen: Box<dyn RewindableSeqStreamGenerator + Sync + Send> = Box::new(VecVecRewindableGen::new(all_seqs));
        let key_kmers = mark_key_kmers(&sbwt, &lcs, sample_distance, &mut gen, n_threads, 1, false);
        let sampled_ids: Vec<usize> = colex_to_id.iter().enumerate().filter(|(i, _)| key_kmers[*i]).map(|(_,x)| *x).collect();
        assert!(key_kmers.count_ones() == sampled_ids.len());

        let mut key_kmers = crate::util::bitvec_to_simple_sds_bitvec(key_kmers);
        key_kmers.enable_rank();

        let colex_map = ColexToColorSetMap{
            sampling: key_kmers,
            color_set_ids: CompactIntVec::from_vec(sampled_ids),
        };

        CompactColexKmers::new(sbwt, lcs, colex_map, storage, Some(names))
    }

    /// Merges the colorings of input_seqs_1 and input_seqs_2 (one color per sequence, all colors
    /// distinct) and checks the result against a coloring built directly from all sequences.
    fn check_merge(k: usize, input_seqs_1: &[Vec<u8>], input_seqs_2: &[Vec<u8>], input_sample_distance: usize, merge_sample_distance: usize, n_threads: usize) {
        let names: Vec<String> = (0..input_seqs_1.len() + input_seqs_2.len()).map(|i| format!("s{}", i)).collect();
        let (names1, names2) = names.split_at(input_seqs_1.len());
        let colors1: Vec<(&str, Vec<u8>)> = names1.iter().zip(input_seqs_1).map(|(name, seq)| (name.as_str(), seq.clone())).collect();
        let colors2: Vec<(&str, Vec<u8>)> = names2.iter().zip(input_seqs_2).map(|(name, seq)| (name.as_str(), seq.clone())).collect();

        // With distinct names, merging shared colors must not change anything
        for merge_shared_colors in [false, true] {
            check_named_merge(k, &colors1, &colors2, merge_shared_colors, input_sample_distance, merge_sample_distance, n_threads);
        }
    }

    /// Merges two colorings given as (color name, sequence) pairs (see assign_color_ids) and checks
    /// that every k-mer gets the same color set as in a coloring built directly from the expected
    /// merged colors. Colors with the same name are shared if merge_shared_colors is true.
    fn check_named_merge(k: usize, colors1: &[(&str, Vec<u8>)], colors2: &[(&str, Vec<u8>)], merge_shared_colors: bool, input_sample_distance: usize, merge_sample_distance: usize, n_threads: usize) {
        let mut names1 = Vec::<String>::new();
        let mut names2 = Vec::<String>::new();
        let input_seqs_1 = assign_color_ids(&mut names1, colors1);
        let input_seqs_2 = assign_color_ids(&mut names2, colors2);

        // Expected merged colors: the colors of the first coloring keep their ids, and the new
        // colors of the second coloring come after them.
        let mut expected_names = names1.clone();
        let expected_seqs_2 = if merge_shared_colors {
            // A color of the second coloring with a name of the first coloring gets that id
            assign_color_ids(&mut expected_names, colors2)
        } else {
            // Every color of the second coloring is new, even if its name is not
            expected_names.extend(names2.iter().cloned());
            input_seqs_2.iter().map(|(color, seq)| (names1.len() + color, seq.clone())).collect()
        };
        let expected_seqs = [input_seqs_1.clone(), expected_seqs_2].concat();

        let ccc1 = build_named_coloring(k, &names1, &input_seqs_1, input_sample_distance, n_threads);
        let ccc2 = build_named_coloring(k, &names2, &input_seqs_2, input_sample_distance, n_threads);
        let ccc_both = build_named_coloring(k, &expected_names, &expected_seqs, input_sample_distance, n_threads);

        let ccc_merged = super::merge_compact_colex_kmers(ccc1, ccc2, merge_shared_colors, true, merge_sample_distance, n_threads);
        assert_eq!(ccc_merged.get_color_names(), &expected_names);
        assert_eq!(ccc_merged.get_set_storage().n_colors(), expected_names.len());
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

    // =====================================================================
    // Tests
    // =====================================================================

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

    #[test]
    fn test_merge_shared_colors() {
        let _ = env_logger::try_init();
        let n_threads = 3;

        for k in [5_usize, 8, 11] {
            for sample_distance in [NO_SAMPLING, 3] {
                let genome = gen_random_dna_string(200, 5000 + k as u64);
                let seed = 6000 + 10 * k as u64;

                // A and D are private to one index. B has unrelated sequences in the two indexes.
                // C has overlapping substrings of the same genome in both indexes, so some k-mers
                // have C on both sides, and it must be reported only once.
                let colors1 = vec![
                    ("A", gen_random_dna_string(60, seed)),
                    ("B", gen_random_dna_string(60, seed + 1)),
                    ("B", genome[150..200].to_vec()),
                    ("C", genome[0..100].to_vec()),
                ];
                let colors2 = vec![
                    ("B", gen_random_dna_string(60, seed + 2)),
                    ("D", genome[40..160].to_vec()),
                    ("C", genome[50..130].to_vec()),
                    ("C", genome[140..190].to_vec()),
                ];

                check_named_merge(k, &colors1, &colors2, true, sample_distance, sample_distance, n_threads);
                check_named_merge(k, &colors2, &colors1, true, sample_distance, sample_distance, n_threads);

                // Without the flag, equally named colors stay distinct
                check_named_merge(k, &colors1, &colors2, false, sample_distance, sample_distance, n_threads);
            }
        }
    }

    #[test]
    fn test_merge_all_colors_shared() {
        let _ = env_logger::try_init();
        let k = 9;
        let genome = gen_random_dna_string(200, 77);
        let colors1 = vec![("X", genome[0..120].to_vec()), ("Y", genome[60..200].to_vec())];
        let colors2 = vec![("Y", genome[0..90].to_vec()), ("X", genome[100..200].to_vec())];
        check_named_merge(k, &colors1, &colors2, true, NO_SAMPLING, NO_SAMPLING, 3);
        check_named_merge(k, &colors1, &colors2, true, 3, 3, 3);
    }

    #[test]
    fn test_merge_shared_colors_rejects_duplicate_names() {
        let a = vec!["A".to_string(), "A".to_string()];
        let b = vec!["B".to_string()];
        assert!(super::merged_color_mapping(&a, &b, true).is_err());
        assert!(super::merged_color_mapping(&b, &a, true).is_err());
        assert!(super::merged_color_mapping(&a, &b, false).is_ok()); // Names are not matched without the flag
    }

    #[test]
    fn test_merged_color_mapping() {
        let to_strings = |v: &[&str]| -> Vec<String> { v.iter().map(|s| s.to_string()).collect() };
        let names1 = to_strings(&["A", "B", "C"]);
        let names2 = to_strings(&["D", "C", "E", "A"]);

        let (map, names) = super::merged_color_mapping(&names1, &names2, true).unwrap();
        assert_eq!(map, vec![3, 2, 4, 0]);
        assert_eq!(names, to_strings(&["A", "B", "C", "D", "E"]));

        let (map, names) = super::merged_color_mapping(&names1, &names2, false).unwrap();
        assert_eq!(map, vec![3, 4, 5, 6]);
        assert_eq!(names, to_strings(&["A", "B", "C", "D", "C", "E", "A"]));
    }
}
