use std::{cmp::max, sync::Arc};

use sbwt::{dbg::Dbg, LcsArray};
use simple_sds_sbwt::ops::{BitVec, Rank, Select};

use crate::{atomic_bitmap::AtomicBitmap, colex_colored_kmers::{ColexToColorSetMap, CompactColexKmers}, coloring_interface::ColorSetStorage, set_operations::result_positions::ToIntersection, set_operations::key_kmers::{mark_key_kmers_for, mark_structural_key_kmers}, parallel_ms_iteration::{ColorCombination, ElementGeneratorFromIntersectionInterleaving}, set_of_sets_construction::{build_color_set_storage, find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates}};

/// How the color set of a k-mer of the intersection is computed from its color sets in the two
/// input indexes.
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntersectColors {
    /// The union of the two color sets, as when merging the indexes
    Union,
    /// The intersection of the two color sets, matching colors by name. The result has only the
    /// colors whose names are in both indexes.
    Intersect,
}

/// Intersects two colored indexes: the result has the k-mers that are in both indexes. The color
/// set of each k-mer is the union or the intersection of its color sets in the two indexes,
/// depending on `colors`. With the union, merge_shared_colors is as in
/// merge::merge_compact_colex_kmers. With the intersection, colors are always matched by name, and
/// a k-mer whose color sets do not intersect has an empty color set.
///
/// Both SBWTs need select support.
#[allow(clippy::too_many_arguments)]
pub fn intersect_compact_colex_kmers<CSS: ColorSetStorage + Send + Sync>(coloring1: CompactColexKmers<CSS>, coloring2: CompactColexKmers<CSS>, colors: IntersectColors, merge_shared_colors: bool, optimize_peak_ram: bool, sample_distance: usize, n_threads: usize) -> CompactColexKmers<CSS> {

    log::info!("Computing the sbwt interleaving");
    let interleaving = sbwt::MergeInterleaving::new(coloring1.sbwt(), coloring2.sbwt(), optimize_peak_ram, n_threads);
    assert_eq!(interleaving.s1.len(), interleaving.s2.len());

    // Color mapping first, so that an error is reported before the expensive work
    let (color1_to_result, color2_to_result, result_color_names) = match colors {
        IntersectColors::Union => crate::set_operations::colors::merged_color_mapping(coloring1.get_color_names(), coloring2.get_color_names(), merge_shared_colors)
            .map(|(color2_to_merged, names)| (vec![], color2_to_merged.into_iter().map(Some).collect::<Vec<_>>(), names)),
        IntersectColors::Intersect => crate::set_operations::colors::intersected_color_mapping(coloring1.get_color_names(), coloring2.get_color_names()),
    }.unwrap_or_else(|e| {
        log::error!("{}", e);
        panic!("{}", e);
    });

    log::info!("Intersecting SBWTs");
    let precalc_len = max(coloring1.sbwt().get_lookup_table().prefix_length, coloring2.sbwt().get_lookup_table().prefix_length);
    // Temporarily destructure the colorings to put the SBWTs into Arcs for sbwt::intersect, as in
    // merge::merge_compact_colex_kmers.
    let (sbwt1, lcs1, map1, sets1, color_names_1) = coloring1.into_parts();
    let (sbwt2, lcs2, map2, sets2, color_names_2) = coloring2.into_parts();
    let sbwt1 = Arc::new(sbwt1);
    let sbwt2 = Arc::new(sbwt2);
    let interleaving = Arc::new(interleaving);

    // TODO: sbwt::intersect repairs the dummy chains of the result by building a three-way
    // interleaving of both inputs and an auxiliary index from scratch, which takes about as long as
    // the interleaving above and is the memory peak of the intersection. Reusing this interleaving
    // and only interleaving the (usually small) auxiliary index into it would save most of that.
    let mut result_sbwt = sbwt::intersect(sbwt1.clone(), sbwt2.clone(), interleaving.clone(), precalc_len, optimize_peak_ram, n_threads)
        .unwrap_or_else(|_| panic!("The SBWT of the first index has no select support"));
    result_sbwt.build_select();

    // Put the coloring structs back together
    let sbwt1 = Arc::try_unwrap(sbwt1).unwrap();
    let sbwt2 = Arc::try_unwrap(sbwt2).unwrap();
    let interleaving = Arc::try_unwrap(interleaving).unwrap();
    let coloring1 = CompactColexKmers::<CSS>::new(sbwt1, lcs1, map1, sets1, Some(&color_names_1));
    let coloring2 = CompactColexKmers::<CSS>::new(sbwt2, lcs2, map2, sets2, Some(&color_names_2));

    log::info!("Building the LCS array for the intersection SBWT");
    let result_lcs = LcsArray::from_sbwt(&result_sbwt, n_threads, true);

    log::info!("Initializing DBG for the intersection SBWT");
    let result_dbg = Dbg::new(&result_sbwt, Some(&result_lcs), n_threads);

    log::info!("Locating the k-mers of the intersection");
    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    let mut result_kmers = crate::util::bitvec_to_simple_sds_bitvec(!thread_pool.install(|| result_sbwt.compute_dummy_node_marks()));
    result_kmers.enable_select();
    assert_eq!(result_kmers.count_ones(), result_sbwt.n_kmers());
    let mut in_result = {
        let words: Vec<usize> = interleaving.s1.as_raw_slice().iter()
            .zip(interleaving.s2.as_raw_slice())
            .zip(interleaving.is_dummy.as_raw_slice())
            .map(|((&a, &b), &d)| (a & b & !d) as usize).collect();
        let mut in_result = bitvec::vec::BitVec::<usize, bitvec::order::Lsb0>::from_vec(words);
        in_result.truncate(interleaving.s1.len());
        crate::util::bitvec_to_simple_sds_bitvec(in_result)
    };
    in_result.enable_rank();
    assert_eq!(in_result.count_ones(), result_sbwt.n_kmers());

    log::info!("=== Phase 1/3: marking key k-mers ===");
    let key_kmer_marks = AtomicBitmap::new(result_sbwt.n_sets());
    mark_structural_key_kmers(&result_sbwt, &result_dbg, &key_kmer_marks, sample_distance, n_threads);
    // A color set changes between consecutive k-mers x and y of a unitig of the intersection only
    // if the color set of x or y changes in one of the inputs, which has the edge from x to y too.
    // Then x is the last k-mer of a colored subunitig of that input, so it is marked here.
    for (coloring, in_input) in [(&coloring1, &interleaving.s1), (&coloring2, &interleaving.s2)] {
        let to_result = ToIntersection::new(in_input, &in_result, &result_kmers);
        mark_key_kmers_for(coloring, |colex| to_result.result_colex(colex), &result_sbwt, &result_lcs, &result_dbg, &key_kmer_marks, None, n_threads);
    }
    drop(in_result);
    let key_kmer_marks = key_kmer_marks.into_bitvec();
    log::info!("Marked {:.2} % of all k-mers", key_kmer_marks.count_ones() as f64 / result_sbwt.n_kmers() as f64 * 100.0);

    let n_colors = u32::try_from(result_color_names.len()).unwrap_or_else( |_| {
        log::error!("Maximum number of colors 2^32 exceeded");
        panic!();
    });
    let color2_to_merged: Vec<usize>; // The union needs the mapping without Options
    let color_combination = match colors {
        IntersectColors::Union => {
            color2_to_merged = color2_to_result.iter().map(|c| c.unwrap()).collect();
            ColorCombination::Union { color2_to_result: &color2_to_merged }
        },
        IntersectColors::Intersect => ColorCombination::Intersection { color1_to_result: &color1_to_result, color2_to_result: &color2_to_result },
    };
    let gen = || ElementGeneratorFromIntersectionInterleaving {
        interleaving: &interleaving,
        coloring1: &coloring1,
        coloring2: &coloring2,
        result_kmers: &result_kmers,
        result_key_kmer_marks: &key_kmer_marks,
        filter: None,
        colors: color_combination,
        n_result_colors: n_colors as usize,
    };

    log::info!("=== PHASE 2/3: Building color set fingerprints for key k-mers ===");
    let random_seed = 123123; // Todo: be more random
    // TODO: do not clone key_kmer_marks
    let (repr_kmer_marks, distinct_set_sizes, key_kmer_idx_to_set_id, _) = find_kmers_that_cover_all_distinct_sets_from_generator_that_does_not_give_duplicates(gen(), key_kmer_marks.clone(), n_colors, n_threads, random_seed);

    log::info!("=== PHASE 3/3: Build the distinct color set storage ===");
    let css = build_color_set_storage(n_colors as usize, repr_kmer_marks, distinct_set_sizes, gen(), n_threads);

    log::info!("Building rank support for key k-mer marks");
    let mut key_kmer_marks = crate::util::bitvec_to_simple_sds_bitvec(key_kmer_marks);
    key_kmer_marks.enable_rank();
    assert!(key_kmer_idx_to_set_id.len() == key_kmer_marks.rank(key_kmer_marks.len()));

    let colex_map = ColexToColorSetMap {
        sampling: key_kmer_marks,
        color_set_ids: key_kmer_idx_to_set_id,
    };

    CompactColexKmers::<CSS>::new(result_sbwt, result_lcs, colex_map, css, Some(&result_color_names))
}

#[cfg(test)]
mod tests {
    use crate::{coloring_interface::{ColorSetStorage, ColorSetView}, merge::tests::{assign_color_ids, build_named_coloring, gen_random_dna_string, NO_SAMPLING}};

    use super::IntersectColors;

    /// Intersects two colorings given as (color name, sequence) pairs (see assign_color_ids) and
    /// checks that the result has exactly the k-mers of both colorings, and that each k-mer has the
    /// union or intersection of its color sets in the inputs, computed here from the color names.
    fn check_intersect(k: usize, colors1: &[(&str, Vec<u8>)], colors2: &[(&str, Vec<u8>)], colors: IntersectColors, merge_shared_colors: bool, input_sample_distance: usize, sample_distance: usize) {
        let n_threads = 3;
        let mut names1 = Vec::<String>::new();
        let mut names2 = Vec::<String>::new();
        let seqs1 = assign_color_ids(&mut names1, colors1);
        let seqs2 = assign_color_ids(&mut names2, colors2);

        // Expected result color names, and the result color id of each color of both inputs
        let (expected_names, color1_to_result, color2_to_result): (Vec<String>, Vec<Option<usize>>, Vec<Option<usize>>) = match (colors, merge_shared_colors) {
            (IntersectColors::Union, false) => {
                let names = [names1.clone(), names2.clone()].concat();
                (names, (0..names1.len()).map(Some).collect(), (names1.len()..names1.len() + names2.len()).map(Some).collect())
            },
            (IntersectColors::Union, true) => {
                let mut names = names1.clone();
                for name in names2.iter() {
                    if !names.contains(name) { names.push(name.clone()); }
                }
                let id = |name: &String| names.iter().position(|x| x == name);
                let (map1, map2) = (names1.iter().map(id).collect(), names2.iter().map(id).collect());
                (names, map1, map2)
            },
            (IntersectColors::Intersect, _) => {
                let names: Vec<String> = names1.iter().filter(|name| names2.contains(name)).cloned().collect();
                let id = |name: &String| names.iter().position(|x| x == name);
                let (map1, map2) = (names1.iter().map(id).collect(), names2.iter().map(id).collect());
                (names, map1, map2)
            },
        };

        let ccc1 = build_named_coloring(k, &names1, &seqs1, input_sample_distance, n_threads);
        let ccc2 = build_named_coloring(k, &names2, &seqs2, input_sample_distance, n_threads);

        // The expected k-mers with their expected color sets, from searching the k-mers of ccc1 in ccc2
        let mut expected = Vec::<(Vec<u8>, Vec<usize>)>::new();
        for colex1 in 0..ccc1.sbwt().n_sets() {
            let kmer = ccc1.sbwt().access_kmer(colex1);
            if kmer.contains(&b'$') { continue; } // Dummy
            let Some(range2) = ccc2.sbwt().search(&kmer) else { continue; };
            if range2.is_empty() { continue; }
            assert_eq!(range2.len(), 1);
            let set1: Vec<usize> = ccc1.colex_to_set(colex1).iter().filter_map(|c| color1_to_result[c]).collect();
            let set2: Vec<usize> = ccc2.colex_to_set(range2.start).iter().filter_map(|c| color2_to_result[c]).collect();
            let mut set: Vec<usize> = match colors {
                IntersectColors::Union => [set1, set2].concat(),
                IntersectColors::Intersect => set1.into_iter().filter(|c| set2.contains(c)).collect(),
            };
            set.sort_unstable();
            set.dedup();
            expected.push((kmer, set));
        }

        let result = super::intersect_compact_colex_kmers(ccc1, ccc2, colors, merge_shared_colors, true, sample_distance, n_threads);
        assert_eq!(result.get_color_names(), &expected_names);
        assert_eq!(result.get_set_storage().n_colors(), expected_names.len());
        assert_eq!(result.sbwt().n_kmers(), expected.len());

        for (kmer, expected_set) in expected {
            let range = result.sbwt().search(&kmer).unwrap();
            assert_eq!(range.len(), 1);
            let mut set: Vec<usize> = result.colex_to_set(range.start).iter().collect();
            set.sort_unstable();
            assert_eq!(set, expected_set, "k = {}, k-mer {}, colors = {:?}, merge_shared_colors = {}", k, String::from_utf8_lossy(&kmer), colors, merge_shared_colors);
        }
    }

    fn check_all_modes(k: usize, colors1: &[(&str, Vec<u8>)], colors2: &[(&str, Vec<u8>)], sample_distance: usize) {
        for (colors, merge_shared_colors) in [(IntersectColors::Union, false), (IntersectColors::Union, true), (IntersectColors::Intersect, false)] {
            check_intersect(k, colors1, colors2, colors, merge_shared_colors, sample_distance, sample_distance);
            check_intersect(k, colors2, colors1, colors, merge_shared_colors, sample_distance, sample_distance);
        }
    }

    #[test]
    fn test_intersect_overlapping_substrings() {
        let _ = env_logger::try_init();
        for k in [5_usize, 8, 11] {
            for sample_distance in [NO_SAMPLING, 3] {
                let genome = gen_random_dna_string(200, 7000 + k as u64);
                let seed = 8000 + 10 * k as u64;

                // A and D are private to one index. B has unrelated sequences in the two indexes
                // apart from one shared stretch. C has overlapping substrings of the genome in both.
                let colors1 = vec![
                    ("A", genome[30..120].to_vec()),
                    ("B", gen_random_dna_string(60, seed)),
                    ("B", genome[150..200].to_vec()),
                    ("C", genome[0..100].to_vec()),
                ];
                let colors2 = vec![
                    ("B", gen_random_dna_string(60, seed + 1)),
                    ("B", genome[160..190].to_vec()),
                    ("D", genome[40..160].to_vec()),
                    ("C", genome[50..130].to_vec()),
                    ("C", genome[140..190].to_vec()),
                ];
                check_all_modes(k, &colors1, &colors2, sample_distance);
            }
        }
    }

    #[test]
    fn test_intersect_random_substrings() {
        use rand_chacha::rand_core::{RngCore, SeedableRng};
        let _ = env_logger::try_init();

        // Both indexes consist of random substrings of the same genome, with color names drawn
        // from a small pool, so the intersection cuts the unitigs and colored subunitigs of the
        // inputs in all kinds of ways.
        for k in [5_usize, 8, 11, 15] {
            for (sample_distance, seed) in [(NO_SAMPLING, 0_u64), (3, 1)] {
                let genome = gen_random_dna_string(300, 9000 + k as u64 + seed);
                let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(9500 + k as u64 * 100 + seed);
                let names = ["P", "Q", "R", "S"];
                let mut random_colors = |n: usize, n_names: usize| -> Vec<(&str, Vec<u8>)> {
                    (0..n).map(|_| {
                        let start = rng.next_u64() as usize % (genome.len() - k);
                        let len = k + rng.next_u64() as usize % 120;
                        (names[rng.next_u64() as usize % n_names], genome[start..(start + len).min(genome.len())].to_vec())
                    }).collect()
                };
                let colors1 = random_colors(8, 3);
                let colors2 = random_colors(8, 4);
                check_all_modes(k, &colors1, &colors2, sample_distance);
            }
        }
    }

    #[test]
    fn test_intersect_no_shared_colors() {
        let _ = env_logger::try_init();
        let k = 9;
        let genome = gen_random_dna_string(150, 11);
        let colors1 = vec![("X", genome[0..100].to_vec())];
        let colors2 = vec![("Y", genome[50..150].to_vec())];
        // With the intersection, every k-mer gets an empty color set
        check_all_modes(k, &colors1, &colors2, NO_SAMPLING);
    }


}
