//! Marking the key k-mers of the result of a set operation: the k-mers whose color set is stored.

use sbwt::{dbg::{Dbg, Node}, LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};

use crate::{atomic_bitmap::AtomicBitmap, colex_colored_kmers::CompactColexKmers, coloring_interface::ColorSetStorage, unitig_export::break_to_colored_subunitigs};

pub(crate) fn mark_kmer(colex: usize, marks: &AtomicBitmap) {
    marks.set(colex, true);
}

pub(crate) fn mark_in_neighbors<'a>(colex: usize, dbg: &Dbg<'a, SubsetMatrix>, marks: &AtomicBitmap) {
    let mut in_neighbor_buf = Vec::<(Node, u8)>::new(); // TODO: avoid this allocation
    dbg.push_in_neighbors(Node{id: colex}, &mut in_neighbor_buf);
    for (in_node, _) in in_neighbor_buf.iter() {
        marks.set(in_node.id, true);
    }
}

/// Marks the key k-mers of the merged (or otherwise combined) index that the colored subunitigs of
/// `coloring` require: the last k-mer of each colored subunitig, and the in-neighbors in the merged
/// DBG of the first k-mer of each colored subunitig. `to_merged` maps a colex position of
/// `coloring` to its colex position in the merged index, or None if the k-mer is not there.
#[allow(clippy::too_many_arguments)]
pub(crate) fn mark_key_kmers_for<'a, CSS: ColorSetStorage + Send + Sync>(coloring: &'a CompactColexKmers<CSS>, to_merged: impl Fn(usize) -> Option<usize> + Sync, merged_sbwt: &SbwtIndex<SubsetMatrix>, merged_lcs: &LcsArray, merged_dbg: &Dbg<'_, SubsetMatrix>, key_kmer_marks: &AtomicBitmap, visited_marks: Option<&AtomicBitmap>, n_threads: usize) -> Dbg<'a, SubsetMatrix> {

    let merged_si = StreamingIndex::new(merged_sbwt, merged_lcs);

    let k = merged_sbwt.k();
    assert_eq!(k, coloring.get_k());

    log::info!("Initializing DBG");
    let dbg = Dbg::new(coloring.sbwt(), Some(coloring.lcs()), n_threads);

    log::info!("Iterating unitigs");
    let bar = indicatif::ProgressBar::new(coloring.sbwt().n_kmers() as u64);
    dbg.iter_unitigs_with_callback(|nodes|{
        let unitig_colex_ranks = nodes.iter().map(|v| v.id).collect::<Vec<usize>>(); // TODO: avoid this allocation
        let (_, subunitig_ranges) = break_to_colored_subunitigs(&unitig_colex_ranks, &[], coloring.get_map(), coloring.sbwt());

        // Mark last k-mer of each colored subunitig, and the in-neighbors of
        // the first k-mer of each colored subunitig.
        for subunitig_range in subunitig_ranges {
            // (s,e) = (start of first k-mer, start of the k-mer after the last k-mer)
            let (s,e) = (subunitig_range.start, subunitig_range.end); 
            assert!(s < e);

            // The k-mer at position j of the unitig is node j. Its position in the merged index
            // comes from the merge interleaving, so neither the unitig string nor a search in
            // the merged index is needed.
            if let Some(last) = to_merged(nodes[e-1].id) {
                mark_kmer(last, key_kmer_marks);
            }
            if let Some(first) = to_merged(nodes[s].id) {
                mark_in_neighbors(first, merged_dbg, key_kmer_marks);
            }
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
            let mut unitig = Vec::<u8>::with_capacity(nodes.len() + k - 1);
            dbg.push_unitig_string(nodes, &mut unitig);
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

/// Marks the key k-mers around branches of the DBG (the last k-mer of each unitig and the
/// in-neighbors of its first k-mer), and samples every sample_distance-th k-mer of each unitig.
/// This part is independent of coloring.
pub(crate) fn mark_structural_key_kmers(sbwt: &SbwtIndex<SubsetMatrix>, dbg: &Dbg<'_, SubsetMatrix>, key_kmer_marks: &AtomicBitmap, sample_distance: usize, n_threads: usize) {
    let bar = indicatif::ProgressBar::new(sbwt.n_kmers() as u64);
    dbg.iter_unitigs_with_callback(|nodes|{
        mark_in_neighbors(nodes.first().unwrap().id, dbg, key_kmer_marks);
        mark_kmer(nodes.last().unwrap().id, key_kmer_marks);

        for v in nodes.iter().rev().step_by(sample_distance) {
            key_kmer_marks.set(v.id,  true);
        }
        bar.inc(nodes.len() as u64);
    }, n_threads);
    bar.finish();
}
