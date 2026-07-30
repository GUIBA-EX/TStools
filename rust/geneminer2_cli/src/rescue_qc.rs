use std::collections::HashMap;

const KMER: usize = 15;
const ANCHOR: usize = 25;
const MAX_KMER_OCCURRENCES: usize = 16;
pub const MINIMUM_NEW_UNSUPPORTED_GAP: usize = 40;

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

fn maximum_bracketed_zero_run(values: &[usize]) -> usize {
    let mut longest = 0;
    let mut index = 0;
    while index < values.len() {
        if values[index] != 0 {
            index += 1;
            continue;
        }
        let start = index;
        while index < values.len() && values[index] == 0 {
            index += 1;
        }
        if start > 0 && index < values.len() {
            longest = longest.max(index - start);
        }
    }
    longest
}

/// Returns the longest internal interval without a coherent read chain.
/// A chain must contain matching 15-mers on one read/candidate diagonal and
/// retain 25 aligned bases on both sides of every reported boundary.
pub fn maximum_unsupported_internal_gap(sequence: &str, reads: &[(String, String)]) -> usize {
    if sequence.len() <= 2 * (KMER + ANCHOR) || sequence.len() < KMER {
        return 0;
    }
    let sequence = sequence.to_ascii_uppercase();
    let mut positions = HashMap::<String, Vec<usize>>::new();
    for start in 0..=sequence.len() - KMER {
        let word = &sequence[start..start + KMER];
        if word
            .bytes()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
        {
            positions.entry(word.to_owned()).or_default().push(start);
        }
    }
    positions.retain(|_, starts| starts.len() <= MAX_KMER_OCCURRENCES);

    let mut difference = vec![0_i64; sequence.len() + 1];
    for (_, read) in reads {
        for oriented in [read.to_ascii_uppercase(), reverse_complement(read)] {
            if oriented.len() < KMER {
                continue;
            }
            let mut diagonals = HashMap::<i64, (usize, usize)>::new();
            for offset in 0..=oriented.len() - KMER {
                let word = &oriented[offset..offset + KMER];
                for start in positions.get(word).into_iter().flatten().copied() {
                    let diagonal = start as i64 - offset as i64;
                    diagonals
                        .entry(diagonal)
                        .and_modify(|range| {
                            range.0 = range.0.min(start);
                            range.1 = range.1.max(start);
                        })
                        .or_insert((start, start));
                }
            }
            for (minimum, maximum) in diagonals.into_values() {
                let start = minimum + KMER + ANCHOR;
                let end = maximum.saturating_sub(ANCHOR).saturating_add(1);
                if start < end && end <= sequence.len() {
                    difference[start] += 1;
                    difference[end] -= 1;
                }
            }
        }
    }

    let mut support = Vec::with_capacity(sequence.len());
    let mut current = 0_i64;
    for delta in difference.into_iter().take(sequence.len()) {
        current += delta;
        support.push(current.max(0) as usize);
    }
    let margin = KMER + ANCHOR;
    maximum_bracketed_zero_run(&support[margin..sequence.len() - margin])
}

pub fn introduces_unsupported_internal_gap(
    before: &str,
    after: &str,
    reads: &[(String, String)],
) -> bool {
    let before_gap = maximum_unsupported_internal_gap(before, reads);
    let after_gap = maximum_unsupported_internal_gap(after, reads);
    after_gap >= MINIMUM_NEW_UNSUPPORTED_GAP
        && after_gap.saturating_sub(before_gap) >= MINIMUM_NEW_UNSUPPORTED_GAP
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dna(length: usize) -> String {
        let mut state = 0x9e37_79b9_u64;
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize] as char
            })
            .collect()
    }

    fn tiled_reads(sequence: &str) -> Vec<(String, String)> {
        (0..=sequence.len() - 120)
            .step_by(30)
            .enumerate()
            .map(|(index, start)| {
                (
                    format!("read{index}"),
                    sequence[start..start + 120].to_owned(),
                )
            })
            .collect()
    }

    #[test]
    fn detects_a_new_internal_sequence_without_a_spanning_read_chain() {
        let before = dna(500);
        let reads = tiled_reads(&before);
        let after = format!("{}{}{}", &before[..250], dna(80), &before[250..]);
        assert_eq!(maximum_unsupported_internal_gap(&before, &reads), 0);
        assert!(maximum_unsupported_internal_gap(&after, &reads) >= 40);
        assert!(introduces_unsupported_internal_gap(&before, &after, &reads));
    }

    #[test]
    fn does_not_reject_a_contig_with_coherent_tiled_read_support() {
        let before = dna(500);
        let after = format!("{}{}", before, dna(80));
        let reads = tiled_reads(&after);
        assert!(!introduces_unsupported_internal_gap(
            &after[..500],
            &after,
            &reads
        ));
    }
}
