use std::collections::{HashSet};
use rand::seq::SliceRandom;

fn reverse_complement(s: &str) -> String {
    s.chars().map(|c| match c {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        _ => unreachable!(),
    }).collect::<String>().chars().rev().collect()
}

fn extend_simplitig_forward(k: usize, mut simplitig: String, k_set: &mut HashSet<String>) -> String {
    let mut extending = true;
    while extending {
        extending = false;
        let q = &simplitig[simplitig.len() - k + 1..];
        for x in &['A', 'C', 'G', 'T'] {
            let kmer = format!("{}{}", q, x);
            if k_set.contains(&kmer) {
                extending = true;
                simplitig.push(*x);
                k_set.remove(&kmer);
                k_set.remove(&reverse_complement(&kmer));
                break;
            }
        }
    }
    simplitig
}

fn get_maximal_simplitig(k_set: &mut HashSet<String>, initial_kmer: String) -> String {
    let mut simplitig = initial_kmer.clone();
    k_set.remove(&initial_kmer);
    k_set.remove(&reverse_complement(&initial_kmer));
    simplitig = extend_simplitig_forward(initial_kmer.len(), simplitig, k_set);
    let mut simplitig_rc = reverse_complement(&simplitig);
    simplitig_rc = extend_simplitig_forward(initial_kmer.len(), simplitig_rc, k_set);
    simplitig_rc
}

fn compute_simplitigs(kmers: Vec<String>) -> HashSet<String> {
    let mut k_set: HashSet<String> = HashSet::new();
    for kmer in &kmers {
        k_set.insert(kmer.clone());
        k_set.insert(reverse_complement(kmer));
    }
    let mut simplitigs: HashSet<String> = HashSet::new();
    while !k_set.is_empty() {
        let initial_kmer = k_set.iter().next().unwrap().clone(); // "Random choice"
        let simplitig = get_maximal_simplitig(&mut k_set, initial_kmer);
        simplitigs.insert(simplitig);
    }
    simplitigs
}

#[cfg(test)]
mod tests {
    use super::compute_simplitigs;


    #[test]
    fn small_example(){
        let k = 6;
        let raw_input: Vec<&[u8]> = vec![b"AAACCC", b"CCCGGG", b"GGGTTT"];
        let kmers: Vec<String> = raw_input.iter().map(|&s| s.windows(k)).fold(Vec::<String>::new(), |mut acc, it| {
            let s_kmers: Vec<String> = it.into_iter().map(|s| String::from_utf8_lossy(s).into_owned()).collect();
            acc.extend(s_kmers);
            acc
        });
        let result = compute_simplitigs(kmers);
        dbg!(result);
    }

}

// Todo:
// String -> Vec[u8]
// Use unitigs 