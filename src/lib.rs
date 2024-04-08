pub mod dbg;

use dbg::*;

/// A stream of ASCII-encoded DNA-sequences. This is not necessarily a standard Rust iterator
/// because we want to support streaming sequences from disk, which is not possible
/// with a regular iterator due to lifetime constraints of the Iterator trait.
pub trait SeqStream<'a>
{
    /// Return the next sequence in the stream, or None if the stream is exhausted.
    /// You may be wondering: why such complicated lifetime annotations?
    /// Lifetime 'a is not used in this signature at all. We just state that 'a is longer than 'b.
    /// The existence of this free lifetime parameter 'a will be useful for implementors
    /// of this trait who now have one more lifetime parameter to play with. See the comments
    /// in `impl<'a, I: Iterator<Item = &'a [u8]>> SeqStream<'a> for I` for an explanation
    /// for why this is required.
    fn stream_next<'b>(&'b mut self) -> Option<&'b [u8]> where 'a : 'b;
}

/// Any iterator of &'a [u8] that is valid for 'a can work as a SeqStream<'a>.
impl<'a, 'c: 'a, I: Iterator<Item=&'c [u8]>> SeqStream<'a> for I
{
    fn stream_next<'b>(&'b mut self) -> Option<&'b [u8]> where 'a : 'b {
        // We need the bound 'a: 'b because otherwise the compiler is unhappy because
        // we are returning &'a but the function signature says we need to return &'b.
        // This is okay only if 'a is a subtype of 'b, i.e. 'a lasts longer than 'b.
        // In that case, 'a is also a 'b so we can do the conversion. But if we don't
        // explicitly require 'a: 'b, then the compiler does not see any relationship
        // between 'a and 'b and refuses to do the conversion. But if we have 'a: 'b
        // here, then we also need it in the trait definition, even though 'a does not
        // appear in the function signature of the trait. But we need to say already in
        // the trait definition that there exists some lifetime 'a that is longer than
        // 'b in order to be able to relate the lifetime of the items to the lifetime
        // of the borrow in this implementation.
        self.next()
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum Direction{
    Forward,
    Backward,
}

impl Direction {
    pub fn flip(&self) -> Direction {
        match self {
            Direction::Forward => Direction::Backward,
            Direction::Backward => Direction::Forward,
        }
    }
}

/// Given a stream of unitigs, returns a vector of orientations, one for each sequence, aiming to
/// minimize the number of unitigs which do not have an incoming edge in the de Bruijn graph
/// of order k.
pub fn optimize_unitig_orientation<'a, SS: SeqStream<'a>>(mut input: SS, k: usize) -> Vec<Orientation>{
    let mut db = jseqio::seq_db::SeqDB::new();
    let mut rc_db = jseqio::seq_db::SeqDB::new();
    let mut rc_buf = Vec::<u8>::new();
    while let Some(seq) = input.stream_next() {
        db.push_seq(seq);

        rc_buf.clear();
        rc_buf.extend(seq);
        jseqio::reverse_complement_in_place(&mut rc_buf);
        rc_db.push_seq(&rc_buf);
    }

    let dbg = DBG::build(db, rc_db, k);

    pick_orientations_with_non_switching_bfs(&dbg)
}

// BFS forward or backward from the root, such that we do not switch direction during the search
// The forward search visits a tree of nodes starting from v such that each root-to-leaf 
// path is a sequence of unitigs that can be merged together left-to-right starting from 
// the root, using the orientations saved during the search. The backward search is the 
// same but in the other direction.
fn bfs(mut root: usize, dir: Direction, oris: &mut [Orientation], visited: &mut [bool], dbg: &DBG ) {

    if dir == Direction::Backward {
        // If we are going backward, it's the same as walking forward but
        // starting from the reverse complement twin node.
        root = dbg.twin(root);
    }

    let mut queue = std::collections::VecDeque::<usize>::new();
    queue.push_back(root);

    while let Some(v) = queue.pop_front() {

        if visited[v % dbg.n_unitigs] { continue }
        visited[v % dbg.n_unitigs] = true;

        oris[v % dbg.n_unitigs] = match (dir, v < dbg.n_unitigs) {
            (Direction::Forward, true) => Orientation::Forward,
            (Direction::Forward, false) => Orientation::Reverse,
            (Direction::Backward, true) => Orientation::Reverse,
            (Direction::Backward, false) => Orientation::Forward,
        };

        for &u in dbg.out_edges[v].iter() {
            queue.push_back(u);
        }
    }

}

// BFS that can switch direction in the middle of the search.
fn direction_switching_bfs(mut root: usize, root_dir: Direction, oris: &mut [Orientation], visited: &mut [bool], dbg: &DBG ) {

    if root_dir == Direction::Backward {
        // If we are going backward, it's the same as walking forward but
        // starting from the reverse complement twin node.
        root = dbg.twin(root);
    }

    let mut queue = std::collections::VecDeque::<(usize, Direction)>::new(); // Nodes from directed graph on [0, 2n), walking in Direction
    queue.push_back((root, root_dir));

    while let Some((v, dir)) = queue.pop_front() {

        if visited[v % dbg.n_unitigs] { continue }
        visited[v % dbg.n_unitigs] = true;

        oris[v % dbg.n_unitigs] = match (dir, v < dbg.n_unitigs) {
            (Direction::Forward, true) => Orientation::Forward,
            (Direction::Forward, false) => Orientation::Reverse,
            (Direction::Backward, true) => Orientation::Reverse,
            (Direction::Backward, false) => Orientation::Forward,
        };

        // Keep walking in the same direction
        for &u in dbg.out_edges[v].iter() {
            queue.push_back((u, dir));
        }

        // Switch direction
        for &u in dbg.out_edges[dbg.twin(v)].iter() {
            queue.push_back((u, dir.flip()));
        }
    }
}

pub fn pick_orientations_with_non_switching_bfs(dbg: &DBG) -> Vec<Orientation>{
    let n = dbg.n_unitigs;
    let mut orientations = Vec::<Orientation>::new();
    orientations.resize(n, Orientation::Forward);

    let mut visited = vec![false; n];

    for v in 0..n {
        if !visited[v] {
            bfs(v, Direction::Forward, &mut orientations, &mut visited, dbg); // Visits v
            visited[v] = false;
            bfs(v, Direction::Backward, &mut orientations, &mut visited, dbg); // Visits v again
        }
    };


    orientations
}

pub fn pick_orientations_with_switching_bfs(dbg: &DBG) -> Vec<Orientation>{
    let n = dbg.n_unitigs;
    let mut orientations = Vec::<Orientation>::new();
    orientations.resize(n, Orientation::Forward);

    let mut visited = vec![false; n];

    for v in 0..n {
        direction_switching_bfs(v, Direction::Forward, &mut orientations, &mut visited, dbg); // Visits v
    };

    orientations
}

// Returns the number of unitigs that do not have a predecessor
pub fn evaluate(choices: &[Orientation], dbg: &DBG) -> usize{
    let n = dbg.n_unitigs;
    let mut has_pred = vec![false; n];

    #[allow(unused_parens)]
    for v in 0..(n*2){
        for &u in dbg.out_edges[v].iter(){
            let v_flipped = (choices[v%n] == Orientation::Reverse);
            let u_flipped = (choices[u%n] == Orientation::Reverse);
            if (v_flipped == (v >= n)) && (u_flipped == (u >= n)){
                has_pred[u % n] = true;
            }
        }
    }
    
    // Return the number of 1-bits in has_pred
    has_pred.iter().fold(0_usize, |sum, &x| sum + x as usize)
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;

    use super::*;
    use rand::{RngCore, SeedableRng};
    use rand::rngs::StdRng;
    
    fn generate_random_dna_string(length: usize, seed: u64) -> Vec<u8> {
        let mut rng = StdRng::seed_from_u64(seed);
        let alphabet: Vec<u8> = vec![b'A', b'C', b'G', b'T'];
    
        let mut s = Vec::<u8>::new();
        for _ in 0..length{
            s.push(alphabet[(rng.next_u64() % 4) as usize]);
        }
        s
    }

    fn change(c: u8) -> u8 {
        match c {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            b'T' => b'A',
            _ => panic!(),
        }
    }
    
    fn get_dbg(seqs: Vec<&[u8]>, k: usize) -> DBG {
        let rc_seqs = seqs.iter().map(|s| jseqio::reverse_complement(s));

        let mut db = jseqio::seq_db::SeqDB::new();
        for seq in seqs.iter() {
            db.push_seq(seq);
        }

        let mut rc_db = jseqio::seq_db::SeqDB::new();
        for rc_seq in rc_seqs {
            rc_db.push_seq(&rc_seq);
        }

        DBG::build(db, rc_db, k)
    }

    #[test]
    #[allow(non_snake_case)]
    fn sanity_check_algorithms(){
        // Make a graph like this
        // o--\          /--o 
        //     >--------< 
        // o--/          \--o

        let seed = 123;
        let S = generate_random_dna_string(100, seed);
        let k = 10;

        let s1 = S[20..30].to_vec();
        let mut s2 = s1.to_owned();
        s2[0] = change(s2[0]);

        let mut s_middle = S[21..31].to_owned();
        let s_middle2 = S[22..50].to_owned();

        let s3 = S[41..51].to_vec();
        let mut s4 = s3.to_owned();
        s4[9] = change(s4[9]);

        jseqio::reverse_complement_in_place(&mut s_middle);

        // Run switching BFS

        // Start with s1 so that we see that the search switches backward at middle
        // Todo: make sure the assert_eq! later would fail if this does not happen.
        let dbg = get_dbg(vec![&s1, &s2, &s_middle, &s_middle2, &s3, &s4], k);

        let orientations = pick_orientations_with_switching_bfs(&dbg);
        for ori in orientations.iter() {
            dbg!(&ori);
        }

        let n_has_pred = evaluate(&orientations, &dbg);
        let n_sources = dbg.n_unitigs - n_has_pred;
        dbg!(n_sources);
        assert_eq!(n_sources, 2);

        // Run non-switching BFS

        // Start with s_middle, which reversed, so we see that the search goes both ways
        let dbg = get_dbg(vec![&s_middle, &s1, &s2, &s_middle2, &s3, &s4], k);

        let orientations = pick_orientations_with_non_switching_bfs(&dbg);
        for ori in orientations.iter() {
            dbg!(&ori);
        }

        let n_has_pred = evaluate(&orientations, &dbg);
        let n_sources = dbg.n_unitigs - n_has_pred;
        dbg!(n_sources);
        assert_eq!(n_sources, 2);

    }    

    use Orientation::*;

    #[test]
    fn straight_line(){
        // Input
        let k = 3;
        let data: Vec<&[u8]> = vec![b"TCG", b"ATC", b"ATG", b"ACC", b"CCG"];

        // TCG
        // ATC
        // ATG
        // ACC
        // CCG

        // Run the test
        let orientations = optimize_unitig_orientation(data.clone().into_iter(), k);

        // Two possible answers
        let ans1 = vec![Reverse, Reverse, Forward, Forward, Forward];
        let ans2 = vec![Forward, Forward, Reverse, Reverse, Reverse];
        assert!(orientations == ans1 || orientations == ans2);

    }
}

