pub mod dbg;

use dbg::*;

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
    

    #[test]
    fn test_rethink(){
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
        let mut s_middle2 = S[22..50].to_owned();

        let s3 = S[41..51].to_vec();
        let mut s4 = s3.to_owned();
        s4[9] = change(s4[9]);


        jseqio::reverse_complement_in_place(&mut s_middle);
        let seqs = vec![s_middle,s1,s2,s_middle2,s3,s4];

        for s in seqs.iter() {
            eprintln!("{}", String::from_utf8_lossy(s));
        }

        let rc_seqs = seqs.iter().map(|s| jseqio::reverse_complement(s));

        let mut db = jseqio::seq_db::SeqDB::new();
        for seq in seqs.iter() {
            db.push_seq(seq);
        }

        let mut rc_db = jseqio::seq_db::SeqDB::new();
        for rc_seq in rc_seqs {
            rc_db.push_seq(&rc_seq);
        }

        let dbg = DBG::build(db, rc_db, k);

        let orientations = pick_orientations_with_switching_bfs(&dbg);
        for ori in orientations.iter() {
            dbg!(&ori);
        }

        let n_sources = evaluate(&orientations, &dbg);
        dbg!(n_sources);

        assert_eq!(n_sources, 3);

    }    
}

