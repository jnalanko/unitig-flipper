use std::collections::HashMap;
use crate::dbg::Orientation;

// If there are n unitigs, we have 2n nodes
// Nodes 0..n-1 are the nodes in their orientations in the input, and
// nodes n..2n-1 are the reverse complemented versions of those, so that
// nodes v and v+n correnspond to each other.
// Now we have just a regular directed graph on 2n nodes. 
// There is an edge from v to u if v[|v|-k..|v|-1] = u[0..k-1]
// We store neighbor lists for both incoming and outgoing edges.
pub struct DBG {
    out_edges: Vec<Vec<usize>>,
    in_edges: Vec<Vec<usize>>,
    n_unitigs: usize,
}

#[derive(Copy, Clone, Debug, PartialEq)]
enum Position{ // Used internally in construction
    Start,
    End,
}

#[derive(Copy, Clone, Debug, PartialEq)]
struct MapValue{ // Used internally in construction
    unitig_id: usize,
    position: Position,
}

impl DBG {

    fn new(n_unitigs : usize) -> Self {
        DBG{out_edges: vec![Vec::new(); n_unitigs * 2], in_edges: vec![Vec::new(); n_unitigs * 2], n_unitigs}
    }

    pub fn twin(&self, v: usize) -> usize {
        (v + self.n_unitigs) % (2*self.n_unitigs)
    }

    pub fn add_edge(&mut self, from_node: usize, to_node: usize, from_orientation: Orientation, to_orientation: Orientation) {
        let v = from_node + ((from_orientation == Orientation::Reverse) as usize) * self.n_unitigs;
        let u = to_node + ((to_orientation == Orientation::Reverse) as usize) * self.n_unitigs;
        self.out_edges[v].push(u);
        self.in_edges[u].push(v);
    }

    fn insert_if_not_present<'key>(map: &mut HashMap<&'key [u8], Vec<MapValue>>, key: &'key [u8]){
        if !map.contains_key(key){
            map.insert(key, Vec::<MapValue>::new());
        }
    }

    pub fn build(unitigs: jseqio::seq_db::SeqDB, unitigs_rc: jseqio::seq_db::SeqDB, k: usize) -> DBG{

        use Orientation::*;

        let mut borders: HashMap<&[u8], Vec<MapValue>> = HashMap::new(); // (k-1)-mer to locations of that k-mer

        let n = unitigs.sequence_count();

        log::info!("Hashing border k-mers");

        // Build borders map
        for i in 0..n{
            let unitig = unitigs.get(i);

            let first = &unitig.seq[..k-1];
            let last = &unitig.seq[unitig.seq.len()-(k-1)..];

            Self::insert_if_not_present(&mut borders, first);
            Self::insert_if_not_present(&mut borders, last);

            borders.get_mut(first).unwrap().push(
                MapValue{
                    unitig_id: i, 
                    position: Position::Start, 
                }
            );

            borders.get_mut(last).unwrap().push(
                MapValue{
                    unitig_id: i, 
                    position: Position::End, 
                }
            );
        }

        log::info!("Building edges");
        let mut dbg = DBG::new(n);
        for i in 0..n{
            // List all outgoing edges from node i or its rev. comp. twin
            let unitig = unitigs.get(i);
            let unitig_rc = unitigs_rc.get(i);

            let first = &unitig.seq[..k-1];
            let last = &unitig.seq[unitig.seq.len()-(k-1)..];

            let first_rc = &unitig_rc.seq[unitig_rc.seq.len()-(k-1)..];
            let last_rc = &unitig_rc.seq[..k-1];


            if let Some(occs) = borders.get(last){
                for right in occs {
                    if right.position == Position::Start {
                        dbg.add_edge(i, right.unitig_id, Forward, Forward)
                    }
                }
            }

            if let Some(occs) = borders.get(last_rc){
                for right in occs {
                    if right.position == Position::End{
                        dbg.add_edge(i, right.unitig_id, Forward, Reverse)
                    }
                }
            }

            if let Some(occs) = borders.get(first_rc){
                for right in occs {
                    if right.position == Position::Start {
                        dbg.add_edge(i, right.unitig_id, Reverse, Forward)
                    }
                }
            }

            if let Some(occs) = borders.get(first){
                for right in occs {
                    if right.position == Position::End{
                        dbg.add_edge(i, right.unitig_id, Reverse, Reverse)
                    }
                }
            }
        }

        dbg

    }
}

fn bfs_from(root: usize, oris: &mut [Orientation], visited: &mut [bool], dbg: &DBG ) {

    let mut queue = std::collections::VecDeque::<usize>::new();
    queue.push_back(root);

    while let Some(v) = queue.pop_front() {

        if visited[v % dbg.n_unitigs] { continue }
        visited[v % dbg.n_unitigs] = true;
        if v < dbg.n_unitigs {
            oris[v % dbg.n_unitigs] = Orientation::Forward;
        }
        else {
            oris[v % dbg.n_unitigs] = Orientation::Reverse;
        }

        for &u in dbg.out_edges[v].iter() {
            queue.push_back(u);
        }
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

fn walk_from(mut v: usize, dir: Direction, oris: &mut [Orientation], visited: &mut [bool], dbg: &DBG) {

    if dir == Direction::Backward {
        // If we are going backward, it's the same as walking forward but
        // starting from the reverse complement twin node.
        v = dbg.twin(v);
    }

    //eprintln!("walking from {} in direction {:?}", v % dbg.n_unitigs, dir);
    loop {

        if visited[v % dbg.n_unitigs] { return }
        visited[v % dbg.n_unitigs] = true;

        oris[v % dbg.n_unitigs] = match (dir, v < dbg.n_unitigs) {
            (Direction::Forward, true) => Orientation::Forward,
            (Direction::Forward, false) => Orientation::Reverse,
            (Direction::Backward, true) => Orientation::Reverse,
            (Direction::Backward, false) => Orientation::Forward,
        };
        //eprintln!("Oriented {} to {:?}", v % dbg.n_unitigs, oris[v % dbg.n_unitigs]);

        let mut have_next = false;
        for &u in dbg.out_edges[v].iter() {
            if !visited[u % dbg.n_unitigs] {
                have_next = true;
                v = u;
                break;
            }
        }

        if !have_next {
            break;
        }
    }

}

pub fn pick_orientations_rethink(dbg: &DBG) -> Vec<Orientation>{
    let n = dbg.n_unitigs;
    let mut orientations = Vec::<Orientation>::new();
    orientations.resize(n, Orientation::Forward);

    let mut visited = vec![false; n];

    // For all source nodes
    (0..2*n).filter(|&v| !dbg.out_edges[v].is_empty()).for_each(|v|{
        bfs_from(v, &mut orientations, &mut visited, dbg);
    });


    orientations
}

pub fn pick_orientations_simplitigs(dbg: &DBG) -> Vec<Orientation>{
    let n = dbg.n_unitigs;
    let mut orientations = Vec::<Orientation>::new();
    orientations.resize(n, Orientation::Forward);

    let mut visited = vec![false; n];

    let mut string_count = 0_usize;
    for v in 0..n {
        if !visited[v] {
            walk_from(v, Direction::Forward, &mut orientations, &mut visited, dbg); // Visits v
            visited[v] = false;
            walk_from(v, Direction::Backward, &mut orientations, &mut visited, dbg); // Visits v again
            string_count += 1;
        }
    }

    eprintln!("Constructed {} strings", string_count);

    orientations
}

// Returns the number of unitigs that do not have a predecessor
pub fn evaluate(choices: &[Orientation], dbg: &DBG) -> usize{
    let n = dbg.n_unitigs;
    let mut has_pred = vec![false; n];

    for v in 0..(n*2){
        for &u in dbg.out_edges[v].iter(){
            let v_flipped = (choices[v%n] == Orientation::Reverse);
            let u_flipped = (choices[u%n] == Orientation::Reverse);
            if (v_flipped == (v >= n)) && (u_flipped == (u >= n)){
                has_pred[u % n] = true;
            }
        }
    }
    
    // Return the number of 1-bit in has_pred
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

        let orientations = pick_orientations_simplitigs(&dbg);
        for ori in orientations.iter() {
            dbg!(&ori);
        }

        let n_sources = evaluate(&orientations, &dbg);
        dbg!(n_sources);

        assert_eq!(n_sources, 2);

    }    
}