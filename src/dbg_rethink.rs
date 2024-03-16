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