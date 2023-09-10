use jseqio::seq_db::SeqDB;
use std::collections::HashMap;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Orientation{
    Forward,
    Reverse,
}

#[derive(Copy, Clone, Debug, PartialEq)]
enum Position{
    Start,
    End,
}

#[derive(Copy, Clone, Debug, PartialEq)]
struct MapValue{
    unitig_id: usize,
    position: Position,
}


impl Orientation {
    pub fn flip(&self) -> Orientation {
        match self {
            Orientation::Forward => Orientation::Reverse,
            Orientation::Reverse => Orientation::Forward,
        }
    }
}


// Edge from 'from' to 'to' with given orientations means that
// with those orientations, the last k-1 bases of 'from' are the 
// same as the first k-1 bases of 'to'.
// This means that if (u, v, u_o, v_o) is a valid edge, then
// (v, u, flip u_o, flip v_o) is also a valid edge.
#[derive(Copy, Clone, Debug)]
pub struct Edge{
    pub from: usize,
    pub to: usize,
    pub from_orientation: Orientation,
    pub to_orientation: Orientation,
}

fn insert_if_not_present<'key, 'borrow>(map: &'borrow mut HashMap<&'key [u8], Vec<MapValue>>, key: &'key [u8]){
    if !map.contains_key(key){
        map.insert(&key, Vec::<MapValue>::new());
    }
}

pub struct DBG{
    pub unitigs: SeqDB, // A sequence database with random access to the i-th unitig
    pub edges: Vec<Vec<Edge>> // edges[i] = outgoing edges from unitig i
}

fn push_edges(from: usize, from_orientation: Orientation, to_orientation: Orientation, to_position: Position, linking_kmer: &[u8], edges: &mut Vec<Vec<Edge>>, borders: &HashMap<&[u8], Vec<MapValue>>){
    if let Some(vec) = borders.get(linking_kmer){
        for x in vec.iter(){
            if x.position == to_position {
                let edge = Edge{from, to: x.unitig_id, from_orientation, to_orientation};
                edges[from].push(edge);
            }
        }
    }
}

pub fn build_dbg(unitigs: SeqDB, unitigs_rc: SeqDB, k: usize) -> DBG{
    let mut borders: HashMap<&[u8], Vec<MapValue>> = HashMap::new(); // (k-1)-mer to locations of that k-mer

    let n = unitigs.sequence_count();

    log::info!("Hashing border k-mers");

    // Build borders map
    for i in 0..n{
        let unitig = unitigs.get(i);

        let first = &unitig.seq[..k-1];
        let last = &unitig.seq[unitig.seq.len()-(k-1)..];

        insert_if_not_present(&mut borders, first);
        insert_if_not_present(&mut borders, last);

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

    log::info!("Finding overlaps");
    let mut edges = Vec::<Vec::<Edge>>::new();
    edges.resize_with(n, || Vec::<Edge>::new()); // Allocate n edge lists

    // Build edges
    for i in 0..n{
        let unitig = unitigs.get(i);
        let unitig_rc = unitigs_rc.get(i);

        let first = &unitig.seq[..k-1];
        let last = &unitig.seq[unitig.seq.len()-(k-1)..];

        let first_rc = &unitig_rc.seq[unitig_rc.seq.len()-(k-1)..];
        let last_rc = &unitig_rc.seq[..k-1];


        push_edges(i, Orientation::Forward, Orientation::Forward, Position::Start, last, &mut edges, &borders);
        push_edges(i, Orientation::Forward, Orientation::Reverse, Position::End, last_rc, &mut edges, &borders);
        push_edges(i, Orientation::Reverse, Orientation::Reverse, Position::End, first, &mut edges, &borders);
        push_edges(i, Orientation::Reverse, Orientation::Forward, Position::Start, first_rc, &mut edges, &borders);

    }

    DBG {unitigs, edges}
}
