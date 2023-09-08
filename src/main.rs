use jseqio::seq_db::SeqDB;
use jseqio::reader::DynamicFastXReader;
use jseqio::record::{OwnedRecord, RefRecord};
use std::collections::HashMap;
use std::hash::Hash;

#[derive(Copy, Clone, Debug)]
enum Orientation{
    Forward,
    Reverse,
}

#[derive(Copy, Clone, Debug)]
enum Position{
    Start,
    End,
}

#[derive(Copy, Clone, Debug)]
struct Edge{
    from: usize,
    to: usize,
    from_orientation: Orientation,
    to_orientation: Orientation,
}

struct DBG{
    unitigs: SeqDB, // A sequence database with random access to the i-th unitig
    edges: Vec<Vec<Edge>> // edges[i] = outgoing edges from unitig i
}

struct MapValue{
    unitig_id: usize,
    position: Position,
}

fn ensure_is_in_map(map: &mut HashMap<Vec<u8>, Vec<MapValue>>, key: &[u8]){
    if !map.contains_key(key){
        map.insert(key.to_owned(), Vec::<MapValue>::new());
    }
}

fn rc(c: u8) -> u8{
    match c{
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        _ => panic!("Invalid character: {}", c),
    }
}

fn build_dbg(unitigs: SeqDB, k: usize) -> DBG{
    let mut borders: HashMap<Vec<u8>, Vec<MapValue>> = HashMap::new(); // (k-1)-mer to locations of that k-mer

    let n = unitigs.sequence_count();

    // Build borders map
    for i in 0..n{
        let unitig = unitigs.get(i).unwrap();

        let first = &unitig.seq[0..k-1];
        let last = &unitig.seq[unitig.seq.len()-(k-1)..];

        ensure_is_in_map(&mut borders, first);
        ensure_is_in_map(&mut borders, last);

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

    let mut edges = Vec::<Vec::<Edge>>::new();
    edges.resize_with(n, || Vec::<Edge>::new()); // Allocate n edge lists

    // Build edges
    for i in 0..n{
        let unitig = unitigs.get(i).unwrap().to_owned();
        let unitig_rc = OwnedRecord{
            head: unitig.head.clone(),
            seq: unitig.seq.iter().rev().map(|&c| rc(c)).collect(),
            qual: unitig.qual.clone(),
        };

        let first = &unitig.seq[0..k-1];
        let last = &unitig.seq[unitig.seq.len()-(k-1)..];

        let first_rc = &unitig_rc.seq[0..k-1];
        let last_rc = &unitig_rc.seq[unitig_rc.seq.len()-(k-1)..];

        // Forward -> Forward edges
        borders.get(last).unwrap().iter().for_each(|&MapValue{unitig_id, position}|{
            match position{
                Position::Start => {
                    let edge = Edge{
                        from: i,
                        to: unitig_id,
                        from_orientation: Orientation::Forward,
                        to_orientation: Orientation::Forward,
                    };
                    edges[i].push(edge); // Does not work as an edge
                }
                Position::End => ()
            }
        });

        // Forward -> Reverse edges
        if borders.contains_key(last_rc){
            borders.get(last_rc).unwrap().iter().for_each(|&MapValue{unitig_id, position}|{
                match position {
                    Position::End => {
                        let edge = Edge{
                            from: i,
                            to: unitig_id,
                            from_orientation: Orientation::Forward,
                            to_orientation: Orientation::Reverse,
                        };
                        edges[i].push(edge);
                    }
                    Position::Start => () // Does not work as an edge
                }
            });
        }

        // Reverse -> Reverse edges
        if borders.contains_key(first){
            borders.get(first).unwrap().iter().for_each(|&MapValue{unitig_id, position}|{
                match position {
                    Position::End => {
                        let edge = Edge{
                            from: i,
                            to: unitig_id,
                            from_orientation: Orientation::Reverse,
                            to_orientation: Orientation::Reverse,
                        };
                        edges[i].push(edge);
                    }
                    Position::Start => () // Does not work as an edge
                }
            });
        }

        // Reverse -> Forward edges
        if borders.contains_key(first_rc){
            borders.get(first_rc).unwrap().iter().for_each(|&MapValue{unitig_id, position}|{
                match position {
                    Position::Start => {
                        let edge = Edge{
                            from: i,
                            to: unitig_id,
                            from_orientation: Orientation::Reverse,
                            to_orientation: Orientation::Forward,
                        };
                        edges[i].push(edge);
                    }
                    Position::End => () // Does not work as an edge
                }
            });
        }
    }

    DBG {unitigs, edges}
}

fn main() {
    let filename = std::env::args().nth(1).unwrap();
    let k = std::env::args().nth(2).unwrap().parse::<usize>().unwrap();
    let reader = DynamicFastXReader::from_file(&filename).unwrap();
    let db = reader.into_db().unwrap();
    let dbg = build_dbg(db, k);
    dbg!(dbg.edges);
}
