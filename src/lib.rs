pub mod dbg;
mod dbg_rethink;

use log::info;
use dbg::Orientation;
use dbg::Orientation::{Forward, Reverse};

pub fn new_algorithm(dbg: &dbg::DBG) -> Vec<Orientation> {
    let terminal_nodes = (0..dbg.unitigs.sequence_count()).filter(|&v| is_terminal(dbg, v));

    let mut orientations = Vec::<Orientation>::new();
    orientations.resize(dbg.unitigs.sequence_count(), Orientation::Forward);
    let mut visited = vec![false; dbg.unitigs.sequence_count()];

    // Orient the terminals
    for v in terminal_nodes {
        let v_orientation = get_terminal_orientation(dbg, v);
        bfs_from(v, v_orientation, dbg, &mut visited, &mut orientations);
    }

    // BFS from the rest.
    for component_root in 0..dbg.unitigs.sequence_count(){

        if visited[component_root]{
            continue;
        }

        bfs_from(component_root, Forward, dbg, &mut visited, &mut orientations);
    }

    orientations
}

pub fn pick_orientations(dbg: &dbg::DBG) -> Vec<Orientation>{
    let mut orientations = Vec::<Orientation>::new();
    orientations.resize(dbg.unitigs.sequence_count(), Orientation::Forward);

    let mut visited = vec![false; dbg.unitigs.sequence_count()];

    
    let mut n_components: usize = 0;    

    // BFS from terminal nodes
    for component_root in 0..dbg.unitigs.sequence_count(){

        if !is_terminal(dbg, component_root) {
            // This node will be visited from some other node
            continue;
        }

        if visited[component_root]{
            continue;
        }

        bfs_from(component_root, Forward, dbg, &mut visited, &mut orientations);
        n_components += 1;
    }

    let n_visited = visited.iter().fold(0_usize, |sum, &x| sum + x as usize);
    info!("{:.2}% of unitigs visited so far... proceeding to clean up the rest", 100.0 * n_visited as f64 / dbg.unitigs.sequence_count() as f64);

    // BFS from the rest.
    for component_root in 0..dbg.unitigs.sequence_count(){

        if visited[component_root]{
            continue;
        }

        bfs_from(component_root, Forward, dbg, &mut visited, &mut orientations);
        n_components += 1;
    }

    info!("Done. Total {} BFS rounds executed", n_components);

    orientations
}

fn is_terminal(dbg: &dbg::DBG, unitig_id: usize) -> bool{
    // Todo: prove that this is correct
    let mut has_fw = false;
    let mut has_bw = false;
    for e in dbg.edges[unitig_id].iter(){ 
        if e.from == e.to{
            continue; // Self-loops don't count
        }
        has_fw |= e.from_orientation == Forward;
        has_bw |= e.from_orientation == Reverse;
    }

    !(has_fw && has_bw)
}

// Assumes unitig_id is a terminal node
fn get_terminal_orientation(dbg: &dbg::DBG, unitig_id: usize) -> Orientation{
    let has_fw = dbg.edges[unitig_id].iter().filter(|e| e.from != e.to).any(|e| e.from_orientation == Forward);
    let has_bw = dbg.edges[unitig_id].iter().filter(|e| e.from != e.to).any(|e| e.from_orientation == Reverse);
    assert!(!has_fw || !has_bw);
    if has_fw {
        // In forward orientation this unitig is always a source.
        // Which means that in reverse orientation it is always a sink
        Reverse 
    } else {
        Forward
    }
}

fn bfs_from(root: usize, root_orientation: Orientation, dbg: &dbg::DBG, visited: &mut [bool], orientations: &mut [Orientation]){
    let mut queue = std::collections::VecDeque::<(usize, Orientation)>::new();

    // Arbitrarily orient the root as forward
    //queue.push_back((root, Orientation::Forward));
    queue.push_back((root, root_orientation));

    // BFS from root and orient all reachable unitigs the same way
    while let Some((unitig_id, orientation)) = queue.pop_front(){
        if visited[unitig_id]{
            continue;
        }

        visited[unitig_id] = true;
        orientations[unitig_id] = orientation;

        for edge in dbg.edges[unitig_id].iter(){
            if edge.from_orientation != orientation {
                // If we came to this node in the forward orientation, we may only
                // leave on edges that leave in the forward orientation. And vice versa.
                continue;
            }

            queue.push_back((edge.to, edge.to_orientation));
        }
    }
}

// Returns the number of unitigs that do not have a predecessor
pub fn evaluate(choices: &[Orientation], dbg: &dbg::DBG) -> usize{
    let mut has_pred = vec![false; dbg.unitigs.sequence_count()];

    for v in 0..dbg.unitigs.sequence_count(){
        for edge in dbg.edges[v].iter(){
            let u = edge.to;
            if choices[v] == edge.from_orientation && choices[u] == edge.to_orientation{
                has_pred[u] = true;
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
    fn test_new_algorithm(){
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

        println!("{}", String::from_utf8_lossy(&s1));
        println!("{}", String::from_utf8_lossy(&s2));
        println!("{}", String::from_utf8_lossy(&s_middle));
        println!("{}", String::from_utf8_lossy(&s_middle2));
        println!("{}", String::from_utf8_lossy(&s3));
        println!("{}", String::from_utf8_lossy(&s4));

        let original_seqs = vec![s1,s2,s_middle,s_middle2,s3,s4];

        // Should get the same number of source nodes no matter the initial orientations
        for subset_mask in 0..2_u64.pow(original_seqs.len() as u32) {
            let mut seqs = original_seqs.clone();
            for i in 0..original_seqs.len() {
                if ((subset_mask >> i) & 1) > 0 {
                    jseqio::reverse_complement_in_place(&mut seqs[i]);
                }
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

            let dbg = crate::dbg::build_dbg(db, rc_db, k);
            dbg!(&dbg.edges);

            let orientations = new_algorithm(&dbg);
            for ori in orientations.iter() {
                dbg!(&ori);
            }

            let n_sources = evaluate(&orientations, &dbg);
            dbg!(n_sources);

            assert_eq!(n_sources, 2);
        }
        

    }    
}