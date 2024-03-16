pub mod dbg;

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
            if edge.from_orientation != orientation{
                // If we came to this node in the forward orientation, we may only
                // leave on edges that leave in the forward orientation. And vice versa.
                continue;
            }

            match (edge.from_orientation, edge.to_orientation, orientation){
                 // Edge leaves from the forward end of the current unitig
                 (Forward, Forward, _) => queue.push_back((edge.to, orientation)),
                 (Forward, Reverse, _) => queue.push_back((edge.to, orientation.flip())),
                 // Edge leaves from the reverse end of the current unitig
                 (Reverse, Forward, _) => queue.push_back((edge.to, orientation.flip())),
                 (Reverse, Reverse, _) => queue.push_back((edge.to, orientation)),
             };
        }
    }
}
