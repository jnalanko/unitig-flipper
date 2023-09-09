mod dbg;

use jseqio::reader::DynamicFastXReader;
use jseqio::record::*;

use dbg::Orientation;

fn pick_orientations(dbg: &dbg::DBG) -> Vec<Orientation>{
    let mut orientations = Vec::<Orientation>::new();
    orientations.resize(dbg.unitigs.sequence_count(), Orientation::Forward);

    let mut visited = vec![false; dbg.unitigs.sequence_count()];

    let mut stack = Vec::<(usize, Orientation)>::new(); // Reused DFS stack between iterations
    let mut n_components: usize = 0;    
    for component_root in 0..dbg.unitigs.sequence_count(){
        if visited[component_root]{
            continue;
        }

        n_components += 1;
        // Arbitrarily orient the root as forward        
        stack.push((component_root, Orientation::Forward));

        let mut component_size: usize = 0;
        // DFS from root and orient all reachable unitigs the same way
        while let Some((unitig_id, orientation)) = stack.pop(){
            if visited[unitig_id]{
                continue;
            }

            component_size += 1;
            visited[unitig_id] = true;
            orientations[unitig_id] = orientation;
    
            for edge in dbg.edges[unitig_id].iter(){
                let next_orientation = match (edge.from_orientation, edge.to_orientation){
                    (Orientation::Forward, Orientation::Forward) => orientation,
                    (Orientation::Forward, Orientation::Reverse) => orientation.flip(),
                    (Orientation::Reverse, Orientation::Forward) => orientation.flip(),
                    (Orientation::Reverse, Orientation::Reverse) => orientation,
                };
                stack.push((edge.to, next_orientation));
            }
        }
        eprintln!("Component size = {}", component_size);
    }

    eprintln!("Found {} component{}", n_components, match n_components > 1 {true => "s", false => ""});

    orientations
}

fn main() {
    let filename = std::env::args().nth(1).unwrap();
    let k = std::env::args().nth(2).unwrap().parse::<usize>().unwrap();
    let reader = DynamicFastXReader::from_file(&filename).unwrap();
    let filetype = reader.filetype();
    let (db, rc_db) = reader.into_db_with_revcomp().unwrap();
    let dbg = dbg::build_dbg(db, rc_db, k);
    let orientations = pick_orientations(&dbg);

    // Todo: gzip
    let mut writer = jseqio::writer::DynamicFastXWriter::new_to_stdout(filetype, false);
    for i in 0..dbg.unitigs.sequence_count(){
        let orientation = orientations[i];
        let rec: OwnedRecord = match orientation{
            Orientation::Forward => dbg.unitigs.get(i).unwrap().to_owned(),
            Orientation::Reverse => {
                let mut unitig = dbg.unitigs.get(i).unwrap().to_owned();
                unitig.reverse_complement();
                unitig
            }
        };
        writer.write(&rec);
    }

}
