mod dbg;

use jseqio::reader::DynamicFastXReader;
use jseqio::record::*;
use clap::{Command, Arg};

use dbg::Orientation;
use dbg::Orientation::{Forward, Reverse};

fn is_terminal(dbg: &dbg::DBG, unitig_id: usize) -> bool{
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

fn dfs_from(root: usize, dbg: &dbg::DBG, visited: &mut Vec<bool>, orientations: &mut Vec<Orientation>){
    let mut stack = Vec::<(usize, Orientation)>::new();

    // Arbitrarily orient the root as forward
    stack.push((root, Orientation::Forward));

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
            if edge.from_orientation != orientation{
                // If we came to this node in the forward orientation, we may only
                // leave on edges that leave in the forward orientation. And vice versa.
                continue;
            }

            match (edge.from_orientation, edge.to_orientation, orientation){
                 // Edge leaves from the forward end of the current unitig
                 (Forward, Forward, _) => stack.push((edge.to, orientation)),
                 (Forward, Reverse, _) => stack.push((edge.to, orientation.flip())),
                 // Edge leaves from the reverse end of the current unitig
                 (Reverse, Forward, _) => stack.push((edge.to, orientation.flip())),
                 (Reverse, Reverse, _) => stack.push((edge.to, orientation)),
             };
        }
    }
    eprintln!("Component size = {}", component_size);
}

fn pick_orientations(dbg: &dbg::DBG) -> Vec<Orientation>{
    let mut orientations = Vec::<Orientation>::new();
    orientations.resize(dbg.unitigs.sequence_count(), Orientation::Forward);

    let mut visited = vec![false; dbg.unitigs.sequence_count()];

    
    let mut n_components: usize = 0;    

    // DFS from terminal nodes
    for component_root in 0..dbg.unitigs.sequence_count(){

        if !is_terminal(&dbg, component_root) {
            // This node will be visited from some other node
            continue;
        }

        if visited[component_root]{
            continue;
        }

        dfs_from(component_root, &dbg, &mut visited, &mut orientations);
        n_components += 1;
    }

    let n_visited = visited.iter().fold(0_usize, |sum, &x| sum + x as usize);
    eprintln!("{:.2}% of unitigs visited so far... proceeding to clean up cycles", 100.0 * n_visited as f64 / dbg.unitigs.sequence_count() as f64);

    // Only cycles remain. DFS from the rest.
    for component_root in 0..dbg.unitigs.sequence_count(){

        if visited[component_root]{
            continue;
        }

        dfs_from(component_root, &dbg, &mut visited, &mut orientations);
        n_components += 1;
    }

    eprintln!("Done. Total {} DFS iterations executed", n_components);

    orientations
}

fn main() {

    let cli = Command::new("unitig-flipper")
        .about("Orients unitigs heuristically to minimize the number of dummy nodes in the SBWT graph.")
        .author("Jarno N. Alanko <alanko.jarno@gmail.com>")
        .arg(Arg::new("input")
            .help("Input FASTA or FASTQ file, possibly gzipped")
            .long("input")
            .short('i')
            .required(true)
            .value_parser(clap::value_parser!(std::path::PathBuf))
        )
        .arg(Arg::new("output")
            .help("Output FASTA or FASTQ file, possibly gzipped")
            .long("output")
            .short('o')
            .required(true)
            .value_parser(clap::value_parser!(std::path::PathBuf))
        )
        .arg(Arg::new("k")
            .help("k-mer length")
            .short('k')
            .required(true)
            .value_parser(clap::value_parser!(usize))
        );

    let cli_matches = cli.get_matches();
    let infile: &std::path::PathBuf = cli_matches.get_one("input").unwrap();
    let outfile: &std::path::PathBuf = cli_matches.get_one("output").unwrap();
    let k: usize = *cli_matches.get_one("k").unwrap();

    let reader = DynamicFastXReader::from_file(infile).unwrap();

    // Let's also open the writer right away so it errors out if the file cannot be opened
    let mut writer = jseqio::writer::DynamicFastXWriter::new_to_file(outfile).unwrap();

    eprintln!("Reading and reverse-complementing sequences from {} into memory", infile.display());
    let (db, rc_db) = reader.into_db_with_revcomp().unwrap();

    eprintln!("Building bidirected DBG edges");
    let dbg = dbg::build_dbg(db, rc_db, k);

    eprintln!("Choosing unitig orientations");
    let orientations = pick_orientations(&dbg);

    eprintln!("Writing output to {}", outfile.display());

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
