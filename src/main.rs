mod dbg;

use jseqio::reader::DynamicFastXReader;
use jseqio::record::*;
use clap::{Command, Arg};

use dbg::Orientation;

fn pick_orientations(dbg: &dbg::DBG) -> Vec<Orientation>{
    let mut orientations = Vec::<Orientation>::new();
    orientations.resize(dbg.unitigs.sequence_count(), Orientation::Forward);

    let mut visited = vec![false; dbg.unitigs.sequence_count()];

    let mut stack = Vec::<(usize, Orientation)>::new(); // Reused DFS stack between iterations
    let mut n_components: usize = 0;    
    for component_root in 0..dbg.unitigs.sequence_count(){
        let is_source = dbg.edges[component_root].iter().all(|edge| edge.from_orientation == Orientation::Forward);
        let is_self_loop = dbg.edges[component_root].iter().all(|edge| edge.from == edge.to);

        if !is_source && !is_self_loop {
            // This node will be visited from some other node
            continue;
        }

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
               match (edge.from_orientation, edge.to_orientation){
                    (Orientation::Forward, Orientation::Forward) => stack.push((edge.to, orientation)),
                    (Orientation::Forward, Orientation::Reverse) => stack.push((edge.to, orientation.flip())),
                    (Orientation::Reverse, _) => () // This edge would move backwards to a predecessor, so we don't follow it
                };
            }
        }
        eprintln!("Component size = {}", component_size);
    }

    eprintln!("Found {} component{}", n_components, match n_components > 1 {true => "s", false => ""});

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
