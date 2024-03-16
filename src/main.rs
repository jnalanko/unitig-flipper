use std::path::PathBuf;

use jseqio::reader::*;
use jseqio::writer::*;
use jseqio::record::*;
use clap::{Command, Arg};

use log::info;

use unitig_flipper::pick_orientations;
use unitig_flipper::dbg;
use unitig_flipper::dbg::Orientation;
use jseqio::seq_db::SeqDB;

fn run(forward_seqs: SeqDB, reverse_seqs: SeqDB, seqs_out: &mut impl SeqRecordWriter, k: usize){

    info!("Building bidirected DBG edges...");
    let dbg = dbg::build_dbg(forward_seqs, reverse_seqs, k);

    info!("Choosing unitig orientations");
    let orientations = pick_orientations(&dbg);

    info!("Evaluating the solution");
    let n_with_predecessor = evaluate(&orientations, &dbg);

    info!("{}/{} unitigs have a predecessor ({:.2}%)", n_with_predecessor, dbg.unitigs.sequence_count(), 100.0 * n_with_predecessor as f64 / dbg.unitigs.sequence_count() as f64);

    info!("Writing output");

    #[allow(clippy::needless_range_loop)] // We need the index for get(i) also
    for i in 0..dbg.unitigs.sequence_count(){
        let orientation = orientations[i];
        let rec: OwnedRecord = match orientation{
            Orientation::Forward => dbg.unitigs.get(i).to_owned(),
            Orientation::Reverse => {
                let mut unitig = dbg.unitigs.get(i).to_owned();
                unitig.reverse_complement();
                unitig
            }
        };

        seqs_out.write_owned_record(&rec).unwrap();
    }    
}


// Returns the number of unitigs that do not have a predecessor
fn evaluate(choices: &[Orientation], dbg: &dbg::DBG) -> usize{
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


fn main() {

    if std::env::var("RUST_LOG").is_err(){
        std::env::set_var("RUST_LOG", "info");
    }

    env_logger::init();

    let cli = Command::new("unitig-flipper")
        .about("Orients unitigs heuristically to minimize the number of dummy nodes in the SBWT graph.")
        .author("Jarno N. Alanko <alanko.jarno@gmail.com>")
        .arg(Arg::new("input")
            .help("Input FASTA or FASTQ file, possibly gzipped")
            .long("input")
            .short('i')
            .required(true)
            .value_parser(clap::value_parser!(PathBuf))
        )
        .arg(Arg::new("output")
            .help("Output FASTA or FASTQ file, possibly gzipped")
            .long("output")
            .short('o')
            .required(true)
            .value_parser(clap::value_parser!(PathBuf))
        )
        .arg(Arg::new("k")
            .help("k-mer length")
            .short('k')
            .required(true)
            .value_parser(clap::value_parser!(usize))
        );

    let cli_matches = cli.get_matches();
    let infile: &PathBuf = cli_matches.get_one("input").unwrap();
    let outfile: &PathBuf = cli_matches.get_one("output").unwrap();
    let k: usize = *cli_matches.get_one("k").unwrap();

    let reader = DynamicFastXReader::from_file(infile).unwrap();
    let mut writer = DynamicFastXWriter::new_to_file(outfile).unwrap();

    info!("Reading and reverse-complementing sequences into memory");
    let (db, rc_db) = reader.into_db_with_revcomp().unwrap();
    run(db, rc_db, &mut writer, k);

}

#[cfg(test)]
mod tests{
    use super::*;
    use std::io::BufReader;

    fn helper_get_seq_dbs<'a>(seqs: impl IntoIterator<Item = &'a [u8]>) -> (SeqDB, SeqDB){
        let mut fw_db = SeqDB::new();
        let mut rc_db = SeqDB::new();

        for seq in seqs{
            let rc = jseqio::reverse_complement(seq);
            fw_db.push_record(RefRecord{head: b"", seq, qual: None});
            rc_db.push_record(RefRecord{head: b"", seq: rc.as_slice(), qual: None});
        }

        (fw_db, rc_db)
    }

    #[test]
    fn straight_line(){
        // Input
        let k = 3;
        let data = [b"TCG", b"ATC", b"ATG", b"ACC", b"CCG"];

        // Two possible answers
        let ans1 = [b"CGA", b"GAT", b"ATG", b"ACC", b"CCG"];
        let ans2: Vec<Vec<u8>> = ans1.iter().map(|s| jseqio::reverse_complement(s.as_slice())).collect();

        // Run the test

        let seq_iter = data.iter().map(|s| s.as_slice());
    
        let (fw_db, rc_db) = helper_get_seq_dbs(seq_iter);

        let out_buf = Vec::<u8>::new();
        let mut writer = FastXWriter::new(out_buf, jseqio::FileType::FASTA);

        run(fw_db, rc_db, &mut writer, k);
        let out_buf = writer.into_inner().unwrap(); // Get back the out buffer

        let br = BufReader::new(std::io::Cursor::new(out_buf));
        let reader = DynamicFastXReader::new(br).unwrap();
        let out_db = reader.into_db().unwrap();
        assert_eq!(data.len(), out_db.sequence_count());

        let ans1_match = (0..ans1.len()).all(|i| ans1[i] == out_db.get(i).seq);
        let ans2_match = (0..ans2.len()).all(|i| ans2[i] == out_db.get(i).seq);

        assert!(ans1_match || ans2_match);

    }
}