use std::path::PathBuf;

use jseqio::reader::*;
use jseqio::writer::*;
use jseqio::record::*;
use clap::{Command, Arg};

use log::info;

use unitig_flipper::dbg::Orientation;
use unitig_flipper::optimize_unitig_orientation;
use unitig_flipper::SeqStream;

struct MyReader {
    inner: jseqio::reader::DynamicFastXReader,
}

impl<'a> SeqStream<'a> for MyReader{

    fn stream_next<'b>(&'b mut self) -> Option<&'b [u8]> where 'a : 'b {
        self.inner.read_next().unwrap().map(|rec| rec.seq)
    }
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

    let reader = MyReader{inner: reader};

    let orientations = optimize_unitig_orientation(reader, k);

    let n_forward = orientations.iter().fold(0_usize, |acc, &x| (acc + (x == Orientation::Forward) as usize));
    info!("{}% Forward", 100.0 * n_forward as f64 / orientations.len() as f64);

    info!("Writing output");

    let mut reader = DynamicFastXReader::from_file(infile).unwrap();
    let mut seq_idx = 0_usize;
    while let Some(rec) = reader.read_next().unwrap(){
        let orientation = orientations[seq_idx];
        let new_rec: OwnedRecord = match orientation{
            Orientation::Forward => rec.to_owned(),
            Orientation::Reverse => {
                let mut unitig = rec.to_owned();
                unitig.reverse_complement();
                unitig
            }
        };

        writer.write_owned_record(&new_rec).unwrap();
        seq_idx += 1;
    }    

}
