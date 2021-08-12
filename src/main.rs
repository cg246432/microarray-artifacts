use rust_bioutils;
use rust_htslib::{bam, bam::Read};
use std::str;
use std::collections::HashMap;
use clap::{App, Arg};

fn main(){
    // Parse CLI
    let matches = App::new("Cluster Identifier")
    .about("Identifying short clipped clusters from alignment files.")
    .author("Chris Greco, chris.greco11@gmail.com")
    .setting(clap::AppSettings::TrailingVarArg)
    .arg(Arg::with_name("alignment_file")
        .help("BAM/CRAM/SAM File")
        .required(true)
        .takes_value(true)
    )
    .arg(Arg::with_name("min_sc_bases")
        .help("Minimum number of bases in Soft Clipped Read")
        .default_value("10")
        .takes_value(true)
        .required(true)
    )   
    .arg(Arg::with_name("min_sc_reads")
        .help("Minimum number of supporting reads to make a cluster")
        .default_value("5")
        .takes_value(true)
        .required(true)
    )
    .arg(Arg::with_name("score")
        .help("Simple Consensus score of a Cluster to be output.")
        .default_value("0.9")
        .takes_value(true)
        .required(true)
    )
    .arg(Arg::with_name("threads")
        .help("Threads to use.")
        .default_value("1")
        .takes_value(true)
        .required(false)
    )
    .arg(Arg::with_name("region")
        .help("A region to subset the cluster search to. Format: chr:start-stop.")
        .takes_value(true)
        .required(false)
    )
    .arg(Arg::with_name("bed_file")
        .help("A bed file that defines regions to subset the cluster search to. Used for multiple regions.")
        .required(false)
        .takes_value(true)
    )
    .get_matches();
    let bam_file = matches.value_of("alignment_file").unwrap();
    let bed_file = matches.value_of("bed_file").unwrap_or("");
    let min_sc_bases: i64 = matches.value_of("min_sc_bases").unwrap().parse().unwrap();
    let min_sc_reads: i64 = matches.value_of("min_sc_reads").unwrap().parse().unwrap();
    let score_cutoff: f32 = matches.value_of("score").unwrap().parse().unwrap(); 
    //let threads: usize = matches.value_of("threads").unwrap().parse().unwrap();
    let region = matches.value_of("region").unwrap_or(".");
    let mut bam = bam::IndexedReader::from_path(bam_file).unwrap();
    let mut name_to_tid = HashMap::<String, i32>::new(); 
    let mut tid_to_name = HashMap::<i32, String>::new();

    for name in bam.header().target_names(){
        let value = bam.header().tid(name).unwrap() as i32;
        let key = str::from_utf8(name).unwrap();
        name_to_tid.insert(key.to_string(), value);
        tid_to_name.insert(value, key.to_string());
    }
    let genome_interval_hash = rust_bioutils::genomic_intervaltree_hash(bed_file.to_string(), &name_to_tid);
    bam.fetch(region).unwrap();
    // Find soft clipped reads
    let sc_reads = rust_bioutils::iterate_reads(&mut bam, genome_interval_hash, min_sc_bases);
    // Find soft clipped clusters
    if sc_reads.len() > 0{
        let sc_clusters = rust_bioutils::generate_clusters(&sc_reads, min_sc_reads);
        for cluster in &sc_clusters{
            if cluster.score >= score_cutoff{
                println!("{:?}", tid_to_name.get(&cluster.tid).unwrap());
                println!("{:?}", cluster);
            }
        }
    }
}