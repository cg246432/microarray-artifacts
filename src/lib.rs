use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::io::{self, BufRead};
use rust_htslib::{bam::Read};
use rust_htslib::htslib;
use bio::data_structures::interval_tree::{IntervalTree};
use std::fmt;
use std::str;
use std::cmp::Ordering;

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>, {
        let file = File::open(filename)?;
        Ok(io::BufReader::new(file).lines())
}

struct BedInterval{
    chrom: i32,
    start: i64,
    stop: i64,
}

#[derive(Eq, Clone)]
pub struct ScRead{
    pub side: String, // Left or right
    pub mapped_pos: i64,  // Position the read mapped.
    pub sc_pos: i64,  // Position where softclipped
    pub tid: i32,     // chromosome
    pub l_qseq: usize, // Length of the read
    pub sc_len: usize, // Length of soft clipped portion
    pub i_len: usize, // Length of insertion (Right Side Only)
    pub d_len: usize, // Length of deletion (Right Side Only)
    pub seq: String,
}

impl fmt::Debug for ScRead {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("Soft Clipped Read")
         .field("tid", &self.tid)
         .field("mapped_pos", &self.mapped_pos)
         .field("side", &self.side)
         .finish()
    }
}

impl Ord for ScRead{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering{
        // Sort by chromosome, then by sc pos, then by read side
        if self.tid != other.tid{
            if self.tid < other.tid {
                return std::cmp::Ordering::Less;
            } else {
                return std::cmp::Ordering::Greater;
            }
        } else { // Same chromosome
            let retval = self.sc_pos - other.sc_pos;
            if retval < 0 {
                return std::cmp::Ordering::Less;
            } else if retval > 0{
                return std::cmp::Ordering::Greater;
            } else {
                if self.side != other.side {
                    if self.side == "Left" {
                        return std::cmp::Ordering::Less;
                    } else {
                        return std::cmp::Ordering::Greater;
                    }
                } else {
                    return std::cmp::Ordering::Equal;
                }
            }
        }
    }
}

impl PartialOrd for ScRead{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering>{
        Some(self.cmp(other))
    }
}

impl PartialEq for ScRead{
    fn eq(&self, other: &Self) -> bool {
        self.tid == other.tid && self.sc_pos == other.sc_pos && self.side == other.side
    }
}

pub struct ScCluster{
    pub clipped_consensus: String, 
    pub anchored_consensus: String,
    pub reads: Vec<ScRead>,
    pub side: String,
    pub tid: i32,
    pub pos: i64,
    pub score: f32,
}

impl ScCluster {
    fn new() -> ScCluster {
        ScCluster {
            clipped_consensus: "".to_string(),
            anchored_consensus: "".to_string(),
            reads: Vec::new(), // Len is proxy for size of cluster 
            side: "".to_string(),
            tid: 0,
            pos: 0,
            score: 0.0,
        }
    }
}

impl fmt::Debug for ScCluster {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("Soft Clipped Cluster")
         .field("tid", &self.tid)
         .field("pos", &self.pos)
         .field("side", &self.side)
         .field("score", &self.score)
         .field("supporting_reads", &self.reads.len())
         .field("clipped_consensus", &self.clipped_consensus)
         .finish()
    }
}

pub fn genomic_intervaltree_hash(bed_file: String, name_to_tid: &std::collections::HashMap<std::string::String, i32>) -> HashMap<i32, IntervalTree<i64, i64>>{
    let mut bed_hash = HashMap::<i32, IntervalTree<i64,i64>>::new(); 
    if let Ok(lines) = read_lines(bed_file) {
        for line in lines {
            if let Ok(ip) = line {
                if ip != "" {
                    let split_ip = ip.split("\t").collect::<Vec<_>>();
                    let bed_line = BedInterval{
                        // Fetch tid from bam header instead of string representation
                        chrom: *name_to_tid.get(&split_ip[0].to_string()).unwrap(),
                        start: split_ip[1].parse::<i64>().unwrap(),
                        stop: split_ip[2].parse::<i64>().unwrap(),
                    };
                    if !bed_hash.contains_key(&bed_line.chrom){
                        bed_hash.insert(bed_line.chrom, IntervalTree::new());
                    }
                    let v = bed_hash.get_mut(&bed_line.chrom).unwrap();
                    v.insert(bed_line.start..bed_line.stop, 1)
                }
            }
        }
    }
    return bed_hash;
}

fn subset_bam(chromosome: &i32, position: i64, tree: &HashMap<i32, IntervalTree<i64, i64>>) -> bool {
    if tree.is_empty(){
        return true; // No bed file supplied
    } else if tree.contains_key(chromosome) {
        return tree.get(chromosome).unwrap().find(position..position+1).next().is_some();
    } else {
        return false; 
    }

}

pub fn iterate_reads(bam: &mut rust_htslib::bam::IndexedReader, genome_interval_hash: HashMap<i32, IntervalTree<i64, i64>>, sc_len_threshold: i64) -> Vec<ScRead> {
    let mut sc_reads: Vec<ScRead> = Vec::new();
    for read in bam.rc_records()
        .map(|x| x.expect("Failure parsing Bam file"))
        // Rust Filter keeps if True
        // Filter Optical Duplicates
        .filter(|read|
            read.flags()
             & (htslib::BAM_FDUP) as u16
             == 0
        )
        // // Filter Quality Scores of 0
        .filter(|read|
            read.inner().core.qual != 0
        )
        // Filter out those softclipped reads on both sides
        .filter(|read|
            !(read.cigar().leading_softclips() != 0 && read.cigar().trailing_softclips() != 0)
        )
        // Filter if not softclipped or if not long enough 
        .filter(|read|
            read.cigar().leading_softclips() > sc_len_threshold|| read.cigar().trailing_softclips() > sc_len_threshold
        )
        // Filter if not in supplied bed file
        .filter(|read|
            subset_bam(&(read.tid()), read.pos(), &genome_interval_hash)
        )
        {   
            let mut side = "Left".to_string();
            let mut sc_len = read.cigar().leading_softclips() as usize;
            let mut sc_pos = read.pos();
            let mut i_len: usize = 0;
            let mut d_len: usize = 0;
            if read.cigar().leading_softclips() == 0{
                side = "Right".to_string();
                // Iterate over cigar to find insertions/deletions
                for cig in read.cigar().iter(){
                    match cig.char(){
                        'D' => d_len += cig.len() as usize,
                        'I' => i_len += cig.len() as usize,
                        _ => continue,
                    }
                }
                sc_len = read.cigar().trailing_softclips() as usize;
                sc_pos = read.pos() + ((read.seq_len() - sc_len - i_len + d_len) as i64);
            }
                let s = match str::from_utf8(&read.seq().as_bytes()){
                    Ok(v) => v.to_string(),
                    Err(e) => panic!("{}", e),
                };
                let sc_read = ScRead{
                    side: side,
                    mapped_pos: read.pos(),
                    sc_pos: sc_pos,
                    tid: read.tid(),
                    l_qseq: read.seq_len(),
                    sc_len: sc_len,
                    i_len: i_len,
                    d_len: d_len,
                    seq: s,
                };
                if sc_read.tid == 51{
                    println!("{:?}", sc_read);
                    std::process::exit(1);
                }
                sc_reads.push(sc_read);
            }
    sc_reads.sort_by(|l,r| l.cmp(r));
    return sc_reads
}

fn handle_cluster(cluster: &mut ScCluster){
    let mut max_sc_len: usize = 0;
    let mut max_an_len: usize = 0;
    let mut clipped_seqs: Vec<String> = Vec::new();
    let mut anchored_seqs: Vec<String> = Vec::new();
    for read in &cluster.reads{
        if read.sc_len > max_sc_len  {
            max_sc_len = read.sc_len;
        }
        if read.l_qseq - read.sc_len > max_an_len  {
            max_an_len = read.l_qseq - read.sc_len;
        }
        let mut clipped_seq: Vec<u8> = Vec::new();
        let mut anchor_seq: Vec<u8> = Vec::new();
        let mut junction = read.sc_len;
        if read.side == "Right" {
            junction = read.l_qseq - read.sc_len;
        }
        for base_index in 0..read.seq.len(){
            let sequence_char: u8 = read.seq.as_bytes()[base_index];
            if base_index <= junction {
                if read.side == "Left" {
                    clipped_seq.insert(base_index, sequence_char);
                } else {
                    anchor_seq.insert(base_index, sequence_char);
                }
            } else {
                let index = base_index - junction as usize - 1;
                if read.side == "Left" {                    
                    anchor_seq.insert(index, sequence_char);
                } else {
                    clipped_seq.insert(index, sequence_char);
                }
            }   
        }
        clipped_seqs.push(String::from_utf8(clipped_seq).ok().unwrap());
        anchored_seqs.push(String::from_utf8(anchor_seq).ok().unwrap());
    }
    let (clipped_consensus, clipped_score) = get_consensus(max_sc_len, clipped_seqs, &cluster.side);
    let (anchored_consensus, _anchored_score) = get_consensus(max_an_len, anchored_seqs, &cluster.side);
    //println!("{:?}\n{:?}", anchored_consensus, anchored_score);
    cluster.anchored_consensus = anchored_consensus;
    cluster.clipped_consensus = clipped_consensus;
    cluster.score = clipped_score;
}

fn get_consensus(max_len: usize, seqs: Vec<String>, align: &str) -> (String, f32){
    let mut consensus_score = 0.0;
    let mut consensus_seq = vec!["".to_string(); max_len];
    for i in 0..max_len - 1{ // Iterate over sequence position
        let (mut max, mut ac, mut cc, mut gc, mut tc, mut nc) = (0.0, 0.0, 0.0, 0.0, 0.0 ,0.0);
        let mut maxval = "N";
        let mut char_ind = i;
        let mut retval_ind = i;
        if align == "Right"{
            retval_ind = (max_len - i) - 1;
        }
        for seq in seqs.iter(){ // Iterate over each sequence 
            if align == "Right"{
                if seq.len() <= i  {
                    continue;
                } 
                char_ind = (seq.len() - i) - 1;
            } else {
                if i >= seq.len(){
                    continue;
                }
            }

            let base_option = seq.to_string().chars().nth(char_ind as usize);
            match base_option {
                Some('A') => {
                    ac += 1.0;
                    if ac > max {
                        max = ac;
                        maxval = "A";
                    }
                },
                Some('C') => {
                    cc += 1.0;
                    if cc > max {
                        max = cc;
                        maxval = "C";
                    }
                },
                Some('G') => {
                    gc += 1.0;
                    if gc > max {
                        max = gc;
                        maxval = "G";
                    }
                },
                Some('T') => {
                    tc += 1.0;
                    if tc > max {
                        max = tc;
                        maxval = "T";
                    }
                },
                Some('N') => {
                    nc += 1.0;
                    if nc > max {
                        max = nc;
                        maxval = "N";
                    }
                },
                _ => continue,
            };
        }

        let pos_seq_count = ac + tc + gc + nc + cc;
        let pos_score = max / pos_seq_count;
        consensus_score = consensus_score + pos_score as f32;
        consensus_seq[retval_ind] = maxval.to_string();
    }
    consensus_score = consensus_score / max_len as f32;
    return (consensus_seq.join(""), consensus_score, );
}

pub fn generate_clusters(sc_reads: &Vec<ScRead>, min_cluster_size: i64) -> Vec<ScCluster> {
    let mut sc_clusters: Vec<ScCluster> = Vec::new();
    let mut current_position = 0;
    let mut current_side = "";
    let mut current_tid = 0;
    let mut cluster = ScCluster::new(); // Current cluster
    for i in 0..sc_reads.len() {
        let current_read = &sc_reads[i];
        // New cluster, not initial instantiated cluster
        if (current_read.sc_pos != current_position || current_read.side != current_side || current_read.tid != current_tid) && cluster.reads.len() > 0 {
            if cluster.reads.len() as i64 > min_cluster_size {
                // Get details of cluster (consensus, score, etc.)
                cluster.side = current_side.to_string();
                cluster.tid = current_tid;
                cluster.pos = current_position;
                handle_cluster(&mut cluster);
                // Add to array of output clusters
                sc_clusters.push(cluster);
            }
            // Revert defaults to current state
            current_position = current_read.sc_pos;
            current_side = &current_read.side;
            current_tid = current_read.tid;
            // Reset cluster
            cluster = ScCluster::new();

        } else { // New member of cluster
            cluster.reads.push(current_read.clone());
        }
    }
    // Handle last read and last cluster (For loop exhaustion)
    if cluster.reads.len() as i64 > min_cluster_size {
        // Get details of cluster (consensus, score, etc.)
        cluster.side = current_side.to_string();
        cluster.tid = current_tid;
        handle_cluster(&mut cluster);
        // Add to array of output clusters
        sc_clusters.push(cluster);
    } 
    return sc_clusters;
}