use std::{
    collections::HashMap,
    fs::{read_to_string, File},
    io::{stdout, BufWriter, Write},
    str,
};

use anyhow::Context;
use clap::{App, Arg};
use itertools::Itertools;
use nanoserde::DeRon;
use rand::{rngs::SmallRng, Rng, SeedableRng};
use rand_distr::{Distribution, WeightedAliasIndex};
use seq_io::{
    fasta::{write_to as write_to_fa, Reader as faReader, Record as faRecord},
    fastq::{write_to, Reader as fqReader, Record as fqRecord},
};

#[derive(Clone, DeRon)]
struct Region {
    name: String,
    start: u64,
    end: u64,
}

trait Kmer {
    fn from_k_with_key(kmer: usize, key: String) -> Self;
    fn normalize(&mut self);
}

impl Kmer for HashMap<String, f32> {
    fn from_k_with_key(k: usize, key: String) -> Self {
        let nucs = ['A', 'C', 'T', 'G'];
        let mut map: HashMap<String, f32> = HashMap::new();
        if k == 1 {
            nucs.iter().for_each(|n| {
                map.insert(n.to_string(), 0f32);
            });
        } else {
            nucs.iter().for_each(|n| {
                let mut s = String::new();
                for _ in 0..k {
                    s.push(*n);
                }
                map.insert(s, 0f32);
            });
            nucs.iter().permutations(k).for_each(|per| {
                let mut s = String::new();
                per.iter().for_each(|p| {
                    s.push(**p);
                });
                map.insert(s, 0f32);
            });
        }

        if let Some(value) = map.get_mut(&key) {
            *value = 1f32;
        }
        map
    }

    fn normalize(&mut self) {
        let sum: f32 = self.values().sum();
        self.values_mut().for_each(|v| {
            *v /= sum;
        });
    }
}

struct Record {
    head: String,
    qual: Option<String>,
    len: usize,
}

trait Dist {
    fn score_seq(&self, seq: &[u8], k: usize) -> f32;
    fn top_seq<'a>(&self, seqs: &'a [(&[u8], Region)], k: usize) -> &(&'a [u8], Region);
}

impl Dist for Vec<HashMap<String, f32>> {
    #[inline]
    fn score_seq(&self, seq: &[u8], k: usize) -> f32 {
        seq.windows(k)
            .zip(self.iter())
            .fold(0f32, |acc, (elem, kmer)| {
                if let Some(value) =
                    kmer.get(str::from_utf8(elem).expect("sequence is not in utf8 format"))
                {
                    acc + value
                } else {
                    acc
                }
            })
    }
    #[inline]
    fn top_seq<'a>(&self, seqs: &'a [(&[u8], Region)], k: usize) -> &(&'a [u8], Region) {
        seqs.iter()
            .max_by(|x, y| {
                self.score_seq(x.0, k)
                    .partial_cmp(&self.score_seq(y.0, k))
                    .unwrap()
            })
            .unwrap()
    }
}

trait DistMap {
    fn update(&mut self, seq: &[u8]);
    fn normalize(&mut self);
}

fn update_distmap(dist_map: &mut HashMap<usize, Vec<HashMap<String, f32>>>, seq: &[u8], k: usize) {
    if let Some(dist) = dist_map.get_mut(&seq.len()) {
        seq.windows(k)
            .zip(dist.iter_mut())
            .for_each(|(elem, kmer)| {
                if let Some(value) =
                    kmer.get_mut(str::from_utf8(elem).expect("sequence is not in utf8 format"))
                {
                    *value += 1f32;
                };
            });
    } else {
        dist_map.insert(
            seq.len(),
            seq.windows(k)
                .map(|elem| {
                    Kmer::from_k_with_key(
                        k,
                        str::from_utf8(elem)
                            .expect("sequence is not in utf8 format")
                            .to_string(),
                    )
                })
                .collect(),
        );
    }
}

fn normalize_distmap(dist_map: &mut HashMap<usize, Vec<HashMap<String, f32>>>) {
    dist_map.values_mut().for_each(|value| {
        value.iter_mut().for_each(|n| {
            n.normalize();
        })
    })
}

enum Fastx<'a> {
    Fasta(&'a str),
    Fastq(&'a str),
}

impl<'a> Fastx<'a> {
    fn parse(
        &self,
        kmer: usize,
    ) -> anyhow::Result<(Vec<Record>, HashMap<usize, Vec<HashMap<String, f32>>>)> {
        let mut records: Vec<Record> = Vec::new();
        let mut ndist: HashMap<usize, Vec<HashMap<String, f32>>> = HashMap::new();

        match *self {
            Fastx::Fasta(fname) => {
                let mut reader = faReader::from_path(fname)
                    .with_context(|| format!("Failed to open file {}", fname))?;

                while let Some(record) = reader.next() {
                    let record = record.context("Invalid record")?;
                    let record_seq = record.full_seq();
                    if record_seq.len() != 0 {
                        records.push(Record {
                            head: String::from_utf8(record.head().to_owned())?,
                            qual: None,
                            len: record_seq.len(),
                        });

                        update_distmap(&mut ndist, &record_seq, kmer);
                    }
                }
                normalize_distmap(&mut ndist);

                Ok((records, ndist))
            }
            Fastx::Fastq(fname) => {
                let mut reader = fqReader::from_path(fname)
                    .with_context(|| format!("Failed to open file {}", fname))?;

                while let Some(record) = reader.next() {
                    let record = record.context("Invalid record")?;
                    let record_seq = record.seq();
                    if !record_seq.is_empty() {
                        records.push(Record {
                            head: String::from_utf8(record.head().to_owned())?,
                            qual: Some(String::from_utf8(record.qual().to_owned())?),
                            len: record_seq.len(),
                        });

                        update_distmap(&mut ndist, &record_seq, kmer);
                    }
                }
                normalize_distmap(&mut ndist);

                Ok((records, ndist))
            }
        }
    }

    fn sim(
        &self,
        genome_path: &str,
        regions_path: &str,
        bed_path: Option<&str>,
        seed: Option<&str>,
        kmer: usize,
    ) -> anyhow::Result<()> {
        let genome = parse_genome(genome_path)?;
        let (records, ndist) = self.parse(kmer)?;
        let regions_str = read_to_string(regions_path)?;
        let regions: Vec<Region> =
            DeRon::deserialize_ron(&regions_str).context("Failed to parse regions file")?;

        let stdout = stdout();

        let mut rng: SmallRng = match seed {
            Some(seed) => SeedableRng::seed_from_u64(seed.parse()?),
            None => SeedableRng::from_entropy(),
        };

        let mut bed_file = match bed_path {
            Some(bed_path) => {
                let bed_file = File::create(bed_path)?;
                Some(BufWriter::new(bed_file))
            }
            None => None,
        };

        let weights: Vec<u64> = regions
            .iter()
            .map(|region| region.end - region.start)
            .collect();

        let rng_dist = WeightedAliasIndex::new(weights).unwrap();

        for (index, record) in records.into_iter().enumerate() {
            if let Some(dist) = ndist.get(&record.len) {
                let mut n_reads: Vec<(&[u8], Region)> = Vec::new();
                let mut n_reads_rev: Vec<(Vec<u8>, Region)> = Vec::new();
                let d: Vec<(Vec<u8>, Region)>;
                let chr = regions
                    .get(rng_dist.sample(&mut rng))
                    .context("This should never fail")?;
                let genome_seq = genome.get(&chr.name).with_context(|| {
                    format!(
                        "Failed to get chromosome {} from reference genome",
                        chr.name
                    )
                })?;

                for _ in 0..10 {
                    let start_pos = rng
                        .gen_range(chr.start as usize + record.len..chr.end as usize - record.len);
                    let end_pos = start_pos + record.len;
                    let seq = genome_seq.get(start_pos..end_pos).with_context(|| {
                        format!(
                            "Failed to get sequence {}:{}:{} from reference genome",
                            chr.name, chr.start, chr.end
                        )
                    })?;

                    if index % 2 == 0 {
                        n_reads_rev.push((
                            seq.to_vec(),
                            Region {
                                name: chr.name.to_owned(),
                                start: start_pos as u64,
                                end: end_pos as u64,
                            },
                        ));
                    } else {
                        n_reads.push((
                            seq,
                            Region {
                                name: chr.name.to_owned(),
                                start: start_pos as u64,
                                end: end_pos as u64,
                            },
                        ));
                    }
                }

                let (seq, coord) = match index % 2 {
                    0 => {
                        d = n_reads_rev
                            .iter()
                            .map(|elem| (reverse_complement(elem.0.as_slice()), elem.1.to_owned()))
                            .collect();
                        n_reads = d
                            .iter()
                            .map(|elem| (elem.0.as_slice(), elem.1.to_owned()))
                            .collect();
                        dist.top_seq(&n_reads, kmer)
                    }
                    _ => dist.top_seq(&n_reads, kmer),
                };

                let strand = match index % 2 {
                    0 => "-",
                    _ => "+",
                };

                if let Some(ref mut bed_file) = bed_file {
                    bed_file.write_all(
                        format!(
                            "{}\t{}\t{}\t.\t0\t{}\n",
                            coord.name.as_str(),
                            coord.start,
                            coord.end,
                            strand
                        )
                        .as_bytes(),
                    )?;
                }

                match *self {
                    Fastx::Fastq(_) => {
                        write_to(
                            &stdout,
                            record.head.as_bytes(),
                            seq,
                            record.qual.unwrap().as_bytes(),
                        )?;
                    }
                    Fastx::Fasta(_) => {
                        write_to_fa(&stdout, record.head.as_bytes(), seq)?;
                    }
                }
            };
        }

        Ok(())
    }
}

fn parse_genome(filename: &str) -> anyhow::Result<HashMap<String, Vec<u8>>> {
    let mut reader = faReader::from_path(filename)
        .with_context(|| format!("Failed to open file {}", filename))?;

    let mut genome = HashMap::new();

    while let Some(record) = reader.next() {
        let record = record.context("Invalid record")?;
        let record_id = record
            .id()
            .with_context(|| format!("Failed to parse record id of record {:?}", &record))?;
        let record_seq = record.full_seq();

        genome.insert(record_id.to_owned(), record_seq.into_owned());
    }

    Ok(genome)
}

#[inline]
fn reverse_complement(input: &[u8]) -> Vec<u8> {
    input
        .iter()
        .rev()
        .map(|elem| match elem {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => *elem,
        })
        .collect()
}

fn main() -> anyhow::Result<()> {
    let matches = App::new("boquila")
        .version("0.4.1")
        .about("Generate NGS reads with same nucleotide distribution as input file\nGenerated reads will be written to stdout\nBy default input and output format is FASTQ")
        .arg(Arg::new("src").about("Model file").index(1).required(true))
        .arg(
            Arg::new("ref")
                .about("Reference FASTA")
                .long("ref")
                .takes_value(true)
                .value_name("FILE")
                .required(true),
        )
        .arg(
            Arg::new("fasta")
                .about("Change input and output format to FASTA")
                .long("fasta")
        )
        .arg(
            Arg::new("regions")
                .about(
                    "RON formatted file containing genomic regions that generated reads will be selected from",
                )
                .long("regions")
                .takes_value(true)
                .value_name("FILE")
                .required(true),
        ).arg(
            Arg::new("outbed")
                .about(
                    "File name in which the simulated reads will be saved in BED format",
                )
                .long("bed")
                .takes_value(true)
                .value_name("FILE")
        ).arg(
            Arg::new("seed")
                .about(
                    "Random number seed. If not provided system's default source of entropy will be used instead.",
                )
                .long("seed")
                .takes_value(true)
                .value_name("INT")
        ).arg(
            Arg::new("kmer")
                .about(
                    "Kmer size to be used while calculating frequency",
                )
                .long("kmer")
                .takes_value(true)
                .value_name("INT")
                .default_value("1")
        )
        .get_matches();

    let regions_file = matches.value_of("regions").unwrap();
    let reference_file = matches.value_of("ref").unwrap();
    let input_file = matches.value_of("src").unwrap();
    let bed_file = matches.value_of("outbed");
    let seed = matches.value_of("seed");
    let kmer: usize = matches.value_of_t("kmer")?;

    if matches.is_present("fasta") {
        let input = Fastx::Fasta(input_file);
        input.sim(reference_file, regions_file, bed_file, seed, kmer)
    } else {
        let input = Fastx::Fastq(input_file);
        input.sim(reference_file, regions_file, bed_file, seed, kmer)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn reverse_complement_works() {
        let mut vec = vec![b'A', b'C', b'G', b'T'];
        assert_eq!(
            vec![b'A', b'C', b'G', b'T'].as_slice(),
            reverse_complement(vec.as_slice())
        );
        vec.push(b'N');
        assert_eq!(
            vec![b'N', b'A', b'C', b'G', b'T'].as_slice(),
            reverse_complement(vec.as_slice())
        );
    }

    //TODO: tests for other functions
}
