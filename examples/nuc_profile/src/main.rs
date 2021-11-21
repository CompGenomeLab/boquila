use std::{collections::HashMap, str};

use anyhow::Context;
use clap::{App, Arg};
use itertools::Itertools;
use seq_io::{
    fasta::Reader as faReader,
    fastq::{Reader as fqReader, Record},
};

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
                .map(|elem| Kmer::from_k_with_key(k, str::from_utf8(elem).unwrap().to_string()))
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
    fn parse(&self, kmer: usize) -> anyhow::Result<HashMap<usize, Vec<HashMap<String, f32>>>> {
        let mut ndist: HashMap<usize, Vec<HashMap<String, f32>>> = HashMap::new();

        match *self {
            Fastx::Fasta(fname) => {
                let mut reader = faReader::from_path(fname)
                    .with_context(|| format!("Failed to open file {}", fname))?;

                while let Some(record) = reader.next() {
                    let record = record.context("Invalid record")?;
                    let record_seq = record.full_seq();
                    if record_seq.len() != 0 {
                        update_distmap(&mut ndist, &record_seq, kmer);
                    }
                }

                normalize_distmap(&mut ndist);

                Ok(ndist)
            }
            Fastx::Fastq(fname) => {
                let mut reader = fqReader::from_path(fname)
                    .with_context(|| format!("Failed to open file {}", fname))?;

                while let Some(record) = reader.next() {
                    let record = record.context("Invalid record")?;
                    let record_seq = record.seq();
                    if !record_seq.is_empty() {
                        update_distmap(&mut ndist, &record_seq, kmer);
                    }
                }

                normalize_distmap(&mut ndist);

                Ok(ndist)
            }
        }
    }

    fn profile(&self, kmer: usize, oligo_len: usize) -> anyhow::Result<()> {
        let dist = self.parse(kmer)?;
        let profile = dist
            .get(&oligo_len)
            .with_context(|| format!("Reads file don't have any reads of length {}", oligo_len))?;
        println!("A\tC\tT\tG");
        for position in profile {
            println!(
                "{}\t{}\t{}\t{}",
                position["A"], position["C"], position["T"], position["G"]
            );
        }
        Ok(())
    }
}

fn main() -> anyhow::Result<()> {
    let matches = App::new("nuc_profile")
        .version("0.1.0")
        .about("Generate nucleotide profile from FASTA or FASTQ")
        .arg(
            Arg::new("reads")
                .about("File containing input reads")
                .long("reads")
                .short('R')
                .takes_value(true)
                .value_name("FILE")
                .required(true),
        )
        .arg(
            Arg::new("format")
                .about("Input file format")
                .long("format")
                .short('F')
                .takes_value(true)
                .possible_values(&["fastq", "fasta"])
                .default_value("fastq"),
        )
        .arg(
            Arg::new("kmer")
                .about("Kmer size to be used while calculating frequency")
                .long("kmer")
                .takes_value(true)
                .value_name("INT")
                .default_value("1"),
        )
        .arg(
            Arg::new("oligo_len")
                .about("Oligomer length which nucleotide frequency is calculated for")
                .long("len")
                .takes_value(true)
                .value_name("INT")
                .default_value("26"),
        )
        .get_matches();

    let reads = matches.value_of("reads").unwrap();
    let format = matches.value_of("format").unwrap();
    let kmer: usize = matches.value_of_t("kmer")?;
    let oligo_len: usize = matches.value_of_t("oligo_len")?;

    match format {
        "fasta" => {
            let fa = Fastx::Fasta(reads);
            fa.profile(kmer, oligo_len)
        }
        "fastq" => {
            let fq = Fastx::Fastq(reads);
            fq.profile(kmer, oligo_len)
        }
        _ => unreachable!(),
    }
}
