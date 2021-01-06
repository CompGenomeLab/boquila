use std::{collections::HashMap, fs::read_to_string, io::stdout};

use anyhow::Context;
use clap::{App, Arg};
use nanorand::{WyRand, RNG};
use nanoserde::DeRon;
use seq_io::{
    fasta::{write_to as write_to_fa, Reader as faReader, Record as faRecord},
    fastq::{write_to, Reader as fqReader, Record as fqRecord},
};

#[derive(DeRon)]
struct Region {
    name: String,
    start: u64,
    end: u64,
}

#[derive(Clone, Default)]
struct Nuc {
    a: f64,
    c: f64,
    g: f64,
    t: f64,
}

impl Nuc {
    fn sum(&self) -> f64 {
        self.a + self.c + self.g + self.t
    }

    fn normalize(&mut self) {
        let sum = self.sum();
        self.a /= sum;
        self.c /= sum;
        self.g /= sum;
        self.t /= sum;
    }
}

struct Record {
    head: String,
    qual: Option<String>,
    len: usize,
}

trait Dist {
    fn score_seq(&self, seq: &[u8]) -> f64;
    fn top_seq<'a>(&self, seqs: &'a [&[u8]]) -> &'a [u8];
}

impl Dist for Vec<Nuc> {
    #[inline]
    fn score_seq(&self, seq: &[u8]) -> f64 {
        seq.iter()
            .enumerate()
            .fold(0f64, |acc, (i, elem)| match elem {
                b'A' => acc + self[i].a,
                b'C' => acc + self[i].c,
                b'G' => acc + self[i].g,
                b'T' => acc + self[i].t,
                _ => acc,
            })
    }

    #[inline]
    fn top_seq<'a>(&self, seqs: &'a [&[u8]]) -> &'a [u8] {
        seqs.iter()
            .max_by(|x, y| self.score_seq(x).partial_cmp(&self.score_seq(y)).unwrap())
            .unwrap()
            .to_owned()
    }
}

trait DistMap {
    fn update(&mut self, seq: &[u8]);
    fn normalize(&mut self);
}

impl DistMap for HashMap<usize, Vec<Nuc>> {
    fn update(&mut self, seq: &[u8]) {
        if let Some(dist) = self.get_mut(&seq.len()) {
            seq.iter().enumerate().for_each(|(i, elem)| match elem {
                b'A' => dist[i].a += 1f64,
                b'C' => dist[i].c += 1f64,
                b'G' => dist[i].g += 1f64,
                b'T' => dist[i].t += 1f64,
                _ => (),
            })
        } else {
            self.insert(
                seq.len(),
                seq.iter()
                    .map(|elem| match elem {
                        b'A' => Nuc {
                            a: 1f64,
                            ..Default::default()
                        },
                        b'C' => Nuc {
                            c: 1f64,
                            ..Default::default()
                        },
                        b'G' => Nuc {
                            g: 1f64,
                            ..Default::default()
                        },
                        b'T' => Nuc {
                            t: 1f64,
                            ..Default::default()
                        },
                        _ => Nuc {
                            ..Default::default()
                        },
                    })
                    .collect(),
            );
        }
    }

    fn normalize(&mut self) {
        self.values_mut().for_each(|value| {
            value.iter_mut().for_each(|n| {
                n.normalize();
            })
        })
    }
}

enum Fastx<'a> {
    Fasta(&'a str),
    Fastq(&'a str),
}

impl<'a> Fastx<'a> {
    fn parse(&self) -> anyhow::Result<(Vec<Record>, HashMap<usize, Vec<Nuc>>)> {
        let mut records: Vec<Record> = Vec::new();
        let mut ndist: HashMap<usize, Vec<Nuc>> = HashMap::new();

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

                        ndist.update(&record_seq);
                    }
                }
                ndist.normalize();

                Ok((records, ndist))
            }
            Fastx::Fastq(fname) => {
                let mut reader = fqReader::from_path(fname)
                    .with_context(|| format!("Failed to open file {}", fname))?;

                while let Some(record) = reader.next() {
                    let record = record.context("Invalid record")?;
                    let record_seq = record.seq();
                    if record_seq.len() != 0 {
                        records.push(Record {
                            head: String::from_utf8(record.head().to_owned())?,
                            qual: Some(String::from_utf8(record.qual().to_owned())?),
                            len: record_seq.len(),
                        });

                        ndist.update(record_seq);
                    }
                }
                ndist.normalize();

                Ok((records, ndist))
            }
        }
    }

    fn sim(&self, genome_path: &str, regions_path: &str) -> anyhow::Result<()> {
        let genome = parse_genome(genome_path)?;
        let (records, ndist) = self.parse()?;
        let regions_str = read_to_string(regions_path)?;
        let regions: Vec<Region> =
            DeRon::deserialize_ron(&regions_str).context("Failed to parse regions file")?;
        let mut rng = WyRand::new();
        let stdout = stdout();

        for record in records {
            if let Some(dist) = ndist.get(&record.len) {
                let mut n_reads: Vec<&[u8]> = Vec::new();
                let chr_index = rng.generate_range(0, regions.len());
                let chr = regions.get(chr_index).context("This should never fail")?;
                let genome_seq = genome.get(&chr.name).with_context(|| {
                    format!(
                        "Failed to get chromosome {} from reference genome",
                        chr.name
                    )
                })?;

                for _ in 0..10 {
                    let start_pos = rng.generate_range(
                        chr.start as usize + record.len,
                        chr.end as usize - record.len,
                    );
                    let end_pos = start_pos + record.len;
                    let seq = genome_seq.get(start_pos..end_pos).with_context(|| {
                        format!(
                            "Failed to get sequence {}:{}:{} from reference genome",
                            chr.name, chr.start, chr.end
                        )
                    })?;
                    n_reads.push(seq);
                }

                let seq = dist.top_seq(&n_reads);

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

fn main() -> anyhow::Result<()> {
    let matches = App::new("boquila")
        .version("0.2.0")
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
        )
        .get_matches();

    let regions_file = matches.value_of("regions").unwrap();
    let reference_file = matches.value_of("ref").unwrap();
    let input_file = matches.value_of("src").unwrap();

    if matches.is_present("fasta") {
        let input = Fastx::Fasta(input_file);
        input.sim(reference_file, regions_file)
    } else {
        let input = Fastx::Fastq(input_file);
        input.sim(reference_file, regions_file)
    }
}
