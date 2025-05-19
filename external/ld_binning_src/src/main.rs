#![feature(portable_simd)]
use anyhow::{Context, Result, bail};
use clap::{Parser, command};
use indicatif::ProgressBar;
use rand::distr::Uniform;
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;
use rust_htslib::bcf::header::HeaderRecord;
use rust_htslib::bcf::{IndexedReader, Read, Record, record};
use std::error::Error;

// Parse the genotypes and compute the minor allele frequency
fn precompute_standarized_genotypes(
    record: Record,
    parameters: &Cli,
    genotypes: &mut [f64],
) -> Result<Option<()>, Box<dyn Error>> {
    // Get the genotype field
    let n_samples = genotypes.len();
    let buffer = record::Buffer::new();
    let raw_genotypes = record
        .genotypes_shared_buffer(buffer)
        .context("Error getting genotypes")?;

    let mut total = 0;
    for (i, val) in genotypes.iter_mut().enumerate().take(n_samples) {
        *val = 0.0;
        let sample = raw_genotypes.get(i);
        for j in 0..2 {
            if let Some(gt) = sample[j].index() {
                *val += gt as f64;
                total += gt;
            } else {
                return Err("Missing genotype".into());
            }
        }
    }
    let allele_freq = total as f64 / (2 * n_samples) as f64;
    let maf = if allele_freq > 0.5 {
        1.0 - allele_freq
    } else {
        allele_freq
    };
    // If the MAF is less than the threshold, skip this record
    if maf < parameters.maf_threshold {
        return Ok(None);
    }
    // Compute the standardized genotypes
    for val in genotypes.iter_mut() {
        *val = (*val - 2.0 * allele_freq) / (2.0 * allele_freq * (1.0 - allele_freq)).sqrt();
    }
    Ok(Some(()))
}

pub fn linkage_disequilibrium(genotypes1: &[f64], genotypes2: &[f64], n_samples: usize) -> f64 {
    assert!(
        genotypes1.len() >= n_samples && genotypes2.len() >= n_samples,
        "Input length mismatch"
    );

    let s = n_samples as f64;
    let mut ld = 0.0;
    let mut ld_square = 0.0;

    for i in 0..n_samples {
        let a = genotypes1[i];
        let b = genotypes2[i];
        let prod = a * b;
        ld += prod;
        ld_square += prod * prod;
    }

    (ld * ld - ld_square) / (s * (s - 1.0))
}

struct Bins {
    nbins: usize,
    left_edges_in_cm: Vec<f64>,
    right_edges_in_cm: Vec<f64>,
    left_edges_in_bp: Vec<f64>,
    right_edges_in_bp: Vec<f64>,
    minimum: i64,
    maximum: i64,
}

impl Bins {
    // From HapNe supplementary material
    fn hapne_default(recombination_rate: f64) -> Self {
        let nbins = 19;
        let mut left_edges_in_cm = Vec::with_capacity(nbins);
        let mut right_edges_in_cm = Vec::with_capacity(nbins);

        for i in 0..nbins {
            let i = i as f64;
            left_edges_in_cm.push(0.5 + 0.5 * i);
            right_edges_in_cm.push(1.0 + 0.5 * i);
        }
        // Transform to base pairs using x / 100 / recombination_rate
        let left_edges_in_bp = left_edges_in_cm
            .iter()
            .map(|&x| x / 100.0 / recombination_rate)
            .collect::<Vec<f64>>();
        let right_edges_in_bp = right_edges_in_cm
            .iter()
            .map(|&x| x / 100.0 / recombination_rate)
            .collect::<Vec<f64>>();
        let minimum = left_edges_in_bp[0].round() as i64;
        let maximum = right_edges_in_bp[nbins - 1].round() as i64;
        Self {
            nbins,
            left_edges_in_cm,
            right_edges_in_cm,
            left_edges_in_bp,
            right_edges_in_bp,
            minimum,
            maximum,
        }
    }
}

#[derive(Debug)]
struct SufficientSummaryStats {
    pub counts: Vec<u32>,
    pub ld: Vec<f64>,
    pub ld_square: Vec<f64>,
}
impl SufficientSummaryStats {
    pub fn new(bins: &Bins) -> Self {
        Self {
            counts: vec![0; bins.nbins],
            ld: vec![0.0; bins.nbins],
            ld_square: vec![0.0; bins.nbins],
        }
    }
    pub fn should_stop(&self, args: &Cli) -> bool {
        // First, check that all bins have at least the minimum number of samples
        for count in self.counts.iter() {
            if *count < args.min_loci as u32 {
                return false;
            }
        }
        for i in 0..self.ld.len() {
            // Compute sample standard deviation
            let std = (self.ld_square[i] / (self.counts[i] - 1) as f64).sqrt();
            // Compute the CI half-width
            let ci_half_width = std / (self.counts[i] as f64).sqrt();
            if ci_half_width * 1.96 > args.epsilon {
                return false;
            }
        }
        true
    }
}

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Compressed and indexed VCF or BCF file
    infile: String,

    /// Recombination rate
    #[arg(long, default_value_t = 1.0e-8)]
    recombination_rate: f64,

    /// Minor allele frequency threshold
    #[arg(long, default_value_t = 0.25)]
    maf_threshold: f64,

    /// Contig index
    #[arg(long, default_value_t = 0)]
    contig_index: u32,

    /// Minimum number of loci per bin
    #[arg(long, default_value_t = 2000)]
    min_loci: usize,

    /// Epsilon value for the CI Half-Width
    #[arg(long, default_value_t = 0.0001)]
    epsilon: f64,

    /// Random seed
    #[arg(long, default_value_t = 1234)]
    seed: u64,
}

fn find_contig_length(records: Vec<rust_htslib::bcf::HeaderRecord>, rid: u32) -> Result<u64> {
    let mut seen = 0;
    for record in records {
        if let HeaderRecord::Contig { values, .. } = record {
            if seen < rid {
                seen += 1;
                continue;
            }
            if seen > rid {
                break;
            }
            let contig_length = values.get("length").context("Contig length not found")?;
            if let Ok(contig_length) = contig_length.parse::<u64>() {
                return Ok(contig_length);
            }
        }
    }
    bail!("Contig not found")
}

fn main() -> Result<()> {
    // Read parameters from command line
    let args = Cli::parse();
    let src = &args.infile;
    let rid: u32 = args.contig_index;
    // Open indexed VCF or BCF (better)
    let mut file =
        IndexedReader::from_path(src).with_context(|| format!("Error opening file {src}"))?;
    // Get contig information
    let header = file.header();
    let records = header.header_records();
    let contig_length = find_contig_length(records, rid)
        .with_context(|| format!("Error finding contig length for the index {rid}"))?;
    // Initialize data structures
    let bins = Bins::hapne_default(args.recombination_rate);
    let mut summary_stats = SufficientSummaryStats::new(&bins);
    let num_samples = header.samples().len();
    let mut genotypes1: Vec<f64> = vec![0.0; num_samples];
    let mut genotypes2: Vec<f64> = vec![0.0; num_samples];
    // Chromosome to analyze
    let pb = ProgressBar::no_length();
    pb.set_message("Starting iterations...");
    let between = Uniform::try_from(0..contig_length).with_context(|| {
        format!(
            "Error creating uniform distribution from 0 to {contig_length}"
        )
    })?;
    // Get an RNG:
    let mut rng = ChaCha8Rng::seed_from_u64(args.seed);
    loop {
        // First, we draw a random position in the chromosome
        let sampled = between.sample(&mut rng);
        // Fetch the region from pos1 to pos1 + bins.maximum
        file.fetch(rid, sampled, None)
            .with_context(|| format!("Error fetching region {rid}:{sampled}:"))?;
        let record1 = file.records().next();
        if record1.is_none() {
            // If there are no records in the region, we skip this iteration
            continue;
        }
        let record1 = record1.unwrap().context("Error while reading record")?;
        let pos1 = record1.pos() as u64;
        match precompute_standarized_genotypes(record1, &args, &mut genotypes1) {
            Ok(Some(())) => {}
            Ok(None) => continue,
            Err(e) => {
                bail!("Error parsing genotypes: {}", e);
            }
        }
        // Now, we fetch the region from pos1 + bins.minimum to pos1 + bins.maximum
        let region = (pos1 + bins.minimum as u64, pos1 + bins.maximum as u64);
        file.fetch(rid, region.0, Some(region.1))
            .with_context(|| format!("Error fetching region {}:{}:{}", rid, region.0, region.1))?;
        // Most of the time, the second record will be in the first bin
        let mut index = 0;
        for record2 in file.records() {
            let record2 = record2.context("Error while reading record")?;
            let pos2 = record2.pos() as u64;
            let distance = (pos2 - pos1) as f64;
            if distance < (bins.minimum as f64) || distance > (bins.maximum as f64) {
                unreachable!("Distance is less than minimum or greater than maximum");
            }
            // Find current bin index
            while distance > bins.right_edges_in_bp[index] {
                index += 1;
            }
            assert!(index < bins.nbins);
            assert!(bins.left_edges_in_bp[index] <= distance);
            assert!(bins.right_edges_in_bp[index] >= distance);

            // Parse the genotypes of the second record and skip if MAF is too low
            match precompute_standarized_genotypes(record2, &args, &mut genotypes2) {
                Ok(Some(())) => {}
                Ok(None) => continue,
                Err(e) => {
                    bail!("Error parsing genotypes: {}", e);
                }
            };
            // Compute the sufficient statistics
            summary_stats.counts[index] += 1;
            let new_value = linkage_disequilibrium(&genotypes1, &genotypes2, num_samples);
            let delta = new_value - summary_stats.ld[index];
            summary_stats.ld[index] += delta / summary_stats.counts[index] as f64;
            let delta2 = new_value - summary_stats.ld[index];
            summary_stats.ld_square[index] += delta * delta2;
        }
        // Increment the progress bar
        pb.inc(1);
        if summary_stats.should_stop(&args) {
            break;
        }
    }
    pb.finish_with_message("Done!");
    // Finalize the summary statistics
    println!("#bin_index\tleft_bin\tright_bin\tN\tmean\tvar");
    for i in 0..bins.nbins {
        summary_stats.ld_square[i] /= summary_stats.counts[i] as f64;
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            i,
            bins.left_edges_in_cm[i] / 100.0,
            bins.right_edges_in_cm[i] / 100.0,
            summary_stats.counts[i],
            summary_stats.ld[i],
            summary_stats.ld_square[i]
        );
    }
    Ok(())
}
