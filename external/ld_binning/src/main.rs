#![feature(portable_simd)]
use anyhow::{bail, Context, Result};
use clap::{command, Parser};
use indicatif::ProgressBar;
use rust_htslib::bcf::{record, IndexedReader, Read, Record};
use std::error::Error;

// Parse the genotypes and compute the minor allele frequency
fn precompute_standarized_genotypes(
    record: Record,
    parameters: &HyperParameters,
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

pub fn linkage_disequilibrium(
    genotypes1: &[f64],
    genotypes2: &[f64],
    n_samples: usize,
) -> f64 {
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
}

struct HyperParameters {
    pub recombination_rate: f64,
    pub maf_threshold: f64,
}

impl Default for HyperParameters {
    fn default() -> Self {
        Self {
            recombination_rate: 1.0e-8,
            maf_threshold: 0.25,
        }
    }
}

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Compressed and indexed VCF file
    infile: String,
}

fn main() -> Result<()> {
    // Read from file so we can have two readers
    let args = Cli::parse();
    let src = &args.infile;
    // Open the VCF (bgzipped and indexed)
    let mut vcf1 = IndexedReader::from_path(src).context("Error opening VCF file")?;
    let mut vcf2 = IndexedReader::from_path(src).context("Error opening VCF file")?;
    // Get the header to resolve contig names to IDs
    let header = vcf1.header().clone();
    let parameters = HyperParameters::default();
    // Initialize data structures
    let bins = Bins::hapne_default(parameters.recombination_rate);
    let mut summary_stats = SufficientSummaryStats::new(&bins);
    let num_samples = header.samples().len();
    let mut genotypes1: Vec<f64> = vec![0.0; num_samples];
    let mut genotypes2: Vec<f64> = vec![0.0; num_samples];
    // For now, assume there is a single chromosome per file
    let rid = 0;
    vcf1.fetch(rid, 0, None)
        .context("Error fetching the chromosome")?;
    let num_records = vcf1.records().count();
    vcf1.fetch(rid, 0, None)
        .context("Error fetching the chromosome")?;
    let pb = ProgressBar::new((num_records) as u64);
    for record1 in vcf1.records() {
        let record1 = record1.context("Error while reading record")?;
        let pos1 = record1.pos();
        // Parse the genotypes of the first record and compute the minor allele frequency
        match precompute_standarized_genotypes(record1, &parameters, &mut genotypes1) {
            Ok(Some(())) => {}
            Ok(None) => continue,
            Err(e) => {
                bail!("Error parsing genotypes: {}", e);
            }
        };
        // Most of the time, the second record will be in the first bin
        let mut index = 0;
        // Fetch the second region of interest
        vcf2.fetch(rid, (pos1 + bins.minimum) as u64, None)
            .context("Error fetching second region")?;
        for record2 in vcf2.records() {
            let record2 = record2.context("Error while reading record")?;
            let pos2 = record2.pos();
            let distance = (pos2 - pos1) as f64;
            if distance < (bins.minimum as f64) {
                unreachable!("Distance is less than minimum");
            }
            if distance > (bins.maximum as f64) {
                break;
            }
            while distance > bins.right_edges_in_bp[index] {
                index += 1;
            }
            assert!(index < bins.nbins);
            assert!(bins.left_edges_in_bp[index] <= distance);
            assert!(bins.right_edges_in_bp[index] >= distance);
            // Parse the genotypes of the second record
            match precompute_standarized_genotypes(record2, &parameters, &mut genotypes2) {
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
        pb.inc(1);
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
