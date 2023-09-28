/// manysketch: massively parallel sketching of sequence files.
use anyhow::{anyhow, Result};
use rayon::prelude::*;

use crate::utils::{sigwriter, Params, ZipMessage};
use needletail::parse_fastx_file;
use sourmash::cmd::ComputeParameters;
use sourmash::signature::Signature;
use std::path::Path;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;
use std::path::PathBuf;

enum CSVType {
    Assembly,
    Reads,
    Unknown,
}

fn detect_csv_type(headers: &csv::StringRecord) -> CSVType {
    if headers.len() == 3
        && headers.get(0).unwrap() == "name"
        && headers.get(1).unwrap() == "genome_filename"
        && headers.get(2).unwrap() == "protein_filename"
    {
        CSVType::Assembly
    } else if headers.len() == 3
        && headers.get(0).unwrap() == "name"
        && headers.get(1).unwrap() == "read1"
        && headers.get(2).unwrap() == "read2"
    {
        CSVType::Reads
    } else {
        CSVType::Unknown
    }
}

pub fn load_fasta_fromfile<P: AsRef<Path>>(
    sketchlist_filename: &P,
) -> Result<Vec<(String, Vec<PathBuf>, String)>> {
    let mut rdr = csv::Reader::from_path(sketchlist_filename)?;
    let headers = rdr.headers()?;

    match detect_csv_type(&headers) {
        CSVType::Assembly => process_assembly_csv(rdr),
        CSVType::Reads => process_reads_csv(rdr),
        CSVType::Unknown => Err(anyhow!(
            "Invalid header. Expected 'name,genome_filename,protein_filename' or 'name,read1,read2', but got '{}'",
            headers.iter().collect::<Vec<_>>().join(",")
        )),
    }
}

fn process_assembly_csv(mut rdr: csv::Reader<std::fs::File>) -> Result<Vec<(String, Vec<PathBuf>, String)>> {
    let mut results = Vec::new();
    let mut row_count = 0;
    let mut genome_count = 0;
    let mut protein_count = 0;
    let mut processed_rows = std::collections::HashSet::new();
    let mut duplicate_count = 0;

    for result in rdr.records() {
        let record = result?;
        let row_string = record.iter().collect::<Vec<_>>().join(",");
        if processed_rows.contains(&row_string) {
            duplicate_count += 1;
            continue;
        }
        processed_rows.insert(row_string.clone());
        row_count += 1;
        let name = record
            .get(0)
            .ok_or_else(|| anyhow!("Missing 'name' field"))?
            .to_string();

        let genome_filename = record
            .get(1)
            .ok_or_else(|| anyhow!("Missing 'genome_filename' field"))?;
        if !genome_filename.is_empty() {
            results.push((
                name.clone(),
                vec![PathBuf::from(genome_filename)],
                "dna".to_string(),
            ));
            genome_count += 1;
        }

        let protein_filename = record
            .get(2)
            .ok_or_else(|| anyhow!("Missing 'protein_filename' field"))?;
        if !protein_filename.is_empty() {
            results.push((
                name,
                vec![PathBuf::from(protein_filename)],
                "protein".to_string(),
            ));
            protein_count += 1;
        }
    }

    if duplicate_count > 0 {
        println!("Warning: {} duplicated rows were skipped.", duplicate_count);
    }
    println!(
        "Loaded {} rows in total ({} genome and {} protein files)",
        row_count, genome_count, protein_count
    );

    Ok(results)
}

fn process_reads_csv(mut rdr: csv::Reader<std::fs::File>) -> Result<Vec<(String, Vec<PathBuf>, String)>> {
    let mut results = Vec::new();
    let mut processed_rows = std::collections::HashSet::new();
    let mut duplicate_count = 0;

    for result in rdr.records() {
        let record = result?;
        let row_string = record.iter().collect::<Vec<_>>().join(",");
        if processed_rows.contains(&row_string) {
            duplicate_count += 1;
            continue;
        }
        processed_rows.insert(row_string.clone());

        let name = record
            .get(0)
            .ok_or_else(|| anyhow!("Missing 'name' field"))?
            .to_string();

        let read1 = record
            .get(1)
            .ok_or_else(|| anyhow!("Missing 'read1' field"))?;
        let read2 = record
            .get(2)
            .ok_or_else(|| anyhow!("Missing 'read2' field"))?;
        let paths = vec![
            PathBuf::from(read1),
            PathBuf::from(read2),
        ];
        results.push((name, paths, "alternate".to_string()));
    }

    if duplicate_count > 0 {
        println!("Warning: {} duplicated rows were skipped.", duplicate_count);
    }
    println!("Loaded alternate CSV variant.");

    Ok(results)
}


fn parse_params_str(params_strs: String) -> Result<Vec<Params>, String> {
    let mut unique_params: std::collections::HashSet<Params> = std::collections::HashSet::new();

    // split params_strs by _ and iterate over each param
    for p_str in params_strs.split('_').collect::<Vec<&str>>().iter() {
        let items: Vec<&str> = p_str.split(',').collect();

        let mut ksizes = Vec::new();
        let mut track_abundance = false;
        let mut num = 0;
        let mut scaled = 1000;
        let mut seed = 42;
        let mut is_protein = false;
        let mut is_dna = true;

        for item in items.iter() {
            match *item {
                _ if item.starts_with("k=") => {
                    let k_value = item[2..]
                        .parse()
                        .map_err(|_| format!("cannot parse k='{}' as a number", &item[2..]))?;
                    ksizes.push(k_value);
                }
                "abund" => track_abundance = true,
                "noabund" => track_abundance = false,
                _ if item.starts_with("num=") => {
                    num = item[4..]
                        .parse()
                        .map_err(|_| format!("cannot parse num='{}' as a number", &item[4..]))?;
                }
                _ if item.starts_with("scaled=") => {
                    scaled = item[7..]
                        .parse()
                        .map_err(|_| format!("cannot parse scaled='{}' as a number", &item[7..]))?;
                }
                _ if item.starts_with("seed=") => {
                    seed = item[5..]
                        .parse()
                        .map_err(|_| format!("cannot parse seed='{}' as a number", &item[5..]))?;
                }
                "protein" => {
                    is_protein = true;
                    is_dna = false;
                }
                "dna" => {
                    is_protein = false;
                    is_dna = true;
                }
                _ => return Err(format!("unknown component '{}' in params string", item)),
            }
        }

        for &k in &ksizes {
            let param = Params {
                ksize: k,
                track_abundance,
                num,
                scaled,
                seed,
                is_protein,
                is_dna,
            };
            unique_params.insert(param);
        }
    }

    Ok(unique_params.into_iter().collect())
}

fn build_siginfo(
    params: &[Params],
    moltype: &str,
    name: &str,
    filename: &Path,
) -> (Vec<Signature>, Vec<Params>) {
    let mut sigs = Vec::new();
    let mut params_vec = Vec::new();

    for param in params.iter().cloned() {
        match moltype {
            // if dna, only build dna sigs. if protein, only build protein sigs
            "dna" if !param.is_dna => continue,
            "protein" if !param.is_protein => continue,
            _ => (),
        }

        // Adjust ksize value based on the is_protein flag
        let adjusted_ksize = if param.is_protein {
            param.ksize * 3
        } else {
            param.ksize
        };

        let cp = ComputeParameters::builder()
            .ksizes(vec![adjusted_ksize])
            .scaled(param.scaled)
            .protein(param.is_protein)
            .dna(param.is_dna)
            .num_hashes(param.num)
            .track_abundance(param.track_abundance)
            .build();

        // let sig = Signature::from_params(&cp); // cant set name with this
        let template = sourmash::cmd::build_template(&cp);
        let sig = Signature::builder()
            .hash_function("0.murmur64")
            .name(Some(name.to_string()))
            .filename(Some(filename.to_string_lossy().into_owned()))
            .signatures(template)
            .build();
        sigs.push(sig);

        params_vec.push(param);
    }

    (sigs, params_vec)
}


fn sketch_fastas(
    filenames: &[PathBuf], 
    moltype: &str, 
    name: &str, 
    params_vec: &[Params]
) -> Result<(Vec<Signature>, Vec<Params>, PathBuf), String> {

    // initialize signatures once for all files
    // in sourmash sketch, if merging multiple files, the `filename` ends up being the last file processed.
    // do we want to mirror that behavior? or is there a better option?
    let (mut sigs, sig_params) = build_siginfo(&params_vec, moltype, name, &filenames[0]);

    // if no sigs to build, skip
    if sigs.is_empty() {
        return Ok((Vec::new(), Vec::new(), filenames[0].clone()));
    }
    
    for filename in filenames.iter() {
        let mut reader = match needletail::parse_fastx_file(filename) {
            Ok(r) => r,
            Err(err) => {
                eprintln!("Error opening file {}: {:?}", filename.display(), err);
                return Err(format!("Error opening file {}: {:?}", filename.display(), err));
            }
        };

        // parse fasta and add to signature
        while let Some(record_result) = reader.next() {
            match record_result {
                Ok(record) => {
                    for sig in &mut sigs {
                        if moltype == "protein" {
                            sig.add_protein(&record.seq()).unwrap();
                        } else {
                            sig.add_sequence(&record.seq(), true).unwrap();
                        }
                    }
                }
                Err(err) => {
                    eprintln!("Error while processing record: {:?}", err);
                }
            }
        }
    }

    Ok((sigs, sig_params, filenames[0].clone()))
}


fn sketch_singleton(
    filename: &PathBuf, 
    moltype: &str, 
    name: &str, 
    params_vec: &[Params]
) -> Result<(Vec<Signature>, Vec<Params>, PathBuf), String> {
    
    let mut reader = match parse_fastx_file(filename) {
        Ok(r) => r,
        Err(err) => {
            eprintln!("Error opening file {}: {:?}", filename.display(), err);
            return Err(format!("Error opening file {}: {:?}", filename.display(), err));
        }
    };
    //build sig templates
    let (template_sigs, sig_params) = build_siginfo(&params_vec, moltype, name, filename);
    // if no sigs to build, skip
    if template_sigs.is_empty() {
        return Ok((Vec::new(), Vec::new(), filename.clone()));
    }
    
    let mut all_sigs: Vec<Signature> = vec![];
    let mut all_sig_params: Vec<Params> = vec![];
    all_sig_params.extend(sig_params.iter().cloned());
    
    
    // parse fasta and add to a new signature for each record
    while let Some(record_result) = reader.next() {
        match record_result {
            Ok(record) => {
                // copy template sigs for each record
                let mut record_sigs = template_sigs.clone();
                for mut sig in record_sigs.iter_mut() {
                    if moltype == "protein" {
                        sig.add_protein(&record.seq()).unwrap();
                    } else {
                        sig.add_sequence(&record.seq(), true).unwrap();
                    }
                }
                all_sigs.extend(record_sigs);
            }
            Err(err) => {
                eprintln!("Error while processing record: {:?}", err);
            }
        }
    }

    Ok((all_sigs, all_sig_params, filename.clone()))
}



pub fn manysketch<P: AsRef<Path> + Sync>(
    filelist: P,
    param_str: String,
    output: String,
    singleton: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let fileinfo = match load_fasta_fromfile(&filelist) {
        Ok(result) => result,
        Err(e) => bail!("Could not load fromfile csv. Underlying error: {}", e),
    };

    // if no files to process, exit with error
    let n_fastas = fileinfo.len();
    if n_fastas == 0 {
        bail!("No files to load, exiting.");
    }

    // if output doesnt end in zip, bail
    if Path::new(&output)
        .extension()
        .map_or(true, |ext| ext != "zip")
    {
        bail!("Output must be a zip file.");
    }

    // set up a multi-producer, single-consumer channel that receives Signature
    let (send, recv) = std::sync::mpsc::sync_channel::<ZipMessage>(rayon::current_num_threads());
    // need to use Arc so we can write the manifest after all sigs have written
    let send = std::sync::Arc::new(send);

    // & spawn a thread that is dedicated to printing to a buffered output
    let thrd = sigwriter::<&str>(recv, output);

    // parse param string into params_vec, print error if fail
    let param_result = parse_params_str(param_str);
    let params_vec = match param_result {
        Ok(params) => params,
        Err(e) => {
            eprintln!("Error parsing params string: {}", e);
            bail!("Failed to parse params string");
        }
    };

    // iterate over filelist_paths
    let processed_fastas = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);
    let skipped_paths: AtomicUsize = AtomicUsize::new(0);

    // set reporting threshold at every 5% or every 1 fasta, whichever is larger)
    let reporting_threshold = std::cmp::max(n_fastas / 20, 1);

    let send_result = fileinfo
    .par_iter()
    .filter_map(|(name, filenames, moltype)| {
        // increment processed_fastas counter; make 1-based for % reporting
        let i = processed_fastas.fetch_add(1, atomic::Ordering::SeqCst);
        
        // progress report at threshold
        if (i + 1) % reporting_threshold == 0 {
            let percent_processed = (((i + 1) as f64 / n_fastas as f64) * 100.0).round();
            eprintln!(
                "Starting file {}/{} ({}%)",
                (i + 1),
                n_fastas,
                percent_processed
            );
        }
        
        let (sigs, sig_params, filename) = if singleton {
            sketch_singleton(filenames, moltype, name, &params_vec).unwrap_or_default()
        } else {
            sketch_fastas(filenames, moltype, name, &params_vec).unwrap_or_default()
        };
        
        // if no sigs to build, skip
        if sigs.is_empty() {
            let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
            return None;
        }
        
        Some((sigs, sig_params, filename))
    })
    .try_for_each_with(
        send.clone(),
        |s: &mut std::sync::Arc<std::sync::mpsc::SyncSender<ZipMessage>>,
         (sigs, sig_params, filename)| {
            if let Err(e) = s.send(ZipMessage::SignatureData(
                sigs,
                sig_params,
                filename.clone(),
            )) {
                Err(format!("Unable to send internal data: {:?}", e))
            } else {
                Ok(())
            }
        },
    );
This assumes:

Both sketch_singleton and sketch_fastas return Result<(Vec<Signature>, Vec<SigParams>, PathBuf), SomeErrorType> where SomeErrorType is the error type of the functions.
SomeErrorType implements the default trait or you provide a suitable default.
Note: I've added .unwrap_or_default() for handling errors returned from the sketch functions. This will panic if there's an error. Depending on your application, you might want to replace this with a proper error-handling strategy.







    // After the parallel work, send the WriteManifest message
    std::sync::Arc::try_unwrap(send)
        .unwrap()
        .send(ZipMessage::WriteManifest)
        .unwrap();

    // do some cleanup and error handling -
    if let Err(e) = send_result {
        eprintln!("Error during parallel processing: {}", e);
    }

    // join the writer thread
    if let Err(e) = thrd
        .join()
        .unwrap_or_else(|e| Err(anyhow!("Thread panicked: {:?}", e)))
    {
        eprintln!("Error in sigwriter thread: {:?}", e);
    }

    // done!
    let i: usize = processed_fastas.load(atomic::Ordering::SeqCst);
    eprintln!("DONE. Processed {} fasta files", i);

    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

    if failed_paths == i {
        bail!("Could not load fasta files: no signatures created.");
    }
    if failed_paths > 0 {
        eprintln!(
            "WARNING: {} fasta files failed to load. See error messages above.",
            failed_paths
        );
    }

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    if skipped_paths == i {
        bail!("No fasta files compatible with provided sketch parameters: no signatures created.");
    }
    if skipped_paths > 0 {
        eprintln!(
            "WARNING: {} fasta files skipped - no compatible signatures.",
            skipped_paths
        );
    }

    Ok(())
}
