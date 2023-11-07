use pyo3::prelude::*;
use rust_htslib::{ bcf::{ Reader, Read }, bcf::header::HeaderRecord };
use bio_types::genome::AbstractLocus;

use std::{ fs::File, io::Write, str::{ from_utf8, FromStr } };

use itertools::{ iproduct, Itertools };

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

#[pyfunction]
fn generate_csv(input: &str, output: &str) -> PyResult<()> {
    let mut bcf = Reader::from_path(input).expect("Error opening file.");
    // iterate through each row of the vcf body.

    let header_records = bcf.header().header_records();

    let sample_count = usize::try_from(bcf.header().sample_count()).unwrap();

    let info_fields: Vec<&HeaderRecord> = header_records
        .iter()
        .filter(|r| matches!(r, HeaderRecord::Info { .. }))
        .collect();

    let format_fields: Vec<&HeaderRecord> = header_records
        .iter()
        .filter(|r| matches!(r, HeaderRecord::Format { .. }))
        .collect();

    let filters: Vec<&HeaderRecord> = header_records
        .iter()
        .filter(|r| matches!(r, HeaderRecord::Filter { .. }))
        .collect();

    let samples: Vec<&str> = bcf
        .header()
        .samples()
        .iter()
        .map(|s| from_utf8(s).unwrap())
        .collect();

    let mut skipped_records: usize = 0;

    let mut output_file = File::options().write(true).open(output)?;

    // WRITE CSV HEADER

    //Write mandatory fields header

    output_file
        .write(&["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"].join("\t").as_bytes())
        .expect("Cannot write");

    output_file.write(b"\t").expect("Cannot write");

    //Write info fields
    output_file
        .write(
            info_fields
                .iter()
                .filter_map(|info| {
                    if let HeaderRecord::Info { values, .. } = info {
                        Some(
                            format!(
                                "info_{}",
                                values.get("ID").unwrap_or(&String::from("MALFORMED_VCF"))
                            )
                        )
                    } else {
                        None
                    }
                })
                .collect::<Vec<String>>()
                .join("\t")
                .as_bytes()
        )
        .expect("Cannot write");

    output_file.write(b"\t").expect("Cannot write");

    //Write format fields
    output_file
        .write(
            iproduct!(bcf.header().samples().iter(), format_fields.iter())
                .map(|(sample, fmt)| {
                    if let HeaderRecord::Format { values, .. } = fmt {
                        format!(
                            "{}_{}",
                            from_utf8(sample).unwrap_or("UNKOWN"),
                            values.get("ID").unwrap_or(&String::from("MALFORMED_VCF"))
                        )
                    } else {
                        "".to_string()
                    }
                })
                .collect::<Vec<String>>()
                .join("\t")
                .as_bytes()
        )
        .expect("Cannot write");

    output_file.write(b"\n").expect("Cannot write");

    for (i, record_result) in bcf.records().enumerate() {
        let Ok(record) = record_result else {
            skipped_records += 1;
            eprintln!("Skipping record {}", i);
            continue;
        };

        for allele_id in 1..record.alleles().len() {
            //Write CHROM
            output_file.write(record.contig().as_bytes()).expect("Cannot write");
            output_file.write(b"\t").expect("Cannot write");

            //Write POS
            output_file.write(format!("{}", record.pos()).as_bytes()).expect("Cannot write");
            output_file.write(b"\t").expect("Cannot write");

            //Write ID
            output_file
                .write(format!("{}", from_utf8(&record.id()).unwrap()).as_bytes())
                .expect("Cannot write");

            //Write REF
            output_file
                .write(record.alleles().get(0).unwrap_or(&"N".as_bytes()))
                .expect("Cannot write");

            //Write ALT
            output_file
                .write(record.alleles().get(allele_id).unwrap_or(&"N".as_bytes()))
                .expect("Cannot write");

            //Write QUAL
            output_file.write(format!("{:}", record.qual()).as_bytes()).expect("Cannot write");

            //Write FILTER...
            output_file
                .write(
                    record
                        .filters()
                        .filter_map(|filter| (
                            if
                                let HeaderRecord::Filter { values, .. } = filters
                                    .get(usize::try_from(filter.0).unwrap())
                                    .unwrap()
                            {
                                Some(values.get("ID").unwrap_or(&"PASS".to_owned()).to_owned())
                            } else {
                                None
                            }
                        ))
                        .collect::<Vec<String>>()
                        .join(",")
                        .as_bytes()
                )
                .expect("Cannot write");
            // for (sample, fmt) in iproduct!(samples.iter(), format_fields.iter()) {
            // }

            // Write INFO
            for info in info_fields.iter() {
                if let HeaderRecord::Info { values, .. } = info {
                    output_file
                        .write(
                            record
                                .info(values.get("ID").unwrap().as_bytes())
                                .string()
                                .map(|tagval_opt| {
                                    match tagval_opt {
                                        Some(tagval) => {
                                            (*tagval)
                                                .iter()
                                                .map(|s|
                                                    String::from_utf8(s.to_vec()).unwrap_or(
                                                        "".to_owned()
                                                    )
                                                )
                                                .collect::<Vec<String>>()
                                                .join(",")
                                        }
                                        None => { "".to_owned() }
                                    }
                                })
                                .unwrap_or("".to_owned())
                                .as_bytes()
                        )
                        .expect("Cannot write");
                }
            }

            for sample in samples.iter() {
                for fmt in format_fields.iter() {
                    if let HeaderRecord::Format { values, .. } = fmt {
                        let format_field = record
                            .format(values.get("ID").unwrap().as_bytes())
                            .string()
                            .map(|tagval| {
                                (*tagval)
                                    .iter()
                                    .map(|s| String::from_utf8(s.to_vec()).unwrap_or("".to_owned()))
                                    .collect::<Vec<String>>()
                                    .join(",")
                            })
                            .unwrap_or("".to_owned());
                        output_file.write(format_field.as_bytes()).expect("Cannot write");
                    }
                }
            }
        }
    }
    if skipped_records > 0 {
        eprintln!("Skipped {} records", skipped_records);
    }
    Ok(())
}

/// A Python module implemented in Rust.
#[pymodule]
fn vcf2csv(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(generate_csv, m)?)?;
    Ok(())
}
