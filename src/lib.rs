use pyo3::prelude::*;
use rust_htslib::{ bcf, bcf::Read, bcf::HeaderRecord };

use std::fs::File;
use std::io::Write;
use derive_new::new;
use itertools::Itertools;
use rust_htslib::bcf::record::Numeric;
use std::io;
use std::str;
use thiserror::Error;
use anyhow::bail;

#[derive(new)]
pub struct Writer {
    inner: io::BufWriter<File>,
    #[new(value = "0")]
    field_count: usize,
}

impl Writer {
    fn write_integer(&mut self, value: i32) -> anyhow::Result<()> {
        let fmt = if value.is_missing() { "".to_owned() } else { format!("{}", value) };
        self.write_field(fmt.as_bytes())
    }

    fn write_float(&mut self, value: f32) -> anyhow::Result<()> {
        let fmt = if value.is_missing() { "".to_owned() } else { format!("{}", value) };
        self.write_field(fmt.as_bytes())
    }

    fn write_flag(&mut self, value: bool) -> anyhow::Result<()> {
        self.write_field(format!("{}", value).as_bytes())
    }

    fn write_field(&mut self, value: &[u8]) -> anyhow::Result<()> {
        if self.field_count > 0 {
            self.inner.write_all(b"\t")?;
        }
        self.inner.write_all(value)?;
        self.field_count += 1;
        Ok(())
    }

    fn newline(&mut self) -> anyhow::Result<()> {
        self.inner.write_all(b"\n")?;
        self.field_count = 0;
        Ok(())
    }
}

const HEADER_COMMON: &[u8] = b"VARIANT";

#[pyfunction]
pub fn to_txt(vcf_path: &str, csv_path: &str) -> anyhow::Result<()> {
    let mut reader = bcf::Reader::from_path(vcf_path)?;
    let mut writer = Writer::new(
        io::BufWriter::new(File::options().write(true).create(true).open(csv_path)?)
    );

    let header_records = reader.header().header_records();

    let info_tags: Vec<&str> = header_records
        .iter()
        .filter_map(|records| (
            if let HeaderRecord::Info { values, .. } = records {
                if let Some(v) = values.get("ID").clone() { Some(v.as_str()) } else { None }
            } else {
                None
            }
        ))
        .collect();

    let format_tags: Vec<&str> = header_records
        .iter()
        .filter_map(|records| (
            if let HeaderRecord::Format { values, .. } = records {
                if let Some(v) = values.get("ID").clone() { Some(v.as_str()) } else { None }
            } else {
                None
            }
        ))
        .collect();

    println!("{:?}", format_tags);

    let common_n = 6 + info_tags.len();

    writer.write_field(HEADER_COMMON)?;
    for _ in 1..common_n {
        writer.write_field(HEADER_COMMON)?;
    }

    for sample in reader.header().samples() {
        writer.write_field(sample)?;
        for _ in 0..format_tags.len() - 1 {
            writer.write_field(sample)?;
        }
    }

    writer.newline()?;
    writer.write_field(b"CHROM")?;
    writer.write_field(b"POS")?;
    writer.write_field(b"REF")?;
    writer.write_field(b"ALT")?;
    writer.write_field(b"QUAL")?;
    writer.write_field(b"FILTER")?;

    for name in &info_tags {
        writer.write_field(name.as_bytes())?;
    }
    for _ in 0..reader.header().sample_count() {
        for name in &format_tags {
            writer.write_field(name.as_bytes())?;
        }
    }

    writer.newline()?;
    let mut rec = reader.empty_record();
    loop {
        match reader.read(&mut rec) {
            Some(Ok(())) => (),
            None => {
                break;
            }
            Some(Err(e)) => bail!(e),
        }
        let alleles = rec
            .alleles()
            .into_iter()
            .map(|a| a.to_owned())
            .collect_vec();
        for (i, allele) in alleles[1..].iter().enumerate() {
            writer.write_field(reader.header().rid2name(rec.rid().unwrap())?)?;
            writer.write_integer((rec.pos() as i32) + 1)?;
            writer.write_field(&alleles[0])?;
            writer.write_field(allele)?;
            match rec.qual() {
                q if q.is_missing() => writer.write_field(b"")?,
                q => writer.write_float(q)?,
            }

            if rec.has_filter(".".as_bytes()) {
                writer.write_field(b"")?;
            } else if rec.has_filter("PASS".as_bytes()) {
                writer.write_field(b"PASS")?;
            } else {
                let mut filters = Vec::new();
                for (i, filter) in rec.filters().enumerate() {
                    if i != 0 {
                        filters.push(b';');
                    }
                    filters.extend_from_slice(&reader.header().id_to_name(filter));
                }
                writer.write_field(&filters)?;
            }

            for name in &info_tags {
                let _name = name.as_bytes();
                if let Ok((tag_type, tag_length)) = rec.header().info_type(_name) {
                    let get_idx = || {
                        match tag_length {
                            bcf::header::TagLength::Fixed(_) => Ok(0),
                            bcf::header::TagLength::AltAlleles => Ok(i),
                            bcf::header::TagLength::Alleles => Ok(i + 1),
                            bcf::header::TagLength::Variable => Ok(0),
                            bcf::header::TagLength::Genotypes => Ok(0),
                            _ => Err(Box::new(ParseError::UnsupportedTagLength)),
                        }
                    };

                    match tag_type {
                        bcf::header::TagType::Flag => {
                            writer.write_flag(rec.info(_name).flag()?)?;
                        }
                        bcf::header::TagType::Integer => {
                            let i = get_idx()?;
                            if let Some(values) = rec.info(_name).integer().unwrap_or(None) {
                                writer.write_integer(values[i])?;
                            } else {
                                writer.write_field(b"")?;
                            }
                        }
                        bcf::header::TagType::Float => {
                            let i = get_idx()?;
                            if let Some(values) = rec.info(_name).float().unwrap_or(None) {
                                writer.write_float(values[i])?;
                            } else {
                                writer.write_field(b"")?;
                            }
                        }
                        bcf::header::TagType::String => {
                            let i = get_idx()?;
                            if let Some(values) = rec.info(_name).string().unwrap_or(None) {
                                writer.write_field(values[i])?;
                            } else {
                                writer.write_field(b"")?;
                            }
                        }
                    }
                } else {
                    // tag undefined, write NA
                    writer.write_field(b"")?;
                }
            }

            let genotypes = {
                let genotypes = rec.genotypes()?;

                Some(
                    (0..reader.header().sample_count() as usize)
                        .map(|s| format!("{}", genotypes.get(s)))
                        .collect_vec()
                )
            };

            for s in 0..reader.header().sample_count() as usize {
                for name in format_tags.iter() {
                    let _name = name.as_bytes();
                    if let Ok((tag_type, tag_length)) = reader.header().format_type(_name) {
                        let i = match tag_length {
                            bcf::header::TagLength::Fixed(_) => 0,
                            bcf::header::TagLength::AltAlleles => i,
                            bcf::header::TagLength::Alleles => i + 1,
                            _ => bail!(ParseError::UnsupportedTagLength),
                        };

                        match tag_type {
                            bcf::header::TagType::Flag => {
                                panic!(
                                    "Unable to find FORMAT \"{0}\" in the input file! Is \"{0}\" an INFO tag?",
                                    name
                                );
                            }
                            bcf::header::TagType::Integer => {
                                let field = if let Ok(v) = rec.format(_name).integer() {
                                    format!("{}", v[s][i])
                                } else {
                                    "".to_owned()
                                };
                                writer.write_field(field.as_bytes())?;
                            }
                            bcf::header::TagType::Float => {
                                let field = if let Ok(v) = rec.format(_name).float() {
                                    format!("{}", v[s][i])
                                } else {
                                    "".to_owned()
                                };
                                writer.write_field(field.as_bytes())?;
                            }
                            bcf::header::TagType::String => {
                                if _name != b"GT" {
                                    if let Ok(v) = rec.format(_name).string() {
                                        writer.write_field(
                                            std::str::from_utf8(v[s]).unwrap_or("").as_bytes()
                                        )?;
                                    } else {
                                        writer.write_field(b"")?;
                                    }
                                } else {
                                    //Treat genotypes different
                                    if let Some(ref gt) = genotypes {
                                        writer.write_field(gt[s].as_bytes())?;
                                    } else {
                                        writer.write_field(b"")?;
                                    }
                                }
                            }
                        }
                    } else {
                        // tag undefined, write NA
                        writer.write_field(b"")?;
                    }
                }
            }
            writer.newline()?;
        }
    }

    Ok(())
}

#[derive(Error, Debug)]
pub enum ParseError {
    #[error("currently, only R, A, and 1 are supported multiplicities of tags")]
    UnsupportedTagLength,
}

// #[pyfunction]
// fn generate_csv(input: &str, output: &str) -> PyResult<()> {
//     let mut bcf = Reader::from_path(input).expect("Error opening file.");
//     // iterate through each row of the vcf body.

//     let header = bcf.header().clone();
//     let header_records = header.header_records();

//     let samples = get_vcf_samples(&header);
//     let info_fields: Vec<&HeaderRecord> = get_info_fields(&header_records);
//     let format_fields: Vec<&HeaderRecord> = get_format_fields(&header_records);
//     let filters: Vec<&HeaderRecord> = get_filters(&header_records);

//     let mut skipped_records: usize = 0;
//     let mut output_file = File::options().write(true).create(true).open(output)?;

//     // WRITE CSV HEADER

//     //Write mandatory fields header

//     output_file
//         .write(&["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"].join("\t").as_bytes())
//         .expect("Cannot write");

//     output_file.write(b"\t").expect("Cannot write");

//     //Write info fields
//     output_file
//         .write(
//             info_fields
//                 .iter()
//                 .filter_map(|info| {
//                     if let HeaderRecord::Info { values, .. } = info {
//                         Some(
//                             format!(
//                                 "info_{}",
//                                 values.get("ID").unwrap_or(&String::from("MALFORMED_VCF"))
//                             )
//                         )
//                     } else {
//                         None
//                     }
//                 })
//                 .collect::<Vec<String>>()
//                 .join("\t")
//                 .as_bytes()
//         )
//         .expect("Cannot write");

//     output_file.write(b"\t").expect("Cannot write");

//     //Write format fields
//     output_file
//         .write(
//             iproduct!(samples.iter(), format_fields.iter())
//                 .map(|(sample, fmt)| {
//                     if let HeaderRecord::Format { values, .. } = fmt {
//                         format!(
//                             "{}_{}",
//                             sample,
//                             values.get("ID").unwrap_or(&String::from("MALFORMED_VCF"))
//                         )
//                     } else {
//                         "".to_string()
//                     }
//                 })
//                 .collect::<Vec<String>>()
//                 .join("\t")
//                 .as_bytes()
//         )
//         .expect("Cannot write");

//     output_file.write(b"\n").expect("Cannot write");

//     for (i, record_result) in bcf.records().enumerate() {
//         let Ok(record) = record_result else {
//             skipped_records += 1;
//             eprintln!("Skipping record {}", i);
//             continue;
//         };

//         for allele_id in 1..record.alleles().len() {
//             //Write CHROM
//             output_file.write(record.contig().as_bytes()).expect("Cannot write");
//             output_file.write(b"\t").expect("Cannot write");

//             //Write POS
//             output_file.write(format!("{}", record.pos() + 1).as_bytes()).expect("Cannot write");
//             output_file.write(b"\t").expect("Cannot write");

//             //Write ID
//             output_file
//                 .write(format!("{}", from_utf8(&record.id()).unwrap_or("INVALID UTF8")).as_bytes())
//                 .expect("Cannot write");
//             output_file.write(b"\t").expect("Cannot write");

//             //Write REF
//             output_file
//                 .write(record.alleles().get(0).unwrap_or(&"N".as_bytes()))
//                 .expect("Cannot write");
//             output_file.write(b"\t").expect("Cannot write");

//             //Write ALT
//             output_file
//                 .write(record.alleles().get(allele_id).unwrap_or(&"N".as_bytes()))
//                 .expect("Cannot write");
//             output_file.write(b"\t").expect("Cannot write");

//             //Write QUAL
//             output_file.write(format!("{:}", record.qual()).as_bytes()).expect("Cannot write");
//             output_file.write(b"\t").expect("Cannot write");

//             //Write FILTER...
//             output_file
//                 .write(
//                     record
//                         .filters()
//                         .filter_map(|filter| (
//                             if
//                                 let Some(HeaderRecord::Filter { values, .. }) = filters.get(
//                                     usize::try_from(filter.0).unwrap()
//                                 )
//                             {
//                                 Some(values.get("ID").unwrap_or(&"PASS".to_owned()).to_owned())
//                             } else {
//                                 None
//                             }
//                         ))
//                         .collect::<Vec<String>>()
//                         .join(",")
//                         .as_bytes()
//                 )
//                 .expect("Cannot write");
//             // for (sample, fmt) in iproduct!(samples.iter(), format_fields.iter()) {
//             // }
//             output_file.write(b"\t").expect("Cannot write");

//             // Write INFO
//             for info in info_fields.iter() {
//                 if let HeaderRecord::Info { values, .. } = info {
//                     output_file
//                         .write(
//                             format!(
//                                 "{}={}\t",
//                                 values.get("ID").unwrap(),
//                                 info_res
//                                     .map(|tagval_opt| {
//                                         match tagval_opt {
//                                             Some(tagval) => {
//                                                 (*tagval)
//                                                     .iter()
//                                                     .map(|s|
//                                                         String::from_utf8(s.to_vec()).unwrap_or(
//                                                             "".to_owned()
//                                                         )
//                                                     )
//                                                     .collect::<Vec<String>>()
//                                                     .join(",")
//                                             }
//                                             None => { "".to_owned() }
//                                         }
//                                     })
//                                     .unwrap_or("?".to_owned())
//                             ).as_bytes()
//                         )
//                         .expect("Cannot write");
//                 }
//             }

//             for sample in samples.iter() {
//                 for fmt in format_fields.iter() {
//                     if let HeaderRecord::Format { values, .. } = fmt {
//                         let format_field = record
//                             .format(values.get("ID").unwrap().as_bytes())
//                             .string()
//                             .map(|tagval| {
//                                 (*tagval)
//                                     .iter()
//                                     .map(|s| String::from_utf8(s.to_vec()).unwrap_or("".to_owned()))
//                                     .collect::<Vec<String>>()
//                                     .join("\t")
//                             })
//                             .unwrap_or("".to_owned());
//                         output_file.write(format_field.as_bytes()).expect("Cannot write");
//                     }
//                 }
//             }
//             output_file.write(b"\n").expect("Cannot write");
//         }
//     }
//     if skipped_records > 0 {
//         eprintln!("Skipped {} records", skipped_records);
//     }
//     Ok(())
// }

/// A Python module implemented in Rust.
#[pymodule]
fn vcf2csv(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(to_txt, m)?)?;
    Ok(())
}
