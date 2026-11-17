//! Exact int32 reconstruction from compressed interleaved int16 and sparse overflows.
use anyhow::{ensure, Result};
use hdf5::Group;
use std::path::Path;

pub fn write(g: &Group, name: &str, shape: &[usize], values: &[i32]) -> Result<()> {
    ensure!(
        shape.iter().product::<usize>() == values.len(),
        "IQ dimensions"
    );
    let compact = values
        .iter()
        .map(|&v| v.clamp(i16::MIN as i32, i16::MAX as i32) as i16)
        .collect::<Vec<_>>();
    let mut index = vec![];
    let mut words = vec![];
    for (i, &v) in values.iter().enumerate() {
        if v < i16::MIN as i32 || v > i16::MAX as i32 {
            index.push(i as u64);
            words.extend([v as u16 as i16, (v >> 16) as i16]);
        }
    }
    let mut chunk = shape.to_vec();
    chunk[0] = chunk[0].min(16);
    g.new_dataset::<i16>()
        .shape(shape)
        .chunk(chunk)
        .shuffle()
        .deflate(4)
        .create(format!("{name}_iq_i16").as_str())?
        .write_raw(&compact)?;
    g.new_dataset::<u64>()
        .shape([index.len()])
        .create(format!("{name}_overflow_index").as_str())?
        .write_raw(&index)?;
    g.new_dataset::<i16>()
        .shape([index.len(), 2])
        .create(format!("{name}_overflow_words_i16").as_str())?
        .write_raw(&words)?;
    let restored = read(g, name)?;
    ensure!(restored == values, "compressed IQ round trip failed");
    Ok(())
}

pub fn read(g: &Group, name: &str) -> Result<Vec<i32>> {
    let compact = format!("{name}_iq_i16");
    if !g.link_exists(&compact) {
        return Ok(g.dataset(&format!("{name}_iq_i32"))?.read_raw::<i32>()?);
    }
    let mut values = g
        .dataset(&compact)?
        .read_raw::<i16>()?
        .iter()
        .map(|&v| v as i32)
        .collect::<Vec<_>>();
    if !g.link_exists(&format!("{name}_overflow_index")) {
        ensure!(g.link_exists(&format!("{name}_overflow_index_empty")),"missing overflow table");
        return Ok(values);
    }
    let index = g
        .dataset(&format!("{name}_overflow_index"))?
        .read_raw::<u64>()?;
    let words = g
        .dataset(&format!("{name}_overflow_words_i16"))?
        .read_raw::<i16>()?;
    ensure!(words.len() == index.len() * 2, "invalid overflow table");
    for (i, &j) in index.iter().enumerate() {
        ensure!((j as usize) < values.len(), "invalid overflow index");
        values[j as usize] = ((words[2 * i + 1] as i32) << 16) | (words[2 * i] as u16 as i32);
    }
    Ok(values)
}

fn copy_group(from: &Group, to: &Group, skip_voltage: bool) -> Result<()> {
    for name in from.member_names()? {
        if skip_voltage && (name == "zrx_echoes_iq_i32" || name == "ztx_history_iq_i32") {
            continue;
        }
        if let Ok(d) = from.dataset(&name) {
            let dtype = d.dtype()?;
            macro_rules! typed {
                ($t:ty) => {
                    if dtype.is::<$t>() {
                        to.new_dataset::<$t>()
                            .shape(d.shape())
                            .create(name.as_str())?
                            .write_raw(&d.read_raw::<$t>()?)?;
                        continue;
                    }
                };
            }
            typed!(u8);
            typed!(u16);
            typed!(u32);
            typed!(u64);
            typed!(i8);
            typed!(i16);
            typed!(i32);
            typed!(i64);
            typed!(f32);
            typed!(f64);
            anyhow::bail!("unsupported metadata type in {name}");
        } else {
            copy_group(&from.group(&name)?, &to.create_group(&name)?, false)?;
        }
    }
    Ok(())
}

pub fn repack(source: &Path, destination: &Path) -> Result<()> {
    ensure!(!destination.exists(), "destination exists");
    let partial = destination.with_extension("partial");
    ensure!(!partial.exists(), "partial destination exists");
    std::fs::create_dir_all(partial.join("cuts"))?;
    let index = hdf5::File::open(source.join("index.h5"))?;
    ensure!(
        index.dataset("complete")?.read_scalar::<u8>()? == 1,
        "incomplete source"
    );
    let new_index = hdf5::File::create(partial.join("index.h5"))?;
    copy_group(&index, &new_index, false)?;
    new_index
        .dataset("schema")?
        .write_raw(b"tibet_meteor_cut_v2")?;
    for entry in std::fs::read_dir(source.join("cuts"))? {
        let path = entry?.path();
        if path.extension().and_then(|s| s.to_str()) != Some("h5") {
            continue;
        }
        let old = hdf5::File::open(&path)?;
        let new = hdf5::File::create(partial.join("cuts").join(path.file_name().unwrap()))?;
        copy_group(&old, &new, true)?;
        for prefix in ["zrx_echoes", "ztx_history"] {
            let ds = old.dataset(&format!("{prefix}_iq_i32"))?;
            write(&new, prefix, &ds.shape(), &ds.read_raw::<i32>()?)?;
        }
        new.dataset("schema")?.write_raw(b"tibet_meteor_cut_v2")?;
        new.flush()?;
    }
    new_index.flush()?;
    drop(new_index);
    std::fs::rename(partial, destination)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_extreme_values_and_filters() {
        let path = std::env::temp_dir().join(format!("tibet-voltage-{}.h5", std::process::id()));
        let f = hdf5::File::create(&path).unwrap();
        let v = vec![i32::MIN, -32769, -32768, -1, 0, 32767, 32768, i32::MAX];
        write(&f, "echo", &[1, 1, 4, 2], &v).unwrap();
        assert_eq!(read(&f, "echo").unwrap(), v);
        assert_eq!(f.dataset("echo_iq_i16").unwrap().filters().len(), 2);
        drop(f);
        std::fs::remove_file(path).unwrap();
    }
}
