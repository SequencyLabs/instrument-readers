//! Pure acquisition-readiness validation for Waters MassLynx `.raw` directories.

use std::collections::BTreeMap;

use crate::binary::file_name;

use super::metadata::parse_func_number;
use super::{FUNCTNS_RECORD_LEN, pair_has_complete_scan};

/// Pure metadata view of one file inside a candidate Waters `.raw` directory.
#[derive(Debug, Clone, Copy)]
pub struct WatersBundleFile<'a> {
    pub relative_path: &'a str,
    pub size_bytes: u64,
    /// Contents of `_FUNC*.IDX` files (only needed for `.IDX` files; `None` for other files).
    pub idx_bytes: Option<&'a [u8]>,
}

#[derive(Default)]
struct WatersDirGroup<'a> {
    header_size: Option<u64>,
    extern_size: Option<u64>,
    functns_size: Option<u64>,
    dats: BTreeMap<u32, u64>,
    idxs: BTreeMap<u32, (u64, Option<&'a [u8]>)>,
}

impl WatersDirGroup<'_> {
    fn has_waters_files(&self) -> bool {
        self.header_size.is_some()
            || self.extern_size.is_some()
            || self.functns_size.is_some()
            || !self.dats.is_empty()
            || !self.idxs.is_empty()
    }

    fn is_ready(&self) -> bool {
        if self.header_size.is_none_or(|sz| sz == 0) {
            return false;
        }
        if self.extern_size.is_some_and(|sz| sz == 0) {
            return false;
        }
        if self.dats.is_empty() || self.dats.len() != self.idxs.len() {
            return false;
        }
        if let Some(functns_size) = self.functns_size {
            let record_len = FUNCTNS_RECORD_LEN as u64;
            if functns_size == 0
                || !functns_size.is_multiple_of(record_len)
                || functns_size / record_len != self.dats.len() as u64
            {
                return false;
            }
        }
        let mut saw_non_empty = false;
        for (&func_num, &dat_size) in &self.dats {
            let Some(&(idx_size, idx_bytes)) = self.idxs.get(&func_num) else {
                return false;
            };
            if dat_size == 0 && idx_size == 0 {
                continue;
            }
            let Ok(dat_len) = usize::try_from(dat_size) else {
                return false;
            };
            let Some(bytes) = idx_bytes else {
                return false;
            };
            if !pair_has_complete_scan(bytes, dat_len) {
                return false;
            }
            saw_non_empty = true;
        }
        saw_non_empty
    }
}

fn parent_dir(path: &str) -> &str {
    path.rsplit_once('/').map_or("", |(parent, _)| parent)
}

/// Pure canonical readiness check for Waters MassLynx `.raw` bundles (including nested `.raw`
/// directories).
///
/// Requires every Waters directory in `files` to have a non-empty `_HEADER.TXT`, a non-empty
/// `_EXTERN.INF` when present, a `_FUNCTNS.INF` whose size is a non-zero multiple of 416 bytes
/// matching the `_FUNC*.DAT` count when present, matched `_FUNC*.DAT`/`_FUNC*.IDX` pairs, and
/// every non-empty pair passing [`pair_has_complete_scan`] (with at least one non-empty pair).
#[must_use]
pub fn is_bundle_ready(files: &[WatersBundleFile<'_>]) -> bool {
    let mut groups: BTreeMap<&str, WatersDirGroup<'_>> = BTreeMap::new();
    for file in files {
        let dir = parent_dir(file.relative_path);
        let upper = file_name(file.relative_path).to_ascii_uppercase();
        if upper == "_HEADER.TXT" {
            groups.entry(dir).or_default().header_size = Some(file.size_bytes);
        } else if upper == "_EXTERN.INF" {
            groups.entry(dir).or_default().extern_size = Some(file.size_bytes);
        } else if upper == "_FUNCTNS.INF" {
            groups.entry(dir).or_default().functns_size = Some(file.size_bytes);
        } else if let Some(num) = parse_func_number(&upper, ".DAT") {
            groups
                .entry(dir)
                .or_default()
                .dats
                .insert(num, file.size_bytes);
        } else if let Some(num) = parse_func_number(&upper, ".IDX") {
            groups
                .entry(dir)
                .or_default()
                .idxs
                .insert(num, (file.size_bytes, file.idx_bytes));
        }
    }
    let waters_groups: Vec<&WatersDirGroup<'_>> =
        groups.values().filter(|g| g.has_waters_files()).collect();
    !waters_groups.is_empty() && waters_groups.iter().all(|g| g.is_ready())
}
