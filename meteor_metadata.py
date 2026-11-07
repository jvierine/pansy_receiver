"""Digital Metadata reads that never suppress I/O failures or delete files."""
import h5py
import digital_rf
import numpy as np

MF_FIELDS=("beam_pos_idx","tx_idxs","tx_pwr","max_snr","max_range","max_dopvel","noise_floor","tx_std")


def write_complete_record(writer,root,key,payload):
    """Caller serializes producers. Keep complete duplicates; replay partial rows."""
    try:
        writer.write([int(key)],payload)
        return True
    except (ValueError,OSError) as exc:
        if "already in data" not in str(exc): raise
    reader=StrictDigitalMetadataReader(str(root))
    required=set(payload)
    for path in reader._get_file_list(key,key):
        with h5py.File(path,"a") as stream:
            name=str(int(key))
            if name not in stream: continue
            if required.issubset(stream[name].keys()): return False
            del stream[name]
            stream.flush()
            print(f"Replaying incomplete metadata sample {key} in {root}",flush=True)
        writer.write([int(key)],payload)
        return True
    raise RuntimeError(f"duplicate metadata sample {key} could not be verified")


def matched_filter_record_exists(reader,key):
    """Caller holds the exclusive MF lock. Recompute interrupted partial rows."""
    def complete(row):
        return all(name in row and np.asarray(row[name]).size==20 for name in MF_FIELDS) and np.array_equal(
            row["tx_idxs"],int(key)+np.arange(20)*1600)
    try:
        rows=reader.read(key,key,MF_FIELDS)
        if not rows: return False
        if complete(rows[key]): return True
    except KeyError:
        pass
    # Remove only the demonstrably incomplete sample, never its HDF5 file or
    # valid neighbours. The producer still has the raw batch and will replay it.
    for path in reader._get_file_list(key,key):
        with h5py.File(path,"a") as stream:
            name=str(int(key))
            if name not in stream: continue
            group=stream[name]
            values={field:group[field][()] for field in MF_FIELDS if field in group}
            if complete(values): return True
            del stream[name]
            stream.flush()
            print(f"Replaying incomplete MF metadata sample {key}",flush=True)
    return False


class StrictDigitalMetadataReader(digital_rf.DigitalMetadataReader):
    def _add_metadata(self, result, path, columns, sample0, sample1, is_edge):
        # The upstream implementation catches IOError, sometimes deletes the
        # file, and returns a partial result. A completion-based pipeline must
        # instead leave its checkpoint unchanged and retry the entire window.
        with h5py.File(path, "r") as stream:
            for sample in sorted(int(key) for key in stream):
                if is_edge and not sample0 <= sample <= sample1:
                    continue
                record=stream[str(sample)]
                if columns is None:
                    self._populate_data(result,record,sample)
                elif isinstance(columns,str):
                    self._populate_data(result,record[columns],sample)
                else:
                    row={}
                    for name in columns:
                        self._populate_data(row,record[name],name)
                    result[sample]=row
