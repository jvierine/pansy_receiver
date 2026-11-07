"""Digital Metadata reads that never suppress I/O failures or delete files."""
import h5py
import digital_rf


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
