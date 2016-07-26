#https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/qc/fastqc.py

import os,sys

import pandas as pd


try:
    from fadapa import Fadapa
except ImportError:
    Fadapa = None

class FastQCParser:
    def __init__(self, base_dir, sample=None):
        self._dir = base_dir
        self.sample = sample

    def get_fastqc_summary(self):
        ignore = set(["Total Sequences", "Filtered Sequences",
                      "Filename", "File type", "Encoding"])
        stats = {}
        for stat_line in self._fastqc_data_section("Basic Statistics")[1:]:
            k, v = stat_line.split("\t")[:2]
            if k not in ignore:
                stats[k] = v
        return stats

    def _fastqc_data_section(self, section_name):
        out = []
        in_section = False
        data_file = os.path.join(self._dir, "fastqc_data.txt")
        if os.path.exists(data_file):
            with open(data_file) as in_handle:
                for line in in_handle:
                    if line.startswith(">>%s" % section_name):
                        in_section = True
                    elif in_section:
                        if line.startswith(">>END"):
                            break
                        out.append(line.rstrip("\r\n"))
        return out

    def save_sections_into_file(self):

        data_file = os.path.join(self._dir, "fastqc_data.txt")
        if os.path.exists(data_file) and Fadapa:
            parser = Fadapa(data_file)
            module = [m[1] for m in parser.summary()][2:9]
            for m in module:
                out_file = os.path.join(self._dir, m.replace(" ", "_") + ".tsv")
                dt = self._get_module(parser, m)
                dt.to_csv(out_file, sep="\t", index=False)

    def _get_module(self, parser, module):
        """
        Get module using fadapa package
        """
        dt = []
        lines = parser.clean_data(module)
        header = lines[0]
        for data in lines[1:]:
            if data[0].startswith("#"):  # some modules have two headers
                header = data
                continue
            if data[0].find("-") > -1:  # expand positions 1-3 to 1, 2, 3
                f, s = map(int, data[0].split("-"))
                for pos in range(f, s):
                    dt.append([str(pos)] + data[1:])
            else:
                dt.append(data)
        dt = pd.DataFrame(dt)
        dt.columns = [h.replace(" ", "_") for h in header]
        dt['sample'] = self.sample
        return dt

if __name__ == "__main__":
	fastqc_outdir = sys.argv[1]
	sample = sys.argv[2]
	
	parser = FastQCParser(fastqc_outdir, sample)
	stats = parser.get_fastqc_summary()
	parser.save_sections_into_file()


