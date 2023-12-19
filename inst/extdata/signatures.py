import sys
import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerExtractor import sigpro as sig
from SigProfilerMatrixGenerator import install as genInstall

def install_genome(genome, custom, bash):
  genInstall.install(genome, rsync=False, bash=False, ftp=True, custom = False)

def cosmic_fit_DupSeqR(samples,
                       output,
                       genome_build="GRCh38",
                       signatures=None,
                       cosmic_version=3.3,
                       verbose=True,
                       exome=False):
  Analyze.cosmic_fit(samples,
                       output,
                       genome_build="GRCh38",
                       signatures=None,
                       cosmic_version=3.3,
                       verbose=True,
                       exome=False)

def sigProfilerExtractorDupSeqR(sample_file, output_dir, reference_genome):
  sig.sigProfilerExtractor(input_data = sample_file,
                       output = output_dir,
                       input_type = "matrix",
                       opportunity_genome = reference_genome,
                       cosmic_version=3.3,
                       minimum_signatures=1,
                       maximum_signatures=10,
                       nmf_replicates=100)

