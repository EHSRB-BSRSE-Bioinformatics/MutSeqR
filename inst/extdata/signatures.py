import sys
import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerExtractor import sigpro as sig

def cosmic_fitR(samples, output):
  Analyze.cosmic_fit(samples,
                       output,
                       genome_build="GRCh38",
                       signatures=None,
                       cosmic_version=3.3,
                       verbose=True,
                       exome=False)

def sigProfilerExtractorR(sample_file, output_dir, reference_genome):
  sig.sigProfilerExtractor(input_data = sample_file,
                       output = output_dir,
                       input_type = "matrix",
                       opportunity_genome = reference_genome,
                       cosmic_version=3.3,
                       minimum_signatures=1,
                       maximum_signatures=10,
                       nmf_replicates =100)

