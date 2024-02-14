import sys
from SigProfilerAssignment import Analyzer as Analyze
 
def cosmic_fit_MutSeqR(samples,
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

