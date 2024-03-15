import sys
from SigProfilerAssignment import Analyzer as Analyze
 
def cosmic_fit_MutSeqR(samples,
                       output,
                       input_type="matrix",
                       context_type="96",
                       cosmic_version=3.3,
                       exome = False,
                       genome_build="GRCh38",
                       signature_database=None,
                       exclude_signature_subgroups=None,
                       export_probabilities=True,
                       export_probabilities_per_mutation=False,
                       make_plots=True,
                       sample_reconstruction_plots = "png",
                       verbose=True):
  Analyze.cosmic_fit(samples,
                       output,
                       input_type=input_type,
                       context_type=context_type,
                       cosmic_version=cosmic_version,
                       exome = exome,
                       genome_build=genome_build,
                       signature_database=signature_database,
                       exclude_signature_subgroups=exclude_signature_subgroups,
                       export_probabilities=export_probabilities,
                       export_probabilities_per_mutation=export_probabilities_per_mutation,
                       make_plots=make_plots,
                       sample_reconstruction_plots = sample_reconstruction_plots,
                       verbose=verbose
                       )

