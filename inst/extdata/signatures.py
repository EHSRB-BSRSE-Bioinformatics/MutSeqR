from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerExtractor import sigpro as sig

 
def cosmic_fit_MutSeqR(samples,
                       output,
                       input_type="matrix",
                       context_type="96",
                       cosmic_version=3.4,
                       exome=False,
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

def denovo_fit_MutSeqR(samples,
                       output,
                       input_type="matrix",
                       signatures="path/to/input/denovo/signatures/file",
                       genome_build="GRCh38"):
    Analyze.denovo_fit(samples,
                        output,
                        input_type=input_type,
                        signatures=signatures,
                        genome_build=genome_build)

def decompose_fit_MutSeqR(samples,
                       output,
                       input_type="matrix",
                       signatures="path/to/input/denovo/signatures/file",
                       signature_database="path/to/optional/reference/signatures/database/file",
                       genome_build="GRCh38"):
    Analyze.decompose_fit(samples,
                        output,
                        input_type=input_type,
                        signatures=signatures,
                        signature_database=signature_database,
                        genome_build=genome_build)

