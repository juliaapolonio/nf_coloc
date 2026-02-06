process COLOC {
  """
  Perform colocalization
  """

  label 'process_medium'

  container "${ workflow.containerEngine == 'singularity' ? 'docker://juliaapolonio/coloc:5.2.3dev':
            'docker.io/juliaapolonio/coloc:5.2.3dev' }"

  input:
    path qtl
    path sumstats
    path windows

  output:
    path "coloc_summary.csv", emit: colocalization
    path "*regional.png", optional: true, emit: regional_plot

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
    #!/bin/bash
    nf_coloc.R ${sumstats} ${qtl} ${windows}
    """
}