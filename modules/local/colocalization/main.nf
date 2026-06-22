process COLOC {
  """
  Perform colocalization
  """

  tag "${meta.id} - ${chr}:${start}-${end}"

  label 'process_medium'

  container "${ workflow.containerEngine == 'singularity' ? 'docker://juliaapolonio/coloc:5.2.3dev':
            'docker.io/juliaapolonio/coloc:5.2.3dev' }"

  input:
    path sumstats
    tuple val(meta), path(qtl), val(chr), val(start), val(end)

  output:
    path "${meta.id}_coloc_summary_*.csv", emit: colocalization
    path "${meta.id}_*regional.png", optional: true, emit: regional_plot

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
    #!/bin/bash
    nf_coloc.R ${qtl} ${sumstats} ${chr} ${start} ${end} ${meta.id}
    """
}
