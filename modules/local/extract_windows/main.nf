process WINDOWS {
  """
  Extract colocalization windows from a given sumstats
  """

  label 'process_medium'

  container "${ workflow.containerEngine == 'singularity' ? 'docker://juliaapolonio/coloc:5.2.3dev':
            'docker.io/juliaapolonio/coloc:5.2.3dev' }"

  input:
    path sumstats

  output:
    path "windows.tsv", emit: windows

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
    #!/bin/bash
    extract_windows_nf.R ${sumstats}
    """
}
