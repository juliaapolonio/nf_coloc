/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { WINDOWS } from "./modules/local/extract_windows/main.nf"
include { COLOC } from "./modules/local/colocalization/main.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    sumstats = file(params.sumstats)
    qtl_file = file(params.qtl_file)

    WINDOWS( sumstats )

    ch_windows = WINDOWS.out.windows
    .splitCsv(header:true, sep:'\t')
    .map { row -> tuple(row.seqnames, row.start, row.end) }

    COLOC( qtl_file, sumstats, ch_windows )

}
