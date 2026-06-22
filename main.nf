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

    clump_file = file(params.clump_file)
    sumstats = file(params.sumstats)

    ch_samplesheet = Channel.fromPath(params.input)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            if( row.sample && row.qtl_file ) {
		tuple( [ id: row.sample ], file(row.qtl_file) )
            }
        }

    WINDOWS( clump_file )

    ch_windows = WINDOWS.out.windows
    .splitCsv(header:true, sep:'\t')
    .map { row -> tuple(row.chr, row.start, row.end) }

    ch_coloc_input = ch_samplesheet.combine(ch_windows)
        .map { meta, qtl, chr, start, end -> 
            tuple( meta, qtl, chr, start, end ) 
        }
    
    COLOC( sumstats, ch_coloc_input )

}
