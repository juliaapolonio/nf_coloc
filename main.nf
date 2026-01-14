/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { WINDOWS } from "./modules/local/extract_windows/main.nf"
//include { COLOC } from "./modules/local/colocalization/main.nf"
//include { REPORT } from "./modules/local/locuszoom_report/main.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    sumstats = file(params.sumstats)
    //qtl_file = file(params.qtl_file)

    WINDOWS( sumstats )

    //WINDOWS.out.windows
    //.flatten()
    //.set { ext_win }

    //COLOC( ext_win, qtl_file )

    //COLOC.out.coloc_result
    //.collect()
    //.set { result }

    //REPORT( result ) 

}
