class NeurobridgeRefs {

    static def hm3Snplist(projectDir) {
        file("${projectDir}/ref/ldsc/w_hm3.snplist")
    }

    static def ldChrDir(projectDir) {
        file("${projectDir}/ref/ldsc/eur_w_ld_chr")
    }

    static def wldDir(projectDir) {
        file("${projectDir}/ref/ldsc/weights_hm3_no_hla")
    }
}