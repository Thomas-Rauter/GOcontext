.go_cc <- NULL
.go_cc <- load_go(ont = "CC")

.go_cc_ecoli <- NULL
if (requireNamespace("org.EcK12.eg.db", quietly = TRUE)) {
    .org_ecoli <- org.EcK12.eg.db::org.EcK12.eg.db

    .go_cc_ecoli <- GOcontext::attach_org(
        go = .go_cc,
        OrgDb = .org_ecoli,
        keytype = "ENTREZID"
    )
}
