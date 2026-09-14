"""L1, the lung identity covariate, fixed in PREREGISTRATION 6y before the
TCGA-LUAD matrix was obtained. Not to be modified."""
L1_EFFECTORS = ["SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "SFTPD", "NAPSA",
                "SLC34A2", "PGC", "LAMP3", "ABCA3", "AGER", "CLDN18",
                "SCGB1A1", "SCGB3A2", "CTSH", "SFTA2", "CLIC5"]
L1_REGULATORS = ["NKX2-1", "FOXJ1", "HOPX", "CEBPA", "ETV5"]
L1 = L1_EFFECTORS + L1_REGULATORS
