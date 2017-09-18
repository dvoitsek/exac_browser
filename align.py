INFO = {
    "AC":"AC",
    "AF":"AF",
    "AN":"AN",
    "BaseQRankSum":"BaseQRankSum",
    "ClippingRankSum":"ClippingRankSum",
    "DB":"DB",
    "DP":"DP",
    "DS":"DS",
    "END":"END",
    "ExcessHet":"ExcessHet",
    "FS":"FS",
    "HaplotypeScore":"HaplotypeScore",
    "InbreedingCoeff":"InbreedingCoeff",
    "MLEAC":"MLEAC",
    "MLEAF":"MLEAF",
    "MQ":"MQ",
    "MQRankSum":"MQRankSum",
    "QD":"QD",
    "RAW_MQ":"RAW_MQ",
    "ReadPosRankSum":"ReadPosRankSum",
    "SOR":"SOR",
    "ANN":"ANN",
    "LOF":"LOF",
    "NMD":"NMD",
    "dbNSFP_ExAC_NFE_AF":"dbNSFP_ExAC_NFE_AF",
    "dbNSFP_ExAC_SAS_AF":"dbNSFP_ExAC_SAS_AF",
    "dbNSFP_GERP___RS":"dbNSFP_GERP___RS",
    "dbNSFP_GERP___NR":"dbNSFP_GERP___NR",
    "dbNSFP_1000Gp1_AMR_AF":"dbNSFP_1000Gp1_AMR_AF",
    "dbNSFP_ExAC_Adj_AC":"dbNSFP_ExAC_Adj_AC",
    "dbNSFP_ExAC_Adj_AF":"dbNSFP_ExAC_Adj_AF",
    "dbNSFP_ExAC_SAS_AC":"dbNSFP_ExAC_SAS_AC",
    "dbNSFP_MetaSVM_pred":"dbNSFP_MetaSVM_pred",
    "dbNSFP_Interpro_domain":"dbNSFP_Interpro_domain",
    "dbNSFP_FATHMM_pred":"dbNSFP_FATHMM_pred",
    "dbNSFP_ExAC_AFR_AF":"dbNSFP_ExAC_AFR_AF",
    "dbNSFP_ExAC_AFR_AC":"dbNSFP_ExAC_AFR_AC",
    "dbNSFP_ExAC_AF":"dbNSFP_ExAC_AF",
    "dbNSFP_1000Gp1_AF":"dbNSFP_1000Gp1_AF",
    "dbNSFP_Uniprot_acc":"dbNSFP_Uniprot_acc",
    "dbNSFP_ExAC_AC":"dbNSFP_ExAC_AC",
    "dbNSFP_LRT_pred":"dbNSFP_LRT_pred",
    "dbNSFP_PROVEAN_pred":"dbNSFP_PROVEAN_pred",
    "dbNSFP_ExAC_FIN_AC":"dbNSFP_ExAC_FIN_AC",
    "dbNSFP_phastCons100way_vertebrate":"dbNSFP_phastCons100way_vertebrate",
    "dbNSFP_ExAC_FIN_AF":"dbNSFP_ExAC_FIN_AF",
    "dbNSFP_CADD_phred":"dbNSFP_CADD_phred",
    "dbNSFP_Polyphen2_HDIV_pred":"dbNSFP_Polyphen2_HDIV_pred",
    "dbNSFP_1000Gp1_ASN_AF":"dbNSFP_1000Gp1_ASN_AF",
    "dbNSFP_1000Gp1_AFR_AF":"dbNSFP_1000Gp1_AFR_AF",
    "dbNSFP_ExAC_AMR_AF":"dbNSFP_ExAC_AMR_AF",
    "dbNSFP_MutationTaster_pred":"dbNSFP_MutationTaster_pred",
    "dbNSFP_1000Gp1_EUR_AF":"dbNSFP_1000Gp1_EUR_AF",
    "dbNSFP_MutationAssessor_pred":"dbNSFP_MutationAssessor_pred",
    "dbNSFP_ESP6500_AA_AF":"dbNSFP_ESP6500_AA_AF",
    "dbNSFP_Polyphen2_HVAR_pred":"dbNSFP_Polyphen2_HVAR_pred",
    "dbNSFP_ExAC_AMR_AC":"dbNSFP_ExAC_AMR_AC",
    "dbNSFP_ExAC_NFE_AC":"dbNSFP_ExAC_NFE_AC",
    "dbNSFP_SIFT_pred":"dbNSFP_SIFT_pred",
    "dbNSFP_ExAC_EAS_AC":"dbNSFP_ExAC_EAS_AC",
    "dbNSFP_ExAC_EAS_AF":"dbNSFP_ExAC_EAS_AF",
    "dbNSFP_ESP6500_EA_AF":"dbNSFP_ESP6500_EA_AF"
}

ANN = {
    "Allele":"Allele",
    "Annotation":"Annotation",
    "Annotation_Impact":"Annotation_Impact",
    "Gene_Name":"Gene_Name",
    "Gene_ID":"Gene_ID",
    "Feature_Type":"Feature_Type",
    "Feature_ID":"Feature_ID",
    "Transcript_BioType":"Transcript_BioType",
    "Rank":"Rank",
    "HGVS.c":"HGVS_c",
    "HGVS.p":"HGVS_p",
    "cDNA.pos / cDNA.length":"cDNA_pos___cDNA_length",
    "CDS.pos / CDS.length":"CDS_pos___CDS_length",
    "AA.pos / AA.length":"AA_pos___AA_length",
    "Distance":"Distance",
    "ERRORS / WARNINGS / INFO":"ERRORS___WARNINGS___INFO"
}

LOF = {
    "Gene_Name":"Gene_Name",
    "Gene_ID":"Gene_ID",
    "Number_of_transcripts_in_gene":"Number_of_transcripts_in_gene",
    "Percent_of_transcripts_affected":"Percent_of_transcripts_affected"
}

NMD = {
    "Gene_Name":"Gene_Name",
    "Gene_ID":"Gene_ID",
    "Number_of_transcripts_in_gene":"Number_of_transcripts_in_gene",
    "Percent_of_transcripts_affected":"Percent_of_transcripts_affected"
}

TYPE_CONVERT = {
    "AC":int,
    "AF":float,
    "AN":int,
    "BaseQRankSum":float,
    "ClippingRankSum":float,
    "DB":None,
    "DP":int,
    "DS":None,
    "END":int,
    "ExcessHet":float,
    "FS":float,
    "HaplotypeScore":float,
    "InbreedingCoeff":float,
    "MLEAC":int,
    "MLEAF":float,
    "MQ":float,
    "MQRankSum":float,
    "QD":float,
    "RAW_MQ":float,
    "ReadPosRankSum":float,
    "SOR":float,
    "ANN":None,
    "LOF":None,
    "NMD":None,
    "dbNSFP_ExAC_NFE_AF":int,
    "dbNSFP_ExAC_SAS_AF":int,
    "dbNSFP_GERP___RS":float,
    "dbNSFP_GERP___NR":float,
    "dbNSFP_1000Gp1_AMR_AF":None,
    "dbNSFP_ExAC_Adj_AC":int,
    "dbNSFP_ExAC_Adj_AF":float,
    "dbNSFP_ExAC_SAS_AC":int,
    "dbNSFP_MetaSVM_pred":None,
    "dbNSFP_Interpro_domain":None,
    "dbNSFP_FATHMM_pred":None,
    "dbNSFP_ExAC_AFR_AF":float,
    "dbNSFP_ExAC_AFR_AC":int,
    "dbNSFP_ExAC_AF":float,
    "dbNSFP_1000Gp1_AF":None,
    "dbNSFP_Uniprot_acc":None,
    "dbNSFP_ExAC_AC":int,
    "dbNSFP_LRT_pred":None,
    "dbNSFP_PROVEAN_pred":None,
    "dbNSFP_ExAC_FIN_AC":int,
    "dbNSFP_phastCons100way_vertebrate":float,
    "dbNSFP_ExAC_FIN_AF":float,
    "dbNSFP_CADD_phred":float,
    "dbNSFP_Polyphen2_HDIV_pred":None,
    "dbNSFP_1000Gp1_ASN_AF":None,
    "dbNSFP_1000Gp1_AFR_AF":None,
    "dbNSFP_ExAC_AMR_AF":int,
    "dbNSFP_MutationTaster_pred":None,
    "dbNSFP_1000Gp1_EUR_AF":None,
    "dbNSFP_MutationAssessor_pred":None,
    "dbNSFP_ESP6500_AA_AF":None,
    "dbNSFP_Polyphen2_HVAR_pred":None,
    "dbNSFP_ExAC_AMR_AC":int,
    "dbNSFP_ExAC_NFE_AC":int,
    "dbNSFP_SIFT_pred":None,
    "dbNSFP_ExAC_EAS_AC":int,
    "dbNSFP_ExAC_EAS_AF":int,
    "dbNSFP_ESP6500_EA_AF":None,
    # Annotations (ANN)
    "Allele":None,
    "Annotation":None,
    "Annotation_Impact":None,
    "Gene_Name":None,
    "Gene_ID":None,
    "Feature_Type":None,
    "Feature_ID":None,
    "Transcript_BioType":None,
    "Rank":None,
    "HGVS.c":None,
    "HGVS.p":None,
    "cDNA.pos / cDNA.length":None,
    "CDS.pos / CDS.length":None,
    "AA.pos / AA.length":None,
    "Distance":None,
    "ERRORS / WARNINGS / INFO":None,
    # LOF & NMD
    "Gene_Name":None,
    "Gene_ID":None,
    "Number_of_transcripts_in_gene":None,
    "Percent_of_transcripts_affected":None
}
