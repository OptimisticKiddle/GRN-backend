from pydantic import BaseModel

'''
分页
'''


# 分页paging
class Paging(BaseModel):
    start: int = 0
    length: int = 10

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


'''
overall_data
'''


# 总体dataset
class OverallData(BaseModel):
    id: int
    pb_gene: str
    pb_ensembl: str
    n_sample: int
    celline: str
    method: str
    conditions: str
    accession: str
    datasource: str

    class Config:
        orm_mode = True


# 无id版本的 overall data
class OverallDataNoId(BaseModel):
    pb_gene: str
    pb_ensembl: str
    n_sample: int
    celline: str
    method: str
    conditions: str
    accession: str
    datasource: str

    class Config:
        orm_mode = True


# browser页filter
class BrowserFilter(BaseModel):
    id: int = None
    pb_gene: str = None
    pb_ensembl: str = None
    n_sample_greater: int = None
    n_sample_less: int = None
    celline: str = None
    method: str = None
    datasource: str = None


# 单个数据 resp
class SingleDataResponse(BaseModel):
    data: OverallDataNoId

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


'''
diff peaks
'''


class DiffPeaks(BaseModel):
    # id: int  # 不需要好像？
    peak_id: int
    chr: str
    start: int
    end: int
    width: int
    conc: float
    conc_ctrl: float
    conc_treat: float
    fold: float
    p_value: float
    FDR: float
    annotation: str
    gene_chr: int
    gene_start: int
    gene_end: int
    gene_len: int
    gene_strand: int
    gene_id: str
    transcript_id: str
    distance_to_tss: int

    class Config:
        orm_mode = True


# 数值范围（最大最小值）
class DiffPeaksDataRange(BaseModel):
    start_min: int
    start_max: int
    end_min: int
    end_max: int
    width_min: int
    width_max: int
    conc_min: float = 0
    conc_max: float
    conc_ctrl_min: float = 0
    conc_ctrl_max: float
    conc_treat_min: float = 0
    conc_treat_max: float
    fold_min: float
    fold_max: float
    p_value_min: float = 0
    p_value_max: float
    FDR_min: float = 0
    FDR_max: float
    gene_start_min: int
    gene_start_max: int
    gene_end_min: int
    gene_end_max: int
    gene_len_min: int
    gene_len_max: int
    distance_to_tss_min: int
    distance_to_tss_max: int

    class Config:
        orm_mode = True


# diff peaks 请求体
class DiffPeaksRequest(BaseModel):
    class Filter(BaseModel):
        chr: str = None
        start_min: int = None
        start_max: int = None
        end_min: int = None
        end_max: int = None
        width_min: int = None
        width_max: int = None
        conc_min: float = None
        conc_max: float = None
        conc_ctrl_min: float = None
        conc_ctrl_max: float = None
        conc_treat_min: float = None
        conc_treat_max: float = None
        fold_min: float = None
        fold_max: float = None
        p_value_min: float = None
        p_value_max: float = None
        FDR_min: float = None
        FDR_max: float = None
        gene_chr: int = None
        gene_start_min: int = None
        gene_start_max: int = None
        gene_end_min: int = None
        gene_end_max: int = None
        gene_len_min: int = None
        gene_len_max: int = None
        gene_id: str = None
        transcript_id: str = None
        distance_to_tss_min: int = None
        distance_to_tss_max: int = None

    class Sequence(BaseModel):  # 1 代表升序，-1代表降序，仅传一个参数
        start_seq: int = None
        end_seq: int = None
        width_seq: int = None
        conc_seq: int = None
        conc_ctrl_seq: int = None
        conc_treat_seq: int = None
        fold_seq: int = None
        p_value_seq: int = None
        FDR_seq: int = None
        gene_start_seq: int = None
        gene_end_seq: int = None
        gene_len_seq: int = None
        distance_to_tss_seq: int = None

    id: int
    filter: Filter = None
    seq: Sequence = None
    paging: Paging = None


'''
diff_go_enrich  、 diff kegg enrich
'''


class DiffGOEnrich(BaseModel):
    GO_id: str
    description: str
    gene_ratio: str
    bg_ratio: str
    p_value: float
    p_adjust: float
    q_value: float
    gene_ids: str
    count: int

    class Config:
        orm_mode = True


class DiffKEGGEnrich(BaseModel):
    # cluster: str
    KEGG_id: str
    description: str
    gene_ratio: str
    bg_ratio: str
    p_value: float
    p_adjust: float
    q_value: float
    gene_ids: str
    count: int

    class Config:
        orm_mode = True


# diff go/kegg enrich 请求体
class DiffEnrichRequest(BaseModel):
    class Sequence(BaseModel):
        p_value_seq: int = None
        p_adjust_seq: int = None
        q_value_seq: int = None
        count_seq: int = None

    id: int
    seq: Sequence = None
    paging: Paging = None


'''
diff motif
'''


class DiffMotif(BaseModel):
    rank: int
    motif: str
    name: str
    consensus: str
    p_value: float
    log_p_value: float
    q_value: float
    target_sequences: float
    target_sequences_percent: float
    bg_sequences: float
    bg_sequences_percent: float

    class Config:
        orm_mode = True


# diff motif  请求体
class DiffMotifRequest(BaseModel):
    class Sequence(BaseModel):
        p_value_seq: int = None
        q_value_seq: int = None
        target_sequences_seq: int = None
        bg_sequences_seq: int = None

    id: int
    seq: Sequence = None
    paging: Paging = None


'''
diff footprint
'''


class DiffFootprint(BaseModel):
    motif: str
    tf: str
    num: int
    protection_score_ctrl: float
    protection_score_treat: float
    tc_ctrl: float
    tc_treat: float
    tf_activity: float
    z_score: float
    p_value: float

    class Config:
        orm_mode = True


# 数值范围（最大最小值）
class DiffFootprintDataRange(BaseModel):
    num_min: int
    num_max: int
    protection_score_ctrl_min: float = 0
    protection_score_ctrl_max: float
    protection_score_treat_min: float = 0
    protection_score_treat_max: float
    tc_ctrl_min: float = 0
    tc_ctrl_max: float
    tc_treat_min: float = 0
    tc_treat_max: float
    tf_activity_min: float
    tf_activity_max: float
    z_score_min: float
    z_score_max: float
    p_value_min: float = 0
    p_value_max: float

    class Config:
        orm_mode = True


#  请求体
class DiffFootprintRequest(BaseModel):
    class Filter(BaseModel):
        motif: str = None
        tf: str = None
        num_min: int = None
        num_max: int = None
        protection_score_ctrl_min: float = None
        protection_score_ctrl_max: float = None
        protection_score_treat_min: float = None
        protection_score_treat_max: float = None
        tc_ctrl_min: float = None
        tc_ctrl_max: float = None
        tc_treat_min: float = None
        tc_treat_max: float = None
        tf_activity_min: float = None
        tf_activity_max: float = None
        z_score_min: float = None
        z_score_max: float = None
        p_value_min: float = None
        p_value_max: float = None

    class Sequence(BaseModel):  # 1 代表升序，-1代表降序，仅传一个参数
        num_seq: int = None
        protection_score_ctrl_seq: int = None
        protection_score_treat_seq: int = None
        tc_ctrl_seq: int = None
        tc_treat_seq: int = None
        tf_activity_seq: int = None
        z_score_seq: int = None
        p_value_seq: int = None

    id: int
    filter: Filter = None
    seq: Sequence = None
    paging: Paging = None


'''
ctrl/treat peaks
'''


class Peaks(BaseModel):
    peak_id: int
    chr: str
    start: int
    end: int
    width: int
    score: int
    signal_value: float
    log_p_value: float
    log_q_value: float
    peak_offset: int
    annotation: str
    gene_chr: int
    gene_start: int
    gene_end: int
    gene_len: int
    gene_strand: int
    gene_id: str
    transcript_id: str
    distance_to_tss: int

    class Config:
        orm_mode = True


# 数值属性列范围（最大最小值）
class PeaksDataRange(BaseModel):
    start_min: int
    start_max: int
    end_min: int
    end_max: int
    width_min: int
    width_max: int
    score_min: int
    score_max: int
    signal_value_min: float
    signal_value_max: float
    log_p_value_min: float
    log_p_value_max: float
    log_q_value_min: float
    log_q_value_max: float
    peak_offset_min: int
    peak_offset_max: int
    gene_start_min: int
    gene_start_max: int
    gene_end_min: int
    gene_end_max: int
    gene_len_min: int
    gene_len_max: int
    distance_to_tss_min: int
    distance_to_tss_max: int

    class Config:
        orm_mode = True


#  请求体
class PeaksRequest(BaseModel):
    class Filter(BaseModel):
        chr: str = None
        start_min: int = None
        start_max: int = None
        end_min: int = None
        end_max: int = None
        width_min: int = None
        width_max: int = None
        score_min: int = None
        score_max: int = None
        signal_value_min: float = None
        signal_value_max: float = None
        log_p_value_min: float = None
        log_p_value_max: float = None
        log_q_value_min: float = None
        log_q_value_max: float = None
        peak_offset_min: int = None
        peak_offset_max: int = None
        gene_chr: int = None
        gene_start_min: int = None
        gene_start_max: int = None
        gene_end_min: int = None
        gene_end_max: int = None
        gene_len_min: int = None
        gene_len_max: int = None
        gene_id: str = None
        transcript_id: str = None
        distance_to_tss_min: int = None
        distance_to_tss_max: int = None

    class Sequence(BaseModel):  # 1 代表升序，-1代表降序，仅传一个参数
        start_seq: int = None
        end_seq: int = None
        width_seq: int = None
        score_seq: int = None
        signal_value_seq: int = None
        log_p_value_seq: int = None
        log_q_value_seq: int = None
        peak_offset_seq: int = None
        gene_start_seq: int = None
        gene_end_seq: int = None
        gene_len_seq: int = None
        distance_to_tss_seq: int = None

    id: int
    filter: Filter = None
    seq: Sequence = None
    paging: Paging = None


'''
通用 resp
'''


# 分页表格 resp
class TableDataResponse(BaseModel):
    data: list[OverallData] | list[DiffPeaks] | list[DiffGOEnrich] | list[DiffKEGGEnrich] | list[DiffMotif] \
          | list[DiffFootprint] | list[Peaks] | list[object]
    records_sum: int

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


# 枚举值 resp
class EnumDataResp(BaseModel):
    data: list[str] | list[object]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


'''
statistics resp
'''


class StatisticsResponse(BaseModel):
    class Statistics(BaseModel):
        perturb_datasets_num: int
        perturb_gene_num: int
        cell_line_num: int
        n_samples_num: int

    class DiffData(BaseModel):
        diff_peaks_num: int
        diff_motif_num: int
        diff_GO_enrich_num: int
        diff_KEGG_enrich_num: int
        diff_footprint_num: int

    class EchartsPieData(BaseModel):
        name: str
        value: int

    # class DiffPeakNumInStruct(BaseModel):
    #     data: list[EchartsPieData]
    #
    #     promoter: int = 0
    #     UTR_5: int = 0
    #     UTR_3: int = 0
    #     exon: int = 0
    #     intron: int = 0
    #     downstream: int = 0
    #     distal_intergenic: int = 0

    class DiffPeaksNumInChr(BaseModel):
        chr1: int = 0
        chr2: int = 0
        chr3: int = 0
        chr4: int = 0
        chr5: int = 0
        chr6: int = 0
        chr7: int = 0
        chr8: int = 0
        chr9: int = 0
        chr10: int = 0
        chr11: int = 0
        chr12: int = 0
        chr13: int = 0
        chr14: int = 0
        chr15: int = 0
        chr16: int = 0
        chr17: int = 0
        chr18: int = 0
        chr19: int = 0
        chr20: int = 0
        chr21: int = 0
        chr22: int = 0
        chrX: int = 0
        chrY: int = 0

    class DiffPeakDistToTSSDistribution(BaseModel):
        less_negative_50000: int = 0
        between_negative_50000_20000: int = 0
        between_negative_20000_10000: int = 0
        between_negative_10000_5000: int = 0
        between_negative_5000_0: int = 0
        dist_zero: int = 0
        between_positive_0_5000: int = 0
        between_positive_5000_10000: int = 0
        between_positive_10000_20000: int = 0
        between_positive_20000_50000: int = 0
        more_positive_50000: int = 0

    statistics: Statistics = None
    diff_data: DiffData = None
    perturbed_gene_frequency: list[EchartsPieData] = None
    cell_line_frequency: list[EchartsPieData] = None
    diff_peaks_num_in_struct: list[EchartsPieData] = None
    diff_peaks_num_in_chrs: DiffPeaksNumInChr = None
    diff_peaks_dist_to_TSS_distribution: DiffPeakDistToTSSDistribution = None


'''
通用 request
'''


# download request
class DownloadRequest(BaseModel):
    id: int
    pb_gene: str
    celline: str
    filter: None = None
    seq: None = None
    paging: Paging = None
