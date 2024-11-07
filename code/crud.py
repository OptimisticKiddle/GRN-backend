from fastapi import Query
from sqlalchemy.orm import Session
from sqlalchemy import func, desc, distinct
import models
from schemas import *

PAGING_START = 0
PAGING_LEN = 10

overall_data = models.OverallData
diff_peaks = models.DiffPeaks
diff_go_enrich = models.DiffGOEnrich
diff_kegg_enrich = models.DiffKEGGEnrich
diff_motif = models.DiffMotif
diff_footprint = models.DiffFootprint
ctrl_peaks = models.CtrlPeaks
treat_peaks = models.TreatPeaks

'''
通用
'''


def get_paging_attrs(paging: Paging) -> (int, int):
    if paging is not None:
        return paging.start, paging.length
    return PAGING_START, PAGING_LEN


def get_order_attr(model_name: str, seq: object):
    order_attr = None
    if seq is not None:
        for key, val in seq.dict().items():
            if val is not None:
                seq_key = eval(f'{model_name}.{key[:-4]}')  # 去除末尾_seq
                order_attr = seq_key if val == 1 else desc(seq_key)
    return order_attr


'''
db查询操作
'''


# 获取总体overall_data
def get_overall_data(db: Session, browser_filter: BrowserFilter = None, paging: Paging = None) -> (Query, int):
    filters = []
    if browser_filter is not None:
        if browser_filter.id is not None:
            filters.append(overall_data.id == browser_filter.id)
        if browser_filter.pb_gene is not None:
            filters.append(overall_data.pb_gene.like(f'%{browser_filter.pb_gene}%'))
        if browser_filter.pb_ensembl is not None:
            filters.append(overall_data.pb_ensembl.like(f'%{browser_filter.pb_ensembl}%'))
        if browser_filter.celline is not None:
            filters.append(overall_data.celline.like(f'%{browser_filter.celline}%'))
        if browser_filter.method is not None:
            filters.append(overall_data.method.like(f'%{browser_filter.method}%'))
        if browser_filter.datasource is not None:
            filters.append(overall_data.datasource.like(f'%{browser_filter.datasource}%'))
        if browser_filter.n_sample_greater is not None:
            filters.append(overall_data.n_sample >= browser_filter.n_sample_greater)
        if browser_filter.n_sample_less is not None:
            filters.append(overall_data.n_sample <= browser_filter.n_sample_less)

    paging_start, paging_len = get_paging_attrs(paging)

    query_data = db.query(overall_data).filter(*filters).offset(paging_start).limit(paging_len).all()
    records_sum = db.query(func.count(overall_data.id)).filter(*filters).scalar()

    return query_data, records_sum


# 获取数据源枚举值
def get_datasource_enum(db: Session) -> Query:
    enum_data = db.query(distinct(overall_data.datasource)).all()
    return enum_data


# 获取敲除方法枚举值
def get_method_enum(db: Session) -> Query:
    enum_data = db.query(distinct(overall_data.method)).all()
    return enum_data


# 获取单个dataset info
def get_single_data(db: Session, id: int) -> Query:
    query_data = db.query(overall_data).filter(overall_data.id == id).first()
    return query_data


# 获取 diff peaks
def get_diff_peaks(db: Session, req: DiffPeaksRequest | DownloadRequest = None) -> (Query, int):
    filters = [diff_peaks.id == req.id]
    if req.filter is not None:
        if req.filter.chr is not None:
            filters.append(diff_peaks.chr.like(f'%{req.filter.chr}%'))
        if req.filter.start_min is not None:
            filters.append(diff_peaks.start >= req.filter.start_min)
        if req.filter.start_max is not None:
            filters.append(diff_peaks.start <= req.filter.start_max)
        if req.filter.end_min is not None:
            filters.append(diff_peaks.end >= req.filter.end_min)
        if req.filter.end_max is not None:
            filters.append(diff_peaks.end <= req.filter.end_max)
        if req.filter.width_min is not None:
            filters.append(diff_peaks.width >= req.filter.width_min)
        if req.filter.width_max is not None:
            filters.append(diff_peaks.width <= req.filter.width_max)
        if req.filter.conc_min is not None:
            filters.append(diff_peaks.conc >= req.filter.conc_min)
        if req.filter.conc_max is not None:
            filters.append(diff_peaks.conc <= req.filter.conc_max)
        if req.filter.conc_ctrl_min is not None:
            filters.append(diff_peaks.conc_ctrl >= req.filter.conc_ctrl_min)
        if req.filter.conc_ctrl_max is not None:
            filters.append(diff_peaks.conc_ctrl <= req.filter.conc_ctrl_max)
        if req.filter.conc_treat_min is not None:
            filters.append(diff_peaks.conc_treat >= req.filter.conc_treat_min)
        if req.filter.conc_treat_max is not None:
            filters.append(diff_peaks.conc_treat <= req.filter.conc_treat_max)
        if req.filter.fold_min is not None:
            filters.append(diff_peaks.fold >= req.filter.fold_min)
        if req.filter.fold_max is not None:
            filters.append(diff_peaks.fold <= req.filter.fold_max)
        if req.filter.p_value_min is not None:
            filters.append(diff_peaks.p_value >= req.filter.p_value_min)
        if req.filter.p_value_max is not None:
            filters.append(diff_peaks.p_value <= req.filter.p_value_max)
        if req.filter.FDR_min is not None:
            filters.append(diff_peaks.FDR >= req.filter.FDR_min)
        if req.filter.FDR_max is not None:
            filters.append(diff_peaks.FDR <= req.filter.FDR_max)
        if req.filter.gene_start_min is not None:
            filters.append(diff_peaks.gene_start >= req.filter.gene_start_min)
        if req.filter.gene_start_max is not None:
            filters.append(diff_peaks.gene_start <= req.filter.gene_start_max)
        if req.filter.gene_end_min is not None:
            filters.append(diff_peaks.gene_end >= req.filter.gene_end_min)
        if req.filter.gene_end_max is not None:
            filters.append(diff_peaks.gene_end <= req.filter.gene_end_max)
        if req.filter.gene_len_min is not None:
            filters.append(diff_peaks.gene_len >= req.filter.gene_len_min)
        if req.filter.gene_len_max is not None:
            filters.append(diff_peaks.gene_len <= req.filter.gene_len_max)
        if req.filter.distance_to_tss_min is not None:
            filters.append(diff_peaks.distance_to_tss >= req.filter.distance_to_tss_min)
        if req.filter.distance_to_tss_max is not None:
            filters.append(diff_peaks.distance_to_tss <= req.filter.distance_to_tss_max)

        if req.filter.gene_chr is not None:
            filters.append(diff_peaks.gene_chr == req.filter.gene_chr)
        if req.filter.gene_id is not None:
            filters.append(diff_peaks.gene_id.like(f'%{req.filter.gene_id}%'))
        if req.filter.transcript_id is not None:
            filters.append(diff_peaks.transcript_id.like(f'%{req.filter.transcript_id}%'))

    # 找到排序key
    order_attr = get_order_attr("diff_peaks", req.seq)

    paging_start, paging_len = get_paging_attrs(req.paging)

    query_data = db.query(diff_peaks).filter(*filters).order_by(order_attr) \
        .offset(paging_start).limit(paging_len).all()
    records_sum = db.query(func.count(diff_peaks.id)).filter(*filters).scalar()
    return query_data, records_sum


# 获取 diff peaks 数值属性列范围
def get_diff_peaks_numerical_range(db: Session, id: int) -> dict:
    filters = [diff_peaks.id == id]
    start_min = db.query(func.min(diff_peaks.start)).filter(*filters).scalar()
    start_max = db.query(func.max(diff_peaks.start)).filter(*filters).scalar()
    end_min = db.query(func.min(diff_peaks.end)).filter(*filters).scalar()
    end_max = db.query(func.max(diff_peaks.end)).filter(*filters).scalar()
    width_min = db.query(func.min(diff_peaks.width)).filter(*filters).scalar()
    width_max = db.query(func.max(diff_peaks.width)).filter(*filters).scalar()
    conc_max = db.query(func.max(diff_peaks.conc)).filter(*filters).scalar()
    conc_ctrl_max = db.query(func.max(diff_peaks.conc_ctrl)).filter(*filters).scalar()
    conc_treat_max = db.query(func.max(diff_peaks.conc_treat)).filter(*filters).scalar()
    fold_min = db.query(func.min(diff_peaks.fold)).filter(*filters).scalar()
    fold_max = db.query(func.max(diff_peaks.fold)).filter(*filters).scalar()
    p_value_max = db.query(func.max(diff_peaks.p_value)).filter(*filters).scalar()
    FDR_max = db.query(func.max(diff_peaks.FDR)).filter(*filters).scalar()
    gene_start_min = db.query(func.min(diff_peaks.gene_start)).filter(*filters).scalar()
    gene_start_max = db.query(func.max(diff_peaks.gene_start)).filter(*filters).scalar()
    gene_end_min = db.query(func.min(diff_peaks.gene_end)).filter(*filters).scalar()
    gene_end_max = db.query(func.max(diff_peaks.gene_end)).filter(*filters).scalar()
    gene_len_min = db.query(func.min(diff_peaks.gene_len)).filter(*filters).scalar()
    gene_len_max = db.query(func.max(diff_peaks.gene_len)).filter(*filters).scalar()
    distance_to_tss_min = db.query(func.min(diff_peaks.distance_to_tss)).filter(*filters).scalar()
    distance_to_tss_max = db.query(func.max(diff_peaks.distance_to_tss)).filter(*filters).scalar()
    data = {
        "start_min": start_min,
        "start_max": start_max,
        "end_min": end_min,
        "end_max": end_max,
        "width_min": width_min,
        "width_max": width_max,
        "conc_max": conc_max,
        "conc_ctrl_max": conc_ctrl_max,
        "conc_treat_max": conc_treat_max,
        "fold_min": fold_min,
        "fold_max": fold_max,
        "p_value_max": p_value_max,
        "FDR_max": FDR_max,
        "gene_start_min": gene_start_min,
        "gene_start_max": gene_start_max,
        "gene_end_min": gene_end_min,
        "gene_end_max": gene_end_max,
        "gene_len_min": gene_len_min,
        "gene_len_max": gene_len_max,
        "distance_to_tss_min": distance_to_tss_min,
        "distance_to_tss_max": distance_to_tss_max
    }
    return data


# 获取diff go enrich
def get_diff_go_enrich(db: Session, req: DiffEnrichRequest | DownloadRequest = None) -> (Query, int):
    filters = [diff_go_enrich.id == req.id]
    order_attr = get_order_attr("diff_go_enrich", req.seq)
    if order_attr is None:
        order_attr = desc(diff_go_enrich.count)

    paging_start, paging_len = get_paging_attrs(req.paging)

    query_data = db.query(diff_go_enrich).filter(*filters).order_by(order_attr) \
        .offset(paging_start).limit(paging_len).all()
    records_sum = db.query(func.count(diff_go_enrich.id)).filter(*filters).scalar()
    return query_data, records_sum


# 获取diff kegg enrich
def get_diff_kegg_enrich(db: Session, req: DiffEnrichRequest | DownloadRequest = None) -> (Query, int):
    filters = [diff_kegg_enrich.id == req.id]
    order_attr = get_order_attr("diff_kegg_enrich", req.seq)
    if order_attr is None:
        order_attr = desc(diff_kegg_enrich.count)

    paging_start, paging_len = get_paging_attrs(req.paging)

    query_data = db.query(diff_kegg_enrich).filter(*filters).order_by(order_attr) \
        .offset(paging_start).limit(paging_len).all()
    records_sum = db.query(func.count(diff_kegg_enrich.id)).filter(*filters).scalar()
    return query_data, records_sum


# 获取diff motif
def get_diff_motif(db: Session, req: DiffMotifRequest | DownloadRequest = None) -> (Query, int):
    filters = [diff_motif.id == req.id]
    order_attr = get_order_attr("diff_motif", req.seq)
    paging_start, paging_len = get_paging_attrs(req.paging)

    query_data = db.query(diff_motif).filter(*filters).order_by(order_attr) \
        .offset(paging_start).limit(paging_len).all()
    records_sum = db.query(func.count(diff_motif.id)).filter(*filters).scalar()
    return query_data, records_sum


# 获取diff footprint
def get_diff_footprint(db: Session, req: DiffFootprintRequest | DownloadRequest = None) -> (Query, int):
    filters = [diff_footprint.id == req.id]
    if req.filter is not None:
        if req.filter.motif is not None:
            filters.append(diff_footprint.motif.like(f'%{req.filter.motif}%'))
        if req.filter.tf is not None:
            filters.append(diff_footprint.tf.like(f'%{req.filter.tf}%'))
        if req.filter.num_min is not None:
            filters.append(diff_footprint.num >= req.filter.num_min)
        if req.filter.num_max is not None:
            filters.append(diff_footprint.num <= req.filter.num_max)
        if req.filter.protection_score_ctrl_min is not None:
            filters.append(diff_footprint.protection_score_ctrl >= req.filter.protection_score_ctrl_min)
        if req.filter.protection_score_ctrl_max is not None:
            filters.append(diff_footprint.protection_score_ctrl <= req.filter.protection_score_ctrl_max)
        if req.filter.protection_score_treat_min is not None:
            filters.append(diff_footprint.protection_score_treat >= req.filter.protection_score_treat_min)
        if req.filter.protection_score_treat_max is not None:
            filters.append(diff_footprint.protection_score_treat <= req.filter.protection_score_treat_max)
        if req.filter.tc_ctrl_min is not None:
            filters.append(diff_footprint.tc_ctrl >= req.filter.tc_ctrl_min)
        if req.filter.tc_ctrl_max is not None:
            filters.append(diff_footprint.tc_ctrl <= req.filter.tc_ctrl_max)
        if req.filter.tc_treat_min is not None:
            filters.append(diff_footprint.tc_treat >= req.filter.tc_treat_min)
        if req.filter.tc_treat_max is not None:
            filters.append(diff_footprint.tc_treat <= req.filter.tc_treat_max)
        if req.filter.tf_activity_min is not None:
            filters.append(diff_footprint.tf_activity >= req.filter.tf_activity_min)
        if req.filter.tf_activity_max is not None:
            filters.append(diff_footprint.tf_activity <= req.filter.tf_activity_max)
        if req.filter.z_score_min is not None:
            filters.append(diff_footprint.z_score >= req.filter.z_score_min)
        if req.filter.z_score_max is not None:
            filters.append(diff_footprint.z_score <= req.filter.z_score_max)
        if req.filter.p_value_min is not None:
            filters.append(diff_footprint.p_value >= req.filter.p_value_min)
        if req.filter.p_value_max is not None:
            filters.append(diff_footprint.p_value <= req.filter.p_value_max)

    order_attr = get_order_attr("diff_footprint", req.seq)
    if order_attr is None:
        order_attr = desc(diff_footprint.z_score)

    paging_start, paging_len = get_paging_attrs(req.paging)

    query_data = db.query(diff_footprint).filter(*filters).order_by(order_attr) \
        .offset(paging_start).limit(paging_len).all()
    records_sum = db.query(func.count(diff_footprint.id)).filter(*filters).scalar()
    return query_data, records_sum


# 获取 diff footprint 数值属性列范围
def get_diff_footprint_numerical_range(db: Session, id: int) -> dict:
    filters = [diff_footprint.id == id]
    num_min = db.query(func.min(diff_footprint.num)).filter(*filters).scalar()
    num_max = db.query(func.max(diff_footprint.num)).filter(*filters).scalar()
    protection_score_ctrl_max = db.query(func.max(diff_footprint.protection_score_ctrl)).filter(*filters).scalar()
    protection_score_treat_max = db.query(func.max(diff_footprint.protection_score_treat)).filter(*filters).scalar()
    tc_ctrl_max = db.query(func.max(diff_footprint.tc_ctrl)).filter(*filters).scalar()
    tc_treat_max = db.query(func.max(diff_footprint.tc_treat)).filter(*filters).scalar()
    tf_activity_min = db.query(func.min(diff_footprint.tf_activity)).filter(*filters).scalar()
    tf_activity_max = db.query(func.max(diff_footprint.tf_activity)).filter(*filters).scalar()
    z_score_min = db.query(func.min(diff_footprint.z_score)).filter(*filters).scalar()
    z_score_max = db.query(func.max(diff_footprint.z_score)).filter(*filters).scalar()
    p_value_max = db.query(func.max(diff_footprint.p_value)).filter(*filters).scalar()
    data = {
        "num_min": num_min,
        "num_max": num_max,
        "protection_score_ctrl_max": protection_score_ctrl_max,
        "protection_score_treat_max": protection_score_treat_max,
        "tc_ctrl_max": tc_ctrl_max,
        "tc_treat_max": tc_treat_max,
        "tf_activity_min": tf_activity_min,
        "tf_activity_max": tf_activity_max,
        "z_score_min": z_score_min,
        "z_score_max": z_score_max,
        "p_value_max": p_value_max
    }
    return data


# 获取 ctrl peaks
def get_ctrl_peaks(db: Session, req: PeaksRequest | DownloadRequest = None) -> (Query, int):
    filters = [ctrl_peaks.id == req.id]
    if req.filter is not None:
        if req.filter.chr is not None:
            filters.append(ctrl_peaks.chr.like(f'%{req.filter.chr}%'))
        if req.filter.gene_chr is not None:
            filters.append(ctrl_peaks.gene_chr == req.filter.gene_chr)
        if req.filter.gene_id is not None:
            filters.append(ctrl_peaks.gene_id.like(f'%{req.filter.gene_id}%'))
        if req.filter.transcript_id is not None:
            filters.append(ctrl_peaks.transcript_id.like(f'%{req.filter.transcript_id}%'))
        if req.filter.start_min is not None:
            filters.append(ctrl_peaks.start >= req.filter.start_min)
        if req.filter.start_max is not None:
            filters.append(ctrl_peaks.start <= req.filter.start_max)
        if req.filter.end_min is not None:
            filters.append(ctrl_peaks.end >= req.filter.end_min)
        if req.filter.end_max is not None:
            filters.append(ctrl_peaks.end <= req.filter.end_max)
        if req.filter.width_min is not None:
            filters.append(ctrl_peaks.width >= req.filter.width_min)
        if req.filter.width_max is not None:
            filters.append(ctrl_peaks.width <= req.filter.width_max)
        if req.filter.score_min is not None:
            filters.append(ctrl_peaks.score >= req.filter.score_min)
        if req.filter.score_max is not None:
            filters.append(ctrl_peaks.score <= req.filter.score_max)
        if req.filter.signal_value_min is not None:
            filters.append(ctrl_peaks.signal_value >= req.filter.signal_value_min)
        if req.filter.signal_value_max is not None:
            filters.append(ctrl_peaks.signal_value <= req.filter.signal_value_max)
        if req.filter.log_p_value_min is not None:
            filters.append(ctrl_peaks.log_p_value >= req.filter.log_p_value_min)
        if req.filter.log_p_value_max is not None:
            filters.append(ctrl_peaks.log_p_value <= req.filter.log_p_value_max)
        if req.filter.log_q_value_min is not None:
            filters.append(ctrl_peaks.log_q_value >= req.filter.log_q_value_min)
        if req.filter.log_q_value_max is not None:
            filters.append(ctrl_peaks.log_q_value <= req.filter.log_q_value_max)
        if req.filter.peak_offset_min is not None:
            filters.append(ctrl_peaks.peak_offset >= req.filter.peak_offset_min)
        if req.filter.peak_offset_max is not None:
            filters.append(ctrl_peaks.peak_offset <= req.filter.peak_offset_max)
        if req.filter.gene_start_min is not None:
            filters.append(ctrl_peaks.gene_start >= req.filter.gene_start_min)
        if req.filter.gene_start_max is not None:
            filters.append(ctrl_peaks.gene_start <= req.filter.gene_start_max)
        if req.filter.gene_end_min is not None:
            filters.append(ctrl_peaks.gene_end >= req.filter.gene_end_min)
        if req.filter.gene_end_max is not None:
            filters.append(ctrl_peaks.gene_end <= req.filter.gene_end_max)
        if req.filter.gene_len_min is not None:
            filters.append(ctrl_peaks.gene_len >= req.filter.gene_len_min)
        if req.filter.gene_len_max is not None:
            filters.append(ctrl_peaks.gene_len <= req.filter.gene_len_max)
        if req.filter.distance_to_tss_min is not None:
            filters.append(ctrl_peaks.distance_to_tss >= req.filter.distance_to_tss_min)
        if req.filter.distance_to_tss_max is not None:
            filters.append(ctrl_peaks.distance_to_tss <= req.filter.distance_to_tss_max)

    order_attr = get_order_attr("ctrl_peaks", req.seq)
    paging_start, paging_len = get_paging_attrs(req.paging)

    query_data = db.query(ctrl_peaks).filter(*filters).order_by(order_attr) \
        .offset(paging_start).limit(paging_len).all()
    records_sum = db.query(func.count(ctrl_peaks.id)).filter(*filters).scalar()
    return query_data, records_sum


# 获取 ctrl peaks 数值属性列范围
def get_ctrl_peaks_numerical_range(db: Session, id: int) -> dict:
    filters = [ctrl_peaks.id == id]
    start_min = db.query(func.min(ctrl_peaks.start)).filter(*filters).scalar()
    start_max = db.query(func.max(ctrl_peaks.start)).filter(*filters).scalar()
    end_min = db.query(func.min(ctrl_peaks.end)).filter(*filters).scalar()
    end_max = db.query(func.max(ctrl_peaks.end)).filter(*filters).scalar()
    width_min = db.query(func.min(ctrl_peaks.width)).filter(*filters).scalar()
    width_max = db.query(func.max(ctrl_peaks.width)).filter(*filters).scalar()
    score_min = db.query(func.min(ctrl_peaks.score)).filter(*filters).scalar()
    score_max = db.query(func.max(ctrl_peaks.score)).filter(*filters).scalar()
    signal_value_min = db.query(func.min(ctrl_peaks.signal_value)).filter(*filters).scalar()
    signal_value_max = db.query(func.max(ctrl_peaks.signal_value)).filter(*filters).scalar()
    log_p_value_min = db.query(func.min(ctrl_peaks.log_p_value)).filter(*filters).scalar()
    log_p_value_max = db.query(func.max(ctrl_peaks.log_p_value)).filter(*filters).scalar()
    log_q_value_min = db.query(func.min(ctrl_peaks.log_q_value)).filter(*filters).scalar()
    log_q_value_max = db.query(func.max(ctrl_peaks.log_q_value)).filter(*filters).scalar()
    peak_offset_min = db.query(func.min(ctrl_peaks.peak_offset)).filter(*filters).scalar()
    peak_offset_max = db.query(func.max(ctrl_peaks.peak_offset)).filter(*filters).scalar()
    gene_start_min = db.query(func.min(ctrl_peaks.gene_start)).filter(*filters).scalar()
    gene_start_max = db.query(func.max(ctrl_peaks.gene_start)).filter(*filters).scalar()
    gene_end_min = db.query(func.min(ctrl_peaks.gene_end)).filter(*filters).scalar()
    gene_end_max = db.query(func.max(ctrl_peaks.gene_end)).filter(*filters).scalar()
    gene_len_min = db.query(func.min(ctrl_peaks.gene_len)).filter(*filters).scalar()
    gene_len_max = db.query(func.max(ctrl_peaks.gene_len)).filter(*filters).scalar()
    distance_to_tss_min = db.query(func.min(ctrl_peaks.distance_to_tss)).filter(*filters).scalar()
    distance_to_tss_max = db.query(func.max(ctrl_peaks.distance_to_tss)).filter(*filters).scalar()

    data = {
        "start_min": start_min,
        "start_max": start_max,
        "end_min": end_min,
        "end_max": end_max,
        "width_min": width_min,
        "width_max": width_max,
        "score_min": score_min,
        "score_max": score_max,
        "signal_value_min": signal_value_min,
        "signal_value_max": signal_value_max,
        "log_p_value_min": log_p_value_min,
        "log_p_value_max": log_p_value_max,
        "log_q_value_min": log_q_value_min,
        "log_q_value_max": log_q_value_max,
        "peak_offset_min": peak_offset_min,
        "peak_offset_max": peak_offset_max,
        "gene_start_min": gene_start_min,
        "gene_start_max": gene_start_max,
        "gene_end_min": gene_end_min,
        "gene_end_max": gene_end_max,
        "gene_len_min": gene_len_min,
        "gene_len_max": gene_len_max,
        "distance_to_tss_min": distance_to_tss_min,
        "distance_to_tss_max": distance_to_tss_max
    }
    return data


# 获取 treat peaks
def get_treat_peaks(db: Session, req: PeaksRequest | DownloadRequest = None) -> (Query, int):
    filters = [treat_peaks.id == req.id]
    if req.filter is not None:
        if req.filter.chr is not None:
            filters.append(treat_peaks.chr.like(f'%{req.filter.chr}%'))
        if req.filter.gene_chr is not None:
            filters.append(treat_peaks.gene_chr == req.filter.gene_chr)
        if req.filter.gene_id is not None:
            filters.append(treat_peaks.gene_id.like(f'%{req.filter.gene_id}%'))
        if req.filter.transcript_id is not None:
            filters.append(treat_peaks.transcript_id.like(f'%{req.filter.transcript_id}%'))
        if req.filter.start_min is not None:
            filters.append(treat_peaks.start >= req.filter.start_min)
        if req.filter.start_max is not None:
            filters.append(treat_peaks.start <= req.filter.start_max)
        if req.filter.end_min is not None:
            filters.append(treat_peaks.end >= req.filter.end_min)
        if req.filter.end_max is not None:
            filters.append(treat_peaks.end <= req.filter.end_max)
        if req.filter.width_min is not None:
            filters.append(treat_peaks.width >= req.filter.width_min)
        if req.filter.width_max is not None:
            filters.append(treat_peaks.width <= req.filter.width_max)
        if req.filter.score_min is not None:
            filters.append(treat_peaks.score >= req.filter.score_min)
        if req.filter.score_max is not None:
            filters.append(treat_peaks.score <= req.filter.score_max)
        if req.filter.signal_value_min is not None:
            filters.append(treat_peaks.signal_value >= req.filter.signal_value_min)
        if req.filter.signal_value_max is not None:
            filters.append(treat_peaks.signal_value <= req.filter.signal_value_max)
        if req.filter.log_p_value_min is not None:
            filters.append(treat_peaks.log_p_value >= req.filter.log_p_value_min)
        if req.filter.log_p_value_max is not None:
            filters.append(treat_peaks.log_p_value <= req.filter.log_p_value_max)
        if req.filter.log_q_value_min is not None:
            filters.append(treat_peaks.log_q_value >= req.filter.log_q_value_min)
        if req.filter.log_q_value_max is not None:
            filters.append(treat_peaks.log_q_value <= req.filter.log_q_value_max)
        if req.filter.peak_offset_min is not None:
            filters.append(treat_peaks.peak_offset >= req.filter.peak_offset_min)
        if req.filter.peak_offset_max is not None:
            filters.append(treat_peaks.peak_offset <= req.filter.peak_offset_max)
        if req.filter.gene_start_min is not None:
            filters.append(treat_peaks.gene_start >= req.filter.gene_start_min)
        if req.filter.gene_start_max is not None:
            filters.append(treat_peaks.gene_start <= req.filter.gene_start_max)
        if req.filter.gene_end_min is not None:
            filters.append(treat_peaks.gene_end >= req.filter.gene_end_min)
        if req.filter.gene_end_max is not None:
            filters.append(treat_peaks.gene_end <= req.filter.gene_end_max)
        if req.filter.gene_len_min is not None:
            filters.append(treat_peaks.gene_len >= req.filter.gene_len_min)
        if req.filter.gene_len_max is not None:
            filters.append(treat_peaks.gene_len <= req.filter.gene_len_max)
        if req.filter.distance_to_tss_min is not None:
            filters.append(treat_peaks.distance_to_tss >= req.filter.distance_to_tss_min)
        if req.filter.distance_to_tss_max is not None:
            filters.append(treat_peaks.distance_to_tss <= req.filter.distance_to_tss_max)

    order_attr = get_order_attr("treat_peaks", req.seq)
    paging_start, paging_len = get_paging_attrs(req.paging)

    query_data = db.query(treat_peaks).filter(*filters).order_by(order_attr) \
        .offset(paging_start).limit(paging_len).all()
    records_sum = db.query(func.count(treat_peaks.id)).filter(*filters).scalar()
    return query_data, records_sum


# 获取 treat peaks 数值属性列范围
def get_treat_peaks_numerical_range(db: Session, id: int) -> dict:
    filters = [treat_peaks.id == id]
    start_min = db.query(func.min(treat_peaks.start)).filter(*filters).scalar()
    start_max = db.query(func.max(treat_peaks.start)).filter(*filters).scalar()
    end_min = db.query(func.min(treat_peaks.end)).filter(*filters).scalar()
    end_max = db.query(func.max(treat_peaks.end)).filter(*filters).scalar()
    width_min = db.query(func.min(treat_peaks.width)).filter(*filters).scalar()
    width_max = db.query(func.max(treat_peaks.width)).filter(*filters).scalar()
    score_min = db.query(func.min(treat_peaks.score)).filter(*filters).scalar()
    score_max = db.query(func.max(treat_peaks.score)).filter(*filters).scalar()
    signal_value_min = db.query(func.min(treat_peaks.signal_value)).filter(*filters).scalar()
    signal_value_max = db.query(func.max(treat_peaks.signal_value)).filter(*filters).scalar()
    log_p_value_min = db.query(func.min(treat_peaks.log_p_value)).filter(*filters).scalar()
    log_p_value_max = db.query(func.max(treat_peaks.log_p_value)).filter(*filters).scalar()
    log_q_value_min = db.query(func.min(treat_peaks.log_q_value)).filter(*filters).scalar()
    log_q_value_max = db.query(func.max(treat_peaks.log_q_value)).filter(*filters).scalar()
    peak_offset_min = db.query(func.min(treat_peaks.peak_offset)).filter(*filters).scalar()
    peak_offset_max = db.query(func.max(treat_peaks.peak_offset)).filter(*filters).scalar()
    gene_start_min = db.query(func.min(treat_peaks.gene_start)).filter(*filters).scalar()
    gene_start_max = db.query(func.max(treat_peaks.gene_start)).filter(*filters).scalar()
    gene_end_min = db.query(func.min(treat_peaks.gene_end)).filter(*filters).scalar()
    gene_end_max = db.query(func.max(treat_peaks.gene_end)).filter(*filters).scalar()
    gene_len_min = db.query(func.min(treat_peaks.gene_len)).filter(*filters).scalar()
    gene_len_max = db.query(func.max(treat_peaks.gene_len)).filter(*filters).scalar()
    distance_to_tss_min = db.query(func.min(treat_peaks.distance_to_tss)).filter(*filters).scalar()
    distance_to_tss_max = db.query(func.max(treat_peaks.distance_to_tss)).filter(*filters).scalar()

    data = {
        "start_min": start_min,
        "start_max": start_max,
        "end_min": end_min,
        "end_max": end_max,
        "width_min": width_min,
        "width_max": width_max,
        "score_min": score_min,
        "score_max": score_max,
        "signal_value_min": signal_value_min,
        "signal_value_max": signal_value_max,
        "log_p_value_min": log_p_value_min,
        "log_p_value_max": log_p_value_max,
        "log_q_value_min": log_q_value_min,
        "log_q_value_max": log_q_value_max,
        "peak_offset_min": peak_offset_min,
        "peak_offset_max": peak_offset_max,
        "gene_start_min": gene_start_min,
        "gene_start_max": gene_start_max,
        "gene_end_min": gene_end_min,
        "gene_end_max": gene_end_max,
        "gene_len_min": gene_len_min,
        "gene_len_max": gene_len_max,
        "distance_to_tss_min": distance_to_tss_min,
        "distance_to_tss_max": distance_to_tss_max
    }
    return data


# 获取statistics info
def get_statistics_info(db: Session) -> dict:
    # statistics
    perturb_datasets_num = db.query(func.count(overall_data.id)).scalar()
    perturb_gene_num = db.query(func.count(distinct(overall_data.pb_gene))).scalar()
    cell_line_num = db.query(func.count(distinct(overall_data.celline))).scalar()
    n_samples_num = db.query(func.sum(overall_data.n_sample)).scalar()
    statistics = {
        "perturb_datasets_num": perturb_datasets_num,
        "perturb_gene_num": perturb_gene_num,
        "cell_line_num": cell_line_num,
        "n_samples_num": n_samples_num
    }

    # diff_data
    diff_peaks_num = db.query(func.count(diff_peaks.id)).scalar()
    diff_motif_num = db.query(func.count(diff_motif.id)).scalar()
    diff_GO_enrich_num = db.query(func.count(diff_go_enrich.id)).scalar()
    diff_KEGG_enrich_num = db.query(func.count(diff_kegg_enrich.id)).scalar()
    diff_footprint_num = db.query(func.count(diff_footprint.id)).scalar()
    diff_data = {
        "diff_peaks_num": diff_peaks_num,
        "diff_motif_num": diff_motif_num,
        "diff_GO_enrich_num": diff_GO_enrich_num,
        "diff_KEGG_enrich_num": diff_KEGG_enrich_num,
        "diff_footprint_num": diff_footprint_num
    }

    # perturbed_gene_frequency
    perturbed_gene_frequency = db.query(overall_data.pb_gene, func.count('*')).group_by(overall_data.pb_gene) \
        .order_by(desc(func.count('*'))).limit(10).all()
    perturbed_gene_frequency = [{"name": item[0], "value": item[1]} for item in perturbed_gene_frequency]

    # cell_line_frequency
    cell_line_frequency = db.query(overall_data.celline, func.count('*')).group_by(overall_data.celline) \
        .order_by(desc(func.count('*'))).limit(10).all()
    cell_line_frequency = [{"name": item[0], "value": item[1]} for item in cell_line_frequency]

    # diff_peaks_num_in_struct
    promoter = db.query(func.count(diff_peaks.annotation)).filter(diff_peaks.annotation.like('%Promoter%')).scalar()
    UTR_5 = db.query(func.count(diff_peaks.annotation)).filter(diff_peaks.annotation.like("%5' UTR%")).scalar()
    UTR_3 = db.query(func.count(diff_peaks.annotation)).filter(diff_peaks.annotation.like("%3' UTR%")).scalar()
    exon = db.query(func.count(diff_peaks.annotation)).filter(diff_peaks.annotation.like("%Exon%")).scalar()
    intron = db.query(func.count(diff_peaks.annotation)).filter(diff_peaks.annotation.like("%Intron%")).scalar()
    downstream = db.query(func.count(diff_peaks.annotation)).filter(diff_peaks.annotation.like("%Downstream%")).scalar()
    distal_intergenic = db.query(func.count(diff_peaks.annotation)).filter(
        diff_peaks.annotation.like("%Distal Intergenic%")).scalar()
    diff_peaks_num_in_struct = [
        {"name": "promoter", "value": promoter},
        {"name": "UTR_5", "value": UTR_5},
        {"name": "UTR_3", "value": UTR_3},
        {"name": "exon", "value": exon},
        {"name": "intron", "value": intron},
        {"name": "downstream", "value": downstream},
        {"name": "distal_intergenic", "value": distal_intergenic}
    ]

    # diff_peaks_num_in_chrs
    diff_peaks_num_in_chrs = db.query(diff_peaks.chr, func.count('*')).group_by(diff_peaks.chr).all()
    diff_peaks_num_in_chrs = {item[0]: item[1] for item in diff_peaks_num_in_chrs}

    # diff_peaks_dist_to_TSS_distribution
    less_negative_50000 = db.query(func.count(diff_peaks.distance_to_tss)).filter(
        diff_peaks.distance_to_tss < -50000).scalar()
    between_negative_50000_20000 = db.query(func.count(diff_peaks.distance_to_tss)).filter(
        diff_peaks.distance_to_tss >= -50000, diff_peaks.distance_to_tss < -20000).scalar()
    between_negative_20000_10000 = db.query(func.count(diff_peaks.distance_to_tss)).filter(
        diff_peaks.distance_to_tss >= -20000, diff_peaks.distance_to_tss < -10000).scalar()
    between_negative_10000_5000 = db.query(func.count(diff_peaks.distance_to_tss)).filter(
        diff_peaks.distance_to_tss >= -10000, diff_peaks.distance_to_tss < -5000).scalar()
    between_negative_5000_0 = db.query(func.count(diff_peaks.distance_to_tss)).filter(
        diff_peaks.distance_to_tss >= -5000, diff_peaks.distance_to_tss < 0).scalar()
    dist_zero = db.query(func.count(diff_peaks.distance_to_tss)).filter(diff_peaks.distance_to_tss == 0).scalar()
    between_positive_0_5000 = db.query(func.count(diff_peaks.distance_to_tss)).filter(
        diff_peaks.distance_to_tss > 0, diff_peaks.distance_to_tss < 5000).scalar()
    between_positive_5000_10000 = db.query(func.count(diff_peaks.distance_to_tss)).filter(
        diff_peaks.distance_to_tss >= 5000, diff_peaks.distance_to_tss < 10000).scalar()
    between_positive_10000_20000 = db.query(func.count(diff_peaks.distance_to_tss)).filter(
        diff_peaks.distance_to_tss >= 10000, diff_peaks.distance_to_tss < 20000).scalar()
    between_positive_20000_50000 = db.query(func.count(diff_peaks.distance_to_tss)).filter(
        diff_peaks.distance_to_tss >= 20000, diff_peaks.distance_to_tss < 50000).scalar()
    more_positive_50000 = db.query(func.count(diff_peaks.distance_to_tss)).filter(
        diff_peaks.distance_to_tss >= 50000).scalar()
    diff_peaks_dist_to_TSS_distribution = {
        "less_negative_50000": less_negative_50000,
        "between_negative_50000_20000": between_negative_50000_20000,
        "between_negative_20000_10000": between_negative_20000_10000,
        "between_negative_10000_5000": between_negative_10000_5000,
        "between_negative_5000_0": between_negative_5000_0,
        "dist_zero": dist_zero,
        "between_positive_0_5000": between_positive_0_5000,
        "between_positive_5000_10000": between_positive_5000_10000,
        "between_positive_10000_20000": between_positive_10000_20000,
        "between_positive_20000_50000": between_positive_20000_50000,
        "more_positive_50000": more_positive_50000
    }

    data = {
        "statistics": statistics,
        "diff_data": diff_data,
        "perturbed_gene_frequency": perturbed_gene_frequency,
        "cell_line_frequency": cell_line_frequency,
        "diff_peaks_num_in_struct": diff_peaks_num_in_struct,
        "diff_peaks_num_in_chrs": diff_peaks_num_in_chrs,
        "diff_peaks_dist_to_TSS_distribution": diff_peaks_dist_to_TSS_distribution
    }
    return data
