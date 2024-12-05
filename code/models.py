from sqlalchemy import Boolean, Column, ForeignKey, Integer, String, Float, DECIMAL, Text
from sqlalchemy.orm import relationship
from database import Base


class OverallData(Base):
    __tablename__ = "overall_data"

    id = Column(Integer, primary_key=True, nullable=False)
    gsm = Column(String)
    gse = Column(String)
    organism = Column(Integer)
    sample_source = Column(String)
    sample_type = Column(String)



class DiffPeaks(Base):
    __tablename__ = "diff_peaks"

    id = Column(Integer, primary_key=True, nullable=False)
    peak_id = Column(Integer, primary_key=True, nullable=False)
    chr = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    width = Column(Integer)
    conc = Column(DECIMAL)
    conc_ctrl = Column(DECIMAL)
    conc_treat = Column(DECIMAL)
    fold = Column(DECIMAL)
    p_value = Column(DECIMAL)
    FDR = Column(DECIMAL)
    annotation = Column(String)
    gene_chr = Column(Integer)
    gene_start = Column(Integer)
    gene_end = Column(Integer)
    gene_len = Column(Integer)
    gene_strand = Column(Integer)
    gene_id = Column(String)
    transcript_id = Column(String)
    distance_to_tss = Column(Integer)


class DiffGOEnrich(Base):
    __tablename__ = "diff_go_enrich"

    id = Column(Integer, primary_key=True, nullable=False)
    GO_id = Column(String, primary_key=True, nullable=False)
    description = Column(String)
    gene_ratio = Column(String)
    bg_ratio = Column(String)
    p_value = Column(DECIMAL)
    p_adjust = Column(DECIMAL)
    q_value = Column(DECIMAL)
    gene_ids = Column(Text)
    count = Column(Integer)


class DiffKEGGEnrich(Base):
    __tablename__ = "diff_kegg_enrich"

    id = Column(Integer, primary_key=True, nullable=False)
    cluster = Column(String)
    KEGG_id = Column(String, primary_key=True, nullable=False)
    description = Column(String)
    gene_ratio = Column(String)
    bg_ratio = Column(String)
    p_value = Column(DECIMAL)
    p_adjust = Column(DECIMAL)
    q_value = Column(DECIMAL)
    gene_ids = Column(Text)
    count = Column(Integer)


class DiffMotif(Base):
    __tablename__ = "diff_motif"

    id = Column(Integer, primary_key=True, nullable=False)
    rank = Column(Integer, primary_key=True, nullable=False)
    motif = Column(Text)
    name = Column(String)
    consensus = Column(String)
    p_value = Column(DECIMAL)
    log_p_value = Column(DECIMAL)
    q_value = Column(DECIMAL)
    target_sequences = Column(DECIMAL)
    target_sequences_percent = Column(DECIMAL)
    bg_sequences = Column(DECIMAL)
    bg_sequences_percent = Column(DECIMAL)


class DiffFootprint(Base):
    __tablename__ = "diff_footprint"

    id = Column(Integer, primary_key=True, nullable=False)
    motif = Column(String, primary_key=True, nullable=False)
    tf = Column(String)
    num = Column(Integer)
    protection_score_ctrl = Column(DECIMAL)
    protection_score_treat = Column(DECIMAL)
    tc_ctrl = Column(DECIMAL)
    tc_treat = Column(DECIMAL)
    tf_activity = Column(DECIMAL)
    z_score = Column(DECIMAL)
    p_value = Column(DECIMAL)


class CtrlPeaks(Base):
    __tablename__ = "ctrl_peaks"

    id = Column(Integer, primary_key=True, nullable=False)
    peak_id = Column(Integer, primary_key=True, nullable=False)
    chr = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    width = Column(Integer)
    score = Column(Integer)
    signal_value = Column(DECIMAL)
    log_p_value = Column(DECIMAL)
    log_q_value = Column(DECIMAL)
    peak_offset = Column(Integer)
    annotation = Column(String)
    gene_chr = Column(Integer)
    gene_start = Column(Integer)
    gene_end = Column(Integer)
    gene_len = Column(Integer)
    gene_strand = Column(Integer)
    gene_id = Column(String)
    transcript_id = Column(String)
    distance_to_tss = Column(Integer)


class TreatPeaks(Base):
    __tablename__ = "treat_peaks"

    id = Column(Integer, primary_key=True, nullable=False)
    peak_id = Column(Integer, primary_key=True, nullable=False)
    chr = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    width = Column(Integer)
    score = Column(Integer)
    signal_value = Column(DECIMAL)
    log_p_value = Column(DECIMAL)
    log_q_value = Column(DECIMAL)
    peak_offset = Column(Integer)
    annotation = Column(String)
    gene_chr = Column(Integer)
    gene_start = Column(Integer)
    gene_end = Column(Integer)
    gene_len = Column(Integer)
    gene_strand = Column(Integer)
    gene_id = Column(String)
    transcript_id = Column(String)
    distance_to_tss = Column(Integer)
