from datetime import datetime
from fastapi import Depends, FastAPI, Body, HTTPException,Response,UploadFile,File
from sqlalchemy.orm import Session
# from fastapi.responses import FileResponse
from starlette.responses import FileResponse
import sys
import crud
import models
from database import SessionLocal, engine
from schemas import *
import pandas
import os, tarfile
import paramiko
from scp import SCPClient, SCPException
import subprocess
import shutil
import os
from fastapi.staticfiles import StaticFiles

models.Base.metadata.create_all(bind=engine)

app = FastAPI()
app.mount("/api/static", StaticFiles(directory="../static"), name="static")

# Dependency
def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()

@app.post("/api/upload/{gse_id}/{gsm_id}/{mark}")
async def upload(gse_id:str,gsm_id:str,mark:str,file: UploadFile = File(...)):
    file_extension = file.filename.split(".")[-1]
    if file_extension.lower()!= "csv":
        return {"message": "Please upload files in CSV format!"}

    fn = file.filename
    path = '../static/GSE' + gse_id + "/GSM"  +  gsm_id + "/" + mark + "/"
    listdir = os.listdir(path)
    for dirname in listdir:
            dirname = path + "//" + dirname
            if os.path.isfile(dirname): # 是文件
                print(dirname)
                os.remove(dirname)      # 删除文件
    file_path = os.path.join(path, fn)
    with open(file_path, "wb") as buffer:
        shutil.copyfileobj(file.file, buffer)
    return {"message": "Upload Successful"}

@app.get("/api/download/{gse_id}/{gsm_id}/{filename}")
async def download_file(gse_id:str,gsm_id:str,filename: str):
    headers = {
               'Content-Disposition': f'attachment; filename="{filename}"'
           }
    path = '../static/GSE' + gse_id + "/GSM" + gsm_id + "/" + filename
    return FileResponse(path=path,headers=headers)

@app.post("/api/get_overall_data", response_model=TableDataResponse)
def get_overall_data(db: Session = Depends(get_db), filter: BrowserFilter = None, paging: Paging = None):
    data, records_sum = crud.get_overall_data(db, filter, paging)
    overall_data = TableDataResponse(data=data, records_sum=records_sum)

    # gfor i in range(len(overall_data.data)):
    #     overall_data.data[i].gse = f"<a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={overall_data.data[i].gse}' target='_blank'>{overall_data.data[i].gse}</a>"
    #     overall_data.data[i].gsm = f"<a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={overall_data.data[i].sm}' target='_blank'>{overall_data.data[i].gsm}</a>"
    return overall_data


@app.post("/get_single_data")
def get_single_data(db: Session = Depends(get_db), id: int = Body(embed=True)):
    data = crud.get_single_data(db, id)

    if data is None:  # 校验防止为空
        return HTTPException(status_code=404, detail="This id data was not found")

    single_data = SingleDataResponse(data=data)
    GSE_id = single_data.data.accession.split('acc=')[1]
    single_data.data.accession = f"<a href='{single_data.data.accession}' target='_blank'>{GSE_id}</a>"
    return single_data


@app.post("/get_diff_peak_data", response_model=TableDataResponse)
def get_diff_peak_data(db: Session = Depends(get_db), req: DiffPeaksRequest = None):
    data, records_sum = crud.get_diff_peaks(db, req)
    diff_peak_data = TableDataResponse(data=data, records_sum=records_sum)
    return diff_peak_data


@app.post("/get_diff_peak_data_range", response_model=DiffPeaksDataRange)
def get_diff_peak_data_range(db: Session = Depends(get_db), id: int = Body(embed=True)):
    data = crud.get_diff_peaks_numerical_range(db, id)
    return data


@app.post("/get_diff_GO_enrichment_data", response_model=TableDataResponse)
def get_diff_GO_enrichment_data(db: Session = Depends(get_db), req: DiffEnrichRequest = None):
    data, records_sum = crud.get_diff_go_enrich(db, req)
    diff_GO_enrichment_data = TableDataResponse(data=data, records_sum=records_sum)
    return diff_GO_enrichment_data


@app.post("/get_diff_KEGG_enrichment_data", response_model=TableDataResponse)
def get_diff_KEGG_enrichment_data(db: Session = Depends(get_db), req: DiffEnrichRequest = None):
    data, records_sum = crud.get_diff_kegg_enrich(db, req)
    diff_KEGG_enrichment_data = TableDataResponse(data=data, records_sum=records_sum)
    return diff_KEGG_enrichment_data


@app.post("/get_diff_motif_data", response_model=TableDataResponse)
def get_diff_motif_data(db: Session = Depends(get_db), req: DiffMotifRequest = None):
    data, records_sum = crud.get_diff_motif(db, req)
    diff_motif_data = TableDataResponse(data=data, records_sum=records_sum)
    return diff_motif_data


@app.post("/get_diff_footprint_data", response_model=TableDataResponse)
def get_diff_footprint_data(db: Session = Depends(get_db), req: DiffFootprintRequest = None):
    data, records_sum = crud.get_diff_footprint(db, req)
    diff_footprint_data = TableDataResponse(data=data, records_sum=records_sum)
    return diff_footprint_data


@app.post("/get_diff_footprint_data_range", response_model=DiffFootprintDataRange)
def get_diff_footprint_data_range(db: Session = Depends(get_db), id: int = Body(embed=True)):
    data = crud.get_diff_footprint_numerical_range(db, id)
    return data


@app.post("/get_ctrl_peak_data", response_model=TableDataResponse)
def get_ctrl_peak_data(db: Session = Depends(get_db), req: PeaksRequest = None):
    data, records_sum = crud.get_ctrl_peaks(db, req)
    ctrl_peak_data = TableDataResponse(data=data, records_sum=records_sum)
    return ctrl_peak_data


@app.post("/get_ctrl_peak_data_range", response_model=PeaksDataRange)
def get_ctrl_peak_data_range(db: Session = Depends(get_db), id: int = Body(embed=True)):
    data = crud.get_ctrl_peaks_numerical_range(db, id)
    return data


@app.post("/get_treat_peak_data", response_model=TableDataResponse)
def get_treat_peak_data(db: Session = Depends(get_db), req: PeaksRequest = None):
    data, records_sum = crud.get_treat_peaks(db, req)
    treat_peak_data = TableDataResponse(data=data, records_sum=records_sum)
    return treat_peak_data


@app.post("/get_treat_peak_data_range", response_model=PeaksDataRange)
def get_treat_peak_data_range(db: Session = Depends(get_db), id: int = Body(embed=True)):
    data = crud.get_treat_peaks_numerical_range(db, id)
    return data


@app.get("/api/get_samplesource_enum", response_model=EnumDataResp)
def get_samplesource_enum(db: Session = Depends(get_db)):
    data = crud.get_datasource_enum(db)
    datasource_enum_list = [d[0] for d in data]
    enum_resp = EnumDataResp(data=datasource_enum_list)
    return enum_resp

@app.get("/api/get_tissue_enum", response_model=EnumDataResp)
def get_tissue_enum(db: Session = Depends(get_db),sample_source:str=''):

    data = crud.get_tissue_enum(db,sample_source)
    datasource_enum_list = [d[0] for d in data]
    enum_resp = EnumDataResp(data=datasource_enum_list)
    return enum_resp

@app.get("/api/get_celltype_enum", response_model=EnumDataResp)
def get_celltype_enum(db: Session = Depends(get_db),sample_source:str='',tissue:str=''):
    print('get_celltype_enum: ',sample_source,tissue)
    data = crud.get_celltype_enum(db,sample_source,tissue)
    datasource_enum_list = [d[0] for d in data]
    enum_resp = EnumDataResp(data=datasource_enum_list)
    return enum_resp


@app.get("/api/get_sampletype_enum", response_model=EnumDataResp)
def get_sampletype_enum(db: Session = Depends(get_db)):
    data = crud.get_method_enum(db)
    enum_list = [d[0] for d in data]
    enum_resp = EnumDataResp(data=enum_list)
    return enum_resp


@app.get("/get_statistics_info", response_model=StatisticsResponse)
def get_statistics_info(db: Session = Depends(get_db)):
    data = crud.get_statistics_info(db)
    return data


'''
download 系列
'''


@app.get("/download_overall_data")
async def download_overall_data(db: Session = Depends(get_db)):
    paging = Paging(start=0, length=sys.maxsize)
    data, records_sum = crud.get_overall_data(db, None, paging)
    overall_data = TableDataResponse(data=data, records_sum=records_sum)
    overall_data = [i.dict() for i in overall_data.data]

    cur_time = datetime.now().replace(microsecond=0)
    file = "../data/overall_data.csv"
    df = pandas.DataFrame(overall_data)
    df.columns = ["id", "pb_gene", "pb_ensembl", "n_sample", "celline", "method", "conditions", "accession",
                  "datasource"]
    df.to_csv(file, index=False)
    return FileResponse(file, filename=f"overall_data_{cur_time}.csv")


@app.post("/download_diff_peak_data")
async def download_diff_peak_data(db: Session = Depends(get_db), req: DownloadRequest = None):
    req.paging = Paging(start=0, length=sys.maxsize)
    data, records_sum = crud.get_diff_peaks(db, req)
    diff_peak_data = TableDataResponse(data=data, records_sum=records_sum)
    diff_peak_data = [item.dict() for item in diff_peak_data.data]

    path = f"../data/{req.id}/csv"
    if not os.path.exists(path):
        os.makedirs(path)

    file = f"{path}/id-{req.id}_{req.pb_gene}_{req.celline}_diff_peak_data.csv"
    df = pandas.DataFrame(diff_peak_data)
    df.columns = ["peak_id", "chr", "start", "end", "width", "conc", "conc_ctrl", "conc_treat", "fold", "p_value",
                  "FDR", "annotation", "gene_chr", "gene_start", "gene_end", "gene_len", "gene_strand", "gene_id",
                  "transcript_id", "distance_to_tss"]
    df.to_csv(file, index=False)


@app.post("/download_diff_GO_enrichment_data")
async def download_diff_GO_enrichment_data(db: Session = Depends(get_db), req: DownloadRequest = None):
    req.paging = Paging(start=0, length=sys.maxsize)
    data, records_sum = crud.get_diff_go_enrich(db, req)
    diff_GO_enrichment_data = TableDataResponse(data=data, records_sum=records_sum)
    diff_GO_enrichment_data = [item.dict() for item in diff_GO_enrichment_data.data]

    path = f"../data/{req.id}/csv"
    if not os.path.exists(path):
        os.makedirs(path)

    file = f"{path}/id-{req.id}_{req.pb_gene}_{req.celline}_diff_peaks_GO_enrich.csv"
    df = pandas.DataFrame(diff_GO_enrichment_data)
    df.columns = ["GO_id", "description", "gene_ratio", "bg_ratio", "p_value", "p_adjust", "q_value", "gene_ids",
                  "count"]
    df.to_csv(file, index=False)


@app.post("/download_diff_KEGG_enrichment_data")
async def download_diff_KEGG_enrichment_data(db: Session = Depends(get_db), req: DownloadRequest = None):
    req.paging = Paging(start=0, length=sys.maxsize)
    data, records_sum = crud.get_diff_kegg_enrich(db, req)
    diff_KEGG_enrichment_data = TableDataResponse(data=data, records_sum=records_sum)
    diff_KEGG_enrichment_data = [item.dict() for item in diff_KEGG_enrichment_data.data]

    path = f"../data/{req.id}/csv"
    if not os.path.exists(path):
        os.makedirs(path)

    file = f"{path}/id-{req.id}_{req.pb_gene}_{req.celline}_diff_peaks_KEGG_enrich.csv"
    df = pandas.DataFrame(diff_KEGG_enrichment_data)
    df.columns = ["KEGG_id", "description", "gene_ratio", "bg_ratio", "p_value", "p_adjust", "q_value",
                  "gene_ids", "count"]
    df.to_csv(file, index=False)


@app.post("/download_diff_motif_data")
async def download_diff_motif_data(db: Session = Depends(get_db), req: DownloadRequest = None):
    req.paging = Paging(start=0, length=sys.maxsize)
    data, records_sum = crud.get_diff_motif(db, req)
    diff_motif_data = TableDataResponse(data=data, records_sum=records_sum)
    diff_motif_data = [item.dict() for item in diff_motif_data.data]

    path = f"../data/{req.id}/csv"
    if not os.path.exists(path):
        os.makedirs(path)

    file = f"{path}/id-{req.id}_{req.pb_gene}_{req.celline}_diff_motif.csv"
    df = pandas.DataFrame(diff_motif_data)
    df.columns = ["rank", "motif", "name", "consensus", "p_value", "log_p_value", "q_value", "target_sequences",
                  "target_sequences_percent", "bg_sequences", "bg_sequences_percent"]
    df = df.drop(columns="motif")
    df.to_csv(file, index=False)


@app.post("/download_diff_footprint_data")
async def download_diff_footprint_data(db: Session = Depends(get_db), req: DownloadRequest = None):
    req.paging = Paging(start=0, length=sys.maxsize)
    data, records_sum = crud.get_diff_footprint(db, req)
    diff_footprint_data = TableDataResponse(data=data, records_sum=records_sum)
    diff_footprint_data = [item.dict() for item in diff_footprint_data.data]

    path = f"../data/{req.id}/csv"
    if not os.path.exists(path):
        os.makedirs(path)

    file = f"{path}/id-{req.id}_{req.pb_gene}_{req.celline}_diff_footprint.csv"
    df = pandas.DataFrame(diff_footprint_data)
    df.columns = ["motif", "tf", "num", "protection_score_ctrl", "protection_score_treat", "tc_ctrl", "tc_treat",
                  "tf_activity", "z_score", "p_value"]
    df.to_csv(file, index=False)


@app.post("/download_ctrl_peak_data")
async def download_ctrl_peak_data(db: Session = Depends(get_db), req: DownloadRequest = None):
    req.paging = Paging(start=0, length=sys.maxsize)
    data, records_sum = crud.get_ctrl_peaks(db, req)
    ctrl_peak_data = TableDataResponse(data=data, records_sum=records_sum)
    ctrl_peak_data = [item.dict() for item in ctrl_peak_data.data]

    path = f"../data/{req.id}/csv"
    if not os.path.exists(path):
        os.makedirs(path)

    file = f"{path}/id-{req.id}_{req.pb_gene}_{req.celline}_ctrl_peaks.csv"
    df = pandas.DataFrame(ctrl_peak_data)
    df.columns = ["peak_id", "chr", "start", "end", "width", "score", "signal_value",
                  "log_p_value", "log_q_value", "peak_offset", "annotation", "gene_chr", "gene_start", "gene_end",
                  "gene_len", "gene_strand", "gene_id", "transcript_id", "distance_to_tss"]
    df.to_csv(file, index=False)


@app.post("/download_treat_peak_data")
async def download_treat_peak_data(db: Session = Depends(get_db), req: DownloadRequest = None):
    req.paging = Paging(start=0, length=sys.maxsize)
    data, records_sum = crud.get_treat_peaks(db, req)
    treat_peak_data = TableDataResponse(data=data, records_sum=records_sum)
    treat_peak_data = [item.dict() for item in treat_peak_data.data]

    path = f"../data/{req.id}/csv"
    if not os.path.exists(path):
        os.makedirs(path)

    file = f"{path}/id-{req.id}_{req.pb_gene}_{req.celline}_treat_peaks.csv"
    df = pandas.DataFrame(treat_peak_data)
    df.columns = ["peak_id", "chr", "start", "end", "width", "score", "signal_value",
                  "log_p_value", "log_q_value", "peak_offset", "annotation", "gene_chr", "gene_start", "gene_end",
                  "gene_len", "gene_strand", "gene_id", "transcript_id", "distance_to_tss"]
    df.to_csv(file, index=False)


@app.post("/transfer_files")
async def transfer_files(db: Session = Depends(get_db), req: DownloadRequest = None):
    # 获取footprint数据
    req.paging = Paging(start=0, length=sys.maxsize)
    data, records_sum = crud.get_diff_footprint(db, req)
    diff_footprint_data = TableDataResponse(data=data, records_sum=records_sum)
    diff_footprint_motifs = [item.dict()["motif"][:8] for item in diff_footprint_data.data]

    # 连接到服务器
    host = '81.70.41.9'
    username = 'wshao'
    password = 'sw6813329'
    port = 2222

    coon = paramiko.SSHClient()
    coon.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    coon.load_system_host_keys()
    coon.connect(host, username=username, password=password, port=port)

    path_plots = f'../data/{req.id}/plots'
    path_footprint_lineplots = f'../data/{req.id}/footprint_lineplots'
    path_narrowPeak = f'../data/{req.id}/narrowPeak'
    path_DBA = f'../data/{req.id}/DBA'
    path_csv = f'../data/{req.id}/csv'
    if not os.path.exists(path_plots):
        os.makedirs(path_plots)
    if not os.path.exists(path_footprint_lineplots):
        os.makedirs(path_footprint_lineplots)
    if not os.path.exists(path_narrowPeak):
        os.makedirs(path_narrowPeak)
    if not os.path.exists(path_DBA):
        os.makedirs(path_DBA)
    if not os.path.exists(path_csv):
        os.makedirs(path_csv)

    # 传输静态资源
    with SCPClient(coon.get_transport()) as scp:
        scp.get(f'/data/cabins/wshao/ATAC-seq/data/{req.id}_{req.pb_gene}_{req.celline}/plots', f'../data/{req.id}',
                recursive=True)
        if len(diff_footprint_motifs) != 0:
            scp.get(
                f'/data/cabins/wshao/ATAC-seq/data/{req.id}_{req.pb_gene}_{req.celline}/footprint/diff/differential_statistics.png',
                path_plots)

            for motif_id in diff_footprint_motifs:
                png_name = f'{motif_id}.png'
                scp.get(
                    f'/data/cabins/wshao/ATAC-seq/data/{req.id}_{req.pb_gene}_{req.celline}/footprint/diff/Lineplots/pngs/{png_name}',
                    path_footprint_lineplots)

        scp.get(f'/data/cabins/wshao/ATAC-seq/data/{req.id}_{req.pb_gene}_{req.celline}/peak/narrowPeaks.tar.gz',
                path_narrowPeak)
        scp.get(f'/data/cabins/wshao/ATAC-seq/data/{req.id}_{req.pb_gene}_{req.celline}/csv/ctrl_peaks.tar.gz',
                path_csv)
        scp.get(f'/data/cabins/wshao/ATAC-seq/data/{req.id}_{req.pb_gene}_{req.celline}/csv/treat_peaks.tar.gz',
                path_csv)
        try:
            scp.get(f'/data/cabins/wshao/ATAC-seq/data/{req.id}_{req.pb_gene}_{req.celline}/DBA/DBA_obj.tar.gz',
                    path_DBA)
        except SCPException:
            print('文件传输错误')

    # 赋权
    cmd = ['chmod', "-R", "775", f"../data/{req.id}"]
    subprocess.run(cmd, universal_newlines=True, stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, shell=False)


@app.post("/transfer_files_fanrui")
async def transfer_files_fanrui(db: Session = Depends(get_db), req: DownloadRequest = None):
    # 获取footprint数据
    req.paging = Paging(start=0, length=sys.maxsize)
    data, records_sum = crud.get_diff_footprint(db, req)
    diff_footprint_data = TableDataResponse(data=data, records_sum=records_sum)
    diff_footprint_motifs = [item.dict()["motif"][:8] for item in diff_footprint_data.data]

    # 连接到服务器
    host = '10.10.1.208'
    username = 'fanrui'
    password = '2Qp-D$b8x$1qaz'
    port = 22

    coon = paramiko.SSHClient()
    coon.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    coon.load_system_host_keys()
    coon.connect(host, username=username, password=password, port=port)

    path_plots = f'../data/{req.id}/plots'
    path_footprint_lineplots = f'../data/{req.id}/footprint_lineplots'
    path_narrowPeak = f'../data/{req.id}/narrowPeak'
    path_DBA = f'../data/{req.id}/DBA'
    path_csv = f'../data/{req.id}/csv'
    if not os.path.exists(path_plots):
        os.makedirs(path_plots)
    if not os.path.exists(path_footprint_lineplots):
        os.makedirs(path_footprint_lineplots)
    if not os.path.exists(path_narrowPeak):
        os.makedirs(path_narrowPeak)
    if not os.path.exists(path_DBA):
        os.makedirs(path_DBA)
    if not os.path.exists(path_csv):
        os.makedirs(path_csv)

    # 传输静态资源
    with SCPClient(coon.get_transport()) as scp:
        scp.get(f'/home/fanrui/wshao_tmp/ATAC-seq/data/{req.id}_{req.pb_gene}_{req.celline}/plots', f'../data/{req.id}',
                recursive=True)
        if len(diff_footprint_motifs) != 0:
            scp.get(
                f'/home/fanrui/wshao_tmp/ATAC-seq/data/{req.id}_{req.pb_gene}_{req.celline}/footprint/diff/differential_statistics.png',
                path_plots)

            for motif_id in diff_footprint_motifs:
                png_name = f'{motif_id}.png'
                scp.get(
                    f'/home/fanrui/wshao_tmp/ATAC-seq/data/{req.id}_{req.pb_gene}_{req.celline}/footprint/diff/Lineplots/pngs/{png_name}',
                    path_footprint_lineplots)

        scp.get(f'/home/fanrui/wshao_tmp/ATAC-seq/data/{req.id}_{req.pb_gene}_{req.celline}/peak/narrowPeaks.tar.gz',
                path_narrowPeak)
        scp.get(f'/home/fanrui/wshao_tmp/ATAC-seq/data/{req.id}_{req.pb_gene}_{req.celline}/csv/ctrl_peaks.tar.gz',
                path_csv)
        scp.get(f'/home/fanrui/wshao_tmp/ATAC-seq/data/{req.id}_{req.pb_gene}_{req.celline}/csv/treat_peaks.tar.gz',
                path_csv)
        try:
            scp.get(f'/home/fanrui/wshao_tmp/ATAC-seq/data/{req.id}_{req.pb_gene}_{req.celline}/DBA/DBA_obj.tar.gz',
                    path_DBA)
        except SCPException:
            print('文件传输错误')

    # 赋权
    cmd = ['chmod', "-R", "775", f"../data/{req.id}"]
    subprocess.run(cmd, universal_newlines=True, stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, shell=False)