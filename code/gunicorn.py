# 设置守护进程
daemon = True
# 监听内网端口8000
bind = '0.0.0.0:8000'
# 设置进程文件目录
chdir = '/data/GRN-backend/code'  # 工作目录
# 工作模式
worker_class = 'uvicorn.workers.UvicornWorker'
# 并行工作进程数 核心数*2+1个
workers = 4  # multiprocessing.cpu_count()+1
# 指定每个工作者的线程数
threads = 2
# 设置最大并发量
worker_connections = 2000
loglevel = 'debug'  # 错误日志的日志级别
access_log_format = '%(t)s %(p)s %(h)s "%(r)s" %(s)s %(L)s %(b)s %(f)s" "%(a)s"'
# 设置访问日志和错误信息日志路径
accesslog = "/data/GRN-backend/code/gunicorn_access.log"
errorlog = "/data/GRN-backend/code/gunicorn_error.log"

timeout = 360000  # 设置超时时间
