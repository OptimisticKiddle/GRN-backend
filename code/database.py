from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import pymysql

user = "root"
password = "sw6813329"
host = "43.143.155.140"
port = 3306
database = "atac_seq_db"

SQLALCHEMY_DATABASE_URL = f'mysql+pymysql://{user}:{password}@{host}:{port}/{database}'

engine = create_engine(SQLALCHEMY_DATABASE_URL)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()
