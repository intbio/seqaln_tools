# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 18:50:45 2013

@author: alexeyshaytan
"""
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import *
from datetime import datetime

Base = declarative_base()

class PDB_UNP(Base):
    """"""
    __tablename__ = "pdb_unp"
 
    id = Column(Integer, primary_key=True)
    pdb_id = Column(String(4))
    unp_id = Column(String(6))
    unp_num = Column(Integer)
 
    def __init__(self,pdb_id, unp_id, unp_num):
        """"""
        self.pdb_id = pdb_id
        self.unp_id = unp_id
        self.unp_num = unp_num 

class History(Base):
    """"""
    __tablename__ = "history"
 
    id = Column(Integer, primary_key=True)
    dbname = Column(String)
    dbdescr=Column(String)
    histrec = Column(String)
    time=Column(DateTime)
 
    def __init__(self,dbname,dbdescr,histrec,time=datetime.now()):
        """"""
        self.dbname =dbname
        self.histrec =histrec
        self.dbdescr=dbdescr
        self.time=time
