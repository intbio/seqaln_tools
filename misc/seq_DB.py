# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 18:50:45 2013

@author: alexeyshaytan
"""
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import *
from datetime import datetime
from sqlalchemy.orm import relationship, backref

Base = declarative_base()

class Prot_Fams(Base):
    """"""
    __tablename__ = "prot_fams"
 
    id = Column(Integer, primary_key=True)
    name = Column(String)
 
    def __init__(self,name):
        """"""
        self.name = name


class Proteins(Base):
    """"""
    __tablename__ = "proteins"
 
    id = Column(Integer, primary_key=True)
    unp_id = Column(String(6))
    descr = Column(String)
    organism = Column(String)
    fam_id = Column(Integer, ForeignKey("prot_fams.id"))
    family = relationship("Prot_Fams", backref=backref("proteins", order_by=id))
 
    
 
    def __init__(self,unp_id, descr, organism):
        """"""
        self.unp_id = unp_id
        self.descr = descr
        self.organism=organism

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
