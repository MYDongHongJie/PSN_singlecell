#!/usr/bin/env python3
# encoding: utf-8
import sys
import pymysql

#######################
# print("传入参数的总长度为：", len(sys.argv))
# print("function name:", sys.argv[0])
#######################

def update_sql(ddh , status = 'finish' ):
    db = pymysql.connect(
    host = '192.168.12.165',
    port = 3306,
    user = 'singlecell',
    password = '440F-da395b[',
    database = 'single_cell',
    charset = 'utf8'
    )
    cur = db.cursor()
    sql = f'''
    UPDATE single_cell.sc_02 
    SET status_01 = '{status}' , update_by = 'scrna'
    WHERE task_code = '{ddh}' ;
    '''
    cur.execute(sql)
    db.commit()
    cur.close()
    db.close()
######################
ddh = sys.argv[0]
if __name__ == '__main__':
    update_sql(ddh)
