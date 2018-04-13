def parsedx(filename) :
    """PARSEDX
    ---------------------------------------------------------------------------
    DESCRIPTION
    Parsedx accepts an Excel file containing strings recorded from a Dillon 
    EDXtreme dynamometer and returns data parsed into separate datetime, 
    elapsed seconds, instantaneous force, and peak force columns. The strings
    are printed by the Dillon EDXtreme in continuous recording mode in #4 
    format with the dynamometer connected to a computer by a serial port 
    connection and WedgeLink software. 
    
    INPUTS
    filename: string - filepath of Excel file containing original strings
    
    OUTPUTS
    B: DataFrame - m x 4 DataFrame containing datetime, elapsed seconds,
    instantaneous force, and peak force records in separate columns
    ---------------------------------------------------------------------------    
    """
    import pandas as pd
    import re
    import matplotlib.pyplot as plt
    df=pd.read_excel(filename,sheet_name=0,header=None)
    df.columns = ['strings']
    dtstrings = [re.search('\d*\s\w*\s\d{4}[\,]\d*\:\d{2}\:\d{2}',string) for string in df.strings]
    dt = [pd.to_datetime(dtstring.group(),format='%d %b %Y,%H:%M:%S') for dtstring in dtstrings]
    t = [(dts - dt[0]).total_seconds() for dts in dt]
    fstrings = [re.findall('\d*\.\d{1}',string) for string in df.strings]
    f_ins = []
    f_pk = []
    for b in fstrings :
        f_ins.append(pd.to_numeric(b[0]))
        f_pk.append(pd.to_numeric(b[1]))
    A = pd.DataFrame([dt,t,f_ins,f_pk]).transpose()
    A.columns = ['Datetime','Seconds','F_ins','F_pk']
    A.to_excel(filename,sheetname=2,index=0)
    return A
    A['F_ins'].plot(x='Seconds')
    plt.ylabel('Force, F (kg)')
    plt.xlabel('Time, t (sec)')
    plt.show()